from e3nn import o3
import torch
from torch import nn
from torch.nn import functional as F
from torch_cluster import radius, radius_graph
from torch_scatter import scatter, scatter_mean
import numpy as np
from e3nn.nn import BatchNorm

from models.layers import GaussianSmearing, OldAtomEncoder, AtomEncoder
from models.tensor_layers import OldTensorProductConvLayer
from utils import so3, torus
from datasets.process_mols import lig_feature_dims, rec_residue_feature_dims, rec_atom_feature_dims

AGGREGATORS = {"mean": lambda x: torch.mean(x, dim=1),
               "max": lambda x: torch.max(x, dim=1)[0],
               "min": lambda x: torch.min(x, dim=1)[0],
               "std": lambda x: torch.std(x, dim=1)}


class AAOldModel(torch.nn.Module):
    def __init__(self, t_to_sigma, device, timestep_emb_func, in_lig_edge_features=4, sigma_embed_dim=32, sh_lmax=2,
                 ns=16, nv=4, num_conv_layers=2, lig_max_radius=5, rec_max_radius=30, cross_max_distance=250,
                 center_max_distance=30, distance_embed_dim=32, cross_distance_embed_dim=32, no_torsion=False,
                 scale_by_sigma=True, norm_by_sigma=True, use_second_order_repr=False, batch_norm=True,
                 dynamic_max_cross=False, dropout=0.0, smooth_edges=False, odd_parity=False,
                 separate_noise_schedule=False, lm_embedding_type=False, confidence_mode=False,
                 confidence_dropout=0, confidence_no_batchnorm = False,
                 asyncronous_noise_schedule=False, affinity_prediction=False, parallel=1,
                 parallel_aggregators="mean max min std", num_confidence_outputs=1, fixed_center_conv=False,
                 no_aminoacid_identities=False, include_miscellaneous_atoms=False, use_old_atom_encoder=False):
        super(AAOldModel, self).__init__()
        assert (not no_aminoacid_identities) or (lm_embedding_type is None), "no language model emb without identities"
        if parallel > 1: assert affinity_prediction
        self.t_to_sigma = t_to_sigma
        self.in_lig_edge_features = in_lig_edge_features
        sigma_embed_dim *= (3 if separate_noise_schedule else 1)
        self.sigma_embed_dim = sigma_embed_dim
        self.lig_max_radius = lig_max_radius
        self.rec_max_radius = rec_max_radius
        self.cross_max_distance = cross_max_distance
        self.dynamic_max_cross = dynamic_max_cross
        self.center_max_distance = center_max_distance
        self.distance_embed_dim = distance_embed_dim
        self.cross_distance_embed_dim = cross_distance_embed_dim
        self.sh_irreps = o3.Irreps.spherical_harmonics(lmax=sh_lmax)
        self.ns, self.nv = ns, nv
        self.scale_by_sigma = scale_by_sigma
        self.norm_by_sigma = norm_by_sigma
        self.device = device
        self.no_torsion = no_torsion
        self.smooth_edges = smooth_edges
        self.odd_parity = odd_parity
        self.num_conv_layers = num_conv_layers
        self.timestep_emb_func = timestep_emb_func
        self.separate_noise_schedule = separate_noise_schedule
        self.confidence_mode = confidence_mode
        self.num_conv_layers = num_conv_layers
        self.asyncronous_noise_schedule = asyncronous_noise_schedule
        self.affinity_prediction = affinity_prediction
        self.parallel, self.parallel_aggregators = parallel, parallel_aggregators.split(' ')
        self.fixed_center_conv = fixed_center_conv
        self.no_aminoacid_identities = no_aminoacid_identities

        # embedding layers
        atom_encoder_class = OldAtomEncoder if use_old_atom_encoder else AtomEncoder
        self.lig_node_embedding = atom_encoder_class(emb_dim=ns, feature_dims=lig_feature_dims, sigma_embed_dim=sigma_embed_dim)
        self.lig_edge_embedding = nn.Sequential(nn.Linear(in_lig_edge_features + sigma_embed_dim + distance_embed_dim, ns),nn.ReLU(),nn.Dropout(dropout),nn.Linear(ns, ns))

        self.rec_node_embedding = atom_encoder_class(emb_dim=ns, feature_dims=rec_residue_feature_dims, sigma_embed_dim=sigma_embed_dim, lm_embedding_type=lm_embedding_type)
        self.rec_edge_embedding = nn.Sequential(nn.Linear(sigma_embed_dim + distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))

        self.atom_node_embedding = atom_encoder_class(emb_dim=ns, feature_dims=rec_atom_feature_dims, sigma_embed_dim=sigma_embed_dim)
        self.atom_edge_embedding = nn.Sequential(nn.Linear(sigma_embed_dim + distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))

        self.lr_edge_embedding = nn.Sequential(nn.Linear(sigma_embed_dim + cross_distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))
        self.ar_edge_embedding = nn.Sequential(nn.Linear(sigma_embed_dim + distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))
        self.la_edge_embedding = nn.Sequential(nn.Linear(sigma_embed_dim + cross_distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))

        self.lig_distance_expansion = GaussianSmearing(0.0, lig_max_radius, distance_embed_dim)
        self.rec_distance_expansion = GaussianSmearing(0.0, rec_max_radius, distance_embed_dim)
        self.cross_distance_expansion = GaussianSmearing(0.0, cross_max_distance, cross_distance_embed_dim)

        if use_second_order_repr:
            irrep_seq = [
                f'{ns}x0e',
                f'{ns}x0e + {nv}x1o + {nv}x2e',
                f'{ns}x0e + {nv}x1o + {nv}x2e + {nv}x1e + {nv}x2o',
                f'{ns}x0e + {nv}x1o + {nv}x2e + {nv}x1e + {nv}x2o + {ns}x0o'
            ]
        else:
            irrep_seq = [
                f'{ns}x0e',
                f'{ns}x0e + {nv}x1o',
                f'{ns}x0e + {nv}x1o + {nv}x1e',
                f'{ns}x0e + {nv}x1o + {nv}x1e + {ns}x0o'
            ]

        # convolutional layers
        conv_layers = []
        for i in range(num_conv_layers):
            in_irreps = irrep_seq[min(i, len(irrep_seq) - 1)]
            out_irreps = irrep_seq[min(i + 1, len(irrep_seq) - 1)]
            parameters = {
                'in_irreps': in_irreps,
                'sh_irreps': self.sh_irreps,
                'out_irreps': out_irreps,
                'n_edge_features': 3 * ns,
                'residual': False,
                'batch_norm': batch_norm,
                'dropout': dropout
            }

            for _ in range(9): # 3 intra & 6 inter per each layer
                conv_layers.append(OldTensorProductConvLayer(**parameters))

        self.conv_layers = nn.ModuleList(conv_layers)

        # confidence and affinity prediction layers
        if self.confidence_mode:
            if self.affinity_prediction:
                if self.parallel > 1:
                    output_confidence_dim = 1 + ns
                else:
                    output_confidence_dim = num_confidence_outputs +1
            else:
                output_confidence_dim = num_confidence_outputs

            self.confidence_predictor = nn.Sequential(
                nn.Linear(2 * self.ns if num_conv_layers >= 3 else self.ns, ns),
                nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                nn.ReLU(),
                nn.Dropout(confidence_dropout),
                nn.Linear(ns, ns),
                nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                nn.ReLU(),
                nn.Dropout(confidence_dropout),
                nn.Linear(ns, output_confidence_dim)
            )

            if self.parallel > 1:
                self.affinity_predictor = nn.Sequential(
                    nn.Linear(len(self.parallel_aggregators) * ns, ns),
                    nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                    nn.ReLU(),
                    nn.Dropout(confidence_dropout),
                    nn.Linear(ns, ns),
                    nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                    nn.ReLU(),
                    nn.Dropout(confidence_dropout),
                    nn.Linear(ns, 1)
                )

        else:
            # convolution for translational and rotational scores
            self.center_distance_expansion = GaussianSmearing(0.0, center_max_distance, distance_embed_dim)
            self.center_edge_embedding = nn.Sequential(
                nn.Linear(distance_embed_dim + sigma_embed_dim, ns),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(ns, ns)
            )

            self.final_conv = OldTensorProductConvLayer(
                in_irreps=self.conv_layers[-1].out_irreps,
                sh_irreps=self.sh_irreps,
                out_irreps=f'2x1o + 2x1e' if not self.odd_parity else '1x1o + 1x1e',
                n_edge_features=2 * ns,
                residual=False,
                dropout=dropout,
                batch_norm=batch_norm
            )

            self.tr_final_layer = nn.Sequential(nn.Linear(1 + sigma_embed_dim, ns),nn.Dropout(dropout), nn.ReLU(), nn.Linear(ns, 1))
            self.rot_final_layer = nn.Sequential(nn.Linear(1 + sigma_embed_dim, ns),nn.Dropout(dropout), nn.ReLU(), nn.Linear(ns, 1))

            if not no_torsion:
                # convolution for torsional score
                self.final_edge_embedding = nn.Sequential(
                    nn.Linear(distance_embed_dim, ns),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                    nn.Linear(ns, ns)
                )
                self.final_tp_tor = o3.FullTensorProduct(self.sh_irreps, "2e")
                self.tor_bond_conv = OldTensorProductConvLayer(
                    in_irreps=self.conv_layers[-1].out_irreps,
                    sh_irreps=self.final_tp_tor.irreps_out,
                    out_irreps=f'{ns}x0o + {ns}x0e' if not self.odd_parity else f'{ns}x0o',
                    n_edge_features=3 * ns,
                    residual=False,
                    dropout=dropout,
                    batch_norm=batch_norm
                )
                self.tor_final_layer = nn.Sequential(
                    nn.Linear(2 * ns if not self.odd_parity else ns, ns, bias=False),
                    nn.Tanh(),
                    nn.Dropout(dropout),
                    nn.Linear(ns, 1, bias=False)
                )

    def forward(self, data):
        if self.no_aminoacid_identities:
            data['receptor'].x = data['receptor'].x * 0

        if not self.confidence_mode:
            tr_sigma, rot_sigma, tor_sigma = self.t_to_sigma(*[data.complex_t[noise_type] for noise_type in ['tr', 'rot', 'tor']])
        else:
            tr_sigma, rot_sigma, tor_sigma = [data.complex_t[noise_type] for noise_type in ['tr', 'rot', 'tor']]

        # build ligand graph
        lig_node_attr, lig_edge_index, lig_edge_attr, lig_edge_sh, lig_edge_weight = self.build_lig_conv_graph(data)
        lig_node_attr = self.lig_node_embedding(lig_node_attr)
        lig_edge_attr = self.lig_edge_embedding(lig_edge_attr)

        # build receptor graph
        rec_node_attr, rec_edge_index, rec_edge_attr, rec_edge_sh, rec_edge_weight = self.build_rec_conv_graph(data)
        rec_node_attr = self.rec_node_embedding(rec_node_attr)
        rec_edge_attr = self.rec_edge_embedding(rec_edge_attr)

        # build atom graph
        atom_node_attr, atom_edge_index, atom_edge_attr, atom_edge_sh, atom_edge_weight = self.build_atom_conv_graph(data)
        atom_node_attr = self.atom_node_embedding(atom_node_attr)
        atom_edge_attr = self.atom_edge_embedding(atom_edge_attr)

        # build cross graph
        cross_cutoff = (tr_sigma * 3 + 20).unsqueeze(1) if self.dynamic_max_cross else self.cross_max_distance
        lr_edge_index, lr_edge_attr, lr_edge_sh, lr_edge_weight, la_edge_index, la_edge_attr, \
            la_edge_sh, la_edge_weight, ar_edge_index, ar_edge_attr, ar_edge_sh, ar_edge_weight = \
            self.build_cross_conv_graph(data, cross_cutoff)
        lr_edge_attr= self.lr_edge_embedding(lr_edge_attr)
        la_edge_attr = self.la_edge_embedding(la_edge_attr)
        ar_edge_attr = self.ar_edge_embedding(ar_edge_attr)

        for l in range(self.num_conv_layers):
            # LIGAND updates
            lig_edge_attr_ = torch.cat([lig_edge_attr, lig_node_attr[lig_edge_index[0], :self.ns], lig_node_attr[lig_edge_index[1], :self.ns]], -1)
            lig_update = self.conv_layers[9*l](lig_node_attr, lig_edge_index, lig_edge_attr_, lig_edge_sh, edge_weight=lig_edge_weight)

            lr_edge_attr_ = torch.cat([lr_edge_attr, lig_node_attr[lr_edge_index[0], :self.ns], rec_node_attr[lr_edge_index[1], :self.ns]], -1)
            lr_update = self.conv_layers[9*l+1](rec_node_attr, lr_edge_index, lr_edge_attr_, lr_edge_sh,
                                                out_nodes=lig_node_attr.shape[0], edge_weight=lr_edge_weight)

            la_edge_attr_ = torch.cat([la_edge_attr, lig_node_attr[la_edge_index[0], :self.ns], atom_node_attr[la_edge_index[1], :self.ns]], -1)
            la_update = self.conv_layers[9*l+2](atom_node_attr, la_edge_index, la_edge_attr_, la_edge_sh,
                                                out_nodes=lig_node_attr.shape[0], edge_weight=la_edge_weight)

            if l != self.num_conv_layers-1:  # last layer optimisation

                # ATOM UPDATES
                atom_edge_attr_ = torch.cat([atom_edge_attr, atom_node_attr[atom_edge_index[0], :self.ns], atom_node_attr[atom_edge_index[1], :self.ns]], -1)
                atom_update = self.conv_layers[9*l+3](atom_node_attr, atom_edge_index, atom_edge_attr_, atom_edge_sh, edge_weight=atom_edge_weight)

                al_edge_attr_ = torch.cat([la_edge_attr, atom_node_attr[la_edge_index[1], :self.ns], lig_node_attr[la_edge_index[0], :self.ns]], -1)
                al_update = self.conv_layers[9*l+4](lig_node_attr, torch.flip(la_edge_index, dims=[0]), al_edge_attr_,
                                                    la_edge_sh, out_nodes=atom_node_attr.shape[0], edge_weight=la_edge_weight)

                ar_edge_attr_ = torch.cat([ar_edge_attr, atom_node_attr[ar_edge_index[0], :self.ns], rec_node_attr[ar_edge_index[1], :self.ns]],-1)
                ar_update = self.conv_layers[9*l+5](rec_node_attr, ar_edge_index, ar_edge_attr_, ar_edge_sh, out_nodes=atom_node_attr.shape[0], edge_weight=ar_edge_weight)

                # RECEPTOR updates
                rec_edge_attr_ = torch.cat([rec_edge_attr, rec_node_attr[rec_edge_index[0], :self.ns], rec_node_attr[rec_edge_index[1], :self.ns]], -1)
                rec_update = self.conv_layers[9*l+6](rec_node_attr, rec_edge_index, rec_edge_attr_, rec_edge_sh, edge_weight=rec_edge_weight)

                rl_edge_attr_ = torch.cat([lr_edge_attr, rec_node_attr[lr_edge_index[1], :self.ns], lig_node_attr[lr_edge_index[0], :self.ns]], -1)
                rl_update = self.conv_layers[9*l+7](lig_node_attr, torch.flip(lr_edge_index, dims=[0]), rl_edge_attr_,
                                                    lr_edge_sh, out_nodes=rec_node_attr.shape[0], edge_weight=lr_edge_weight)

                ra_edge_attr_ = torch.cat([ar_edge_attr, rec_node_attr[ar_edge_index[1], :self.ns], atom_node_attr[ar_edge_index[0], :self.ns]], -1)
                ra_update = self.conv_layers[9*l+8](atom_node_attr, torch.flip(ar_edge_index, dims=[0]), ra_edge_attr_,
                                                    ar_edge_sh, out_nodes=rec_node_attr.shape[0], edge_weight=ar_edge_weight)

            # padding original features and update features with residual updates
            lig_node_attr = F.pad(lig_node_attr, (0, lig_update.shape[-1] - lig_node_attr.shape[-1]))
            lig_node_attr = lig_node_attr + lig_update + la_update + lr_update

            if l != self.num_conv_layers - 1:  # last layer optimisation
                atom_node_attr = F.pad(atom_node_attr, (0, atom_update.shape[-1] - atom_node_attr.shape[-1]))
                atom_node_attr = atom_node_attr + atom_update + al_update + ar_update
                rec_node_attr = F.pad(rec_node_attr, (0, rec_update.shape[-1] - rec_node_attr.shape[-1]))
                rec_node_attr = rec_node_attr + rec_update + ra_update + rl_update

        # confidence and affinity prediction
        if self.confidence_mode:
            scalar_lig_attr = torch.cat([lig_node_attr[:,:self.ns],lig_node_attr[:,-self.ns:]], dim=1) if self.num_conv_layers >= 3 else lig_node_attr[:,:self.ns]
            confidence = self.confidence_predictor(scatter_mean(scalar_lig_attr, data['ligand'].batch if self.parallel == 1 else data['ligand'].batch_parallel, dim=0)).squeeze(dim=-1)

            if self.parallel > 1:
                confidence, affinity = confidence[:, 0], confidence[:, 1:]
                confidence = confidence.reshape(data.num_graphs, self.parallel)
                affinity = affinity.reshape(data.num_graphs, self.parallel, -1)
                affinity = torch.cat([AGGREGATORS[agg](affinity) for agg in self.parallel_aggregators], dim=-1)
                affinity = self.affinity_predictor(affinity).squeeze(dim=-1)
                confidence = confidence, affinity
            return confidence
        assert self.parallel == 1

        # compute translational and rotational score vectors
        center_edge_index, center_edge_attr, center_edge_sh = self.build_center_conv_graph(data)
        center_edge_attr = self.center_edge_embedding(center_edge_attr)
        if self.fixed_center_conv:
            center_edge_attr = torch.cat([center_edge_attr, lig_node_attr[center_edge_index[1], :self.ns]], -1)
        else:
            center_edge_attr = torch.cat([center_edge_attr, lig_node_attr[center_edge_index[0], :self.ns]], -1)
        global_pred = self.final_conv(lig_node_attr, center_edge_index, center_edge_attr, center_edge_sh, out_nodes=data.num_graphs)

        tr_pred = global_pred[:, :3] + (global_pred[:, 6:9] if not self.odd_parity else 0)
        rot_pred = global_pred[:, 3:6] + (global_pred[:, 9:] if not self.odd_parity else 0)

        if self.separate_noise_schedule:
            data.graph_sigma_emb = torch.cat([self.timestep_emb_func(data.complex_t[noise_type]) for noise_type in ['tr', 'rot', 'tor']], dim=1)
        elif self.asyncronous_noise_schedule:
            data.graph_sigma_emb = self.timestep_emb_func(data.complex_t['t'])
        else:  # tr rot and tor noise is all the same in this case
            data.graph_sigma_emb = self.timestep_emb_func(data.complex_t['tr'])

        # adjust the magniture of the score vectors
        tr_norm = torch.linalg.vector_norm(tr_pred, dim=1).unsqueeze(1)
        tr_pred = tr_pred / tr_norm * self.tr_final_layer(torch.cat([tr_norm, data.graph_sigma_emb], dim=1))

        rot_norm = torch.linalg.vector_norm(rot_pred, dim=1).unsqueeze(1)
        rot_pred = rot_pred / rot_norm * self.rot_final_layer(torch.cat([rot_norm, data.graph_sigma_emb], dim=1))

        if self.scale_by_sigma:
            tr_pred = tr_pred / tr_sigma.unsqueeze(1)
            rot_pred = rot_pred * so3.score_norm(rot_sigma.cpu()).unsqueeze(1).to(data['ligand'].x.device)

        if self.no_torsion or data['ligand'].edge_mask.sum() == 0: return tr_pred, rot_pred, torch.empty(0,device=self.device)

        # torsional components
        tor_bonds, tor_edge_index, tor_edge_attr, tor_edge_sh, tor_edge_weight = self.build_bond_conv_graph(data)
        tor_bond_vec = data['ligand'].pos[tor_bonds[1]] - data['ligand'].pos[tor_bonds[0]]
        tor_bond_attr = lig_node_attr[tor_bonds[0]] + lig_node_attr[tor_bonds[1]]

        tor_bonds_sh = o3.spherical_harmonics("2e", tor_bond_vec, normalize=True, normalization='component')
        tor_edge_sh = self.final_tp_tor(tor_edge_sh, tor_bonds_sh[tor_edge_index[0]])

        tor_edge_attr = torch.cat([tor_edge_attr, lig_node_attr[tor_edge_index[1], :self.ns],
                                   tor_bond_attr[tor_edge_index[0], :self.ns]], -1)
        tor_pred = self.tor_bond_conv(lig_node_attr, tor_edge_index, tor_edge_attr, tor_edge_sh,
                                  out_nodes=data['ligand'].edge_mask.sum(), reduce='mean', edge_weight=tor_edge_weight)
        tor_pred = self.tor_final_layer(tor_pred).squeeze(1)
        edge_sigma = tor_sigma[data['ligand'].batch][data['ligand', 'ligand'].edge_index[0]][data['ligand'].edge_mask]

        if self.scale_by_sigma:
            tor_pred = tor_pred * torch.sqrt(torch.tensor(torus.score_norm(edge_sigma.cpu().numpy())).float()
                                             .to(data['ligand'].x.device))
        return tr_pred, rot_pred, tor_pred

    def get_edge_weight(self, edge_vec, max_norm):
        if self.smooth_edges:
            normalised_norm = torch.clip(edge_vec.norm(dim=-1) * np.pi / max_norm, max=np.pi)
            return 0.5 * (torch.cos(normalised_norm) + 1.0).unsqueeze(-1)
        return 1.0

    def build_lig_conv_graph(self, data):
        # build the graph between ligand atoms
        if self.separate_noise_schedule:
            data['ligand'].node_sigma_emb = torch.cat(
                [self.timestep_emb_func(data['ligand'].node_t[noise_type]) for noise_type in ['tr', 'rot', 'tor']],
                dim=1)
        elif self.asyncronous_noise_schedule:
            data['ligand'].node_sigma_emb = self.timestep_emb_func(data['ligand'].node_t['t'])
        else:
            data['ligand'].node_sigma_emb = self.timestep_emb_func(
                data['ligand'].node_t['tr'])  # tr rot and tor noise is all the same

        if self.parallel == 1:
            radius_edges = radius_graph(data['ligand'].pos, self.lig_max_radius, data['ligand'].batch)
        else:
            batches = torch.zeros(data.num_graphs, device=data['ligand'].x.device).long()
            batches = batches.index_add(0, data['ligand'].batch, torch.ones(len(data['ligand'].batch), device=data['ligand'].x.device).long())
            outer_batches = data.num_graphs
            b = [torch.ones(batches[i].item()//self.parallel, device=data['ligand'].x.device).long() * (self.parallel * i + j)
                 for i in range(outer_batches) for j in range(self.parallel)]
            data['ligand'].batch_parallel = torch.cat(b)
            radius_edges = radius_graph(data['ligand'].pos, self.lig_max_radius, data['ligand'].batch_parallel)
        edge_index = torch.cat([data['ligand', 'ligand'].edge_index, radius_edges], 1).long()
        edge_attr = torch.cat([
            data['ligand', 'ligand'].edge_attr,
            torch.zeros(radius_edges.shape[-1], self.in_lig_edge_features, device=data['ligand'].x.device)
        ], 0)

        edge_sigma_emb = data['ligand'].node_sigma_emb[edge_index[0].long()]
        edge_attr = torch.cat([edge_attr, edge_sigma_emb], 1)
        node_attr = torch.cat([data['ligand'].x, data['ligand'].node_sigma_emb], 1)

        src, dst = edge_index
        edge_vec = data['ligand'].pos[dst.long()] - data['ligand'].pos[src.long()]
        edge_length_emb = self.lig_distance_expansion(edge_vec.norm(dim=-1))

        edge_attr = torch.cat([edge_attr, edge_length_emb], 1)
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')
        edge_weight = self.get_edge_weight(edge_vec, self.lig_max_radius)

        return node_attr, edge_index, edge_attr, edge_sh, edge_weight

    def build_rec_conv_graph(self, data):
        # build the graph between receptor residues
        if self.separate_noise_schedule:
            data['receptor'].node_sigma_emb = torch.cat(
                [self.timestep_emb_func(data['receptor'].node_t[noise_type]) for noise_type in ['tr', 'rot', 'tor']],
                dim=1)
        elif self.asyncronous_noise_schedule:
            data['receptor'].node_sigma_emb = self.timestep_emb_func(data['receptor'].node_t['t'])
        else:
            data['receptor'].node_sigma_emb = self.timestep_emb_func(
                data['receptor'].node_t['tr'])  # tr rot and tor noise is all the same
        node_attr = torch.cat([data['receptor'].x, data['receptor'].node_sigma_emb], 1)

        # this assumes the edges were already created in preprocessing since protein's structure is fixed
        edge_index = data['receptor', 'receptor'].edge_index
        src, dst = edge_index
        edge_vec = data['receptor'].pos[dst.long()] - data['receptor'].pos[src.long()]
        #assert torch.all(edge_vec.norm(dim=-1) < self.rec_max_radius)

        edge_length_emb = self.rec_distance_expansion(edge_vec.norm(dim=-1))
        edge_sigma_emb = data['receptor'].node_sigma_emb[edge_index[0].long()]
        edge_attr = torch.cat([edge_sigma_emb, edge_length_emb], 1)
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')
        edge_weight = self.get_edge_weight(edge_vec, self.rec_max_radius)

        return node_attr, edge_index, edge_attr, edge_sh, edge_weight

    def build_atom_conv_graph(self, data):
        # build the graph between receptor atoms
        if self.separate_noise_schedule:
            data['atom'].node_sigma_emb = torch.cat([self.timestep_emb_func(data['atom'].node_t[noise_type]) for noise_type in ['tr', 'rot', 'tor']],dim=1)
        elif self.asyncronous_noise_schedule:
            data['atom'].node_sigma_emb = self.timestep_emb_func(data['atom'].node_t['t'])
        else:
            data['atom'].node_sigma_emb = self.timestep_emb_func(data['atom'].node_t['tr'])  # tr rot and tor noise is all the same
        node_attr = torch.cat([data['atom'].x, data['atom'].node_sigma_emb], 1)

        # this assumes the edges were already created in preprocessing since protein's structure is fixed
        edge_index = data['atom', 'atom'].edge_index
        src, dst = edge_index
        edge_vec = data['atom'].pos[dst.long()] - data['atom'].pos[src.long()]

        edge_length_emb = self.lig_distance_expansion(edge_vec.norm(dim=-1))
        edge_sigma_emb = data['atom'].node_sigma_emb[edge_index[0].long()]
        edge_attr = torch.cat([edge_sigma_emb, edge_length_emb], 1)
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')
        edge_weight = self.get_edge_weight(edge_vec, self.lig_max_radius)

        return node_attr, edge_index, edge_attr, edge_sh, edge_weight

    def build_cross_conv_graph(self, data, lr_cross_distance_cutoff):
        # build the cross edges between ligan atoms, receptor residues and receptor atoms

        # LIGAND to RECEPTOR
        if torch.is_tensor(lr_cross_distance_cutoff):
            # different cutoff for every graph
            lr_edge_index = radius(data['receptor'].pos / lr_cross_distance_cutoff[data['receptor'].batch],
                                data['ligand'].pos / lr_cross_distance_cutoff[data['ligand'].batch], 1,
                                data['receptor'].batch, data['ligand'].batch, max_num_neighbors=10000)
        else:
            lr_edge_index = radius(data['receptor'].pos, data['ligand'].pos, lr_cross_distance_cutoff,
                            data['receptor'].batch, data['ligand'].batch, max_num_neighbors=10000)

        lr_edge_vec = data['receptor'].pos[lr_edge_index[1].long()] - data['ligand'].pos[lr_edge_index[0].long()]
        lr_edge_length_emb = self.cross_distance_expansion(lr_edge_vec.norm(dim=-1))
        lr_edge_sigma_emb = data['ligand'].node_sigma_emb[lr_edge_index[0].long()]
        lr_edge_attr = torch.cat([lr_edge_sigma_emb, lr_edge_length_emb], 1)
        lr_edge_sh = o3.spherical_harmonics(self.sh_irreps, lr_edge_vec, normalize=True, normalization='component')

        cutoff_d = lr_cross_distance_cutoff[data['ligand'].batch[lr_edge_index[0]]].squeeze() \
            if torch.is_tensor(lr_cross_distance_cutoff) else lr_cross_distance_cutoff
        lr_edge_weight = self.get_edge_weight(lr_edge_vec, cutoff_d)

        # LIGAND to ATOM
        la_edge_index = radius(data['atom'].pos, data['ligand'].pos, self.lig_max_radius,
                               data['atom'].batch, data['ligand'].batch, max_num_neighbors=10000)

        la_edge_vec = data['atom'].pos[la_edge_index[1].long()] - data['ligand'].pos[la_edge_index[0].long()]
        la_edge_length_emb = self.cross_distance_expansion(la_edge_vec.norm(dim=-1))
        la_edge_sigma_emb = data['ligand'].node_sigma_emb[la_edge_index[0].long()]
        la_edge_attr = torch.cat([la_edge_sigma_emb, la_edge_length_emb], 1)
        la_edge_sh = o3.spherical_harmonics(self.sh_irreps, la_edge_vec, normalize=True, normalization='component')
        la_edge_weight = self.get_edge_weight(la_edge_vec, self.lig_max_radius)

        # ATOM to RECEPTOR
        ar_edge_index = data['atom', 'receptor'].edge_index
        ar_edge_vec = data['receptor'].pos[ar_edge_index[1].long()] - data['atom'].pos[ar_edge_index[0].long()]
        ar_edge_length_emb = self.rec_distance_expansion(ar_edge_vec.norm(dim=-1))
        ar_edge_sigma_emb = data['atom'].node_sigma_emb[ar_edge_index[0].long()]
        ar_edge_attr = torch.cat([ar_edge_sigma_emb, ar_edge_length_emb], 1)
        ar_edge_sh = o3.spherical_harmonics(self.sh_irreps, ar_edge_vec, normalize=True, normalization='component')
        ar_edge_weight = 1

        return lr_edge_index, lr_edge_attr, lr_edge_sh, lr_edge_weight, la_edge_index, la_edge_attr, \
               la_edge_sh, la_edge_weight, ar_edge_index, ar_edge_attr, ar_edge_sh, ar_edge_weight

    def build_center_conv_graph(self, data):
        # build the filter for the convolution of the center with the ligand atoms
        # for translational and rotational score
        edge_index = torch.cat([data['ligand'].batch.unsqueeze(0), torch.arange(len(data['ligand'].batch)).to(data['ligand'].x.device).unsqueeze(0)], dim=0)

        center_pos, count = torch.zeros((data.num_graphs, 3)).to(data['ligand'].x.device), torch.zeros((data.num_graphs, 3)).to(data['ligand'].x.device)
        center_pos.index_add_(0, index=data['ligand'].batch, source=data['ligand'].pos)
        center_pos = center_pos / torch.bincount(data['ligand'].batch).unsqueeze(1)

        edge_vec = data['ligand'].pos[edge_index[1]] - center_pos[edge_index[0]]
        edge_attr = self.center_distance_expansion(edge_vec.norm(dim=-1))
        edge_sigma_emb = data['ligand'].node_sigma_emb[edge_index[1].long()]
        edge_attr = torch.cat([edge_attr, edge_sigma_emb], 1)
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')
        return edge_index, edge_attr, edge_sh

    def build_bond_conv_graph(self, data):
        # build graph for the pseudotorque layer
        bonds = data['ligand', 'ligand'].edge_index[:, data['ligand'].edge_mask].long()
        bond_pos = (data['ligand'].pos[bonds[0]] + data['ligand'].pos[bonds[1]]) / 2
        bond_batch = data['ligand'].batch[bonds[0]]
        edge_index = radius(data['ligand'].pos, bond_pos, self.lig_max_radius, batch_x=data['ligand'].batch, batch_y=bond_batch)

        edge_vec = data['ligand'].pos[edge_index[1]] - bond_pos[edge_index[0]]
        edge_attr = self.lig_distance_expansion(edge_vec.norm(dim=-1))

        edge_attr = self.final_edge_embedding(edge_attr)
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')
        edge_weight = self.get_edge_weight(edge_vec, self.lig_max_radius)

        return bonds, edge_index, edge_attr, edge_sh, edge_weight

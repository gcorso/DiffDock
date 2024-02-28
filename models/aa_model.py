from e3nn import o3
import torch
from esm.pretrained import load_model_and_alphabet
from torch import nn
from torch.nn import functional as F
from torch_cluster import radius, radius_graph
from torch_geometric.utils import subgraph
from torch_scatter import scatter_mean
import numpy as np

from models.layers import GaussianSmearing, AtomEncoder
from models.tensor_layers import get_irrep_seq, TensorProductConvLayer
from utils import so3, torus
from datasets.process_mols import lig_feature_dims, rec_residue_feature_dims, rec_atom_feature_dims

AGGREGATORS = {"mean": lambda x: torch.mean(x, dim=1),
               "max": lambda x: torch.max(x, dim=1)[0],
               "min": lambda x: torch.min(x, dim=1)[0],
               "std": lambda x: torch.std(x, dim=1)}


class AAModel(torch.nn.Module):
    def __init__(self, t_to_sigma, device, timestep_emb_func, in_lig_edge_features=4, sigma_embed_dim=32, sh_lmax=2,
                 ns=16, nv=4, num_conv_layers=2, lig_max_radius=5, rec_max_radius=30, cross_max_distance=250,
                 center_max_distance=30, distance_embed_dim=32, cross_distance_embed_dim=32, no_torsion=False,
                 scale_by_sigma=True, norm_by_sigma=True, use_second_order_repr=False, batch_norm=True,
                 dynamic_max_cross=False, dropout=0.0, smooth_edges=False, odd_parity=False,
                 separate_noise_schedule=False, lm_embedding_type=False, confidence_mode=False,
                 confidence_dropout=0, confidence_no_batchnorm = False,
                 asyncronous_noise_schedule=False, affinity_prediction=False, parallel=1,
                 parallel_aggregators="mean max min std", num_confidence_outputs=1, atom_num_confidence_outputs=1, fixed_center_conv=False,
                 no_aminoacid_identities=False, include_miscellaneous_atoms=False,
                 differentiate_convolutions=True, tp_weights_layers=2, num_prot_emb_layers=0,
                 reduce_pseudoscalars=False, embed_also_ligand=False, atom_confidence=False, sidechain_pred=False,
                 depthwise_convolution=False, crop_beyond=None):
        super(AAModel, self).__init__()
        assert (not no_aminoacid_identities) or (lm_embedding_type is None), "no language model emb without identities"
        assert not sidechain_pred, "sidechain prediction not implemented/makes sense for all atom model"
        assert not depthwise_convolution, "depthwise convolution not implemented for all atom model"
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
        self.num_prot_emb_layers = num_prot_emb_layers
        self.asyncronous_noise_schedule = asyncronous_noise_schedule
        self.affinity_prediction = affinity_prediction
        self.parallel, self.parallel_aggregators = parallel, parallel_aggregators.split(' ')
        self.fixed_center_conv = fixed_center_conv
        self.no_aminoacid_identities = no_aminoacid_identities
        self.differentiate_convolutions = differentiate_convolutions
        self.reduce_pseudoscalars = reduce_pseudoscalars
        self.atom_confidence = atom_confidence
        self.atom_num_confidence_outputs = atom_num_confidence_outputs
        self.crop_beyond = crop_beyond

        self.lm_embedding_type = lm_embedding_type
        if lm_embedding_type is None:
            lm_embedding_dim = 0
        elif lm_embedding_type == "precomputed":
            lm_embedding_dim=1280
        else:
            lm, alphabet = load_model_and_alphabet(lm_embedding_type)
            self.batch_converter = alphabet.get_batch_converter()
            lm.lm_head = torch.nn.Identity()
            lm.contact_head = torch.nn.Identity()
            lm_embedding_dim = lm.embed_dim
            self.lm = lm

        # embedding layers
        atom_encoder_class = AtomEncoder
        self.lig_node_embedding = atom_encoder_class(emb_dim=ns, feature_dims=lig_feature_dims, sigma_embed_dim=sigma_embed_dim)
        self.lig_edge_embedding = nn.Sequential(nn.Linear(in_lig_edge_features + sigma_embed_dim + distance_embed_dim, ns),nn.ReLU(),nn.Dropout(dropout),nn.Linear(ns, ns))

        self.rec_sigma_embedding = nn.Sequential(nn.Linear(sigma_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout), nn.Linear(ns, ns))
        self.rec_node_embedding = atom_encoder_class(emb_dim=ns, feature_dims=rec_residue_feature_dims, sigma_embed_dim=0, lm_embedding_dim=lm_embedding_dim)
        self.rec_edge_embedding = nn.Sequential(nn.Linear(distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout), nn.Linear(ns, ns))
        self.atom_node_embedding = atom_encoder_class(emb_dim=ns, feature_dims=rec_atom_feature_dims, sigma_embed_dim=0)
        self.atom_edge_embedding = nn.Sequential(nn.Linear(distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout), nn.Linear(ns, ns))

        self.lr_edge_embedding = nn.Sequential(nn.Linear(sigma_embed_dim + cross_distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))
        self.ar_edge_embedding = nn.Sequential(nn.Linear(distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))
        self.la_edge_embedding = nn.Sequential(nn.Linear(sigma_embed_dim + cross_distance_embed_dim, ns), nn.ReLU(), nn.Dropout(dropout),nn.Linear(ns, ns))

        self.lig_distance_expansion = GaussianSmearing(0.0, lig_max_radius, distance_embed_dim)
        self.rec_distance_expansion = GaussianSmearing(0.0, rec_max_radius, distance_embed_dim)
        self.cross_distance_expansion = GaussianSmearing(0.0, cross_max_distance, cross_distance_embed_dim)

        irrep_seq = get_irrep_seq(ns, nv, use_second_order_repr, reduce_pseudoscalars)
        assert not include_miscellaneous_atoms, "currently not supported"

        rec_emb_layers = []
        for i in range(num_prot_emb_layers):
            in_irreps = irrep_seq[min(i, len(irrep_seq) - 1)]
            out_irreps = irrep_seq[min(i + 1, len(irrep_seq) - 1)]
            layer = TensorProductConvLayer(
                in_irreps=in_irreps,
                sh_irreps=self.sh_irreps,
                out_irreps=out_irreps,
                n_edge_features=3 * ns,
                hidden_features=3 * ns,
                residual=True,
                batch_norm=batch_norm,
                dropout=dropout,
                faster=sh_lmax == 1 and not use_second_order_repr,
                tp_weights_layers=tp_weights_layers,
                edge_groups=1 if not differentiate_convolutions else 4,
            )
            rec_emb_layers.append(layer)
        self.rec_emb_layers = nn.ModuleList(rec_emb_layers)

        self.embed_also_ligand = embed_also_ligand
        if embed_also_ligand:
            lig_emb_layers = []
            for i in range(num_prot_emb_layers):
                in_irreps = irrep_seq[min(i, len(irrep_seq) - 1)]
                out_irreps = irrep_seq[min(i + 1, len(irrep_seq) - 1)]
                layer = TensorProductConvLayer(
                    in_irreps=in_irreps,
                    sh_irreps=self.sh_irreps,
                    out_irreps=out_irreps,
                    n_edge_features=3 * ns,
                    hidden_features=3 * ns,
                    residual=True,
                    batch_norm=batch_norm,
                    dropout=dropout,
                    faster=sh_lmax == 1 and not use_second_order_repr,
                    tp_weights_layers=tp_weights_layers,
                    edge_groups=1,
                )
                lig_emb_layers.append(layer)
            self.lig_emb_layers = nn.ModuleList(lig_emb_layers)

        # convolutional layers
        conv_layers = []
        for i in range(num_prot_emb_layers, num_prot_emb_layers + num_conv_layers):
            in_irreps = irrep_seq[min(i, len(irrep_seq) - 1)]
            out_irreps = irrep_seq[min(i + 1, len(irrep_seq) - 1)]
            layer = TensorProductConvLayer(
                in_irreps=in_irreps,
                sh_irreps=self.sh_irreps,
                out_irreps=out_irreps,
                n_edge_features=3 * ns,
                hidden_features=3 * ns,
                residual=True,
                batch_norm=batch_norm,
                dropout=dropout,
                faster=sh_lmax == 1 and not use_second_order_repr,
                tp_weights_layers=tp_weights_layers,
                edge_groups=1 if not differentiate_convolutions else (3 if i == num_prot_emb_layers + num_conv_layers - 1 else 9),
            )
            conv_layers.append(layer)
        self.conv_layers = nn.ModuleList(conv_layers)

        # confidence and affinity prediction layers
        if self.confidence_mode:
            if self.affinity_prediction:
                if self.parallel > 1:
                    output_confidence_dim = 1 + ns
                else:
                    output_confidence_dim = num_confidence_outputs + 1
            else:
                output_confidence_dim = num_confidence_outputs

            input_size = ns + (nv if reduce_pseudoscalars else ns) if num_conv_layers + num_prot_emb_layers >= 3 else ns

            if self.atom_confidence:
                self.atom_confidence_predictor = nn.Sequential(
                    nn.Linear(input_size, ns),
                    nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                    nn.ReLU(),
                    nn.Dropout(confidence_dropout),
                    nn.Linear(ns, ns),
                    nn.BatchNorm1d(ns) if not confidence_no_batchnorm else nn.Identity(),
                    nn.ReLU(),
                    nn.Dropout(confidence_dropout),
                    nn.Linear(ns, atom_num_confidence_outputs + ns)
                )
                input_size = ns

            self.confidence_predictor = nn.Sequential(
                nn.Linear(input_size, ns),
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

            self.final_conv = TensorProductConvLayer(
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
                self.tor_bond_conv = TensorProductConvLayer(
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

    def embedding(self, data):
        if not hasattr(data['receptor'], "rec_node_attr"):
            if self.lm_embedding_type not in [None, 'precomputed']:
                sequences = [s for l in data['receptor'].sequence for s in l]
                if isinstance(sequences[0], list):
                    sequences = [s for l in sequences for s in l]
                sequences = [(i, s) for i, s in enumerate(sequences)]
                batch_labels, batch_strs, batch_tokens = self.batch_converter(sequences)
                out = self.lm(batch_tokens.to(data['receptor'].x.device), repr_layers=[self.lm.num_layers], return_contacts=False)
                rec_lm_emb = torch.cat([t[:len(sequences[i][1])] for i, t in enumerate(out['representations'][self.lm.num_layers])], dim=0)
                data['receptor'].x = torch.cat([data['receptor'].x, rec_lm_emb], dim=-1)

            rec_node_attr, rec_edge_attr, rec_edge_sh, rec_edge_weight = self.build_rec_conv_graph(data)
            rec_node_attr = self.rec_node_embedding(rec_node_attr)
            rec_edge_attr = self.rec_edge_embedding(rec_edge_attr)

            atom_node_attr, atom_edge_attr, atom_edge_sh, atom_edge_weight = self.build_atom_conv_graph(data)
            atom_node_attr = self.atom_node_embedding(atom_node_attr)
            atom_edge_attr = self.atom_edge_embedding(atom_edge_attr)

            ar_edge_attr, ar_edge_sh, ar_edge_weight = self.build_cross_rec_conv_graph(data)
            ar_edge_attr = self.ar_edge_embedding(ar_edge_attr)

            rec_edge_index = data['receptor', 'receptor'].edge_index.clone()
            atom_edge_index = data['atom', 'atom'].edge_index.clone()
            ar_edge_index = data['atom', 'receptor'].edge_index.clone()

            node_attr = torch.cat([rec_node_attr, atom_node_attr], dim=0)
            ar_edge_index[0] = ar_edge_index[0] + len(rec_node_attr)
            edge_index = torch.cat([rec_edge_index, ar_edge_index, atom_edge_index + len(rec_node_attr), torch.flip(ar_edge_index, dims=[0])], dim=1)
            edge_attr = torch.cat([rec_edge_attr, ar_edge_attr, atom_edge_attr, ar_edge_attr], dim=0)
            edge_sh = torch.cat([rec_edge_sh, ar_edge_sh, atom_edge_sh, ar_edge_sh], dim=0)
            edge_weight = torch.cat([rec_edge_weight, ar_edge_weight, atom_edge_weight, ar_edge_weight], dim=0) \
                if torch.is_tensor(rec_edge_weight) else torch.ones((len(edge_index[0]), 1), device=edge_index.device)
            s1, s2, s3 = len(rec_edge_index[0]), len(rec_edge_index[0]) + len(ar_edge_index[0]), len(rec_edge_index[0]) + len(ar_edge_index[0]) + len(atom_edge_index[0])

            for l in range(len(self.rec_emb_layers)):
                edge_attr_ = torch.cat(
                    [edge_attr, node_attr[edge_index[0], :self.ns], node_attr[edge_index[1], :self.ns]], -1)
                if self.differentiate_convolutions: edge_attr_ = [edge_attr_[:s1], edge_attr_[s1:s2], edge_attr_[s2:s3], edge_attr_[s3:]]
                node_attr = self.rec_emb_layers[l](node_attr, edge_index, edge_attr_, edge_sh, edge_weight=edge_weight)


            data['receptor'].rec_node_attr = node_attr[:len(rec_node_attr)]
            data['receptor', 'receptor'].rec_edge_attr = rec_edge_attr
            data['receptor', 'receptor'].edge_sh = rec_edge_sh
            data['receptor', 'receptor'].edge_weight = rec_edge_weight

            data['atom'].atom_node_attr = node_attr[len(rec_node_attr):]
            data['atom', 'atom'].atom_edge_attr = atom_edge_attr
            data['atom', 'atom'].edge_sh = atom_edge_sh
            data['atom', 'atom'].edge_weight = atom_edge_weight

            data['atom', 'receptor'].edge_attr = ar_edge_attr
            data['atom', 'receptor'].edge_sh = ar_edge_sh
            data['atom', 'receptor'].edge_weight = ar_edge_weight

        # receptor embedding
        rec_sigma_emb = self.rec_sigma_embedding(self.timestep_emb_func(data.complex_t['tr']))
        rec_node_attr = data['receptor'].rec_node_attr + 0
        rec_node_attr[:, :self.ns] = rec_node_attr[:, :self.ns] + rec_sigma_emb[data['receptor'].batch]
        rec_edge_attr = data['receptor', 'receptor'].rec_edge_attr + rec_sigma_emb[data['receptor'].batch[data['receptor', 'receptor'].edge_index[0]]]

        # atom embedding
        atom_node_attr = data['atom'].atom_node_attr + 0
        atom_node_attr[:, :self.ns] = atom_node_attr[:, :self.ns] + rec_sigma_emb[data['atom'].batch]
        atom_edge_attr = data['atom', 'atom'].atom_edge_attr + rec_sigma_emb[data['atom'].batch[data['atom', 'atom'].edge_index[0]]]

        # atom-receptor embedding
        ar_edge_attr = data['atom', 'receptor'].edge_attr + rec_sigma_emb[data['atom'].batch[data['atom', 'receptor'].edge_index[0]]]

        # ligand embedding
        lig_node_attr, lig_edge_index, lig_edge_attr, lig_edge_sh, lig_edge_weight = self.build_lig_conv_graph(data)
        lig_node_attr = self.lig_node_embedding(lig_node_attr)
        lig_edge_attr = self.lig_edge_embedding(lig_edge_attr)

        if self.embed_also_ligand:
            for l in range(len(self.lig_emb_layers)):
                edge_attr_ = torch.cat([lig_edge_attr, lig_node_attr[lig_edge_index[0], :self.ns], lig_node_attr[lig_edge_index[1], :self.ns]], -1)
                lig_node_attr = self.lig_emb_layers[l](lig_node_attr, lig_edge_index, edge_attr_, lig_edge_sh, edge_weight=lig_edge_weight)

        else:
            lig_node_attr = F.pad(lig_node_attr, (0, rec_node_attr.shape[-1] - lig_node_attr.shape[-1]))

        return lig_node_attr, lig_edge_index, lig_edge_attr, lig_edge_sh, lig_edge_weight, \
               rec_node_attr, data['receptor', 'receptor'].edge_index, rec_edge_attr, data['receptor', 'receptor'].edge_sh, data['receptor', 'receptor'].edge_weight, \
               atom_node_attr, data['atom', 'atom'].edge_index, atom_edge_attr, data['atom', 'atom'].edge_sh, data['atom', 'atom'].edge_weight, \
               data['atom', 'receptor'].edge_index, ar_edge_attr, data['atom', 'receptor'].edge_sh, data['atom', 'receptor'].edge_weight

    def forward(self, data):
        if self.crop_beyond is not None:
            # TODO missing filtering atoms
            raise NotImplementedError
            ligand_pos = data['ligand'].pos
            receptor_pos = data['receptor'].pos
            residues_to_keep = torch.any(torch.sum((ligand_pos.unsqueeze(0) - receptor_pos.unsqueeze(1)) ** 2, -1) < self.crop_beyond ** 2, dim=1)

            data['receptor'].pos = data['receptor'].pos[residues_to_keep]
            data['receptor'].x = data['receptor'].x[residues_to_keep]
            data['receptor'].side_chain_vecs = data['receptor'].side_chain_vecs[residues_to_keep]
            data['receptor', 'rec_contact', 'receptor'].edge_index = subgraph(residues_to_keep, data['receptor', 'rec_contact', 'receptor'].edge_index, relabel_nodes=True)[0]

        if self.no_aminoacid_identities:
            data['receptor'].x = data['receptor'].x * 0

        if not self.confidence_mode:
            tr_sigma, rot_sigma, tor_sigma = self.t_to_sigma(*[data.complex_t[noise_type] for noise_type in ['tr', 'rot', 'tor']])
        else:
            tr_sigma, rot_sigma, tor_sigma = [data.complex_t[noise_type] for noise_type in ['tr', 'rot', 'tor']]

        lig_node_attr, lig_edge_index, lig_edge_attr, lig_edge_sh, lig_edge_weight, rec_node_attr, \
            rec_edge_index, rec_edge_attr, rec_edge_sh, rec_edge_weight,\
            atom_node_attr, atom_edge_index, atom_edge_attr, atom_edge_sh, atom_edge_weight, \
            ar_edge_index, ar_edge_attr, ar_edge_sh, ar_edge_weight = self.embedding(data)

        # build lig cross graph
        cross_cutoff = (tr_sigma * 3 + 20).unsqueeze(1) if self.dynamic_max_cross else self.cross_max_distance
        lr_edge_index, lr_edge_attr, lr_edge_sh, lr_edge_weight, la_edge_index, la_edge_attr, \
            la_edge_sh, la_edge_weight = self.build_cross_lig_conv_graph(data, cross_cutoff)
        lr_edge_attr= self.lr_edge_embedding(lr_edge_attr)
        la_edge_attr = self.la_edge_embedding(la_edge_attr)

        n_lig, n_rec = len(lig_node_attr), len(rec_node_attr)

        node_attr = torch.cat([lig_node_attr, rec_node_attr, atom_node_attr], dim=0)
        rec_edge_index, atom_edge_index, lr_edge_index, la_edge_index, ar_edge_index = rec_edge_index.clone(), atom_edge_index.clone(), lr_edge_index.clone(), la_edge_index.clone(), ar_edge_index.clone()
        rec_edge_index[0], rec_edge_index[1] = rec_edge_index[0] + n_lig, rec_edge_index[1] + n_lig
        atom_edge_index[0], atom_edge_index[1] = atom_edge_index[0] + n_lig + n_rec, atom_edge_index[1] + n_lig + n_rec
        lr_edge_index[1] = lr_edge_index[1] + n_lig
        la_edge_index[1] = la_edge_index[1] + n_lig + n_rec
        ar_edge_index[0], ar_edge_index[1] = ar_edge_index[0] + n_lig + n_rec, ar_edge_index[1] + n_lig

        edge_index = torch.cat([lig_edge_index, lr_edge_index, la_edge_index, rec_edge_index,
                                torch.flip(lr_edge_index, dims=[0]), torch.flip(ar_edge_index, dims=[0]),
                                atom_edge_index, torch.flip(la_edge_index, dims=[0]), ar_edge_index], dim=1)
        edge_attr = torch.cat([lig_edge_attr, lr_edge_attr, la_edge_attr, rec_edge_attr, lr_edge_attr,
                               ar_edge_attr, atom_edge_attr, la_edge_attr, ar_edge_attr], dim=0)
        edge_sh = torch.cat([lig_edge_sh, lr_edge_sh, la_edge_sh, rec_edge_sh, lr_edge_sh, ar_edge_sh,
                             atom_edge_sh, la_edge_sh, ar_edge_sh], dim=0)
        edge_weight = torch.cat([lig_edge_weight, lr_edge_weight, la_edge_weight, rec_edge_weight, lr_edge_weight,
                                 ar_edge_weight, atom_edge_weight, la_edge_weight, ar_edge_weight],
                                dim=0) if torch.is_tensor(lig_edge_weight) else torch.ones((len(edge_index[0]), 1),
                                                                                           device=edge_index.device)
        s1, s2, s3, s4, s5, s6, s7, s8, _ = tuple(np.cumsum(list(map(len, [lig_edge_attr, lr_edge_attr, la_edge_attr,
                                                                           rec_edge_attr, lr_edge_attr, ar_edge_attr, atom_edge_attr, la_edge_attr, ar_edge_attr]))).tolist())

        for l in range(len(self.conv_layers)):
            if l < len(self.conv_layers) - 1:
                edge_attr_ = torch.cat([edge_attr, node_attr[edge_index[0], :self.ns], node_attr[edge_index[1], :self.ns]], -1)
                if self.differentiate_convolutions: edge_attr_ = [edge_attr_[:s1], edge_attr_[s1:s2], edge_attr_[s2:s3], edge_attr_[s3:s4],
                                                                  edge_attr_[s4:s5], edge_attr_[s5:s6], edge_attr_[s6:s7], edge_attr_[s7:s8], edge_attr_[s8:]]
                node_attr = self.conv_layers[l](node_attr, edge_index, edge_attr_, edge_sh, edge_weight=edge_weight)
            else:
                edge_attr_ = torch.cat([edge_attr[:s3], node_attr[edge_index[0, :s3], :self.ns], node_attr[edge_index[1, :s3], :self.ns]], -1)
                if self.differentiate_convolutions: edge_attr_ = [edge_attr_[:s1], edge_attr_[s1:s2], edge_attr_[s2:s3]]
                node_attr = self.conv_layers[l](node_attr, edge_index[:, :s3], edge_attr_, edge_sh[:s3], edge_weight=edge_weight[:s3])

        lig_node_attr = node_attr[:len(lig_node_attr)]

        # confidence and affinity prediction
        if self.confidence_mode:
            scalar_lig_attr = torch.cat([lig_node_attr[:,:self.ns], lig_node_attr[:,-(self.nv if self.reduce_pseudoscalars else self.ns):] ], dim=1) \
                if self.num_conv_layers + self.num_prot_emb_layers >= 3 else lig_node_attr[:,:self.ns]

            if self.atom_confidence:
                scalar_lig_attr = self.atom_confidence_predictor(scalar_lig_attr)
                atom_confidence = scalar_lig_attr[:, :self.atom_num_confidence_outputs]
                scalar_lig_attr = scalar_lig_attr[:, self.atom_num_confidence_outputs:]
            else:
                atom_confidence = torch.zeros((len(lig_node_attr),), device=lig_node_attr.device)

            confidence = self.confidence_predictor(scatter_mean(scalar_lig_attr, data['ligand'].batch, dim=0)).squeeze(dim=-1)

            if self.parallel > 1:
                confidence, affinity = confidence[:, 0], confidence[:, 1:]
                confidence = confidence.reshape(data.num_graphs, self.parallel)
                affinity = affinity.reshape(data.num_graphs, self.parallel, -1)
                affinity = torch.cat([AGGREGATORS[agg](affinity) for agg in self.parallel_aggregators], dim=-1)
                affinity = self.affinity_predictor(affinity).squeeze(dim=-1)
                confidence = confidence, affinity
            return confidence, atom_confidence
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

        if self.no_torsion or data['ligand'].edge_mask.sum() == 0: return tr_pred, rot_pred, torch.empty(0,device=self.device), None

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
        return tr_pred, rot_pred, tor_pred, None

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
        node_attr = data['receptor'].x

        # this assumes the edges were already created in preprocessing since protein's structure is fixed
        edge_index = data['receptor', 'receptor'].edge_index
        src, dst = edge_index
        edge_vec = data['receptor'].pos[dst.long()] - data['receptor'].pos[src.long()]

        edge_attr = self.rec_distance_expansion(edge_vec.norm(dim=-1))
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')
        edge_weight = self.get_edge_weight(edge_vec, self.rec_max_radius)

        return node_attr, edge_attr, edge_sh, edge_weight

    def build_atom_conv_graph(self, data):
        # build the graph between receptor atoms
        node_attr = data['atom'].x

        # this assumes the edges were already created in preprocessing since protein's structure is fixed
        edge_index = data['atom', 'atom'].edge_index
        src, dst = edge_index
        edge_vec = data['atom'].pos[dst.long()] - data['atom'].pos[src.long()]

        edge_attr = self.lig_distance_expansion(edge_vec.norm(dim=-1))
        edge_sh = o3.spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization='component')
        edge_weight = self.get_edge_weight(edge_vec, self.lig_max_radius)

        return node_attr, edge_attr, edge_sh, edge_weight

    def build_cross_lig_conv_graph(self, data, lr_cross_distance_cutoff):
        # build the cross edges between ligand atoms and receptor residues + atoms

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
        la_edge_length_emb = self.lig_distance_expansion(la_edge_vec.norm(dim=-1))
        la_edge_sigma_emb = data['ligand'].node_sigma_emb[la_edge_index[0].long()]
        la_edge_attr = torch.cat([la_edge_sigma_emb, la_edge_length_emb], 1)
        la_edge_sh = o3.spherical_harmonics(self.sh_irreps, la_edge_vec, normalize=True, normalization='component')
        la_edge_weight = self.get_edge_weight(la_edge_vec, self.lig_max_radius)

        return lr_edge_index, lr_edge_attr, lr_edge_sh, lr_edge_weight, la_edge_index, la_edge_attr, \
               la_edge_sh, la_edge_weight

    def build_cross_rec_conv_graph(self, data):
        # build the cross edges between ligan atoms, receptor residues and receptor atoms

        # ATOM to RECEPTOR
        ar_edge_index = data['atom', 'receptor'].edge_index
        ar_edge_vec = data['receptor'].pos[ar_edge_index[1].long()] - data['atom'].pos[ar_edge_index[0].long()]
        ar_edge_attr = self.rec_distance_expansion(ar_edge_vec.norm(dim=-1))
        ar_edge_sh = o3.spherical_harmonics(self.sh_irreps, ar_edge_vec, normalize=True, normalization='component')
        ar_edge_weight = 1

        return ar_edge_attr, ar_edge_sh, ar_edge_weight

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

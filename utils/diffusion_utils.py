import math
import numpy as np
import torch
import torch.nn.functional as F
from torch import nn
from scipy.stats import beta

from utils.geometry import axis_angle_to_matrix, rigid_transform_Kabsch_3D_torch, rigid_transform_Kabsch_3D_torch_batch
from utils.torsion import modify_conformer_torsion_angles, modify_conformer_torsion_angles_batch


def sigmoid(t):
    return 1 / (1 + np.e**(-t))


def sigmoid_schedule(t, k=10, m=0.5):
    s = lambda t: sigmoid(k*(t-m))
    return (s(t)-s(0))/(s(1)-s(0))


def t_to_sigma_individual(t, schedule_type, sigma_min, sigma_max, schedule_k=10, schedule_m=0.4):
    if schedule_type == "exponential":
        return sigma_min ** (1 - t) * sigma_max ** t
    elif schedule_type == 'sigmoid':
        return sigmoid_schedule(t, k=schedule_k, m=schedule_m) * (sigma_max - sigma_min) + sigma_min


def t_to_sigma(t_tr, t_rot, t_tor, args):
    tr_sigma = args.tr_sigma_min ** (1-t_tr) * args.tr_sigma_max ** t_tr
    rot_sigma = args.rot_sigma_min ** (1-t_rot) * args.rot_sigma_max ** t_rot
    tor_sigma = args.tor_sigma_min ** (1-t_tor) * args.tor_sigma_max ** t_tor
    return tr_sigma, rot_sigma, tor_sigma


def modify_conformer(data, tr_update, rot_update, torsion_updates, pivot=None):
    lig_center = torch.mean(data['ligand'].pos, dim=0, keepdim=True)
    rot_mat = axis_angle_to_matrix(rot_update.squeeze())
    rigid_new_pos = (data['ligand'].pos - lig_center) @ rot_mat.T + tr_update + lig_center

    if torsion_updates is not None:
        flexible_new_pos = modify_conformer_torsion_angles(rigid_new_pos,
                                                           data['ligand', 'ligand'].edge_index.T[data['ligand'].edge_mask],
                                                           data['ligand'].mask_rotate if isinstance(data['ligand'].mask_rotate, np.ndarray) else data['ligand'].mask_rotate[0],
                                                           torsion_updates).to(rigid_new_pos.device)
        if pivot is None:
            R, t = rigid_transform_Kabsch_3D_torch(flexible_new_pos.T, rigid_new_pos.T)
            aligned_flexible_pos = flexible_new_pos @ R.T + t.T
        else:
            R1, t1 = rigid_transform_Kabsch_3D_torch(pivot.T, rigid_new_pos.T)
            R2, t2 = rigid_transform_Kabsch_3D_torch(flexible_new_pos.T, pivot.T)

            aligned_flexible_pos = (flexible_new_pos @ R2.T + t2.T) @ R1.T + t1.T

        data['ligand'].pos = aligned_flexible_pos
    else:
        data['ligand'].pos = rigid_new_pos
    return data


def modify_conformer_batch(orig_pos, data, tr_update, rot_update, torsion_updates, mask_rotate):
    B = data.num_graphs
    N, M, R = data['ligand'].num_nodes // B, data['ligand', 'ligand'].num_edges // B, data['ligand'].edge_mask.sum().item() // B

    pos, edge_index, edge_mask = orig_pos.reshape(B, N, 3) + 0, data['ligand', 'ligand'].edge_index[:, :M], data['ligand'].edge_mask[:M]
    torsion_updates = torsion_updates.reshape(B, -1) if torsion_updates is not None else None

    lig_center = torch.mean(pos, dim=1, keepdim=True)
    rot_mat = axis_angle_to_matrix(rot_update)
    rigid_new_pos = torch.bmm(pos - lig_center, rot_mat.permute(0, 2, 1)) + tr_update.unsqueeze(1) + lig_center

    if torsion_updates is not None:
        flexible_new_pos = modify_conformer_torsion_angles_batch(rigid_new_pos, edge_index.T[edge_mask], mask_rotate, torsion_updates)
        R, t = rigid_transform_Kabsch_3D_torch_batch(flexible_new_pos, rigid_new_pos)
        aligned_flexible_pos = torch.bmm(flexible_new_pos, R.transpose(1, 2)) + t.transpose(1, 2)
        final_pos = aligned_flexible_pos.reshape(-1, 3)
    else:
        final_pos = rigid_new_pos.reshape(-1, 3)
    return final_pos


def modify_conformer_coordinates(pos, tr_update, rot_update, torsion_updates, edge_mask, mask_rotate, edge_index):
    # Made this function which does the same as modify_conformer because passing a graph would require
    # creating a new heterograph for reach graph when unbatching a batch of graphs
    lig_center = torch.mean(pos, dim=0, keepdim=True)
    rot_mat = axis_angle_to_matrix(rot_update.squeeze())
    rigid_new_pos = (pos - lig_center) @ rot_mat.T + tr_update + lig_center

    if torsion_updates is not None:
        flexible_new_pos = modify_conformer_torsion_angles(rigid_new_pos,edge_index.T[edge_mask],mask_rotate \
            if isinstance(mask_rotate, np.ndarray) else mask_rotate[0], torsion_updates).to(rigid_new_pos.device)

        R, t = rigid_transform_Kabsch_3D_torch(flexible_new_pos.T, rigid_new_pos.T)
        aligned_flexible_pos = flexible_new_pos @ R.T + t.T
        return aligned_flexible_pos
    else:
        return rigid_new_pos


def sinusoidal_embedding(timesteps, embedding_dim, max_positions=10000):
    """ from https://github.com/hojonathanho/diffusion/blob/master/diffusion_tf/nn.py   """
    assert len(timesteps.shape) == 1
    half_dim = embedding_dim // 2
    emb = math.log(max_positions) / (half_dim - 1)
    emb = torch.exp(torch.arange(half_dim, dtype=torch.float32, device=timesteps.device) * -emb)
    emb = timesteps.float()[:, None] * emb[None, :]
    emb = torch.cat([torch.sin(emb), torch.cos(emb)], dim=1)
    if embedding_dim % 2 == 1:  # zero pad
        emb = F.pad(emb, (0, 1), mode='constant')
    assert emb.shape == (timesteps.shape[0], embedding_dim)
    return emb


class GaussianFourierProjection(nn.Module):
    """Gaussian Fourier embeddings for noise levels.
    from https://github.com/yang-song/score_sde_pytorch/blob/1618ddea340f3e4a2ed7852a0694a809775cf8d0/models/layerspp.py#L32
    """

    def __init__(self, embedding_size=256, scale=1.0):
        super().__init__()
        self.W = nn.Parameter(torch.randn(embedding_size//2) * scale, requires_grad=False)

    def forward(self, x):
        x_proj = x[:, None] * self.W[None, :] * 2 * np.pi
        emb = torch.cat([torch.sin(x_proj), torch.cos(x_proj)], dim=-1)
        return emb


def get_timestep_embedding(embedding_type, embedding_dim, embedding_scale=10000):
    if embedding_type == 'sinusoidal':
        emb_func = (lambda x : sinusoidal_embedding(embedding_scale * x, embedding_dim))
    elif embedding_type == 'fourier':
        emb_func = GaussianFourierProjection(embedding_size=embedding_dim, scale=embedding_scale)
    else:
        raise NotImplemented
    return emb_func


def get_t_schedule(sigma_schedule, inference_steps, inf_sched_alpha=1, inf_sched_beta=1, t_max=1):
    if sigma_schedule == 'expbeta':
        lin_max = beta.cdf(t_max, a=inf_sched_alpha, b=inf_sched_beta)
        c = np.linspace(lin_max, 0, inference_steps + 1)[:-1]
        return beta.ppf(c, a=inf_sched_alpha, b=inf_sched_beta)
    raise Exception()


def set_time(complex_graphs, t, t_tr, t_rot, t_tor, batchsize, all_atoms, device, include_miscellaneous_atoms=False):
    complex_graphs['ligand'].node_t = {
        'tr': t_tr * torch.ones(complex_graphs['ligand'].num_nodes).to(device),
        'rot': t_rot * torch.ones(complex_graphs['ligand'].num_nodes).to(device),
        'tor': t_tor * torch.ones(complex_graphs['ligand'].num_nodes).to(device)}
    complex_graphs['receptor'].node_t = {
        'tr': t_tr * torch.ones(complex_graphs['receptor'].num_nodes).to(device),
        'rot': t_rot * torch.ones(complex_graphs['receptor'].num_nodes).to(device),
        'tor': t_tor * torch.ones(complex_graphs['receptor'].num_nodes).to(device)}
    complex_graphs.complex_t = {'tr': t_tr * torch.ones(batchsize).to(device),
                               'rot': t_rot * torch.ones(batchsize).to(device),
                               'tor': t_tor * torch.ones(batchsize).to(device)}
    if all_atoms:
        complex_graphs['atom'].node_t = {
            'tr': t_tr * torch.ones(complex_graphs['atom'].num_nodes).to(device),
            'rot': t_rot * torch.ones(complex_graphs['atom'].num_nodes).to(device),
            'tor': t_tor * torch.ones(complex_graphs['atom'].num_nodes).to(device)}

    if include_miscellaneous_atoms and not all_atoms:
        complex_graphs['misc_atom'].node_t = {
            'tr': t_tr * torch.ones(complex_graphs['misc_atom'].num_nodes).to(device),
            'rot': t_rot * torch.ones(complex_graphs['misc_atom'].num_nodes).to(device),
            'tor': t_tor * torch.ones(complex_graphs['misc_atom'].num_nodes).to(device)}

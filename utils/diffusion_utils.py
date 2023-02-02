import math
import numpy as np
import torch
import torch.nn.functional as F
from torch import nn
from scipy.stats import beta

from utils.geometry import axis_angle_to_matrix, rigid_transform_Kabsch_3D_torch
from utils.torsion import modify_conformer_torsion_angles

"""
    A set of functions and classes for performing various transformations on a 3D molecular structure. 
    The code is written in Python and uses the PyTorch library for tensor operations.

    t_to_sigma: This function takes in three inputs, t_tr, t_rot, and t_tor, as well as an object args 
    that contains parameters, and returns three outputs tr_sigma, rot_sigma, and tor_sigma. It uses 
    the inputs t_tr, t_rot, and t_tor to compute the values of the sigmas, which are the standard 
    deviations of three types of transformations applied to the molecular structure: translation, 
    rotation, and torsion angle.

    modify_conformer: This function takes in a dictionary data that contains information about the 
    molecular structure, as well as three transformations, tr_update, rot_update, and torsion_updates, 
    and returns the modified molecular structure. The function first calculates the center of mass 
    of the ligand and uses it to translate the ligand. Then, it performs rotation using rot_update. 
    If there is a torsion_update, it modifies the torsion angles of the molecular structure using 
    the modify_conformer_torsion_angles function. Finally, it aligns the flexible conformer with 
    the rigid conformer using the rigid_transform_Kabsch_3D_torch function.

    sinusoidal_embedding: This function takes in a tensor timesteps and two parameters, embedding_dim 
    and max_positions, and returns a tensor of sinusoidal embeddings. This function is used to embed 
    the time steps of a molecular simulation in a continuous space, which is useful for learning 
    a smooth representation of the simulation over time.

    GaussianFourierProjection: This class extends the nn.Module class from PyTorch and implements 
    Gaussian Fourier embeddings for noise levels. The class has two parameters, embedding_size and 
    scale, that determine the size and scale of the embeddings. The class has a forward method that 
    computes the embeddings based on a given input.

    get_timestep_embedding: This function takes in three parameters, embedding_type, embedding_dim,
    and embedding_scale, and returns an instance of either sinusoidal_embedding or GaussianFourierProjection, 
    depending on the value of embedding_type.

    A simple example of using the code could be to modify a molecular structure by rotating it 
    around its center of mass. To do this, you would first create an instance of the molecular 
    structure as a dictionary with the necessary information and pass it to the modify_conformer 
    function along with a rotation update. The rotation update can be a tensor with the rotation 
    angles, and torsion_updates could be set to None.
"""

def t_to_sigma(t_tr, t_rot, t_tor, args):
    tr_sigma = args.tr_sigma_min ** (1-t_tr) * args.tr_sigma_max ** t_tr
    rot_sigma = args.rot_sigma_min ** (1-t_rot) * args.rot_sigma_max ** t_rot
    tor_sigma = args.tor_sigma_min ** (1-t_tor) * args.tor_sigma_max ** t_tor
    return tr_sigma, rot_sigma, tor_sigma


def modify_conformer(data, tr_update, rot_update, torsion_updates):
    lig_center = torch.mean(data['ligand'].pos, dim=0, keepdim=True)
    rot_mat = axis_angle_to_matrix(rot_update.squeeze())
    rigid_new_pos = (data['ligand'].pos - lig_center) @ rot_mat.T + tr_update + lig_center

    if torsion_updates is not None:
        flexible_new_pos = modify_conformer_torsion_angles(rigid_new_pos,
                                                           data['ligand', 'ligand'].edge_index.T[data['ligand'].edge_mask],
                                                           data['ligand'].mask_rotate if isinstance(data['ligand'].mask_rotate, np.ndarray) else data['ligand'].mask_rotate[0],
                                                           torsion_updates).to(rigid_new_pos.device)
        R, t = rigid_transform_Kabsch_3D_torch(flexible_new_pos.T, rigid_new_pos.T)
        aligned_flexible_pos = flexible_new_pos @ R.T + t.T
        data['ligand'].pos = aligned_flexible_pos
    else:
        data['ligand'].pos = rigid_new_pos
    return data


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


def get_t_schedule(inference_steps):
    return np.linspace(1, 0, inference_steps + 1)[:-1]


def set_time(complex_graphs, t_tr, t_rot, t_tor, batchsize, all_atoms, device):
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
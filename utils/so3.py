import os
import numpy as np
import torch
from scipy.spatial.transform import Rotation

MIN_EPS, MAX_EPS, N_EPS = 0.0005, 4, 2000
X_N = 2000

"""
    Preprocessing for the SO(3) sampling and score computations, truncated infinite series are computed and then
    cached to memory, therefore the precomputation is only run the first time the repository is run on a machine
"""

omegas = np.linspace(0, np.pi, X_N + 1)[1:]


def _compose(r1, r2):  # R1 @ R2 but for Euler vecs
    return Rotation.from_matrix(Rotation.from_rotvec(r1).as_matrix() @ Rotation.from_rotvec(r2).as_matrix()).as_rotvec()


def _expansion(omega, eps, L=2000):  # the summation term only
    l_vec = np.arange(L).reshape(-1, 1)
    p = ((2 * l_vec + 1) * np.exp(-l_vec * (l_vec + 1) * eps ** 2 / 2)
         * np.sin(omega * (l_vec + 1 / 2)) / np.sin(omega / 2)).sum(0)
    return p


def _density(expansion, omega, marginal=True):  # if marginal, density over [0, pi], else over SO(3)
    if marginal:
        return expansion * (1 - np.cos(omega)) / np.pi
    else:
        return expansion / 8 / np.pi ** 2  # the constant factor doesn't affect any actual calculations though


def _score(exp, omega, eps, L=2000):  # score of density over SO(3)
    l_vec = np.arange(L).reshape(-1, 1)
    hi = np.sin((l_vec + 1 / 2) * omega)
    dhi = (l_vec + 1 / 2) * np.cos((l_vec + 1 / 2) * omega)
    lo = np.sin(omega / 2)
    dlo = 1 / 2 * np.cos(omega / 2)
    dSigma = ((2 * l_vec + 1) * np.exp(-l_vec * (l_vec + 1) * eps**2 / 2) * (lo * dhi - hi * dlo) / lo ** 2).sum(0)
    return dSigma / exp


if os.path.exists('.so3_omegas_array4.npy'):
    _omegas_array = np.load('.so3_omegas_array4.npy')
    _cdf_vals = np.load('.so3_cdf_vals4.npy')
    _score_norms = np.load('.so3_score_norms4.npy')
    _exp_score_norms = np.load('.so3_exp_score_norms4.npy')
else:
    _eps_array = 10 ** np.linspace(np.log10(MIN_EPS), np.log10(MAX_EPS), N_EPS)
    _omegas_array = np.linspace(0, np.pi, X_N + 1)[1:]

    _exp_vals = np.asarray([_expansion(_omegas_array, eps) for eps in _eps_array])
    _pdf_vals = np.asarray([_density(_exp, _omegas_array, marginal=True) for _exp in _exp_vals])
    _cdf_vals = np.asarray([_pdf.cumsum() / X_N * np.pi for _pdf in _pdf_vals])
    _score_norms = np.asarray([_score(_exp_vals[i], _omegas_array, _eps_array[i]) for i in range(len(_eps_array))])

    _exp_score_norms = np.sqrt(np.sum(_score_norms**2 * _pdf_vals, axis=1) / np.sum(_pdf_vals, axis=1) / np.pi)

    np.save('.so3_omegas_array4.npy', _omegas_array)
    np.save('.so3_cdf_vals4.npy', _cdf_vals)
    np.save('.so3_score_norms4.npy', _score_norms)
    np.save('.so3_exp_score_norms4.npy', _exp_score_norms)


def sample(eps):
    eps_idx = (np.log10(eps) - np.log10(MIN_EPS)) / (np.log10(MAX_EPS) - np.log10(MIN_EPS)) * N_EPS
    eps_idx = np.clip(np.around(eps_idx).astype(int), a_min=0, a_max=N_EPS - 1)

    x = np.random.rand()
    return np.interp(x, _cdf_vals[eps_idx], _omegas_array)


def sample_vec(eps):
    x = np.random.randn(3)
    x /= np.linalg.norm(x)
    return x * sample(eps)


def score_vec(eps, vec):
    eps_idx = (np.log10(eps) - np.log10(MIN_EPS)) / (np.log10(MAX_EPS) - np.log10(MIN_EPS)) * N_EPS
    eps_idx = np.clip(np.around(eps_idx).astype(int), a_min=0, a_max=N_EPS - 1)

    om = np.linalg.norm(vec)
    return np.interp(om, _omegas_array, _score_norms[eps_idx]) * vec / om


def score_norm(eps):
    eps = eps.numpy()
    eps_idx = (np.log10(eps) - np.log10(MIN_EPS)) / (np.log10(MAX_EPS) - np.log10(MIN_EPS)) * N_EPS
    eps_idx = np.clip(np.around(eps_idx).astype(int), a_min=0, a_max=N_EPS-1)
    return torch.from_numpy(_exp_score_norms[eps_idx]).float()

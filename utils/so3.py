import os
import numpy as np
import torch
from scipy.spatial.transform import Rotation

MIN_EPS, MAX_EPS, N_EPS = 0.01, 2, 1000
X_N = 2000

"""
    Preprocessing for the SO(3) sampling and score computations, truncated infinite series are computed and then
    cached to memory, therefore the precomputation is only run the first time the repository is run on a machine
"""

"""
General

    This code is preprocessing some data for working with rotations in 3-dimensional space, which can be represented by the group SO(3). 
    The preprocessing calculates and caches some data that is used in later computations to generate random rotations with a particular distribution.

    The distribution of rotations is determined by a value called eps, which sets the scale of the distribution. 
    The code precomputes and caches various pieces of information for different values of eps ranging from a minimum value to a maximum value. 
    This precomputation is only run once, so that it can be used later without having to redo the calculation every time.

    The precomputed information includes the values of the distribution's density and cumulative distribution functions, as well as some 
    information about the score (derivative) of the density and its normalization. All of these values are related to the way that rotations 
    are distributed for different values of eps.

    The code uses several helper functions to perform these computations. The _compose function composes two rotations represented by Euler vectors. 
    The _expansion function calculates the sum term in a truncated infinite series representation of the density of the distribution. 
    The _density function calculates the density of the distribution. The _score function calculates the score (derivative) of the density.

    Finally, the code saves all the precomputed data to disk, so that it can be loaded and used later without having to perform the calculations again. 
    The sample function uses the precomputed data to generate random rotations with the specified distribution.

Specific

    The code implements precomputation for a distribution over the special orthogonal group SO(3), a group of 3D rotations. 
    It provides a method to sample from this distribution and calculate relevant statistics.

    The distribution is parameterized by a scalar eps and it is computed using a truncated infinite series. The code computes a
    nd caches to disk the series terms, CDF values, score norms, and expected score norms for a range of values for eps. If these 
    arrays have already been computed and saved, they are loaded from disk to memory, otherwise they are computed and saved.

    The _compose method calculates the composition of two rotations represented as Euler vectors. The _expansion method calculates 
    the series term for a given omega and eps. The _density method calculates the density over either [0, pi] or SO(3), depending on 
    the value of the marginal parameter. The _score method calculates the score of the density over SO(3).

    The sample method samples from the distribution, given a value of eps, by finding the closest precomputed eps_idx and using 
    linear interpolation on the corresponding precomputed CDF values. The pdf method calculates the density at a given omega and 
    eps by linear interpolation on the precomputed _exp_vals. The score method calculates the score of the density at a given omega 
    and eps by linear interpolation on the precomputed _score_norms. The score_norm method calculates the expected score norm for a 
    given eps by linear interpolation on the precomputed _exp_score_norms.



"""

omegas = np.linspace(0, np.pi, X_N + 1)[1:]


def _compose(r1, r2):  # R1 @ R2 but for Euler vecs
    return Rotation.from_matrix(Rotation.from_rotvec(r1).as_matrix() @ Rotation.from_rotvec(r2).as_matrix()).as_rotvec()


def _expansion(omega, eps, L=2000):  # the summation term only
    p = 0
    for l in range(L):
        p += (2 * l + 1) * np.exp(-l * (l + 1) * eps**2) * np.sin(omega * (l + 1 / 2)) / np.sin(omega / 2)
    return p


def _density(expansion, omega, marginal=True):  # if marginal, density over [0, pi], else over SO(3)
    if marginal:
        return expansion * (1 - np.cos(omega)) / np.pi
    else:
        return expansion / 8 / np.pi ** 2  # the constant factor doesn't affect any actual calculations though


def _score(exp, omega, eps, L=2000):  # score of density over SO(3)
    dSigma = 0
    for l in range(L):
        hi = np.sin(omega * (l + 1 / 2))
        dhi = (l + 1 / 2) * np.cos(omega * (l + 1 / 2))
        lo = np.sin(omega / 2)
        dlo = 1 / 2 * np.cos(omega / 2)
        dSigma += (2 * l + 1) * np.exp(-l * (l + 1) * eps**2) * (lo * dhi - hi * dlo) / lo ** 2
    return dSigma / exp


if os.path.exists('.so3_omegas_array2.npy'):
    _omegas_array = np.load('.so3_omegas_array2.npy')
    _cdf_vals = np.load('.so3_cdf_vals2.npy')
    _score_norms = np.load('.so3_score_norms2.npy')
    _exp_score_norms = np.load('.so3_exp_score_norms2.npy')
else:
    print("Precomputing and saving to cache SO(3) distribution table")
    _eps_array = 10 ** np.linspace(np.log10(MIN_EPS), np.log10(MAX_EPS), N_EPS)
    _omegas_array = np.linspace(0, np.pi, X_N + 1)[1:]

    _exp_vals = np.asarray([_expansion(_omegas_array, eps) for eps in _eps_array])
    _pdf_vals = np.asarray([_density(_exp, _omegas_array, marginal=True) for _exp in _exp_vals])
    _cdf_vals = np.asarray([_pdf.cumsum() / X_N * np.pi for _pdf in _pdf_vals])
    _score_norms = np.asarray([_score(_exp_vals[i], _omegas_array, _eps_array[i]) for i in range(len(_eps_array))])

    _exp_score_norms = np.sqrt(np.sum(_score_norms**2 * _pdf_vals, axis=1) / np.sum(_pdf_vals, axis=1) / np.pi)

    np.save('.so3_omegas_array2.npy', _omegas_array)
    np.save('.so3_cdf_vals2.npy', _cdf_vals)
    np.save('.so3_score_norms2.npy', _score_norms)
    np.save('.so3_exp_score_norms2.npy', _exp_score_norms)


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


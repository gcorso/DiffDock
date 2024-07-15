# From Nick Polizzi
import numpy as np
from collections import defaultdict
import prody as pr
import os

from datasets.constants import chi, atom_order, aa_long2short, aa_short2aa_idx, aa_idx2aa_short


def get_dihedral_indices(resname, chi_num):
    """Return the atom indices for the specified dihedral angle.
    """
    if resname not in chi:
        return np.array([np.nan]*4)
    if chi_num not in chi[resname]:
        return np.array([np.nan]*4)
    return np.array([atom_order[resname].index(x) for x in chi[resname][chi_num]])


dihedral_indices = defaultdict(list)
for aa in atom_order.keys():
    for i in range(1, 5):
        inds = get_dihedral_indices(aa, i)
        dihedral_indices[aa].append(inds)
    dihedral_indices[aa] = np.array(dihedral_indices[aa])


def vector_batch(a, b):
    return a - b


def unit_vector_batch(v):
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def dihedral_angle_batch(p):
    b0 = vector_batch(p[:, 0], p[:, 1])
    b1 = vector_batch(p[:, 1], p[:, 2])
    b2 = vector_batch(p[:, 2], p[:, 3])
    
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    
    m1 = np.cross(n1, b1 / np.linalg.norm(b1, axis=1, keepdims=True))
    
    x = np.sum(n1 * n2, axis=1)
    y = np.sum(m1 * n2, axis=1)
    
    deg = np.degrees(np.arctan2(y, x))

    deg[deg < 0] += 360

    return deg


def batch_compute_dihedral_angles(sidechains):
    sidechains_np = np.array(sidechains)
    dihedral_angles = dihedral_angle_batch(sidechains_np)
    return dihedral_angles


def get_coords(prody_pdb):
    resindices = sorted(set(prody_pdb.ca.getResindices()))
    coords = np.full((len(resindices), 14, 3), np.nan)
    for i, resind in enumerate(resindices):
        sel = prody_pdb.select(f'resindex {resind}')
        resname = sel.getResnames()[0]
        for j, name in enumerate(atom_order[aa_long2short[resname] if resname in aa_long2short else 'X']):
            sel_resnum_name = sel.select(f'name {name}')
            if sel_resnum_name is not None:
                coords[i, j, :] = sel_resnum_name.getCoords()[0]
            else:
                coords[i, j, :] = [np.nan, np.nan, np.nan]
    return coords


def get_onehot_sequence(seq):
    onehot = np.zeros((len(seq), 20))
    for i, aa in enumerate(seq):
        idx = aa_short2aa_idx[aa] if aa in aa_short2aa_idx else 7 # 7 is the index for GLY
        onehot[i, idx] = 1
    return onehot


def get_dihedral_indices(onehot_sequence):
    return np.array([dihedral_indices[aa_idx2aa_short[aa_idx]] for aa_idx in np.where(onehot_sequence)[1]])


def _get_chi_angles(coords, indices):
    X = coords
    Y = indices.astype(int)
    N = coords.shape[0]
    mask = np.isnan(indices)
    Y[mask] = 0
    Z = X[np.arange(N)[:, None, None], Y, :]
    Z[mask] = np.nan
    chi_angles = batch_compute_dihedral_angles(Z.reshape(-1, 4, 3)).reshape(N, 4)
    return chi_angles


def get_chi_angles(coords, seq, return_onehot=False):
    """

    Parameters
    ----------
    prody_pdb : prody.AtomGroup
        prody pdb object or selection
    return_coords : bool, optional
        return coordinates of prody_pdb in (N, 14, 3) array format, by default False
    return_onehot : bool, optional
        return one-hot sequence of prody_pdb, by default False

    Returns
    -------
    numpy array of shape (N, 4)
        Array contains chi angles of sidechains in row-order of residue indices in prody_pdb.
        If a chi angle is not defined for a residue, due to missing atoms or GLY / ALA, it is set to np.nan.
    """
    onehot = get_onehot_sequence(seq)
    dihedral_indices = get_dihedral_indices(onehot)
    if return_onehot:
        return _get_chi_angles(coords, dihedral_indices), onehot
    return _get_chi_angles(coords, dihedral_indices)


def test_get_chi_angles(print_chi_angles=False):
    # need internet connection of '6w70.pdb' in working directory
    pdb = pr.parsePDB('6w70')
    prody_pdb = pdb.select('chain A')
    chi_angles = get_chi_angles(prody_pdb)
    assert chi_angles.shape == (prody_pdb.ca.numAtoms(), 4)
    assert chi_angles[0,0] < 56.0 and chi_angles[0,0] > 55.0
    print('test_get_chi_angles passed')
    try:
        os.remove('6w70.pdb.gz')
    except:
        pass
    if print_chi_angles:
        print(chi_angles)
    return True


if __name__ == '__main__':
    test_get_chi_angles(print_chi_angles=True)



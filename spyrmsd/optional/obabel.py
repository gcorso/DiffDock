from typing import List, Optional, Tuple

import numpy as np
from openbabel import openbabel as ob
from openbabel import pybel

from spyrmsd import molecule, utils


def load(fname: str):
    """
    Load molecule from file.

    Parameters
    ----------
    fname: str
        File name

    Returns
    -------
    Molecule
    """

    fmt = utils.molformat(fname)

    obmol = next(pybel.readfile(fmt, fname))

    return obmol


def loadall(fname: str):
    """
    Load molecules from file.

    Parameters
    ----------
    fname: str
        File name

    Returns
    -------
    List of molecules
    """

    fmt = utils.molformat(fname)

    obmols = [obmol for obmol in pybel.readfile(fmt, fname)]

    # FIXME: Special handling for multi-model PDB files
    # See OpenBabel Issue #2097
    if fmt == "pdb":
        if len(obmols) > 1:  # Multi-model PDB file
            obmols = obmols[:-1]

    return obmols


def adjacency_matrix(mol) -> np.ndarray:
    """
    Adjacency matrix from OpenBabel molecule.

    Parameters
    ----------
    mol:
        Molecule

    Returns
    -------
    np.ndarray
        Adjacency matrix of the molecule
    """

    n = len(mol.atoms)

    # Pre-allocate memory for  the adjacency matrix
    A = np.zeros((n, n), dtype=int)

    # Loop over molecular bonds
    for bond in ob.OBMolBondIter(mol.OBMol):
        # Bonds are 1-indexed
        i: int = bond.GetBeginAtomIdx() - 1
        j: int = bond.GetEndAtomIdx() - 1

        # A molecular graph is undirected
        A[i, j] = A[j, i] = 1

    return A


def to_molecule(mol, adjacency: bool = True):
    """
    Transform molecule to `pyrmsd` molecule.

    Parameters
    ----------
    mol:
        Molecule
    adjacency: boolean, optional
        Flag to decide wether to build the adjacency matrix from molecule

    Returns
    -------
    pyrmsd.molecule.Molecule
        `pyrmsd` molecule
    """

    n = len(mol.atoms)

    atomicnums = np.zeros(n, dtype=int)
    coordinates = np.zeros((n, 3))

    for i, atom in enumerate(mol.atoms):
        atomicnums[i] = atom.atomicnum
        coordinates[i] = atom.coords

    A: Optional[np.ndarray] = adjacency_matrix(mol) if adjacency else None

    return molecule.Molecule(atomicnums, coordinates, A)


def numatoms(mol) -> int:
    """
    Number of atoms.

    Parameters
    ----------
    mol:
        Molecule

    Returns
    -------
    int
        Number of atoms
    """
    return mol.OBMol.NumAtoms()


def numbonds(mol) -> int:
    """
    Number of bonds.

    Parameters
    ----------
    mol:
        Molecule

    Returns
    -------
    int
        Number of bonds
    """
    return mol.OBMol.NumBonds()


def bonds(mol) -> List[Tuple[int, int]]:
    """
    List of bonds

    Parameters
    ----------
    mol:
        Molecule

    Returns
    -------
    List[Tuple[int, int]]
        List of bonds

    Notes
    -----
    A bond is defined by a tuple of (0-based) indices of two atoms.
    """
    b = []

    for bond in ob.OBMolBondIter(mol.OBMol):
        i = bond.GetBeginAtomIdx() - 1
        j = bond.GetEndAtomIdx() - 1

        b.append((i, j))

    return b

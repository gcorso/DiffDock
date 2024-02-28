import gzip
import os
from typing import List, Optional, Tuple

import numpy as np
import rdkit.Chem as Chem

from spyrmsd import molecule, utils


def _load_block_gzipped(loader, fname: str):
    """
    Load gzipped files using MolBlocks.

    Parameters
    ----------
    loader:
        RDKit MolBlock loader (MolFromMol2Block, MolFromPDBBlock, ...)
    fname: str
        File name

    Returns
    -------
    Molecule
    """
    with gzip.open(fname, "r") as fgz:
        content = fgz.read()
        rdmol = loader(content, removeHs=False)

    return rdmol


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

    gzipped = os.path.splitext(fname)[-1] == ".gz"
    fmt = utils.molformat(fname)

    if fmt == "mol2":
        if not gzipped:
            rdmol = Chem.MolFromMol2File(fname, removeHs=False)
        else:
            rdmol = _load_block_gzipped(Chem.MolFromMol2Block, fname)
    elif fmt == "sdf":
        if not gzipped:
            rdmol = next(Chem.SDMolSupplier(fname, removeHs=False))
        else:
            with gzip.open(fname, "r") as fgz:
                rdmol = next(Chem.ForwardSDMolSupplier(fgz, removeHs=False))
    elif fmt == "pdb":
        if not gzipped:
            rdmol = Chem.MolFromPDBFile(fname, removeHs=False)
        else:
            rdmol = _load_block_gzipped(Chem.MolFromPDBBlock, fname)
    else:
        raise NotImplementedError

    return rdmol


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

    gzipped = os.path.splitext(fname)[-1] == ".gz"
    fmt = utils.molformat(fname)

    if fmt == "mol2":
        raise NotImplementedError  # See RDKit Issue #415
    elif fmt == "sdf":
        if not gzipped:
            rdmols = Chem.SDMolSupplier(fname, removeHs=False)
            mols = [rdmol for rdmol in rdmols]
        else:
            with gzip.open(fname, "r") as fgz:
                rdmols = Chem.ForwardSDMolSupplier(fgz, removeHs=False)

                # Load all molecules before closing file
                mols = [rdmol for rdmol in rdmols]
    elif fmt == "pdb":
        # TODO: Implement
        raise NotImplementedError
    else:
        raise NotImplementedError

    return mols


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

    return Chem.rdmolops.GetAdjacencyMatrix(mol)


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
    spyrmsd.molecule.Molecule
        `spyrmsd` molecule
    """

    if mol is None:
        # Propagate RDKit parsing failure
        return None

    atoms = mol.GetAtoms()

    n = len(atoms)

    atomicnums = np.zeros(n, dtype=int)
    coordinates = np.zeros((n, 3))

    conformer = mol.GetConformer()

    for i, atom in enumerate(atoms):
        atomicnums[i] = atom.GetAtomicNum()

        pos = conformer.GetAtomPosition(i)

        coordinates[i] = np.array([pos.x, pos.y, pos.z])

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
    return mol.GetNumAtoms()


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
    return mol.GetNumBonds()


def bonds(mol) -> List[Tuple[int, int]]:
    """
    List of bonds.

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

    for bond in mol.GetBonds():
        b.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    return b

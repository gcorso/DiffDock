try:
    from spyrmsd.optional.obabel import (
        adjacency_matrix,
        bonds,
        load,
        loadall,
        numatoms,
        numbonds,
        to_molecule,
    )

except ImportError:
    try:
        from spyrmsd.optional.rdkit import (
            adjacency_matrix,
            bonds,
            load,
            loadall,
            numatoms,
            numbonds,
            to_molecule,
        )
    except ImportError:
        # Use sPyRMSD as standalone library
        __all__ = []
    else:
        # Avoid flake8 complaint "imported but unused"
        __all__ = [
            "load",
            "loadall",
            "adjacency_matrix",
            "to_molecule",
            "numatoms",
            "numbonds",
            "bonds",
        ]
else:
    # Avoid flake8 complaint "imported but unused"
    __all__ = [
        "load",
        "loadall",
        "adjacency_matrix",
        "to_molecule",
        "numatoms",
        "numbonds",
        "bonds",
    ]


def loadmol(fname: str, adjacency: bool = True):
    """
    Load molecule from file.

    Parameters
    ----------
    fname: str
        File name

    Returns
    -------
    molecule.Molecule
        `spyrmsd` molecule
    """

    mol = load(fname)

    return to_molecule(mol, adjacency=adjacency)


def loadallmols(fname: str, adjacency: bool = True):
    """
    Load molecules from file.

    Parameters
    ----------
    fname: str
        File name

    Returns
    -------
    List[molecule.Molecule]
        `spyrmsd` molecules
    """

    mols = loadall(fname)

    return [to_molecule(mol, adjacency=adjacency) for mol in mols]

from typing import Any, List, Optional, Tuple, Union

import numpy as np

from spyrmsd import graph, hungarian, molecule, qcp, utils


def rmsd(
    coords1: np.ndarray,
    coords2: np.ndarray,
    atomicn1: np.ndarray,
    atomicn2: np.ndarray,
    center: bool = False,
    minimize: bool = False,
    atol: float = 1e-9,
) -> float:
    """
    Compute RMSD

    Parameters
    ----------
    coords1: np.ndarray
        Coordinate of molecule 1
    coords2: np.ndarray
        Coordinates of molecule 2
    atomicn1: np.ndarray
        Atomic numbers for molecule 1
    atomicn2: np.ndarray
        Atomic numbers for molecule 2
    center: bool
        Center molecules at origin
    minimize: bool
        Compute minimum RMSD (with QCP method)
    atol: float
        Absolute tolerance parameter for QCP method (see :func:`qcp_rmsd`)

    Returns
    -------
    float
        RMSD

    Notes
    -----
    When `minimize=True`, the QCP method is used. [1]_ The molecules are
    centred at the origin according to the center of geometry and superimposed
    in order to minimize the RMSD.

    .. [1] D. L. Theobald, *Rapid calculation of RMSDs using a quaternion-based
       characteristic polynomial*, Acta Crys. A **61**, 478-480 (2005).
    """

    assert np.all(atomicn1 == atomicn2)
    assert coords1.shape == coords2.shape

    # Center coordinates if required
    c1 = utils.center(coords1) if center or minimize else coords1
    c2 = utils.center(coords2) if center or minimize else coords2

    if minimize:
        rmsd = qcp.qcp_rmsd(c1, c2, atol)
    else:
        n = coords1.shape[0]

        rmsd = np.sqrt(np.sum((c1 - c2) ** 2) / n)

    return rmsd


def hrmsd(
    coords1: np.ndarray,
    coords2: np.ndarray,
    atomicn1: np.ndarray,
    atomicn2: np.ndarray,
    center=False,
):
    """
    Compute minimum RMSD using the Hungarian method.

    Parameters
    ----------
    coords1: np.ndarray
        Coordinate of molecule 1
    coords2: np.ndarray
        Coordinates of molecule 2
    atomicn1: np.ndarray
        Atomic numbers for molecule 1
    atomicn2: np.ndarray
        Atomic numbers for molecule 2

    Returns
    -------
    float
        Minimum RMSD (after assignment)

    Notes
    -----
    The Hungarian algorithm is used to solve the linear assignment problem, which is
    a minimum weight matching of the molecular graphs (bipartite). [2]_

    The linear assignment problem is solved for every element separately.

    .. [2] W. J. Allen and R. C. Rizzo, *Implementation of the Hungarian Algorithm to
        Account for Ligand Symmetry and Similarity in Structure-Based Design*,
        J. Chem. Inf. Model. **54**, 518-529 (2014)
    """

    assert atomicn1.shape == atomicn2.shape
    assert coords1.shape == coords2.shape

    # Center coordinates if required
    c1 = utils.center(coords1) if center else coords1
    c2 = utils.center(coords2) if center else coords2

    return hungarian.hungarian_rmsd(c1, c2, atomicn1, atomicn2)


def _rmsd_isomorphic_core(
    coords1: np.ndarray,
    coords2: np.ndarray,
    aprops1: np.ndarray,
    aprops2: np.ndarray,
    am1: np.ndarray,
    am2: np.ndarray,
    center: bool = False,
    minimize: bool = False,
    isomorphisms: Optional[List[Tuple[List[int], List[int]]]] = None,
    atol: float = 1e-9,
) -> Tuple[float, List[Tuple[List[int], List[int]]], Tuple[List[int], List[int]]]:
    """
    Compute RMSD using graph isomorphism.

    Parameters
    ----------
    coords1: np.ndarray
        Coordinate of molecule 1
    coords2: np.ndarray
        Coordinates of molecule 2
    aprops1: np.ndarray
        Atomic properties for molecule 1
    aprops2: np.ndarray
        Atomic properties for molecule 2
    am1: np.ndarray
        Adjacency matrix for molecule 1
    am2: np.ndarray
        Adjacency matrix for molecule 2
    center: bool
        Centering flag
    minimize: bool
        Compute minized RMSD
    isomorphisms: Optional[List[Dict[int,int]]]
        Previously computed graph isomorphism
    atol: float
        Absolute tolerance parameter for QCP (see :func:`qcp_rmsd`)

    Returns
    -------
    Tuple[float, List[Dict[int, int]]]
        RMSD (after graph matching) and graph isomorphisms
    """

    assert coords1.shape == coords2.shape

    n = coords1.shape[0]

    # Center coordinates if required
    c1 = utils.center(coords1) if center or minimize else coords1
    c2 = utils.center(coords2) if center or minimize else coords2

    # No cached isomorphisms
    if isomorphisms is None:
        # Convert molecules to graphs
        G1 = graph.graph_from_adjacency_matrix(am1, aprops1)
        G2 = graph.graph_from_adjacency_matrix(am2, aprops2)

        # Get all the possible graph isomorphisms
        isomorphisms = graph.match_graphs(G1, G2)

    # Minimum result
    # Squared displacement (not minimize) or RMSD (minimize)
    min_result = np.inf
    min_isomorphisms = None

    # Loop over all graph isomorphisms to find the lowest RMSD
    for idx1, idx2 in isomorphisms:
        # Use the isomorphism to shuffle coordinates around (from original order)
        c1i = c1[idx1, :]
        c2i = c2[idx2, :]

        if not minimize:
            # Compute square displacement
            # Avoid dividing by n and an expensive sqrt() operation
            result = np.sum((c1i - c2i) ** 2)
        else:
            # Compute minimized RMSD using QCP
            result = qcp.qcp_rmsd(c1i, c2i, atol)

        if result < min_result:
            min_result = result
            min_isomorphisms = (idx1, idx2)

    if not minimize:
        # Compute actual RMSD from square displacement
        min_result = np.sqrt(min_result / n)

    # Return the actual RMSD
    return min_result, isomorphisms, min_isomorphisms


def symmrmsd(
    coordsref: np.ndarray,
    coords: Union[np.ndarray, List[np.ndarray]],
    apropsref: np.ndarray,
    aprops: np.ndarray,
    amref: np.ndarray,
    am: np.ndarray,
    center: bool = False,
    minimize: bool = False,
    cache: bool = True,
    atol: float = 1e-9,
    return_permutation: bool = False,
) -> Any:
    """
    Compute RMSD using graph isomorphism for multiple coordinates.

    Parameters
    ----------
    coordsref: np.ndarray
        Coordinate of reference molecule
    coords: List[np.ndarray]
        Coordinates of other molecule
    apropsref: np.ndarray
        Atomic properties for reference
    aprops: np.ndarray
        Atomic properties for other molecule
    amref: np.ndarray
        Adjacency matrix for reference molecule
    am: np.ndarray
        Adjacency matrix for other molecule
    center: bool
        Centering flag
    minimize: bool
        Minimum RMSD
    cache: bool
        Cache graph isomorphisms
    atol: float
        Absolute tolerance parameter for QCP (see :func:`qcp_rmsd`)

    Returns
    -------
    float: Union[float, List[float]]
        Symmetry-corrected RMSD(s) and graph isomorphisms

    Notes
    -----

    Graph isomorphism is introduced for symmetry corrections. However, it is also
    useful when two molecules do not have the atoms in the same order since atom
    matching according to atomic numbers and the molecular connectivity is
    performed. If atoms are in the same order and there is no symmetry, use the
    `rmsd` function.
    """

    if isinstance(coords, list):  # Multiple RMSD calculations
        RMSD: Any = []
        isomorphism = None
        min_iso = []

        for c in coords:
            if not cache:
                # Reset isomorphism
                isomorphism = None

            srmsd, isomorphism, min_i = _rmsd_isomorphic_core(
                coordsref,
                c,
                apropsref,
                aprops,
                amref,
                am,
                center=center,
                minimize=minimize,
                isomorphisms=isomorphism,
                atol=atol,
            )
            min_iso.append(min_i)
            RMSD.append(srmsd)

    else:  # Single RMSD calculation
        RMSD, isomorphism, min_iso = _rmsd_isomorphic_core(
            coordsref,
            coords,
            apropsref,
            aprops,
            amref,
            am,
            center=center,
            minimize=minimize,
            isomorphisms=None,
            atol=atol,
        )

    if return_permutation:
        return RMSD, min_iso
    return RMSD


def rmsdwrapper(
    molref: molecule.Molecule,
    mols: Union[molecule.Molecule, List[molecule.Molecule]],
    symmetry: bool = True,
    center: bool = False,
    minimize: bool = False,
    strip: bool = True,
    cache: bool = True,
) -> Any:
    """
    Compute RMSD between two molecule.

    Parameters
    ----------
    molref: molecule.Molecule
        Reference molecule
    mols: Union[molecule.Molecule, List[molecule.Molecule]]
        Molecules to compare to reference molecule
    symmetry: bool, optional
        Symmetry-corrected RMSD (using graph isomorphism)
    center: bool, optional
        Center molecules at origin
    minimize: bool, optional
        Minimised RMSD (using the quaternion polynomial method)
    strip: bool, optional
        Strip hydrogen atoms

    Returns
    -------
    List[float]
        RMSDs
    """

    if not isinstance(mols, list):
        mols = [mols]

    if strip:
        molref.strip()

        for mol in mols:
            mol.strip()

    if minimize:
        center = True

    cref = molecule.coords_from_molecule(molref, center)
    cmols = [molecule.coords_from_molecule(mol, center) for mol in mols]

    RMSDlist = []

    if symmetry:
        RMSDlist = symmrmsd(
            cref,
            cmols,
            molref.atomicnums,
            mols[0].atomicnums,
            molref.adjacency_matrix,
            mols[0].adjacency_matrix,
            center=center,
            minimize=minimize,
            cache=cache,
        )
    else:  # No symmetry
        for c in cmols:
            RMSDlist.append(
                rmsd(
                    cref,
                    c,
                    molref.atomicnums,
                    mols[0].atomicnums,
                    center=center,
                    minimize=minimize,
                )
            )

    return RMSDlist

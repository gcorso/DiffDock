try:
    from spyrmsd.graphs.gt import (
        cycle,
        graph_from_adjacency_matrix,
        lattice,
        match_graphs,
        num_edges,
        num_vertices,
        vertex_property,
    )

except ImportError:
    try:
        from spyrmsd.graphs.nx import (
            cycle,
            graph_from_adjacency_matrix,
            lattice,
            match_graphs,
            num_edges,
            num_vertices,
            vertex_property,
        )
    except ImportError:
        raise ImportError("graph_tool or NetworkX libraries not found.")

__all__ = [
    "graph_from_adjacency_matrix",
    "match_graphs",
    "vertex_property",
    "num_vertices",
    "num_edges",
    "lattice",
    "cycle",
    "adjacency_matrix_from_atomic_coordinates",
]

import numpy as np

from spyrmsd import constants


def adjacency_matrix_from_atomic_coordinates(
    aprops: np.ndarray, coordinates: np.ndarray
) -> np.ndarray:
    """
    Compute adjacency matrix from atomic coordinates.

    Parameters
    ----------
    aprops: numpy.ndarray
        Atomic properties
    coordinates: numpy.ndarray
        Atomic coordinates

    Returns
    -------
    numpy.ndarray
        Adjacency matrix

    Notes
    -----

    This function is based on an automatic bond perception algorithm: two
    atoms are considered to be bonded when their distance is smaller than
    the sum of their covalent radii plus a tolerance value. [3]_

    .. warning::
        The automatic bond perceptron rule implemented in this functions
        is very simple and only depends on atomic coordinates. Use
        with care!

    .. [3] E. C. Meng and R. A. Lewis, *Determination of molecular topology and atomic
       hybridization states from heavy atom coordinates*, J. Comp. Chem. **12**, 891-898
       (1991).
    """

    n = len(aprops)

    assert coordinates.shape == (n, 3)

    A = np.zeros((n, n))

    for i in range(n):
        r_i = constants.anum_to_covalentradius[aprops[i]]

        for j in range(i + 1, n):
            r_j = constants.anum_to_covalentradius[aprops[j]]

            distance = np.sqrt(np.sum((coordinates[i] - coordinates[j]) ** 2))

            if distance < (r_i + r_j + constants.connectivity_tolerance):
                A[i, j] = A[j, i] = 1

    return A

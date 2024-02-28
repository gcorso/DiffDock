import warnings
from typing import List, Optional, Union

import numpy as np

from spyrmsd import constants, graph, utils


class Molecule:
    def __init__(
        self,
        atomicnums: Union[np.ndarray, List[int]],
        coordinates: Union[np.ndarray, List[List[float]]],
        adjacency_matrix: Optional[Union[np.ndarray, List[List[int]]]] = None,
    ) -> None:
        """
        Molecule initialisation.

        Parameters
        ----------
        atomicnums: Union[np.ndarray, List[int]]
            Atomic numbers
        coordinates: Union[np.ndarray, List[List[float]]]
            Atomic coordinates
        adjacency_matrix: Union[np.ndarray, List[List[int]]], optional
            Molecular graph adjacency matrix

        Notes
        -----

        A molecule is built from atomic numbers and atomic coordinates only.
        Optionally, a good representation of the molecular graph (obtained with
        OpenBabel or RDKit) can be stored as adjacency matrix.
        """

        atomicnums = np.asarray(atomicnums, dtype=int)
        coordinates = np.asarray(coordinates, dtype=float)

        self.natoms: int = len(atomicnums)

        assert atomicnums.shape == (self.natoms,)
        assert coordinates.shape == (self.natoms, 3)

        self.atomicnums = atomicnums
        self.coordinates = coordinates

        self.stripped: bool = bool(np.all(atomicnums != 1))

        if adjacency_matrix is not None:
            self.adjacency_matrix: np.ndarray = np.asarray(adjacency_matrix, dtype=int)

        # Molecular graph
        self.G = None

        self.masses: Optional[List[float]] = None

    @classmethod
    def from_obabel(cls, obmol, adjacency: bool = True):
        """
        Constructor from OpenBabel molecule.

        Parameters
        ----------
        obmol:
            OpenBabel molecule
        adjacency:
            Flag to compute the adjacency matrix

        Returns
        -------
        spyrmsd.molecule.Molecule
            :code:`spyrmsd` Molecule
        """
        # TODO: Check if OpenBabel is available?
        from spyrmsd.optional import obabel as ob

        return ob.to_molecule(obmol, adjacency=adjacency)

    @classmethod
    def from_rdkit(cls, rdmol, adjacency: bool = True):
        """
        Constructor from RDKit molecule.

        Parameters
        ----------
        rdmol:
            RDKit molecule
        adjacency:
            Flag to compute the adjacency matrix

        Returns
        -------
        spyrmsd.molecule.Molecule
            :code:`spyrmsd` Molecule
        """

        # TODO: Check if RDKit is available?
        from spyrmsd.optional import rdkit as rd

        return rd.to_molecule(rdmol, adjacency=adjacency)

    def translate(self, vector: Union[np.ndarray, List[float]]) -> None:
        """
        Translate molecule.

        Parameters
        ----------
        vector: np.ndarray
            Translation vector (in 3D)
        """
        assert len(vector) == 3
        vector = np.asarray(vector)
        self.coordinates += vector

    def rotate(
        self, angle: float, axis: Union[np.ndarray, List[float]], units: str = "rad"
    ) -> None:
        """
        Rotate molecule.

        Parameters
        ----------
        angle: float
            Rotation angle
        axis: np.ndarray
            Axis of rotation (in 3D)
        units: {"rad", "deg"}
            Units of the angle (radians `rad` or degrees `deg`)
        """
        axis = np.asarray(axis)
        self.coordinates = utils.rotate(self.coordinates, angle, axis, units)

    def center_of_mass(self) -> np.ndarray:
        """
        Center of mass.

        Returns
        -------
        np.ndarray
            Center of mass

        Notes
        -----
        Atomic masses are cached.
        """

        # Get masses and cache them
        if self.masses is None:
            self.masses = [constants.anum_to_mass[anum] for anum in self.atomicnums]

        return np.average(self.coordinates, axis=0, weights=self.masses)

    def center_of_geometry(self) -> np.ndarray:
        """
        Center of geometry.

        Returns
        -------
        np.ndarray
            Center of geometry
        """
        return utils.center_of_geometry(self.coordinates)

    # TODO: Change name (to stripH)
    def strip(self) -> None:
        """
        Strip hydrogen atoms.
        """

        if not self.stripped:
            idx = self.atomicnums != 1  # Non-hydrogen atoms

            # Strip
            self.atomicnums = self.atomicnums[idx]
            self.coordinates = self.coordinates[idx, :]

            # Update number of atoms
            self.natoms = len(self.atomicnums)

            # Update adjacency matrix
            if self.adjacency_matrix is not None:
                self.adjacency_matrix = self.adjacency_matrix[np.ix_(idx, idx)]

            # Reset molecular graph when stripping
            self.G = None

            self.stripped = True

    def to_graph(self):
        """
        Convert molecule to graph.

        Returns
        -------
        Graph
            Molecular graph.

        Notes
        -----
        If the molecule does not have an associated adjacency matrix, a simple
        bond perception is used.

        The molecular graph is cached.
        """
        if self.G is None:
            try:
                self.G = graph.graph_from_adjacency_matrix(
                    self.adjacency_matrix, self.atomicnums
                )
            except AttributeError:
                warnings.warn(
                    "Molecule was not initialized with an adjacency matrix. "
                    + "Using bond perception..."
                )

                # Automatic bond perception (with very simple rule)
                self.adjacency_matrix = graph.adjacency_matrix_from_atomic_coordinates(
                    self.atomicnums, self.coordinates
                )

                self.G = graph.graph_from_adjacency_matrix(
                    self.adjacency_matrix, self.atomicnums
                )

        return self.G

    def __len__(self) -> int:
        """
        Molecule size.

        Returns
        -------
        int
            Number of atoms within the molecule
        """
        return self.natoms


def coords_from_molecule(mol: Molecule, center: bool = False) -> np.ndarray:
    """
    Atomic coordinates from molecule.

    Parameters
    ----------
    mol: molecule.Molecule
        Molecule
    center: bool
        Center flag

    Returns
    -------
    np.ndarray
        Atomic coordinates (possibly centred)

    Notes
    -----
    Atomic coordinates are centred according to the center of geometry, not the center
    of mass.
    """

    if center:
        coords = mol.coordinates - mol.center_of_geometry()
    else:
        coords = mol.coordinates

    return coords

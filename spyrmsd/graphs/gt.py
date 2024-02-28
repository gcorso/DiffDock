import warnings
from typing import Any, List, Optional, Tuple, Union

import graph_tool as gt
import numpy as np
from graph_tool import generation, topology

from spyrmsd.exceptions import NonIsomorphicGraphs
from spyrmsd.graphs._common import (
    error_non_isomorphic_graphs,
    warn_disconnected_graph,
    warn_no_atomic_properties,
)


# TODO: Implement all graph-tool supported types
def _c_type(numpy_dtype):
    """
    Get C type compatible with graph-tool from numpy dtype

    Parameters
    ----------
    numpy_dtype: np.dtype
        Numpy dtype

    Returns
    -------
    str
        C type

    Raises
    ------
    ValueError
        If the data type is not supported

    Notes
    -----
    https://graph-tool.skewed.de/static/doc/quickstart.html#sec-property-maps
    """
    name: str = numpy_dtype.name

    if "int" in name:
        return "int"
    elif "float" in name:
        return "double"
    elif "str" in name:
        return "string"
    else:
        raise ValueError(f"Unsupported property type: {name}")


def graph_from_adjacency_matrix(
    adjacency_matrix: Union[np.ndarray, List[List[int]]],
    aprops: Optional[Union[np.ndarray, List[Any]]] = None,
):
    """
    Graph from adjacency matrix.

    Parameters
    ----------
    adjacency_matrix: Union[np.ndarray, List[List[int]]]
        Adjacency matrix
    aprops: Union[np.ndarray, List[Any]], optional
        Atomic properties

    Returns
    -------
    Graph
        Molecular graph

    Notes
    -----
    It the atomic numbers are passed, they are used as node attributes.
    """

    # Get upper triangular adjacency matrix
    adj = np.triu(adjacency_matrix)

    assert adj.shape[0] == adj.shape[1]
    num_vertices = adj.shape[0]

    G = gt.Graph(directed=False)
    G.add_vertex(n=num_vertices)
    G.add_edge_list(np.transpose(adj.nonzero()))

    # Check if graph is connected, for warning
    cc, _ = topology.label_components(G)
    if set(cc.a) != {0}:
        warnings.warn(warn_disconnected_graph)

    if aprops is not None:
        if not isinstance(aprops, np.ndarray):
            aprops = np.array(aprops)

        assert aprops.shape[0] == num_vertices

        ptype: str = _c_type(aprops.dtype)  # Get C type
        vprop = G.new_vertex_property(ptype, vals=aprops)  # Create property map
        G.vertex_properties["aprops"] = vprop  # Set property map

    return G


def match_graphs(G1, G2) -> List[Tuple[List[int], List[int]]]:
    """
    Compute graph isomorphisms.

    Parameters
    ----------
    G1:
        Graph 1
    G2:
        Graph 2

    Returns
    -------
    List[Tuple[List[int],List[int]]]
        All possible mappings between nodes of graph 1 and graph 2 (isomorphisms)

    Raises
    ------
    NonIsomorphicGraphs
        If the graphs `G1` and `G2` are not isomorphic
    """

    try:
        maps = topology.subgraph_isomorphism(
            G1,
            G2,
            vertex_label=(
                G1.vertex_properties["aprops"],
                G2.vertex_properties["aprops"],
            ),
            subgraph=False,
        )
    except KeyError:
        warnings.warn(warn_no_atomic_properties)

        maps = topology.subgraph_isomorphism(G1, G2, subgraph=False)

    # Check if graphs are actually isomorphic
    if len(maps) == 0:
        raise NonIsomorphicGraphs(error_non_isomorphic_graphs)

    n = num_vertices(G1)

    # Extract all isomorphisms in a list
    return [(np.arange(0, n, dtype=int), m.a) for m in maps]


def vertex_property(G, vproperty: str, idx: int) -> Any:
    """
    Get vertex (node) property from graph

    Parameters
    ----------
    G:
        Graph
    vproperty: str
        Vertex property name
    idx: int
        Vertex index

    Returns
    -------
    Any
        Vertex property value
    """
    return G.vertex_properties[vproperty][idx]


def num_vertices(G) -> int:
    """
    Number of vertices

    Parameters
    ----------
    G:
        Graph

    Returns
    -------
    int
        Number of vertices (nodes)
    """
    return G.num_vertices()


def num_edges(G) -> int:
    """
    Number of edges

    Parameters
    ----------
    G:
        Graph

    Returns
    -------
    int
        Number of edges
    """
    return G.num_edges()


def lattice(n1: int, n2: int):
    """
    Build 2D lattice graph

    Parameters
    ----------
    n1: int
        Number of nodes in dimension 1
    n2: int
        Number of nodes in dimension 2

    Returns
    -------
    Graph
        Lattice graph
    """
    return generation.lattice((n1, n2))


def cycle(n):
    """
    Build cycle graph

    Parameters
    ----------
    n: int
        Number of nodes

    Returns
    -------
    Graph
        Cycle graph
    """
    return generation.circular_graph(n)

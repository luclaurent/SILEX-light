"""Unit tests for the Gmsh IO helper module."""

from pathlib import Path

import numpy as np

from SILEXlight import silex_lib_gmsh as sx_gmsh


def _build_sample_mesh(mesh_stem: Path) -> Path:
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
    elements = np.array([[1, 2, 3]], dtype=int)
    field = np.array([[0.0], [0.1], [0.2]], dtype=float)

    sx_gmsh.WriteResults(mesh_stem, nodes, elements, 2, [[field, 'nodal', 1, 'sample']])
    return Path(str(mesh_stem) + '.msh')


def test_WriteResult(tmp_path):
    mesh_file = _build_sample_mesh(tmp_path / 'gmsh_roundtrip')
    assert mesh_file.exists()


def test_ReadGmshNodes(tmp_path):
    mesh_file = _build_sample_mesh(tmp_path / 'gmsh_nodes')
    nodes = sx_gmsh.ReadGmshNodes(mesh_file, nbcoord=2)

    assert nodes.shape == (3, 2)
    np.testing.assert_allclose(nodes[1], np.array([1.0, 0.0]))


def test_ReadGmshElements(tmp_path):
    mesh_file = _build_sample_mesh(tmp_path / 'gmsh_elements')
    elements, node_ids = sx_gmsh.ReadGmshElements(mesh_file, elmttype=2, prop=1)

    assert elements.shape == (1, 3)
    np.testing.assert_array_equal(elements[0], np.array([0, 1, 2]))
    np.testing.assert_array_equal(node_ids, np.array([0, 1, 2]))
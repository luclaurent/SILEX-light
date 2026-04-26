"""Unit tests for the TRI3 pure Python implementation."""

from SILEXlight import silex_lib_tri3_python as sx
import numpy as np


def test_det33_ligne_de_un():
    mat = np.array([[1.0, 1.0, 1.0], [0.0, 2.0, 1.0], [1.0, 1.0, 3.0]], dtype=float)
    np.testing.assert_allclose(sx.det33_ligne_de_un(mat[1:3, :]), np.linalg.det(mat))


def test_stiffnessmatrix():
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
    elements = np.array([[1, 2, 3]], dtype=int)
    material = np.array([210000.0, 0.3, 1.0], dtype=float)

    ik, jk, vk = sx.stiffnessmatrix(nodes, elements, material)

    assert ik.shape == (36,)
    assert jk.shape == (36,)
    assert vk.shape == (36,)
    assert np.isfinite(vk).all()


def test_compute_stress_strain_error():
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
    elements = np.array([[1, 2, 3]], dtype=int)
    material = np.array([210000.0, 0.3, 1.0], dtype=float)
    qq = np.array([0.0, 0.0, 1.0e-5, 0.0, 0.0, 1.0e-5], dtype=float)

    sigma, sig_smooth, eps_elem, eps_nodes, err_elem, err_glob = sx.compute_stress_strain_error(
        nodes, elements, material, qq
    )

    assert sigma.shape == (1, 4)
    assert sig_smooth.shape == (3, 4)
    assert eps_elem.shape == (1, 3)
    assert eps_nodes.shape == (3, 3)
    assert err_elem.shape == (1,)
    assert np.isfinite(err_glob)


def test_forceonline():
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
    edge_elements = np.array([[1, 2]], dtype=int)
    fs = np.array([0.0, -1000.0, 0.0, -1000.0], dtype=float)
    pts = np.array([0.0, 0.0, 1.0, 0.0], dtype=float)

    forces = sx.forceonline(nodes, edge_elements, fs, pts)

    assert forces.shape == (6,)
    assert np.isclose(forces[1] + forces[3], -1000.0)
    assert np.isclose(forces[0] + forces[2], 0.0)


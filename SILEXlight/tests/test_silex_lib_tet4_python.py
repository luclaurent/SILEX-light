"""Unit tests for the TET4 pure Python implementation."""

from SILEXlight import silex_lib_tet4_python as sx_tet4
import numpy as np


def _tet4_fixture():
    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ], dtype=float)
    elements = np.array([[1, 2, 3, 4]], dtype=int)
    material = np.array([210000.0, 0.3], dtype=float)
    return nodes, elements, material


def test_compute_stress_strain_error():
    nodes, elements, material = _tet4_fixture()
    qq = np.array([
        0.0, 0.0, 0.0,
        1.0e-5, 0.0, 0.0,
        0.0, 1.0e-5, 0.0,
        0.0, 0.0, 1.0e-5,
    ], dtype=float)

    sigma, sig_smooth, eps_elem, eps_nodes, err_elem, err_glob = sx_tet4.compute_stress_strain_error(
        nodes, elements, material, qq
    )

    assert sigma.shape == (1, 7)
    assert sig_smooth.shape == (4, 7)
    assert eps_elem.shape == (1, 6)
    assert eps_nodes.shape == (4, 6)
    assert err_elem.shape == (1,)
    assert np.isfinite(err_glob)


def test_det33_ligne_de_un():
    mat = np.random.randn(3, 3)
    mat[0, :] = 1
    np.testing.assert_almost_equal(
        sx_tet4.det33_ligne_de_un(mat[1:3, :]), np.linalg.det(mat))


def test_det44_ligne_de_un():
    mat = np.random.randn(4, 4)
    mat[0, :] = 1
    np.testing.assert_almost_equal(
        sx_tet4.det44_ligne_de_un(mat[1:4, :]), np.linalg.det(mat))


def test_forceonsurface():
    nodes, _, _ = _tet4_fixture()
    surface_elements = np.array([[1, 2, 3]], dtype=int)
    press = 12.0

    fp = sx_tet4.forceonsurface(nodes, surface_elements, press, [0.0, 0.0, 0.0])

    assert fp.shape == (12,)
    assert np.isclose(fp[2] + fp[5] + fp[8], 0.5 * press)
    assert np.isclose(fp[0] + fp[3] + fp[6], 0.0)
    assert np.isclose(fp[1] + fp[4] + fp[7], 0.0)


def test_stiffnessmatrix():
    nodes, elements, material = _tet4_fixture()
    ik, jk, vk = sx_tet4.stiffnessmatrix(nodes, elements, material)

    assert ik.shape == (144,)
    assert jk.shape == (144,)
    assert vk.shape == (144,)
    assert np.isfinite(vk).all()

def test_massmatrix():
    nodes, elements, _ = _tet4_fixture()
    ik, jk, vk = sx_tet4.massmatrix(nodes, elements, rho=7800.0)

    assert ik.shape == (144,)
    assert jk.shape == (144,)
    assert vk.shape == (144,)
    assert np.isfinite(vk).all()
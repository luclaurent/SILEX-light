"""Unit tests for the TET4 Fortran extension module."""

import numpy as np
import pytest

sx_tet4 = pytest.importorskip('SILEXlight.silex_lib_tet4_fortran')


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


def test_crossproduct():
    vA = np.random.randn(3)
    vB = np.random.randn(3)
    np.testing.assert_almost_equal(
        sx_tet4.crossproduct(vA, vB), np.cross(vA, vB))


def test_det33():
    mat = np.random.randn(3, 3)
    np.testing.assert_almost_equal(sx_tet4.det33(mat), np.linalg.det(mat))


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
    fp = sx_tet4.forceonsurface(nodes, surface_elements, 12.0, [0.0, 0.0, 0.0])

    assert fp.shape == (12,)
    assert np.isclose(fp[2] + fp[5] + fp[8], 6.0)


def test_normvector():
    vA = np.random.randn(3, 1)
    np.testing.assert_almost_equal(sx_tet4.normvector(vA,), np.linalg.norm(vA))


def test_stiffnessmatrix():
    nodes, elements, material = _tet4_fixture()
    ik, jk, vk = sx_tet4.stiffnessmatrix(nodes, elements, material)

    assert ik.shape == (144,)
    assert jk.shape == (144,)
    assert vk.shape == (144,)
    assert np.isfinite(vk).all()
    
def test_massmatrix():
    nodes, elements, _ = _tet4_fixture()
    ik, jk, vk = sx_tet4.massmatrix(nodes, elements, 7800.0)

    assert ik.shape == (144,)
    assert jk.shape == (144,)
    assert vk.shape == (144,)
    assert np.isfinite(vk).all()

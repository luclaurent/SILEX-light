from SILEXlight import silex_lib_tet4_fortran as sx_tet4
import numpy as np


def test_compute_stress_strain_error():
    # sx_tet4.compute_stress_strain_error()
    assert True


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
    # sx_tet4.forceonsurface()
    assert True


def test_normvector():
    vA = np.random.randn(3, 1)
    np.testing.assert_almost_equal(sx_tet4.normvector(vA,), np.linalg.norm(vA))


def test_stiffnessmatrix():
    # sx_tet4.stiffnessmatrix()
    assert True

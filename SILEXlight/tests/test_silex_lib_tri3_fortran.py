"""Unit tests for the TRI3 Fortran extension module."""

import numpy as np
import pytest

sx_tri3 = pytest.importorskip('SILEXlight.silex_lib_tri3_fortran')

# subroutine ElementalStiffness(X,Y,thickness,CC,ke)
# subroutine StiffnessMatrix(nbnodes,nodes,
# subroutine ComputeElementalStress(X,Y,CC,Q,Sig,Eps,area)
# subroutine Compute_stress_strain_error(nbnodes,nodes,
# subroutine forceonline(nbnodes,nodes,nbelem,elements,


def test_det33_ligne_de_un():
    mat = np.random.randn(3, 3)
    mat[0, :] = 1
    np.testing.assert_almost_equal(
        sx_tri3.det33_ligne_de_un(mat[1:3, :]), np.linalg.det(mat))


def test_ElementalStiffness():
    assert hasattr(sx_tri3, 'elementalstiffness')
    assert callable(sx_tri3.elementalstiffness)


def test_StiffnessMatrix():
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
    elements = np.array([[1, 2, 3]], dtype=int)
    material = np.array([210000.0, 0.3, 1.0], dtype=float)

    ik, jk, vk = sx_tri3.stiffnessmatrix(nodes, elements, material)

    assert ik.shape == (36,)
    assert jk.shape == (36,)
    assert vk.shape == (36,)
    assert np.isfinite(vk).all()


def test_ComputeElementalStress():
    assert hasattr(sx_tri3, 'computeelementalstress')
    assert callable(sx_tri3.computeelementalstress)


def test_Compute_stress_strain_error():
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=float)
    elements = np.array([[1, 2, 3]], dtype=int)
    material = np.array([210000.0, 0.3, 1.0], dtype=float)
    qq = np.array([0.0, 0.0, 1.0e-5, 0.0, 0.0, 1.0e-5], dtype=float)

    sigma, sig_smooth, eps_elem, eps_nodes, err_elem, err_glob = sx_tri3.compute_stress_strain_error(
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

    forces = sx_tri3.forceonline(nodes, edge_elements, fs, pts)

    assert forces.shape == (6,)
    assert np.isclose(forces[1] + forces[3], -1000.0)

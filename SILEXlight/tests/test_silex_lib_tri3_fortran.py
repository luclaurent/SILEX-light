from SILEXlight import silex_lib_tri3_fortran as sx_tri3
import numpy as np

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
    assert True


def test_StiffnessMatrix():
    assert True


def test_ComputeElementalStress():
    assert True


def test_Compute_stress_strain_error():
    assert True


def test_forceonline():
    assert True

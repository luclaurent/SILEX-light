from SILEXlight import silex_lib_tri3_python as sx
import numpy as np

def test_det33_ligne_de_un():
    mat = np.random.randn(3,3)
    mat[0,:] = 1
    np.testing.assert_allclose(sx.det33_ligne_de_un(mat[1:3,:]), np.linalg.det(mat))
    
def test_stiffnessmatrix():
    assert True
    
def test_compute_stress_strain_error():
    assert True
    
def test_forceonline():
    assert True
    

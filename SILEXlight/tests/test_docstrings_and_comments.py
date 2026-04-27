"""Tests for API docstrings and Fortran routine-purpose comments."""

from pathlib import Path

from SILEXlight import silex_lib_gmsh
from SILEXlight import silex_lib_tet4_python
from SILEXlight import silex_lib_tri3_python


def _has_meaningful_docstring(obj) -> bool:
    doc = getattr(obj, '__doc__', None)
    return bool(doc and len(doc.strip()) >= 20)


def test_tri3_python_function_docstrings():
    assert _has_meaningful_docstring(silex_lib_tri3_python.det33_ligne_de_un)
    assert _has_meaningful_docstring(silex_lib_tri3_python.stiffnessmatrix)
    assert _has_meaningful_docstring(silex_lib_tri3_python.compute_stress_strain_error)
    assert _has_meaningful_docstring(silex_lib_tri3_python.forceonline)


def test_tet4_python_function_docstrings():
    assert _has_meaningful_docstring(silex_lib_tet4_python.det33_ligne_de_un)
    assert _has_meaningful_docstring(silex_lib_tet4_python.det44_ligne_de_un)
    assert _has_meaningful_docstring(silex_lib_tet4_python.stiffnessmatrix)
    assert _has_meaningful_docstring(silex_lib_tet4_python.massmatrix)
    assert _has_meaningful_docstring(silex_lib_tet4_python.forceonsurface)
    assert _has_meaningful_docstring(silex_lib_tet4_python.compute_stress_strain_error)


def test_gmsh_python_function_docstrings():
    assert _has_meaningful_docstring(silex_lib_gmsh.WriteResults)
    assert _has_meaningful_docstring(silex_lib_gmsh.ReadGmshNodes)
    assert _has_meaningful_docstring(silex_lib_gmsh.ReadGmshElements)



def test_fortran_files_have_purpose_comments():
    repo_root = Path(__file__).resolve().parents[2]
    tri3_fortran = repo_root / 'SILEXlight' / 'silex_lib_tri3_fortran.f'
    tet4_fortran = repo_root / 'SILEXlight' / 'silex_lib_tet4_fortran.f'

    tri3_text = tri3_fortran.read_text(encoding='utf-8', errors='ignore')
    tet4_text = tet4_fortran.read_text(encoding='utf-8', errors='ignore')

    assert 'Purpose:' in tri3_text
    assert 'Purpose:' in tet4_text

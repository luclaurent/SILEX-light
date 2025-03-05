# SILEXlight

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0) [![pypi release](https://img.shields.io/pypi/v/SILEXlight.svg)](https://test.pypi.org/project/SILEXlight/) ![PyPI - Downloads](https://img.shields.io/pypi/dm/SILEXlight)

[![CI-ubuntu](https://github.com/luclaurent/SILEX-light/actions/workflows/CI-ubuntu.yml/badge.svg)](https://github.com/luclaurent/SILEX-light/actions/workflows/CI-ubuntu.yml) [![CI-mac-os](https://github.com/luclaurent/SILEX-light/actions/workflows/CI-macos.yml/badge.svg)](https://github.com/luclaurent/SILEX-light/actions/workflows/CI-macos.yml) [![Wheels and sdist](https://github.com/luclaurent/SILEX-light/actions/workflows/CI-build-release.yml/badge.svg)](https://github.com/luclaurent/SILEX-light/actions/workflows/CI-build-release.yml)

SILEX is a finite element code written in Python language, eventually with a Fortran part in order to speed up the computations.
    
*   The Python language is used to define parameters, to read the mesh, to solve the system, to write the results.  
*   The Fortran language is eventually used for elemental computations as well as to build the stiffness matrix.  
*   The open source software [Gmsh](http://www.geuz.org/gmsh/) is used to create the meshes as well as to show the results.  
*   The only free routines available on-line concern the 4-node tetrahedral element and the 3-node triangle element in the case of linear static analysis. They can be adapted to other elements.  
*   The following proposed applications are available for education purpose. They allow to understand the code, in order to perform other computations for other mechanical systems. Later on, the user can develop new elements or new method, and thus extend the possibilities of the code.  
*   The following course document available on-line ([here](http://antoinelegay.free.fr/Cours-programmation-english.pdf)) gives a programming introduction.
## Usage

* A complete example of the use of SILEX on a piston is available [here](calculs/piston_readme.md) in english ([french version](calculs/piston_readme_fr.md))
* A example pratical work for education is provided on a landing gear fork [here](calculs/fork_readme.md) in english ([french version](calculs/fork_readme_fr.md))

## Installation
Classical installation by executing

    pip install --user .

Editable installation for developing

    export SETUPTOOLS_ENABLE_FEATURES="legacy-editable"
    pip install --user -e .

## Run tests

The unit tests can be ran by using the following command

    pytest --pyargs SILEXlight.tests  

## License

`SILEXlight` is available under the LGPLv3 license. See the LICENSE file for more info.
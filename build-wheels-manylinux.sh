#!/bin/bash
set -e -x

PKGNAME='SILEXlight'
PATHSRC='/io/'
PATHWHEEL='/io/wheelhouse/'
PATHTEST='/io/'
PYTHONVERSION=("cp37-cp37m" "cp38-cp38" "cp39-cp39" "cp310-cp310" "cp311-cp311")

PYHOME=/home
cd ${PYHOME}

# Compile wheels
for V in ${PYTHONVERSION[@]}; do
    PYBIN=/opt/python/${V}/bin
    "${PYBIN}/pip" wheel ${PATHSRC} -w ${PATHWHEEL}
    # "${PYBIN}/python" /io/setup.py sdist -d /io/wheelhouse/
done

# Bundle external shared libraries into the wheels and fix naming
for whl in ${PATHWHEEL}/${PKGNAME}*-linux*.whl; do
    auditwheel repair "$whl" -w ${PATHWHEEL}
done

# Test
for PYBIN in ${PYTHONVERSION[@]}; do
PYBIN=/opt/python/${V}/bin
    "${PYBIN}/pip" install pytest pytest-cov
    "${PYBIN}/pip" install --no-index -f ${PATHWHEEL} ${PKGNAME}
    (cd "$PYHOME"; "${PYBIN}/pytest" --pyargs ${PKGNAME}.tests)
done


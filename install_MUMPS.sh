#!/bin/bash

# install MUMPS from source
BASE_DATA_PYTHON=$(python -c "from sysconfig import get_paths;print(get_paths()['data'])")
git clone https://github.com/scivision/mumps.git $1
cd $1 
cmake -Bbuild \
    -DBUILD_SINGLE=on \
    -DBUILD_DOUBLE=on \
    -DBUILD_COMPLEX=on \
    -DBUILD_COMPLEX16=on \
    -DMUMPS_parallel=no \
    -DBUILD_SHARED_LIBS=on  \
    -DMUMPS_UPSTREAM_VERSION=5.8.2 \
    --install-prefix=${BASE_DATA_PYTHON}
cmake --build build -j2
cmake --install build  
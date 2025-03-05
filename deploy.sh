#!/bin/bash
# Deploy SILEXlight weel to pypi server.

set +x 

PKGNAME='SILEXlight'
PATHSRC=$(pwd)
PATHWHEELLINUX=$(pwd)/wheelhouse
PATHWHEEL=$(pwd)/wheelhouse-current


function clean()
{
    [[ -d ${PATHWHEEL} ]] && rm ${PATHWHEEL}/*.whl && rm ${PATHWHEEL}/*.tar.gz
    [[ -d ${PATHWHEELLINUX} ]] && rm ${PATHWHEELLINUX}/*.whl
}

function wheellinux()
{
    docker run --rm -ti -v $PATHSRC:/io quay.io/pypa/manylinux_2_28_x86_64:latest /io/build-wheels-manylinux.sh
    SUCCESS=$?
    if [ $SUCCESS -eq 0 ]
    then
        echo "Ready for uploading :" $(ls ${PATHWHEELLINUX}/${PKGNAME}*.whl)
        return 0
    else
        echo "Test failed. Stop deployment." >&2
        return 1
    fi
}

function sdist()
{
    python -m build --sdist --outdir ${PATHWHEEL}
}

function wheel()
{
    # Build and test wheel in an venv to check if all required files are present in
    # the wheel.
    # Clean-up previous version
    python -m build --outdir ${PATHWHEEL}
    # Setup venv
    TEMP_DIR=$(mktemp -d)
    python -m venv $TEMP_DIR
    source $TEMP_DIR/bin/activate
    # Change dir to test really wheel file and not the repo files
    cd $TEMP_DIR
    # Install the wheel
    pip install ${PATHWHEEL}/${PKGNAME}*.whl
    # Run all tests
    pip install pytest pytest-cov
    python -m pytest --pyargs ${PKGNAME}.tests
    SUCCESS=$?
    # Clean-up tmp directory
    cd ${PATHWHEEL}
    trap 'rm -rf "$TEMP_DIR"' EXIT
    if [ $SUCCESS -eq 0 ]
    then
        echo "Ready for uploading :" $(ls ${PATHWHEEL}/${PKGNAME}*.whl)
        return 0
    else
        echo "Test failed. Stop deployment." >&2
        return 1
    fi
}

function listfiles()
{   
    listfiles=()
    if [[ -d ${PATHWHEEL} ]]; then
        for file in $(ls ${PATHWHEEL}/${PKGNAME}*.whl 2> /dev/null); do
            listfiles[${#listfiles[@]}]=$file
        done
    fi    
    if [[ -d ${PATHWHEEL} ]]; then
        for file in $(ls ${PATHWHEEL}/${PKGNAME}*.tar.gz 2> /dev/null); do
            listfiles[${#listfiles[@]}]=$file
        done
    fi
    if [[ -d ${PATHWHEELLINUX} ]]; then
        for file in $(ls ${PATHWHEELLINUX}/${PKGNAME}*manylinux*.whl 2> /dev/null); do
            listfiles[${#listfiles[@]}]=$file
        done
    fi
    echo ${listfiles[@]}
    if [[ ${#listfiles[@]} -gt 0 ]]; then
        return 0
    else
        return 1
    fi
}

function checkfiles()
{
    listfiles > /dev/null 2>&1
    SUCCESS=$?
    if [ ${SUCCESS} -eq 0 ]; then
        echo "Ready to upload:" $(listfiles)
        return 0
    else
        echo "Nothing to upload" $(listfiles)
        return 1
    fi
}
function machine()
{
    unameOut=$(uname -s)
    case "${unameOut}" in
        Linux*)     m="Linux";;
        Darwin*)    m="Darwin";;
        CYGWIN*)    m="Cygwin";;
        MINGW*)     m="MinGw";;
        *)          m="UNKNOWN:${unameOut}";;
    esac
    echo "${m}"
}

function build()
{
    echo "Clean directories"
    clean
    echo "Build wheels linux"
    wheellinux
    SUCCESSA=$?

    echo "Build wheels and/or source for current platform"
    case $(machine) in
        Linux)
            machine
            sdist
            SUCCESSB=$?
            ;;
        Darwin|Cygwin|MinGw)
            machine
            wheel
            SUCCESSB=$?
            ;;
        *)  SUCCESSB=1;;
    esac    
    if [ $((SUCCESSA+SUCCESSB)) -eq 0 ]; then
        echo "Build done"
        return 0
    else
        echo "Error on build"
        return 1
    fi
}

# Main deployment script
for i in "$@"
do
    case "$i" in
    -i|--install) pip install --user twine build virtualenv
    ;;
    -l|--list)
    checkfiles
    ;;
    -b|--build)
    build
    ;;
    -n|--nobuildtest)
    checkfiles
    twine upload --skip-existing --verbose -r testpypi $(listfiles)
    ;;
    -u|--nobuild)
    checkfiles
    twine upload --skip-existing --verbose $(listfiles)
    ;;
    -t|--test)
    build && twine upload --skip-existing --verbose -r testpypi $(listfiles)
    ;;
    -d|--deploy)
    build && twine upload --skip-existing --verbose $(listfiles)
    ;;
    esac
done

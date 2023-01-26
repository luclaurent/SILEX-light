#!/bin/bash
# Deploy SILEXlight weel to pypi server.

# wheel path
WHEEL_DIR=$(pwd)/dist

function wheel()
{
    # Build and test wheel in an venv to check if all required files are present in
    # the wheel.
    # Clean-up previous version
    rm  $WHEEL_DIR/SILEXlight*.whl
    pip3 wheel . -w dist
    # Setup venv
    TEMP_DIR=$(mktemp -d)
    python3 -m venv $TEMP_DIR
    source $TEMP_DIR/bin/activate
    # Change dir to test really wheel file and not the repo files
    cd $TEMP_DIR
    # Install the wheel
    pip3 install $WHEEL_DIR/SILEXlight*.whl
    # Run all tests
    pytest --pyargs SILEXlight.tests
    SUCCESS=$?
    # Clean-up tmp directory
    cd $WHEEL_DIR
    trap 'rm -rf "$TEMP_DIR"' EXIT
    if [ $SUCCESS -eq 0 ]
    then
        echo "Ready for uploading :" $(ls $WHEEL_DIR/SILEXlight*.whl)
        return 0
    else
        echo "Test failed. Stop deployment." >&2
        return 1
    fi
}

# Main deployment script
for i in "$@"
do
    case "$i" in
    -i|--install) pip install --user twine
    ;;
    -t|--test)
    wheel && twine upload --verbose -r testpypi $WHEEL_DIR/SILEXlight*.whl
    ;;
    -d|--deploy)
    wheel && twine upload --verbose $WHEEL_DIR/SILEXlight*.whl
    ;;
    esac
done

# Development

## Run tests on installed version

    python -m pytest --pyargs SILEXlight.tests

## Build and upload wheels

### Build

Build of wheels for Linux require docker. Execute the following command:

    docker run --rm -ti -v $PWD:/io quay.io/pypa/manylinux_2_28_x86_64:latest /io/build-wheels-manylinux.sh
    

### Build and upload

Run `./deploys.sh -i` to install useful packages for deploying and then run `./deploy.sh -t` to deploy on `testpypi`or `./deploy.sh -d` to deploy on `pypi`.

Check correct deployment on `testpypi` with the following command 

    pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ --no-cache-dir SILEXlight

Keyring can be used to store credentials:

    keyring set https://test.pypi.org/legacy/ <username>
    keyring set https://pypi.org/legacy/ <username>
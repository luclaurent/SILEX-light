FROM quay.io/jupyter/scipy-notebook:latest

# install 
USER root
RUN sudo apt-get update && sudo apt-get install -yy \
          build-essential \
          wget \
          cmake \
          gcc \
          g++ \
          gfortran \
          libblas-dev \
          liblapack-dev \
          libopenmpi-dev
          

          
COPY install_MUMPS.sh /tmp/install_MUMPS.sh
RUN mkdir /tmp/mumps &&\
    chmod +x /tmp/install_MUMPS.sh &&\
    /tmp/install_MUMPS.sh /tmp/mumps
USER $NB_UID
    

RUN BASE_DATA_PYTHON=$(python -c "from sysconfig import get_paths;print(get_paths()['data'])") &&\
    pip install pytest &&\
    pip install -v -Ccmake.define.MUMPS_ROOT=${BASE_DATA_PYTHON}/ https://github.com/luclaurent/pymumps/archive/refs/heads/main.zip
RUN pip install git+https://github.com/luclaurent/SILEX-light
RUN pytest --pyargs mumps
COPY calculs work/calculs
USER root
RUN chown -R $NB_USER:$NB_GROUP work/calculs
USER $NB_UID
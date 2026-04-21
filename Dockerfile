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
          libopenmpi-dev \
          gmsh

          
COPY install_MUMPS.sh /tmp/install_MUMPS.sh
RUN mkdir /tmp/mumps &&\
    chmod +x /tmp/install_MUMPS.sh &&\
    /tmp/install_MUMPS.sh /tmp/mumps
    

RUN BASE_DATA_PYTHON=$(python -c "from sysconfig import get_paths;print(get_paths()['data'])") &&\
    pip install pytest gmsh&&\
    pip install -v -Ccmake.define.MUMPS_ROOT=${BASE_DATA_PYTHON}/ https://github.com/luclaurent/pymumps/archive/refs/heads/main.zip
RUN pip install git+https://github.com/luclaurent/SILEX-light
RUN pytest --pyargs mumps
COPY calculs work/calculs

RUN chown -R jovyan:users work/calculs
RUN mkdir /home/jovyan/external && \
    mkdir -p /home/jovyan/install

VOLUME /home/jovyan/external

USER $NB_UID

EXPOSE 8888
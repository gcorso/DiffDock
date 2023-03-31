FROM nvidia/cuda:11.6.0-cudnn8-devel-ubuntu20.04


RUN chsh -s /bin/bash

SHELL ["/bin/bash", "-c"]

WORKDIR /root/

RUN apt-get update
RUN apt-get install -y wget bzip2 apt-utils git
RUN apt-get install -y \
    g++ \
    python3-dev \
    libeigen3-dev \
    wget \
    cmake \
    less \
    unzip
RUN apt-get clean

RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh -O ~/anaconda.sh && \
        /bin/bash ~/anaconda.sh -b -p /opt/conda && \
        rm ~/anaconda.sh && \
        ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
        echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
        find /opt/conda/ -follow -type f -name '*.a' -delete && \
        find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
        /opt/conda/bin/conda clean -afy

# Installing OpenBabel
RUN wget https://github.com/openbabel/openbabel/releases/download/openbabel-3-1-1/openbabel-3.1.1-source.tar.bz2
RUN tar -xjf openbabel-3.1.1-source.tar.bz2
RUN cd openbabel-3.1.1
RUN cmake ../openbabel-3.1.1 \
    -DPYTHON_BINDINGS=ON \
    -DRUN_SWIG=ON \
    -DCMAKE_INSTALL_PREFIX=/opt/conda \
    -DPYTHON_INCLUDE_DIR=/opt/conda/include/python3.10 \
    -DCMAKE_LIBRARY_PATH=/opt/conda/lib \
    -DSWIG_DIR=/opt/conda/share/swig/4.0.2 \
    -DSWIG_EXECUTABLE=/opt/conda/bin/swig \
    -DPYTHON_LIBRARY=/opt/conda/lib/libpython3.10.so \
    -DCMAKE_BUILD_TYPE=DEBUG
RUN make -j4
RUN make install

ENV PATH /opt/conda/bin:$PATH

# setup conda virtual environment
COPY ./environment.yml ./environment.yml

RUN conda update conda \
	&& conda env create --name diffdock -f ./environment.yml

RUN echo "conda activate diffdock" >> ~/.bashrc
ENV PATH /opt/conda/envs/diffdock/bin:$PATH
ENV CONDA_DEFAULT_ENV $diffdock

#install torch specific packages
RUN pip install --upgrade pip 
COPY ./requirements_docker_GPU_oddt.txt ./

RUN pip install --no-cache-dir torch==1.13.0 torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu116
RUN pip install --no-cache-dir -r ./requirements_docker_GPU_oddt.txt

RUN pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-1.13.0+cu116.html
RUN pip install torch-geometric torch-cluster -f https://data.pyg.org/whl/torch-1.13.0+cu116.html

#install oddt
RUN pip3 install oddt

COPY . .
RUN git submodule init
RUN git submodule update
RUN pip install -e ./esm/.

RUN chmod 777 test.sh
RUN oddt_cli test/test.sdf \
    --receptor test/test.pdb \
    #--score autodock_vina \
    #--score rfscore \
    --score rfscore_v1 \
    --score rfscore_v2 \
    --score rfscore_v3 \
    --score nnscore \
    #--score pleclinear \
    #--score plecnn \
    #--score plecrf \
    -O test/test_scored.sdf
CMD ["bash"]

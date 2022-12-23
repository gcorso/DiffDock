FROM jupyter/scipy-notebook:latest
# FROM cschranz/gpu-jupyter

LABEL author="Tuple, LLC <contact@tuple.xyz>"

## Environment Settings
ENV DEBIAN_FRONTEND=noninteractive \
    JUPYTER_ENABLE_LAB=yes \
    JUPYTER_TOKEN=diffdock \
    NB_USER=diffdock \
    CHOWN_HOME=yes \
    GRANT_SUDO=yes \
    RESTARTABLE=yes \
    SETUPTOOLS_USE_DISTUTILS=stdlib


ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility

USER root



## Install Basic Dependencies
RUN apt-get clean && \
    apt-get update && \
    apt-get -y install git curl wget python3.9
    # nvidia-cuda-toolkit

## Install PyTorch with dependencies
# RUN conda install --quiet --yes \
#     pytorch torchvision cudatoolkit=10.1 rdkit \
#     pyyaml mkl mkl-include setuptools cmake cffi typing pot && \
#     conda clean --all -f -y && \
#     fix-permissions $CONDA_DIR 

## Install CUDA Toolkit (to fix libcusparse.so.11 issue - /opt/conda/lib/libcusparse.so.11)
RUN conda install --yes cudatoolkit=11.3 -c pytorch && \
    conda clean --all -f -y && \
    export LD_LIBRARY_PATH="/opt/conda/lib/:$LD_LIBRARY_PATH" && \
    export PATH="/opt/conda/lib/:$PATH"

## Install DiffDock Dependencies
RUN pip install ipython-autotime && \
    pip install pyg==0.7.1 --quiet && \
    pip install pyyaml==6.0 --quiet && \
    pip install scipy==1.7.3 --quiet && \
    pip install networkx==2.6.3 --quiet && \
    pip install biopython==1.79 --quiet && \
    pip install rdkit-pypi==2022.03.5 --quiet && \
    pip install e3nn==0.5.0 --quiet && \
    pip install spyrmsd==0.5.2 --quiet && \
    pip install pandas==1.3.5 --quiet && \
    pip install biopandas==0.4.1 --quiet && \
    pip install torch==1.12.1+cu113 --extra-index-url https://download.pytorch.org/whl/cu113
    # pip install torch==1.12.1+cu113 --quiet
    # pip install torch==1.12.1 --quiet

### Git DiffDock
RUN mkdir /content && \
    cd /content && \
    git clone https://github.com/gcorso/DiffDock.git && \
    cd /content/DiffDock && \
    git checkout 0f9c419

### Install DiffDock-Specific Torch Dependencies
RUN pip uninstall torch-scatter torch-sparse torch-geometric torch-cluster --y && \
    pip install torch-scatter -f https://data.pyg.org/whl/torch-1.12.0%2Bcu113.html --quiet && \
    pip install torch-sparse -f https://data.pyg.org/whl/torch-1.12.0%2Bcu113.html --quiet && \
    pip install torch-cluster -f https://data.pyg.org/whl/torch-1.12.0%2Bcu113.html --quiet && \
    # pip install torch-scatter -f https://data.pyg.org/whl/torch-1.12.1.html --quiet && \
    # pip install torch-sparse -f https://data.pyg.org/whl/torch-1.12.1.html --quiet && \
    # pip install torch-cluster -f https://data.pyg.org/whl/torch-1.12.1.html --quiet && \
    pip install git+https://github.com/pyg-team/pytorch_geometric.git --quiet

### Install Facebook ESM
RUN cd /content/DiffDock && \
    git clone https://github.com/facebookresearch/esm && \
    cd /content/DiffDock/esm && \
    git checkout f07aed6 && \
    pip install -e .  && \
    cd /content/DiffDock

COPY ./DiffDock_sample_notebook.ipynb /content/DiffDock

RUN chmod -R 777 /content

WORKDIR /content/DiffDock

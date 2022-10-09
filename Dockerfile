# Modified from github user cwlchka

FROM pytorch/pytorch:1.12.1-cuda11.3-cudnn8-devel

RUN apt-get update && \
    apt-get install -y sudo \
                       python3-pip \
                       git  \
                       curl


RUN pip install torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric -f https://data.pyg.org/whl/torch-1.12.0+cu113.html
RUN pip install --no-cache-dir torch-spline-conv==1.2.1 -f https://data.pyg.org/whl/torch-1.12.0+cu113.html
#RUN pip install --no-cache-dir torch-cluster==1.6.0 -f https://data.pyg.org/whl/torch-1.9.0+cu111.html
#RUN pip install --no-cache-dir torch-geometric


RUN pip install --no-cache-dir pyyaml==6.0 \
                               scipy==1.7.3 \
                               networkx==2.6.3 \
                               biopython==1.79 \
                               rdkit-pypi==2022.03.5 \
                               e3nn==0.5.0 \
                               spyrmsd==0.5.2 \
                               pandas==1.3.5 \
                               biopandas==0.4.1 --quiet

RUN echo "alias python=python3" >> ~/.bashrc
RUN alias python=python3
RUN python -m pip install pytest 
# Install Jupyter kernel and jupyternotebook
RUN pip install --no-cache-dir notebook
RUN python -m pip install --no-cache-dir ipykernel 
RUN python -m ipykernel install --user --name=DiffDock


WORKDIR /app

RUN cd /app && git clone https://github.com/gcorso/DiffDock.git && cd /app/DiffDock && git checkout 0f9c419

WORKDIR /app/DiffDock

RUN cd /app/DiffDock && \
    git clone https://github.com/facebookresearch/esm && \
    cd esm && \
    git checkout f07aed6 && \
    pip install -e . && \
    cd /app/DiffDock



ENTRYPOINT bash

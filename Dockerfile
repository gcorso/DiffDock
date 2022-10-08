# Modified from github user cwlchka

FROM nvcr.io/nvidia/cuda:11.1.1-cudnn8-devel-ubuntu20.04

RUN apt-get update && \
    apt-get install -y sudo \
                       python-is-python3 \
                       python3-pip \
                       git

RUN pip3 install --no-cache-dir torch==1.9.0+cu111 torchvision==0.10.0+cu111 torchaudio==0.9.0 -f https://download.pytorch.org/whl/torch_stable.html

RUN pip install --no-cache-dir torch-scatter -f https://pytorch-geometric.com/whl/torch-1.9.0+cu111.html
RUN pip install --no-cache-dir torch-sparse -f https://pytorch-geometric.com/whl/torch-1.9.0+cu111.html
RUN pip install --no-cache-dir torch-cluster -f https://pytorch-geometric.com/whl/torch-1.9.0+cu111.html
RUN pip install --no-cache-dir torch-spline-conv -f https://pytorch-geometric.com/whl/torch-1.9.0+cu111.html
RUN pip install --no-cache-dir torch-geometric


RUN pip install --no-cache-dir pyg==0.7.1 \
                               pyyaml==6.0 \
                               scipy==1.7.3 \
                               networkx==2.6.3 \
                               biopython==1.79 \
                               rdkit-pypi==2022.03.5 \
                               e3nn==0.5.0 \
                               spyrmsd==0.5.2 \
                               pandas==1.3.5 \
                               biopandas==0.4.1 \
                               torch==1.12.1+cu113 --quiet

RUN echo "alias python=python3" >> ~/.bashrc
RUN alias python=python3
RUN python -m pip install pytest

WORKDIR /app
ENTRYPOINT bash

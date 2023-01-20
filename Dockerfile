FROM nvidia/cuda:11.6.0-cudnn8-devel-ubuntu20.04


RUN chsh -s /bin/bash

SHELL ["/bin/bash", "-c"]

WORKDIR /root/

RUN apt-get update
RUN apt-get install -y wget bzip2 apt-utils 
RUN apt-get clean



RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh -O ~/anaconda.sh && \
        /bin/bash ~/anaconda.sh -b -p /opt/conda && \
        rm ~/anaconda.sh && \
        ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
        echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
        find /opt/conda/ -follow -type f -name '*.a' -delete && \
        find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
        /opt/conda/bin/conda clean -afy

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
COPY ./requirements_docker_GPU.txt ./

RUN pip install --no-cache-dir torch==1.13.0 torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu116
RUN pip install --no-cache-dir -r ./requirements_docker_GPU.txt

RUN pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-1.13.0+cu116.html
RUN pip install torch-geometric torch-cluster -f https://data.pyg.org/whl/torch-1.13.0+cu116.html


COPY . .
RUN git submodule init
RUN git submodule update
RUN pip install -e ./esm/.

CMD ["bash"]
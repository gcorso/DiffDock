FROM nvidia/cuda:11.6.0-cudnn8-devel-ubuntu20.04


RUN chsh -s /bin/bash
WORKDIR /root/

# basic pod dependencies
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ENV DEBIAN_FRONTEND noninteractive\
    SHELL=/bin/bash
RUN apt-get update --yes && \
    # - apt-get upgrade is run to patch known vulnerabilities in apt-get packages as
    #   the ubuntu base image is rebuilt too seldom sometimes (less than once a month)
    apt-get upgrade --yes && \
    apt install --yes --no-install-recommends\
    git\
    wget\
    curl\
    git\
    bash\
    bzip2\
    openssh-server &&\
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    echo "en_US.UTF-8 UTF-8" > /etc/locale.gen
RUN apt-get clean

# installing model dependencies
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh -O ~/anaconda.sh && \
        /bin/bash ~/anaconda.sh -b -p /opt/conda && \
        rm ~/anaconda.sh && \
        ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
        echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
        find /opt/conda/ -follow -type f -name '*.a' -delete && \
        find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
        /opt/conda/bin/conda clean -afy

ENV PATH /opt/conda/bin:$PATH

## setup conda virtual environment
COPY ./environment.yml ./environment.yml

RUN conda update conda \
	&& conda env create --name DiffDock -f ./environment.yml

RUN echo "conda activate DiffDock" >> ~/.bashrc
ENV PATH /opt/conda/envs/DiffDock/bin:$PATH
ENV CONDA_DEFAULT_ENV $DiffDock

## install torch specific packages
RUN pip install --upgrade pip 
COPY ./requirements_docker_GPU.txt ./
RUN pip install --no-cache-dir torch==1.12.1 torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu116
RUN pip install --no-cache-dir -r ./requirements_docker_GPU.txt
RUN pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-1.12.1+cu116.html
RUN pip install torch-geometric torch-cluster -f https://data.pyg.org/whl/torch-1.12.1+cu116.html
COPY . .
RUN git submodule init
RUN git submodule update
RUN pip install -e ./esm/.

# install jupyter lab extensions
RUN pip install jupyterlab
RUN pip install ipywidgets
RUN pip install jupyter-archive
RUN jupyter nbextension enable --py widgetsnbextension

ADD start.sh /
RUN chmod +x /start.sh
#TODO #68 create proper test
ADD test.sh /
RUN chmod +x /test.sh
CMD [ "./start.sh" ]

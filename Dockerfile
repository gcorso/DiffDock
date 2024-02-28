# Stage 1: Build Environment Setup
FROM nvidia/cuda:11.7.1-devel-ubuntu22.04 as builder

RUN apt-get update -y && apt-get install -y wget curl git tar bzip2 unzip && rm -rf /var/lib/apt/lists/*

# Create a user
ENV APPUSER="appuser"
ENV HOME=/home/$APPUSER
RUN useradd -m -u 1000 $APPUSER
USER $APPUSER
WORKDIR $HOME

ENV ENV_NAME="diffdock"
ENV DIR_NAME="ligbind"

# Install micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
ENV PATH=$HOME/bin:$HOME/.local/bin:$PATH

# Copy and create Conda environment
ENV ENV_FILE_NAME=environment.yml
COPY --chown=$APPUSER:$APPUSER ./$ENV_FILE_NAME .
RUN ~/bin/micromamba env create --file $ENV_FILE_NAME && ~/bin/micromamba clean -afy --quiet

# Copy application code
COPY --chown=$APPUSER:$APPUSER cuda-ubuntu22.04 $HOME/$DIR_NAME

# Download models
# These should download automatically on first inference
# RUN curl -L -o diffdock_models_v1.1.zip "https://www.dropbox.com/scl/fi/drg90rst8uhd2633tyou0/diffdock_models.zip?rlkey=afzq4kuqor2jb8adah41ro2lz&dl=1" \
#     && mkdir -p $HOME/$DIR_NAME/workdir \
#     && unzip diffdock_models_v1.1.zip -d $HOME/$DIR_NAME/workdir


# Stage 2: Runtime Environment
FROM nvidia/cuda:11.7.1-runtime-ubuntu22.04

# Create user and setup environment
ENV APPUSER="appuser"
ENV HOME=/home/$APPUSER
RUN useradd -m -u 1000 $APPUSER
USER $APPUSER
WORKDIR $HOME

ENV ENV_NAME="diffdock"
ENV DIR_NAME="ligbind"

# Copy the Conda environment and application code from the builder stage
COPY --from=builder --chown=$APPUSER:$APPUSER $HOME/micromamba $HOME/micromamba
COPY --from=builder --chown=$APPUSER:$APPUSER $HOME/bin $HOME/bin
COPY --from=builder --chown=$APPUSER:$APPUSER $HOME/$DIR_NAME $HOME/$DIR_NAME
WORKDIR $HOME/$DIR_NAME

# Set the environment variables
ENV MAMBA_ROOT_PREFIX=$HOME/micromamba
ENV PATH=$HOME/bin:$HOME/.local/bin:$PATH
RUN micromamba shell init -s bash --root-prefix $MAMBA_ROOT_PREFIX

# Precompute series for SO(2) and SO(3) groups
RUN micromamba run -n ${ENV_NAME} python utils/precompute_series.py

# Expose ports for streamlit and gradio
EXPOSE 7860 8501

# Default command
CMD ["sh", "-c", "micromamba run -n ${ENV_NAME} python utils/print_device.py"]

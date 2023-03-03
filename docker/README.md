# Running DiffDock in Docker

## Running Docker

### Build
```sh
# with JupyterLab
docker build -t diffdock -f JupterLab_Dockerfile .

# with NVIDIA Base Image and PyTorch (CLI)
docker build -t diffdock -f CLI_Dockerfile .
```

### Run
```sh
# with CPU (and mounting a volume)
docker run -v C:\data:/data --name diffdock --rm -p 8888:8888 diffdock

# with GPUs
docker run --gpus all --name diffdock --rm -p 8888:8888 diffdock
```


## Jupyter Notes

- Password: `diffdock`
- CUDA Device: You may need to set your CUDA device for torch using `torch.cuda.set_device(<DEVICE_ID>)`
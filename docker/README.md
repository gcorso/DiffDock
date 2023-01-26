# Running DiffDock in Docker

## Running Docker

### Build
```sh
docker build -t diffdock .
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
# Running DiffDock in Docker


## Running Docker

### Build and Run (w/ CPU)
```sh
docker build -t diffdock .

# with CPU (and mounting a volume)
docker run -v C:\data:/data --name diffdock --rm -p 8888:8888 diffdock
```

### Build and Run (w/ GPU)
```sh
docker build -t diffdock .

# with GPUs
docker run --gpus all --name diffdock --rm -p 8888:8888 diffdock
```
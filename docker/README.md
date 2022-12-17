# Running DiffDock in Docker


## Running Docker

### Build
```sh
docker build -t diffdock .
```

### Run
```sh
# with GPUs
docker run --gpus all --name diffdock --rm -p 8888:8888 diffdock
# with CPU (and mounting a volume)
docker run -v T:\diffdock_tests:/data --name diffdock --rm -p 8888:8888 diffdock
```
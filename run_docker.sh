#!/bin/bash

# Build the Docker image without cache
# sudo docker build --no-cache -f ./docker/Dockerfile -t hyper_spec-m .

# Run the Docker container with GPU support and dataset volume mounted
sudo docker run --gpus all -v ~/PXD001468/:/dataset/ -it hyper_spec-m /bin/bash

apt install vim

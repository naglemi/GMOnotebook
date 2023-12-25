#!/bin/bash

# Ask user or use a command-line argument to decide on CUDA
read -p "Include CUDA? (y/n): " include_cuda
export INCLUDE_CUDA=$include_cuda

# Build the Singularity container
sudo SINGULARITY_CACHEDIR=/home/gmobot/GMOGUI/singularity_cache TMPDIR=/home/gmobot/GMOGUI/singularity_tmp singularity -d build gmodetector.sif gmodetector-setup-cuda.def &> build.log


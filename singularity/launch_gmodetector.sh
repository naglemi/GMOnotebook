#!/bin/bash

# Prompt the user to ask if they want to use the GPU version
echo "Overview of this launcher script and what to do next:"
echo ""
echo "1) This script will provide a link, which you can copy and paste into your web browser. Remember, copying from a Linux terminal is ctrl+shift+c, whereas ctrl+c will terminate the script."
echo ""
echo "2) Next, use the file browser on the left side of Jupyter Lab to navigate to /home/GMOnotebook/."
echo ""
echo "3) Select the notebook templates you wish to run, then make a copy and Save As to /mnt/archived/notebooks/."
echo ""
echo "4) Follow instructions for filling out and running your chosen notebook template."
echo ""
echo ""
echo "First, a question. Do you want to use the GPU-enabled version? (yes/no)"
read use_gpu

if [[ "$use_gpu" == "yes" ]]; then
    # User chose GPU version
    SIF_IMAGE="gmodetector_cuda.sif"
else
    # User chose non-GPU version or provided an unrecognized response
    SIF_IMAGE="gmodetector.sif"
fi

# Launch the appropriate Singularity image
singularity exec --contain --bind /media/gmobot:/mnt/drives \
                 --bind output:/mnt/output \
                 --bind models:/mnt/models \
                 --bind archived_notebooks:/mnt/archived_notebooks \
                 $SIF_IMAGE jupyter-lab --notebook-dir=/


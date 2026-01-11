#!/bin/bash
# to be used by service

source /APSshare/miniconda/x86_64/etc/profile.d/conda.sh
CONDA_ENV=matilda
conda activate "${CONDA_ENV}"

# Run the Python script
exec python /home/beams/USAXS/Apps/Matilda/matilda/matilda.py

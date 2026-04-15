#!/bin/bash
# to be used by service

source /APSshare/miniconda/x86_64/etc/profile.d/conda.sh
CONDA_ENV=matilda
conda activate "${CONDA_ENV}"

# Log directory — overrides the default ~/.local/share/matilda/log used on dev machines.
export MATILDA_LOG_DIR=/share1/log/matilda

# Run the Python script.
# Must be invoked from the repo root so that 'import matilda' resolves correctly.
# exec python /home/beams/USAXS/Apps/Matilda/matilda/matilda.py  # old bare-script launch
cd /home/beams/USAXS/Apps/Matilda
exec python -m matilda.matilda

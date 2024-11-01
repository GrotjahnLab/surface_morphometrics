#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate morphometrics || {
    echo "Error: Failed to activate morphometrics environment"
    exit 1
}

# Execute command or start bash
if [ $# -eq 0 ]; then
    exec bash
else
    exec "$@"
fi
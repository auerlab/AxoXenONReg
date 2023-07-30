#!/bin/sh -e

#   Dependencies:
#       Requires trimmed reads.

sbatch \
    --output=Logs/07-reference/slurm-%A.out \
    --error=Logs/07-reference/slurm-%A.err \
    SLURM/07-reference.sh

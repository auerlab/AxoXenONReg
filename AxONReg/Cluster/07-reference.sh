#!/bin/sh -e

#   Dependencies:
#       Requires trimmed reads.

# This script is single threaded, so there is no difference between
# running under srun or on a workstation.
cmd="SLURM/07-reference.sh"

# Run interactively under SLURM if srun is found, otherwise run directly
if which srun > /dev/null; then
    sbatch $cmd
else
    $cmd 2>&1 | tee Logs/07-reference/07-reference.out
fi

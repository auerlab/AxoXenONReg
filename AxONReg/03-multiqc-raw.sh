#!/bin/sh -e

#   Dependencies:
#       Requires raw FastQC results.  Run after *-qc-raw.sbatch.

# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
# export LC_ALL=en_US.UTF-8

cd Results/03-multiqc-raw
rm -rf *
cmd="env LC_ALL=en_US.UTF-8 LANG=en_US.utf-8 multiqc ../02-qc-raw"

# Run interactively under SLURM if srun is found, otherwise run directly
if which srun > /dev/null; then
    set -x
    srun $cmd 2>&1 | tee ../../Logs/03-multiqc-raw/multiqc.out
else
    set -x
    $cmd 2>&1 | tee ../../Logs/03-multiqc-raw/multiqc.out
fi


#!/bin/sh -e

#
#   Dependencies:
#       Requires raw FastQC results.  Run after *-qc-raw.sbatch.

# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
# export LC_ALL=en_US.UTF-8

cd Results/06-multiqc-trimmed
rm -rf *
cmd="env LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 \
    multiqc ../05-qc-trimmed 2>&1 | tee ../../Logs/06-multiqc-trimmed/multiqc.out"

# Run interactively under SLURM if srun is found, otherwise run directly
if which srun > /dev/null; then
    srun $cmd
else
    $cmd
fi


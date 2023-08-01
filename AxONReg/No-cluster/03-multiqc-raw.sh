#!/bin/sh -e

#   Dependencies:
#       Requires raw FastQC results.  Run after *-qc-raw.sbatch.

# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
# export LC_ALL=en_US.UTF8
LC_ALL=en_US.UTF-8
LANG=en_US.utf-8
export LC_ALL LANG

# Run interactively under SLURM if srun is found, otherwise run directly
set -x
multiqc --outdir Results/03-multiqc-raw Results/02-qc-raw \
    2>&1 | tee Logs/03-multiqc-raw/multiqc.out

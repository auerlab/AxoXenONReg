#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw and trimmed data for comparison
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   Dependencies:
#       Requires directory structure and reference transcriptome.
#       Run after 01-organize.sh and 07-reference.sh.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

# This job is single-threaded, so there is no difference between
# the SLURM script and running directly (except the #SBATCH comments
# to set SLURM parameters).

if which sbatch; then
    sbatch SLURM/08-kallisto-index.sbatch
else
    SLURM/08-kallisto-index.sbatch \
	2>&1 | tee Logs/08-kallisto-index/08-kallisto-index.out
fi

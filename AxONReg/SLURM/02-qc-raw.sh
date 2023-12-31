#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw data for comparison
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   Dependencies:
#       Requires directory structure.  Run after 01-organize.sh.
##########################################################################

sbatch SLURM/02-qc-raw.sbatch

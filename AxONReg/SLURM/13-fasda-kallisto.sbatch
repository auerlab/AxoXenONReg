#!/bin/sh -e

#
#   Dependencies:
#       Requires kallisto abundances.  Run after *-kallisto-quant.sbatch.

#SBATCH --output=Logs/13-fasda-kallisto/slurm-%A_%a.out
#SBATCH --error=Logs/13-fasda-kallisto/slurm-%A_%a.err

# Document software versions used for publication
uname -a
fasda --version
pwd

kallisto_dir=Results/09-kallisto-quant
fasda_dir=Results/13-fasda-kallisto

# Ex. sample16-time000-rep1
for condition in time000 time026 time050 time055 time120; do
    printf "Normalizing $condition...\n"
    set -x
    time fasda normalize \
	--output $fasda_dir/$condition-all-norm.tsv \
	$kallisto_dir/sample*-$condition-rep*/abundance.tsv
    set +x
done

printf "Computing fold-change...\n"
time fasda fold-change \
    --output $fasda_dir/time000-time026-FC.txt \
    $fasda_dir/time000-all-norm.tsv \
    $fasda_dir/time026-all-norm.tsv

time fasda fold-change \
    --output $fasda_dir/time000-time050-FC.txt \
    $fasda_dir/time000-all-norm.tsv \
    $fasda_dir/time050-all-norm.tsv

time fasda fold-change \
    --output $fasda_dir/time000-time055-FC.txt \
    $fasda_dir/time000-all-norm.tsv \
    $fasda_dir/time055-all-norm.tsv

time fasda fold-change \
    --output $fasda_dir/time000-time120-FC.txt \
    $fasda_dir/time000-all-norm.tsv \
    $fasda_dir/time120-all-norm.tsv

time fasda fold-change \
    --output $fasda_dir/time026-time050-FC.txt \
    $fasda_dir/time026-all-norm.tsv \
    $fasda_dir/time050-all-norm.tsv

time fasda fold-change \
    --output $fasda_dir/time026-time055-FC.txt \
    $fasda_dir/time026-all-norm.tsv \
    $fasda_dir/time055-all-norm.tsv

time fasda fold-change \
    --output $fasda_dir/time026-time120-FC.txt \
    $fasda_dir/time026-all-norm.tsv \
    $fasda_dir/time120-all-norm.tsv

time fasda fold-change \
    --output $fasda_dir/time050-time055-FC.txt \
    $fasda_dir/time050-all-norm.tsv \
    $fasda_dir/time055-all-norm.tsv

time fasda fold-change \
    --output $fasda_dir/time050-time120-FC.txt \
    $fasda_dir/time050-all-norm.tsv \
    $fasda_dir/time120-all-norm.tsv

time fasda fold-change \
    --output $fasda_dir/time055-time120-FC.txt \
    $fasda_dir/time055-all-norm.tsv \
    $fasda_dir/time120-all-norm.tsv


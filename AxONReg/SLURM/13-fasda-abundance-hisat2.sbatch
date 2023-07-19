#!/bin/sh -e

#
#   Dependencies:
#       Requires hisat2 alignments.  Run after *-hisat2-align.sbatch.

# Computing abundances takes a little while, so we do this in parallel
# One core for fasda abundance, one for samtools
#SBATCH --array=1-15
#SBATCH --cpus-per-task=2
#SBATCH --output=Logs/13-fasda-abundance-hisat2/slurm-%A_%a.out
#SBATCH --error=Logs/13-fasda-abundance-hisat2/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
fasda --version
samtools --version
pwd

hisat2_dir=Results/12-hisat2-align
abundance_dir=Results/13-fasda-abundance-hisat2

# Need to fetch gff3 for computing abundances
gff_filename=$(Reference/gff-filename.sh)
reference_dir=Results/07-reference

##########################################################################
#   Compute abundances
##########################################################################

sample=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
bam_file=`ls $hisat2_dir/sample$sample-*.bam`
abundance_file=$abundance_dir/$(basename ${bam_file%.bam}-abundance.tsv)

fasda abundance \
    --feature-type transcript \
    --output-dir $abundance_dir \
    --exact \
    150 \
    $reference_dir/$gff_filename \
    $bam_file

# Rename and then generate a stringtie abundance file with the same name
head $abundance_file
mv $abundance_file ${abundance_file%.tsv}-exact.tsv

fasda abundance \
    --feature-type transcript \
    --output-dir $abundance_dir \
    150 \
    $reference_dir/$gff_filename \
    $bam_file

head $abundance_file

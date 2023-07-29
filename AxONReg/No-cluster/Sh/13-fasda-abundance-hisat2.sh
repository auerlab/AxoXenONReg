#!/bin/sh -e

#
#   Dependencies:
#       Requires hisat2 alignments.  Run after *-hisat2-align.sbatch.

if [ $# != 1 ]; then
    printf "Usage: $0 bam-file\n"
    exit 1
fi

# Document software versions used for publication
uname -a
fasda --version
samtools --version
pwd

hisat2_dir=Results/12-hisat2-align
abundance_dir=Results/13-fasda-abundance-hisat2

# Need to fetch gff3 for computing abundances
reference_dir=Results/07-reference

##########################################################################
#   Compute abundances
##########################################################################

bam_file=$1
abundance_file=$abundance_dir/$(basename ${bam_file%.bam}-abundance.tsv)

fasda abundance \
    --feature-type transcript \
    --output-dir $abundance_dir \
    150 \
    $reference_dir/reference.gff3 \
    $bam_file

column -t $abundance_file | head

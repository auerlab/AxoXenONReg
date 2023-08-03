#!/bin/sh -e

##########################################################################
#   Description:
#       Compute abundances from hisat2 alignments for one condition.
#       
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0\n"
    exit 1
}

##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    printf "Usage: $0 bam-file\n"
    exit 1
fi
bam_file=$1

# Document software versions used for publication
hostname
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

abundance_file=$abundance_dir/$(basename ${bam_file%.bam}-abundance.tsv)

fasda abundance \
    --feature-type transcript \
    --output-dir $abundance_dir \
    150 \
    $reference_dir/reference.gff3 \
    $bam_file

column -t $abundance_file | head

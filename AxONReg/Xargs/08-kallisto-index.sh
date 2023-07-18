#!/bin/sh -e

# Document software versions used for publication
uname -a
kallisto version
samtools --version
pwd

ref_dir=Results/07-reference
transcriptome=$(../Common/transcriptome-filename.sh)
index=${transcriptome%.fa}.index

printf "Using reference $transcriptome...\n"

# Needed for kallisto --genomebam
if [ ! -e $ref_dir/$transcriptome.fai ]; then
    printf "Building samtools index...\n"
    samtools faidx $ref_dir/$transcriptome
fi
printf "Building kallisto index...\n"
set -x
kallisto index --index=Results/08-kallisto-index/$index $ref_dir/$transcriptome

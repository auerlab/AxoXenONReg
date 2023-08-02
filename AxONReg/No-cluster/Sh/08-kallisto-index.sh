#!/bin/sh -e

##########################################################################
#   Description:
#       Build kallisto index for reference transcriptome.
#
#   Dependencies:
#       Requires trimmed reads and reference transcriptome.
#       Run after *-trim.sh and *-reference.sh.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

# Document software versions used for publication
hostname
uname -a
kallisto version
samtools --version
pwd

ref_dir=Results/07-reference
index_dir=Results/08-kallisto-index
transcriptome=transcriptome-reference.fa
printf "Using reference $transcriptome...\n"

# Needed for kallisto --genomebam
if [ ! -e $ref_dir/$transcriptome.fai ]; then
    printf "Building samtools index...\n"
    samtools faidx $ref_dir/$transcriptome
fi

index=$index_dir/${transcriptome%.fa}.index
if [ ! -e $index ]; then
    threads=$(../../Common/get-hw-threads.sh)
    printf "Building kallisto index...\n"
    set -x
    kallisto index --threads=$threads --index=$index $ref_dir/$transcriptome
else
    printf "$index already exists.  Remove it and run again to recreate.\n"
fi

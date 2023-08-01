#!/bin/sh -e

##########################################################################
#   Description:
#       Build kallisto index for reference transcriptome.
#
#   Dependencies:
#       Requires reference transriptome.  Run after *-reference.sbatch.
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

# Document software versions used for publication
hostname
uname -a
kallisto version
samtools --version
pwd

ref_dir=Results/07-reference
transcriptome=transcriptome-reference.fa
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

#!/bin/sh -e

##########################################################################
#   Description:
#       Run STAR aligner on one trimmed sample
#
#   Dependencies:
#       Requires trimmed reads.  Run after *-trim.sh.
#
#   History:
#   Date        Name        Modification
#   2023-08     Jason Bacon Begin
##########################################################################

if [ $# != 1 ]; then
    printf "Usage: $0 file-R1.fastq\n"
    exit 1
fi
reads_file1=$(echo $1 | sed -e 's|Results/|../../|')
reads_file2=${reads_file1%-R1.fastq.zst}-R2.fastq.zst
base=$(basename $reads_file1)
sample=$(echo $base | cut -d - -f 1)

# Document software versions used for publication
hostname
uname -a
printf "STAR "
STAR --version
samtools --version
pwd

# FIXME: Do this in *-reference.sbatch
printf "Running STAR aligner...\n"
mkdir -p Results/12-star-align/$sample
cd Results/12-star-align/$sample
set -x

# Perhaps after configuring shared memory settings on nodes
# --genomeLoad LoadAndKeep

STAR \
    --runThreadN 1 \
    --genomeDir ../../11-star-index \
    --readFilesIn $reads_file1 $reads_file2 \
    --readFilesCommand zstdcat
ls -l

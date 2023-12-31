#!/bin/sh -e

##########################################################################
#   Description:
#       Run hisat2 aligner on each RNA sample.
#
#   Dependencies:
#       Requires hisat2 index.  Run after *-hisat2-index.sbatch.
##########################################################################

if [ $# -ne 2 ]; then
    printf "Usage: $0 threads input-file-R1.fastq.zst\n"
    exit 1
fi

printf "===\n"

threads=$1
input_file=$2

base=$(basename $input_file | cut -d - -f 1-3)

date=$(date +%Y-%m-%d-%H:%M)

output_log=Logs/12-hisat2-align/xargs-$base-$date.out
error_log=Logs/12-hisat2-align/xargs-$base-$date.err

##############################################################################
# Align with hisat2, which can handle splice junctions in RNA reads

# Document software versions used for publication
hostname
uname -a
hisat2 --version
pwd

# Hisat2 index has the same name as the genome reference
index=genome-reference.fa

# samtools sort dumps temp files in CWD
cd Results/12-hisat2-align
output_log="../../$output_log"
error_log="../../$error_log"

zst1="../../$input_file"
zst2=$(echo $zst1 | sed -e 's|R1|R2|g')
ls $zst1 $zst2

base=$(basename $input_file)
bam=${base%-R*}.bam

# Hisat2 performs seeks on the input, so we can't use a pipe
gzip1=${zst1%.zst}.gz
gzip2=${zst2%.zst}.gz
# Use cheap compression level here, since these are removed
# right after the alignment
if [ ! -e $gzip1 ]; then
    printf "Recompressing $zst1...\n"
    zstdcat $zst1 | gzip -1 > $gzip1 &
fi
if [ ! -e $gzip2 ]; then
    printf "Recompressing $zst2...\n"
    zstdcat $zst2 | gzip -1 > $gzip2
fi
# Wait for zst1 recompress to complete
wait

set -x
hisat2 --threads $threads \
    --time \
    --met-stderr \
    -x ../11-hisat2-index/$index \
    -1 $gzip1 -2 $gzip2 | samtools sort -o $bam

# No further need for the non-zstd
# rm -f $gzip1 $gzip2

# Not sure how helpful multithreading is here, but since we allocated
# the cores for hisat2, might as well use them
samtools index -c -@ $threads $bam

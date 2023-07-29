#!/bin/sh -e

##########################################################################
#   Description:
#       Run hisat2 aligner on each RNA sample.
#
#   Dependencies:
#       Requires hisat2 index.  Run after *-hisat2-index.sbatch.
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

# Set job array to number of samples.
#SBATCH --array=1-15
#SBATCH --cpus-per-task=16
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or run under a tool that
# reports maximum memory use.
# xenopus only needs about 3g, but we're using all the cores on a node anyway
#SBATCH --mem=50g
#SBATCH --output=Logs/12-hisat2-align/slurm-%A_%a.out
#SBATCH --error=Logs/12-hisat2-align/slurm-%A_%a.err

##############################################################################
# Align with hisat2, which can handle splice junctions in RNA reads

hw_threads=$(../Common/get-hw-threads.sh)

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_CPUS_PER_TASK:=$hw_threads}

# Document software versions used for publication
uname -a
hisat2 --version
pwd

# Hisat2 index has the same name as the genome reference
index=genome-reference.fa

# samtools sort dumps temp files in CWD
cd Results/12-hisat2-align

sample=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
zst1=$(echo ../04-trim/sample$sample-*-R1.fastq.zst)
zst2=$(echo ../04-trim/sample$sample-*-R2.fastq.zst)

# Hisat2 performs seeks on the input, so we can't use a pipe
gzip1=${zst1%.zst}.gz
gzip2=${zst2%.zst}.gz
printf "Recompressing FASTQs in gzip format for hisat2...\n"
# Use cheap compression level here, since these are removed
# right after the alignment
if [ ! -e $gzip1 ]; then
    zstdcat $zst1 | gzip -1 > $gzip1 &
fi
if [ ! -e $gzip2 ]; then
    zstdcat $zst2 | gzip -1 > $gzip2
fi
wait

zstbase=$(basename $zst1)
bam=${zstbase%-R*}.bam
threads=$(($SLURM_CPUS_PER_TASK - 1))

set -x
hisat2 --threads $threads \
    --time \
    --met-stderr \
    -x ../11-hisat2-index/$index \
    -1 $gzip1 -2 $gzip2 | samtools sort -T /tmp -o $bam

# No further need for the uncompressed FASTQs
rm -f $gzip1 $gzip2

# Not sure how helpful multithreading is here, but since we allocated
# the cores for hisat2, might as well use them
samtools index -c -@ $SLURM_CPUS_PER_TASK $bam

#!/bin/sh -e

##########################################################################
#   Description:
#       Run hisat2 aligner on each RNA sample.
#
#   Dependencies:
#       Requires hisat2 index.  Run after *-hisat2-index.sbatch.
##########################################################################

# Set job array to number of samples.
#SBATCH --array=1-15
#SBATCH --cpus-per-task=16
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or run under a tool that
# reports maximum memory use.
# xenopus only needs about 3g, but we're using all the cores on a node anyway
# #SBATCH --mem=50g
#SBATCH --output=Logs/12-hisat2-align/slurm-%A_%a.out
#SBATCH --error=Logs/12-hisat2-align/slurm-%A_%a.err

if [ -z "$SLURM_MEM_PER_NODE" ]; then
    printf "Usage: sbatch --mem=X $0\n"
    printf "Use ../12-hisat2-align.sh rather than sbatch.\n"
    exit 1
fi

##############################################################################
# Align with hisat2, which can handle splice junctions in RNA reads

hw_threads=$(../../Common/get-hw-threads.sh)

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_CPUS_PER_TASK:=$hw_threads}

# Document software versions used for publication
hostname
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

# Set this to a local partition with plenty of space or a high-speed
# networked filesystem
tmp_dir=/usr/tmp

# Hisat2 performs seeks on the input, so we can't use a pipe
# Place gzip copies on local drives for better parallism
base1=$(basename $zst1)
base2=$(basename $zst2)
gzip1=$tmp_dir/${base1%.zst}.gz
gzip2=$tmp_dir/${base2%.zst}.gz
bam=${base1%-R*}.bam

if [ ! -e $bam ]; then
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
    
    threads=$(($SLURM_CPUS_PER_TASK - 1))
    set -x
    hisat2 --threads $threads \
	--time \
	--met-stderr \
	-x ../11-hisat2-index/$index \
	-1 $gzip1 -2 $gzip2 | samtools sort -T $tmp_dir -o $bam
    
    # No further need for the gzipped fastas
    rm -f $gzip1 $gzip2
else
    printf "$bam already exists.  Remove it and run again to regenerate it.\n"
fi

# Not sure how helpful multithreading is here, but since we allocated
# the cores for hisat2, might as well use them
if [ ! -e $bam.csi ]; then
    samtools index -c -@ $SLURM_CPUS_PER_TASK $bam
fi

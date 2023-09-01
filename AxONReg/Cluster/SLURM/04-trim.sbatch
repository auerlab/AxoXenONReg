#!/bin/sh -e

##########################################################################
#   Script description:
#       Trim adapters and low quality ends from Lumina reads
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   Dependencies:
#       Requires directory structure.  Run after *-organize.sh.
##########################################################################

##########################################################################
# Cutadapt:
# Each job in the array will run a cutadapt (python) process and a
# compression process for part of the time.  If you don't want to
# oversubscribe compute nodes even for a little while, add --cpus-per-task=2
# There may be 2 pigz processes per job, but --cpus-per-task=3 doesn't help

##########################################################################
# Fastq-trim:
# Limit concurrent jobs to avoid becoming I/O-bound.
# Fastq-trim is so fast it ends up using only about 40% CPU while waiting
# for NFS on albacore (only gigabit Ethernet).  Clusters with higher
# speed networks and file servers might handle more jobs gracefully.
# Processes per job include 2 xzcat, 2 gzip, and 1 fastq-trim, but xzcat
# and gzip use less than half a core each
 
#SBATCH --array=1-15%5
#SBATCH --cpus-per-task=3
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem-per-cpu=250
#SBATCH --output=Logs/04-trim/slurm-%A_%a.out
#SBATCH --error=Logs/04-trim/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
fastq-trim --version
pwd

sample=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
input1=$(ls Results/01-organize/Raw-renamed/sample$sample-*R1*)
input2=$(ls Results/01-organize/Raw-renamed/sample$sample-*R2*)
base=$(basename $input1)
stem=${base%%-R*.fastq.xz}

# https://www.rootusers.com/gzip-vs-bzip2-vs-xz-performance-comparison/
# xz offers the best compression by far, but is slow at mid (-5) to high (-9)
# compression levels.  At -1, xz is faster than bzip2 while providing
# comparable compression.  If you want even faster compression and are
# willing to sacrifice compression ratio, use .zst or no compression for
# outputs instead.  A ZFS filesystem with lz4 or zstd compression will provide
# enough compression for intermediate files without gzip, bzip2, or xz.
# However, this may cause a network bottleneck as all processes write
# uncompressed data over NFS to be compressed by the file server.
# Using compression on the compute nodes distributes the compression
# workload and reduces the volume of data transfered over NFS.

#suffix=.gz
suffix=.zst
output1=Results/04-trim/${stem}-R1.fastq$suffix
output2=Results/04-trim/${stem}-R2.fastq$suffix

# Maximize compression throughput so gzip is not a bottleneck.  Can determine
# performance from fastq-trim CPU % in job-top.
# Don't use -1 if job is I/O bound.  Use the idle CPU to get better compression
# and reduce I/O.  The goal is to maximum CPU utilization of the fastq-trim
# process and the gzip processes, as shown by job-top.
# -4 seems to work best on albacore with its
# gigabit network.  Lower values are probably better with a high-speed network.
# Use the highest value that does not result in a lower CPU utilization
# for fastq-trim.
# export GZIP=-4

adapter=AGATCGGAAGAG

# fastq-trim is 2.5x faster with 1 core than cutadapt with 2 cores
# https://github.com/auerlab/fastq-trim

set -x      # Show full command in log file
time fastq-trim --3p-adapter1 $adapter --3p-adapter2 $adapter \
    --min-qual 24 --polya-min-length 4 $input1 $output1 $input2 $output2

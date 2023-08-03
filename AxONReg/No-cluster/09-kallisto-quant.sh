#!/bin/sh -e

##########################################################################
#   Script description:
#       Perform kallisto pseudo-alignments.
#
#   Dependencies:
#       Requires trimmed reads and kallisto index.  Run after
#       *-trim.sh and *-kallisto-index.sh.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

# Kallisto multithreading scales fairly well, so more jobs with
# fewer cores each won't lead to much better performance.  If might
# reduce performance by increasing competition for disk bandwidth.
# Using minimum number of jobs with at least 2 GB each to ensure
# plenty of RAM per job, or 4 cores per job, whichever is greater.
hw_threads=$(../../Common/get-hw-threads.sh)
hw_mem=$(../../Common/get-hw-mem.sh)

# Find optimal # cores per job with at least 2 GB per job
jobs=$hw_threads
mem_per_job=$(($hw_mem / $jobs ))
while [ $mem_per_job -lt 2000000000 ]; do
    jobs=$(( $jobs / 2 ))
    mem_per_job=$(($hw_mem / $jobs ))
done
threads_per_job=$(( $hw_threads / $jobs ))

# Use at least 4 thread per job, since the jobs will run almost 4 times
# as fast, and fewer jobs means less disk contention.  If CPU utilization
# is less than 90% per core (360% for a 4-thread job), you probably have
# too much disk contention.
if [ $threads_per_job -lt 4 ]; then
    threads_per_job=4
    jobs=$(( $hw_threads / $threads_per_job ))
fi

# Tried GNU parallel and ran into bugs.  Xargs just works.
printf "Running $jobs jobs with $threads_per_job threads each...\n"
ls Results/04-trim/*-R1.fastq.zst | \
    xargs -n 1 -P $jobs ../../Common/redirect.sh \
    Sh/09-kallisto-quant.sh $threads_per_job

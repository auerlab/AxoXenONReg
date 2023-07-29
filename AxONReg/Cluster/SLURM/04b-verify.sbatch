#!/bin/sh -e

#SBATCH --array=1-15
#SBATCH --cpus-per-task=2
#SBATCH --output=Logs/04b-verify/fastx-stats-%A_%a.out
#SBATCH --error=Logs/04b-verify/fastx-stats-%A_%a.err

# Default job array index for testing
: ${SLURM_ARRAY_TASK_ID:=1}

sample=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
raw1=$(ls Results/01-organize/Raw-renamed/sample$sample*-R1.fastq.xz)
raw2=$(ls Results/01-organize/Raw-renamed/sample$sample*-R2.fastq.xz)
trimmed1=$(ls Results/04-trim/sample$sample*-R1.fastq.zst)
trimmed2=$(ls Results/04-trim/sample$sample*-R2.fastq.zst)
    
printf "Scanning $raw1...\n"
blt fastx-stats $raw1

printf "Scanning $trimmed1...\n"
blt fastx-stats $trimmed1

printf "Scanning $raw2...\n"
blt fastx-stats $raw2

printf "Scanning $trimmed2...\n"
blt fastx-stats $trimmed2

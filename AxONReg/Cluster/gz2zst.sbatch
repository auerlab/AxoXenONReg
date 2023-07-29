#!/bin/sh -e

#SBATCH --array=16-30

: ${SLURM_ARRAY_TASK_ID:=16}

cd Results/04-trim
for read in 1 2; do
    gzip=$(ls sample${SLURM_ARRAY_TASK_ID}-*-R$read.fastq.gz)
    zst=${gzip%.gz}.zst
    echo $gzip "->" $zst
    gunzip -c $gzip | zstd > $zst
done

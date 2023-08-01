#!/bin/sh -e

# Raw files from sequencing center
fastq=$1

# Filename stems for fastqc output
stem_fastq=$(basename ${fastq%.fastq.xz})
results=Results/02-qc-raw/$(basename ${stem_fastq})_fastqc.html

if [ -e $results ]; then
    printf "$fastq already processed.\n"
else
    printf "Processing $fastq with fastqc...\n"
    
    # Document software versions used for publication
    hostname
    uname -a
    fastqc --version
    pwd

    xzcat $fastq | fastqc -o Results/02-qc-raw stdin:$stem_fastq
fi

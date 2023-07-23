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
    
    log_stem=Logs/02-qc-raw/$stem_fastq
    # Document software versions used for publication
    uname -a > $log_stem.out
    fastqc --version >> $log_stem.out
    pwd >> $log_stem.out

    xzcat $fastq | fastqc -o Results/02-qc-raw stdin:$stem_fastq \
	>> $log_stem.out 2>> $log_stem.err
fi

#!/bin/sh -e

# Trimmed files
fastq=$1

# Filename stems for fastqc output
stem_fastq=$(basename ${fastq%.fastq.xz})
results=Results/05-qc-trimmed/$(basename ${stem_fastq})_fastqc.html

if [ -e $results ]; then
    printf "$fastq already processed.\n"
else
    printf "Processing $fastq with fastqc...\n"
    
    log_stem=Logs/05-qc-trimmed/$stem_fastq
    # Document software versions used for publication
    hostname
    uname -a > $log_stem.out
    fastqc --version >> $log_stem.out
    pwd >> $log_stem.out

    zstdcat $fastq | fastqc -o Results/05-qc-trimmed stdin:$stem_fastq \
	>> $log_stem.out 2>> $log_stem.err
fi

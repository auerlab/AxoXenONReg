#!/bin/sh -e

# Trimmed files
trimmed=$1

# Filename stems for fastqc output
stem_trimmed=$(basename ${trimmed%.fastq.xz})
results=Results/05-qc-trimmed/$(basename ${stem_trimmed})_fastqc.html

if [ -e $results ]; then
    printf "$trimmed already processed.\n"
else
    printf "Processing $trimmed with fastqc...\n"
    
    log_stem=Logs/05-qc-trimmed/$stem_trimmed
    # Document software versions used for publication
    uname -a > $log_stem.out
    fastqc --version >> $log_stem.out
    pwd >> $log_stem.out

    zstdcat $trimmed | fastqc -o Results/05-qc-trimmed stdin:$stem_trimmed \
	>> $log_stem.out 2>> $log_stem.err
fi

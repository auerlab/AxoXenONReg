#!/bin/sh -e

#fastq-scum Raw/A16Na1/A16Na1_1.fq.gz 
#fastq-scum Raw/A16Na1/A16Na1_2.fq.gz

export GZIP=-5
fastq-trim \
    --3p-adapter1 AGATCGGAAGAG \
    --3p-adapter2 AGATCGGAAGAG \
    Raw/A16Na1/A16Na1_1.fq.gz out_1.fq.gz \
    Raw/A16Na1/A16Na1_2.fq.gz out_2.fq.gz


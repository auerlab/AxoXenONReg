#!/bin/sh -e

##########################################################################
#   Description:
#       Perform quality checks on one trimmed fastq file.
#       
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 trimmed-file.fastq\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi
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
    uname -a
    fastqc --version
    pwd

    zstdcat $fastq | fastqc -o Results/05-qc-trimmed stdin:$stem_fastq
fi

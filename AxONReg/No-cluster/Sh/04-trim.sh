#!/bin/sh -e

raw1=$1
base=$(basename $raw1)
stem=${base%%-R1.fastq.xz}
trimmed1=Results/04-trim/$stem-R1.fastq.zst

raw2=${raw1%%-R1.fastq.xz}-R2.fastq.xz
trimmed2=Results/04-trim/$stem-R2.fastq.zst

if [ -e $trimmed1 ]; then
    printf "$raw1 + $raw2 already processed.\n"
else
    # Document software versions used for publication
    hostname
    uname -a
    fastq-trim --version
    pwd

    adapter=AGATCGGAAGAG
    
    # fastq-trim is 2.5x faster with 1 core than cutadapt with 2 cores
    # Our reads use the default Illumina universal adapter,
    # but we'll state it explicitly anyway
    # export GZIP=-1
    set -x
    time fastq-trim --3p-adapter1 $adapter --3p-adapter2 $adapter \
	--min-qual 24 --polya-min-length 4 $raw1 $trimmed1 $raw2 $trimmed2
fi

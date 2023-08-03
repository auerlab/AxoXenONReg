#!/bin/sh -e

##########################################################################
#   Description:
#       Run kallisto quantification for one trimmed reads file.
#
#   Dependencies:
#       Requires kallisto index.  Run after *-kallisto-index.sbatch.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

if [ $# != 2 ]; then
    printf "Usage: $0 threads forward-sample-name.fastq.zst\n"
    exit 1
fi
threads=$1
zst1=$2

zst2=${zst1%R1.fastq.zst}R1.fastq.zst
base1=$(basename $zst1)
base2=$(basename $zst2)

# Kallisto requires an output subdirectory for each sample
stem=$(basename ${zst1%-R1*})
out_dir=Results/09-kallisto-quant/$stem
mkdir -p $out_dir

log_stem=Logs/09-kallisto-quant/$stem

# Document software versions used for publication
hostname
uname -a
kallisto version
pwd

# kallisto only supports gzip compression as of 0.48.0, so use FIFOs
# feed it raw fastq
pipe1=/tmp/pipe1-kallisto-$base1
pipe2=/tmp/pipe2-kallisto-$base2
rm -f $pipe1 $pipe2
mkfifo $pipe1 $pipe2
zstdcat $zst1 > $pipe1 &
zstdcat $zst2 > $pipe2 &

set -x
time kallisto quant \
    --threads=$threads \
    --index=Results/08-kallisto-index/transcriptome-reference.index \
    --output-dir=$out_dir $pipe1 $pipe2
rm -f $pipe1 $pipe2

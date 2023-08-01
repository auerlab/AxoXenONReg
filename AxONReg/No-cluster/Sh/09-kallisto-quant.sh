#!/bin/sh -e

##########################################################################
#   Description:
#       Run kallisto quantification for each RNA sample.
#
#   Dependencies:
#       Requires kallisto index.  Run after *-kallisto-index.sbatch.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
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

# kallisto 0.46.1 can't handle xz and will simply seg fault rather than
# issue an error message.  If your trimmed fastq files are in xz format,
# this will convert to zstd format.
# Convert xz to zst rather than raw to reduce NFS load from compute nodes
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
# feed it raw text
pipe1=/tmp/pipe1-kallisto-$base1
pipe2=/tmp/pipe2-kallisto-$base2
printf "Pipes: $pipe1 $pipe2\n"

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

#!/bin/sh -e

##########################################################################
#   Script description:
#       Organize files and create directories needed for analysis
#
#       Link raw files to standardized names for that clearly describe
#       conditions and replicates, and are easily parsed by subsequent
#       scripts.
#
#       Researchers involved in sample-prep generally don't think about
#       how filename conventions impact bioinformatics analysis, so this
#       simple step can avoid confusion throughout the pipeline.  In
#       addition, linking this way can correct for sample mixups, etc.
#
#       Use links to preserve the original files and document the mapping.
#
#   Dependencies:
#       Creates required directory structure.
#       Run this before all other scripts.
#       
#   History:
#   Date        Name        Modification
#   2021-09-25  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi

mkdir -p Results Logs
scripts=$(ls 0[2-9]-*)
for script in $scripts; do
    stage=${script%.*}
    mkdir -p Results/$stage Logs/$stage
done

##############################################################################
# RNA-Seq:
#
#   Explain sequence file naming here
##############################################################################

# CE1A_S1_L002-R1.fastq.xz
cd Results
rm -rf 01-organize/Raw-renamed
mkdir -p 01-organize/Raw-renamed
cd 01-organize/Raw-renamed
pwd
for path in ../../../Raw/*/*.fq.gz; do
    file=$(basename $path)
    # FIXME: Generate a non-cryptic name for each file
    sample=$(echo $file | cut -c 2-3)
    sample=$(($sample - 15))
    if echo $file | fgrep -q Na; then
	time=0
	rep=$(echo $file | cut -c 7)
	strand=$(echo $file | cut -c 9)
    elif echo $file | fgrep -q h; then
	time=$(echo $file | cut -c 4-5)
	rep=$(echo $file | cut -c 7)
	strand=$(echo $file | cut -c 9)
    else
	time=$(echo $file | cut -c 4)
	time=$(($time * 24))
	rep=$(echo $file | cut -c 6)
	strand=$(echo $file | cut -c 8)
    fi
    printf "$path = sample-$sample-cond-$time-rep-$rep-$strand.fq.xz\n"
done

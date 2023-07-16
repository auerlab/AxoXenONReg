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
#   2023-06     Jason Bacon Begin
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

cd Results
rm -rf 01-organize/Raw-renamed
mkdir -p 01-organize/Raw-renamed
cd 01-organize/Raw-renamed
pwd
for path in ../../../Raw/*/*.fq.gz; do
    file=$(basename $path)
    # FIXME: Generate a non-cryptic name for each file
    if [ $(echo $file | cut -c 3-4) = Na ]; then
	sample=$(echo $file | cut -c 2)
	day=0
	rep=$(echo $file | cut -c 5)
	strand=$(echo $file | cut -c 7)
    elif [ $(echo $file | cut -c 3) = d ]; then
	sample=$(echo $file | cut -c 2)
	if [ $(echo $file | cut -c 6) = '_' ]; then
	    day=$(echo $file | cut -c 4)
	    rep=$(echo $file | cut -c 5)
	    strand=$(echo $file | cut -c 7)
	else
	    day=$(echo $file | cut -c 4-5)
	    rep=$(echo $file | cut -c 6)
	    strand=$(echo $file | cut -c 8)
	fi
    elif [ $(echo $file | cut -c 4) = d ]; then
	sample=$(echo $file | cut -c 2-3)
	if [ $(echo $file | cut -c 7) = '_' ]; then
	    day=$(echo $file | cut -c 4-5)
	    rep=$(echo $file | cut -c 6)
	    strand=$(echo $file | cut -c 8)
	else
	    day=$(echo $file | cut -c 5-6)
	    rep=$(echo $file | cut -c 7)
	    strand=$(echo $file | cut -c 9)
	fi
    fi
    day=$(printf "%02d" $day)
    sample=$(printf "%02d" $sample)
    rep=$(echo $rep | tr "ABC" "123")
    link=sample$sample-day$day-rep$rep-R$strand.fastq.gz
    ln -sf $path $link
    printf "$path = $link\n"
done

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
scripts=$(ls 0[2-9]-* 1[0-9]-* 2[0-9]-*)
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
for path in ../../../Raw/*/*.fq.xz; do
    file=$(basename $path)
    # FIXME: Generate a non-cryptic name for each file
    sample=$(echo $file | cut -c 2-3)
    # Normalize to match xenopus, so we can use the same scripts
    sample=`printf "%02d" $(($sample - 15))`
    if echo $file | fgrep -q Na; then
	time=0
	rep=$(echo $file | cut -c 6)
	strand=$(echo $file | cut -c 8)
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
    
    # Reduce time points to ranks to simplify analysis and allow
    # sharing scripts between axolotl and xenopus.  It can easily
    # be converted back at any time.
    # 0 26 50 55 120
    case $time in
    0)
	time=1
	;;
    26)
	time=2
	;;
    50)
	time=3
	;;
    55)
	time=4
	;;
    120)
	time=5
	;;
    *)
	printf "Invalid time point.\n"
	exit 1
    esac
    link=sample$sample-time$time-rep$rep-R$strand.fastq.xz
    printf "$path = $link\n"
    ln -sf $path $link
done
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
    printf "Usage: $0 full|test\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
mode=$1
if [ $mode != full ] && [ $mode != test ]; then
    usage
fi

mkdir -p Results Logs
scripts=$(ls 0[2-9]*.sh 1[0-9]*.sh)
for script in $scripts; do
    stage=${script%.*}
    mkdir -p Results/$stage Logs/$stage
done

##############################################################################
# RNA-Seq:
#
#   Examples of raw file naming:
#
#   A16Na1  = Axolotl sample 16 naive       replicate 1
#   A2450h3 = Axolotl sample 24 50 hours    replicate 3
##############################################################################

cd Results
rm -rf 01-organize/Raw-renamed
mkdir -p 01-organize/Raw-renamed
cd 01-organize/Raw-renamed
pwd
if [ ! -e ../../../../Raw ]; then
    printf "Missing Raw directory.\n"
    exit 1
fi

for path in ../../../../Raw/*/*.fq.xz; do
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
	cond=1
	;;
    26)
	cond=2
	;;
    50)
	cond=3
	;;
    55)
	cond=4
	;;
    120)
	cond=5
	;;
    *)
	printf "Invalid time point.\n"
	exit 1
    esac
    link=sample$sample-cond$cond-rep$rep-R$strand.fastq.xz
    printf "$path = $link\n"
    if [ $mode = full ]; then
	ln -sf $path $link
    else
	# Just take the first 1,000,000 reads for testing
	xzcat $path | head -n 4000000 | xz -1 > $link
    fi
done

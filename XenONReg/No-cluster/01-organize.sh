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
scripts=$(ls 0[2-9]-* 1[0-9]-*)
for script in $scripts; do
    stage=${script%.*}
    mkdir -p Results/$stage Logs/$stage
done

##############################################################################
# RNA-Seq:
#
#   Examples of raw file naming:
#
#   X13d27A = Xenopus sample 13 27 days     replicate 1
#   X14d27B = Xenopus sample 14 27 days     replicate 2
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
    
    # Reduce time points to ranks to simplify analysis and allow
    # sharing scripts between axolotl and xenopus.  It can easily
    # be converted back at any time.
    # Time offsets: 0 (naive) 26 50 55 120
    case $day in
    0)
	cond=1
	;;
    7)
	cond=2
	;;
    12)
	cond=3
	;;
    18)
	cond=4
	;;
    27)
	cond=5
	;;
    *)
	printf "Invalid time point.\n"
	exit 1
    esac
    sample=$(printf "%02d" $sample)
    rep=$(echo $rep | tr "ABC" "123")
    link=sample$sample-cond$cond-rep$rep-R$strand.fastq.xz
    printf "$path = $link\n"
    if [ $mode = full ]; then
	ln -sf $path $link
    else
	# Just take the first 1,000,000 reads for testing
	xzcat $path | head -n 4000000 | xz -1 > $link
    fi
done

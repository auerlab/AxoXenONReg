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
for rep in 1 2 3; do
done

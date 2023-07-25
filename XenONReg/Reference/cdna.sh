#!/bin/sh -e

if [ $0 != Reference/cdna.sh ]; then 
    printf "$0 must be run\n"
    printf "from inside the Reference directory.\n"
    exit 1
fi

site=https://download.xenbase.org/xenbase/Genomics/Sequences
cdna=xlaevisMRNA.fasta

# Can't guarantee this file will always be available.
# You may need to edit this.
cd Results/07-reference

curl -O --continue-at - $site/$cdna.gz
if [ ! -e $cdna ]; then
    gunzip --keep $cdna.gz
fi
ln -s $cdna transcriptome-reference.fa

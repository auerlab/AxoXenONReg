#!/bin/sh -e

cd Results/07-reference
rm -f $gff
site=https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1
gff=XENLA_10.1_Xenbase.gff3

printf "Fetching $gff...\n"
curl -O --continue-at - $site/$gff.gz
if [ ! -e $gff ]; then
    printf "Uncompressing...\n"
    gunzip --keep $gff.gz
fi
ln -sf $gff reference.gff3

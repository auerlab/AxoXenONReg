#!/bin/sh -e

cd Results/07-reference

#   https://www.xenbase.org/xenbase/static-xenbase/ftpDatafiles.jsp
version=10.1

# GFF3
site=https://download.xenbase.org/xenbase/Genomics/JGI/Xenla$version/
gff=XENLA_${version}_Xenbase.gff3
curl -O --continue-at - $site/$gff.gz
if [ ! -e $gff ]; then
    gunzip --keep $gff.gz
fi

# GTF
gtf=XENLA_${version}_Xenbase.gtf
curl -O --continue-at - $site/$gtf.gz
if [ ! -e $gtf ]; then
    gunzip --keep $gtf.gz
fi

# Genome
genome=XENLA_10.1_genome.fa
curl -O --continue-at - $site/$genome.gz
if [ ! -e $genome ]; then
    gunzip --keep $genome.gz
fi
ln -f $genome genome-reference.fa

# Transcriptome
if [ ! -e transcriptome-reference.fa ]; then
    printf "Building transcriptome from $gtf...\n"
    gffread -w transcriptome-reference.fa -g $genome $gtf
else
    printf "Using existing transcriptome-reference.fa.\n"
fi

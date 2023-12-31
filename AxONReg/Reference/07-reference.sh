#!/bin/sh -e

##########################################################################
#   Description:
#       Download/build genome and transcriptome references.
@s
cd Results/07-reference

##########################################################################
#   Transcriptomes (CDS and cDNA)
##########################################################################

cds=AmexT_v47_cds.fa
if [ ! -e $cds ]; then
    curl --continue-at - --location --output $cds.gz \
	'http://www.axolotl-omics.org/api?method=Assembly.getSequences&assembly=47&type=cds'
    printf "Deflating...\n"
    gunzip $cds.gz
else
    printf "$cds already exists.\n"
fi

dna=AmexT_v47_dna.fa
if [ ! -e $dna ]; then
    curl --continue-at - --location --output $dna.gz \
	'http://www.axolotl-omics.org/api?method=Assembly.getSequences&assembly=47&type=dna'
    printf "Deflating...\n"
    gunzip $dna.gz
else
    printf "$dna already exists.\n"
fi

# Use DNA so that transcript lengths match the GTF/GFF3
ln -sf $dna transcriptome-reference.fa

##########################################################################
#   GTF
##########################################################################

gtf=AmexT_v47-AmexG_v6.0-DD.gtf
if [ ! -e $gtf ]; then
    curl --continue-at - --remote-name \
	https://www.axolotl-omics.org/dl/$gtf.gz
    printf "Deflating...\n"
    gunzip $gtf.gz
else
    printf "$gtf already exists.\n"
fi
ln -sf $gtf reference.gtf

##########################################################################
#   GFF
##########################################################################

gff=AmexT_v47-AmexG_v6.0-DD.gff3
if [ ! -e $gff ]; then
    whole_gff=${gtf%.gtf}.gff3
    printf "Converting $gtf to $whole_gff...\n"
    gffread -o $whole_gff $gtf
    printf "Removing loose contigs -> $gff...\n"
    awk '$1 ~ "chr"' $whole_gff > $gff
else
    printf "$gff already exists.\n"
fi
ln -sf $gff reference.gff3

##########################################################################
#   Genome
#   Do this last, because it could take hours
##########################################################################

fasta=AmexG_v6.0-DD.fa
if [ ! -e $fasta ]; then
    curl --continue-at - --remote-name \
	https://www.axolotl-omics.org/dl/$fasta.gz
    printf "Deflating...\n"
    gunzip $fasta.gz
else
    printf "$fasta already exists.\n"
fi

if [ ! -e chromosome-sizes.tsv ]; then
    printf "Getting chromosome sizes...\n"
    blt chrom-lens < $fasta > chromosome-sizes.tsv
else
    printf "chromosome-sizes.tsv already exists.\n"
fi
ln -sf $fasta genome-reference.fa

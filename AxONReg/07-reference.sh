#!/bin/sh -e

# From Ava email: https://www.axolotl-omics.org/
# "file" command revealed name within gzip data to use with --output
# Firefox "save link as" comes up with .zip, which is wrong.  Should be .gz.

transcriptome=$(../Common/transcriptome-filename.sh)
genome=$(../Common/genome-filename.sh)
gtf=$(Reference/gtf-filename.sh)
gff=$(Reference/gff-filename.sh)

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
ln -sf $cds $transcriptome

dna=AmexT_v47_dna.fa
if [ ! -e $dna ]; then
    curl --continue-at - --location --output $dna.gz \
	'http://www.axolotl-omics.org/api?method=Assembly.getSequences&assembly=47&type=dna'
    printf "Deflating...\n"
    gunzip $dna.gz
else
    printf "$dna already exists.\n"
fi

##########################################################################
#   GTF
##########################################################################

if [ ! -e $gtf ]; then
    curl --continue-at - --remote-name \
	https://www.axolotl-omics.org/dl/$gtf.gz
    printf "Deflating...\n"
    gunzip $gtf.gz
else
    printf "$gtf already exists.\n"
fi

##########################################################################
#   GFF
##########################################################################

if [ ! -e $gff ]; then
    whole_gff=${gtf%.gtf}.gff3
    printf "Converting $gtf to $whole_gff...\n"
    gffread -o $whole_gff $gtf
    printf "Removing loose contigs -> $gff...\n"
    awk '$1 ~ "chr"' $whole_gff > $gff
else
    printf "$gff already exists.\n"
fi

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

if [ ! -e ${fasta%.fa}-chr-only.fa ]; then
    printf "Generating ${fasta%.fa}-chr-only.fa...\n"
    # FIXME: More than 1 line in each sequence
    awk '{ if ( $1 ~ "^>C" ) { exit } else { print }}' $fasta > ${fasta%.fa}-chr-only.fa
fi

if [ ! -e chromosome-sizes.tsv ]; then
    printf "Getting chromosome sizes...\n"
    blt chrom-lens < $fasta > chromosome-sizes.tsv
else
    printf "chromosome-sizes.tsv already exists.\n"
fi
ln -sf $fasta $genome
ln -sf ${fasta%.fa}-chr-only.fa ${genome%.fa}-chr-only.fa

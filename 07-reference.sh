#!/bin/sh -e

# From Ava email: https://www.axolotl-omics.org/
# "file" command revealed name within gzip data to use with --output
# Firefox "save link as" comes up with .zip, which is wrong.  Should be .gz.

cd Results/07-reference

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

##########################################################################
#   Genome
#   Do this last, because it could take hours
##########################################################################

genome=AmexG_v6.0-DD.fa
if [ ! -e $genome ]; then
    curl --continue-at - --remote-name \
	https://www.axolotl-omics.org/dl/$genome.gz
    printf "Deflating...\n"
    gunzip $genome
else
    printf "$genome already exists.\n"
fi

if [ ! -e chromosome-sizes.tsv ]; then
    printf "Getting chromosome sizes...\n"
    blt chrom-lens < $genome > chromosome-sizes.tsv
else
    printf "chromosome-sizes.tsv already exists.\n"
fi

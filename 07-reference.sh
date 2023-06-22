#!/bin/sh -e

# From Ava email: https://www.axolotl-omics.org/
# "file" command revealed name within gzip data to use with --output
# Firefox "save link as" comes up with .zip, which is wrong.  Should be .gz.

cd Results/07-reference

if [ ! -e AmexT_v47_cds.fa ]; then
    curl --continue-at - --location --output AmexT_v47_cds.fa.gz \
	'http://www.axolotl-omics.org/api?method=Assembly.getSequences&assembly=47&type=cds'
    printf "Deflating...\n"
    gunzip AmexT_v47_cds.fa.gz
fi

if [ ! -e AmexT_v47_dna.fa ]; then
    curl --continue-at - --location --output AmexT_v47_dna.fa.gz \
	'http://www.axolotl-omics.org/api?method=Assembly.getSequences&assembly=47&type=dna'
    printf "Deflating...\n"
    gunzip AmexT_v47_dna.fa.gz
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
    gunzip $genome.gz
fi


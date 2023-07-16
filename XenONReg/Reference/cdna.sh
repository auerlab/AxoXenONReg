#!/bin/sh -e

proper_name=Reference/cdna.sh
if [ $0 != "$proper_name" ]; then
    printf "$0 must be run as $proper_name\n"
    printf "from inside the Reference directory.\n"
    exit 1
fi

# Need GTF for kallisto quant --genomebam in any case
Reference/fetch-gtf.sh

fetch=$(Reference/find-fetch.sh)
build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)
awk=$(Reference/find-awk.sh)
transcriptome=$(../Common/transcriptome-filename.sh)

# Can't guarantee this file will always be available.
# You may need to edit this.
cd Results/07-reference
cdna=Xenopus_tropicalis.UCB_Xtro_$build.cdna.all.fa.gz
if [ ! -e $cdna ]; then
    $fetch https://ftp.ensembl.org/pub/release-$release/fasta/xenopus_tropicalis/cdna/$cdna
else
    printf "$cdna already exists.  Remove and rerun to replace.\n"
fi


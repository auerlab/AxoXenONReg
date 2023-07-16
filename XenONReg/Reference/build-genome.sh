#!/bin/sh -e

if [ $0 != "Reference/build-genome.sh" ]; then
    cat << EOM

$0 must be run as Reference/build-genome.sh.

EOM
    exit 1
fi

printf "Building genome...\n"
fetch=$(../Common/find-fetch.sh)
build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)
genome=$(../Common/genome-filename.sh)

# Chromosome files
mkdir -p Results/07-reference
cd Results/07-reference
for chrom in $(seq 1 10); do
    file=Xenopus_tropicalis.UCB_Xtro_$build.dna.primary_assembly.$chrom.fa.gz
    if [ ! -e $file ]; then
	printf "Fetching $file...\n"
	$fetch http://ftp.ensembl.org/pub/release-$release/fasta/xenopus_tropicalis/dna/$file
    else
	printf "Already have $file...\n"
    fi
    chrom=$((chrom + 1))
done

if [ ! -e $genome ]; then
    printf "Concatenating chromosome FASTAs...\n"
    for chrom in $(seq 1 10); do
	printf "$chrom "
	zcat Xenopus_tropicalis.UCB_Xtro_10.0.dna.primary_assembly.$chrom.fa.gz >> $genome
    done
    printf "\n"
else
    printf "Using existing $genome...\n"
fi

if [ ! -e $genome.fai ]; then
    printf "Creating index $genome.fai...\n"
    samtools faidx $genome      # Speed up gffread
fi

#!/bin/sh -e

##########################################################################
#   Description:
#       Build star index for reference genome.
#
#   Dependencies:
#       Requires a reference genome.  Run after *-reference.sbatch.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

# Document software versions used for publication
hostname
uname -a
printf "STAR "
STAR --version
samtools --version
pwd

# Run STAR genomeGenerate on a copy in 11-star-index
genome=genome-reference.fa
ln -f Results/07-reference/$genome Results/11-star-index
printf "Using reference $genome...\n"

# FIXME: Do this in *-reference.sbatch
if [ ! -e $genome.8.ht2 ]; then
    printf "Building STAR index...\n"
    cd Results/11-star-index
    set -x
    
    STAR --runThreadN $SLURM_CPUS_PER_TASK \
	--limitGenomeGenerateRAM 120000000000 \
	--runMode genomeGenerate \
	--genomeDir . \
	--genomeFastaFiles ../../Results/07-reference/genome-reference.fa \
	--sjdbGTFfile ../../Results/07-reference/reference.gtf \
	--sjdbOverhang 149
fi
ls -l

#!/bin/sh -e

##########################################################################
#   Description:
#       Run kallisto quantification for each RNA sample.
#
#   Dependencies:
#       Requires kallisto index.  Run after *-kallisto-index.sbatch.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#       
#   History:
#   Date        Name        Modification
#   2019-11-??  Jason Bacon Begin
##########################################################################

# Set job array to number of samples.
#SBATCH --array=16-30
#SBATCH --cpus-per-task=4
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=6000
#SBATCH --output=Logs/09-kallisto-quant/slurm-%A_%a.out
#SBATCH --error=Logs/09-kallisto-quant/slurm-%A_%a.err

##############################################################################
# Run kallisto with 500 bootstraps for Sleuth
#
# --genomebam is needed to generate a genome-mapped BAM file for browsing with
# IGV.  It requires --gtf and --chromosomes. --chromosomes requires a TSV file
# with chromosome name and length on each line.  The chromosome names in the
# TSV must exactly match the names in the GTF.
# https://github.com/pachterlab/kallisto/issues/155
#
# The format and source of the chromosomes TSV is not clearly documented.
# I generated one using an Ensemble GFF with Reference/create-chrom-sizes.sh.
# GTF does not contain chromosome features.

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=16}

# Document software versions used for publication
uname -a
kallisto version
pwd

# gtf=$(Reference/gtf-filename.sh)

# If using hdf5, you may need this:
# https://github.com/pachterlab/kallisto/issues/197
# export HDF5_USE_FILE_LOCKING=FALSE

# 6-merge-bams.sbatch relies on sample N being in Results/09-kallisto-quant/N
# The sample number comes after -sample in the filename, e.g.
# chondro-sample4-rep2-time1-R1.fastq.xz is sample 4

# kallisto 0.46.1 can't handle xz and will simply seg fault rather than
# issue an error message.  If your trimmed fastq files are in xz format,
# this will convert to zstd format.
# Convert xz to zst rather than raw to reduce NFS load from compute nodes
zst1=$(echo Results/04-trim/sample${SLURM_ARRAY_TASK_ID}-*-R1.fastq.zst)
zst2=$(echo Results/04-trim/sample${SLURM_ARRAY_TASK_ID}-*-R2.fastq.zst)

printf "$zst1\n$zst2\n"

# Kallisto requires an output subdirectory for each sample
stem=$(basename ${zst1%-R1*})
out_dir=Results/09-kallisto-quant/$stem
echo $out_dir
mkdir -p $out_dir

# kallisto only support gzip compression as of 0.48.0, so use FIFOs
# feed it raw text
pipe1=/tmp/pipe1-$SLURM_ARRAY_TASK_ID
pipe2=/tmp/pipe2-$SLURM_ARRAY_TASK_ID
rm -f $pipe1 $pipe2
mkfifo $pipe1 $pipe2
zstdcat $zst1 > $pipe1 &
zstdcat $zst2 > $pipe2 &

# Manual says a GTF is needed.  Kallisto aborts using GFF3.
# Needed for --genomebam
# We'll get BAMs from hisat2 instead
# gtf=AmexT_v47-AmexG_v6.0-DD.gtf

set -x
kallisto quant \
    --threads=$SLURM_CPUS_PER_TASK \
    --index=Results/08-kallisto-index/AmexT_v47_cds.index \
    --output-dir=$out_dir $pipe1 $pipe2

rm -f $pipe1 $pipe2

# AxoXenONReg

RNA-Seq analysis of Axolotl and Xenopus

## PIs

Ava Udvadia, Fiona Watson

## Bioinformatics Analysis

Jason W. Bacon

Based on the work of Dr. Andrea Rau:

[https://github.com/andreamrau/OpticRegen_2019](https://github.com/andreamrau/OpticRegen_2019)

[https://www.nature.com/articles/s41598-019-50485-6](https://www.nature.com/articles/s41598-019-50485-6)

## Using this pipeline

The scripts in this pipeline are meant to serve as an example of a
basic RNA-Seq analysis, requiring only basic Unix skills to complete.
The goal is to allow budding bioinformaticians to achieve success in
an RNA-Seq analysis as early as possible in their training.  This will
hopefully develop confidence and allow them to see the big picture before
moving on to refine their computer skills and bioinformatics analysis skills.

To run the analysis, cd into Species/No-cluster if you are not using an HPC
cluster, or into Species/Cluster if you are using a SLURM HPC cluster.
Then simple run the scripts in order, e.g.

./01-organize.sh
./02-qc-raw.sh
./03-multiqc-raw.sh
./03b-verify.sh
./04-trim.sh
./04b-verify.sh

Under No-cluster, these scripts will run as many parallel jobs as possible
using the Unix xargs utility.

Under Cluster, these scripts will submit a SLURM job.

Note that there are example checkpoint scripts provided in *b-verify.sh.
These are meant mainly to remind you where you should stop and examine
the results of the stages so far.  They do not perform exhaustive
verification, but only provide a quick example summary of the results.
Closer examination of output files may be necessary to ensure that you
are ready to proceed to the next stage.

See Doc/pipeline.pdf for a detailed description of each stage.

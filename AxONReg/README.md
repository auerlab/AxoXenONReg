# Using these scripts

Simply run the scripts in this directory in order.

These scripts will automatically run the scripts in the SLURM
subdirectory if an "srun" command is present, or the scripts
in the Xargs directory if no "srun" command is found (which
presumably means we are not using a cluster).

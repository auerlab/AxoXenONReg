#!/bin/sh -e

##########################################################################
#   Function description:
#       Pause until user presses return
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

cmd="more Logs/04-trim/*.out"
printf "Running '$cmd'...\n"
pause
$cmd

cmd="more Logs/04-trim/*.err"
printf "\n\nRunning '$cmd'...\n"
pause
$cmd

cat << EOM

Read counts in the trimmed files should be nearly the same as the
raw files, minus reads removed for being too short (shown in the
fastq-trim output above).  This comparison can take a long time.

EOM
printf "Compare stats between raw and trimmed? [y]/n "
read compare
if [ 0"$compare" != 0n ]; then
    sbatch SLURM/04b-verify.sbatch
    cat << EOM

Output will be in Logs/04b-verify.  When the job is complete, you can
view it by running:

    more Logs/04b-verify/fastx-stats*.out
    more Logs/04b-verify/fastx-stats*.err

EOM
fi

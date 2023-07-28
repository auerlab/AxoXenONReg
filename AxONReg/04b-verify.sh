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

cmd="cat Logs/04-trim/*.out"
printf "Running '$cmd'...\n"
pause
$cmd | more

cmd="cat Logs/04-trim/*.err"
printf "\n\nRunning '$cmd'...\n"
pause
$cmd | more

cat << EOM

Read counts in the trimmed files should be nearly the same as the
raw files, minus reads removed for being too short (shown in the
fastq-trim output above).  This comparison can take a long time.

EOM
printf "Compare stats between raw and trimmed? [y]/n "
read compare
if [ 0"$compare" != 0n ]; then
    for raw in Results/01-organize/Raw-renamed/*; do
	printf "Scanning $raw...\n"
	blt fastx-stats $raw
	trimmed=$(basename $raw)
	trimmed=${trimmed%.xz}.zst
	printf "Scanning $trimmed...\n"
	blt fastx-stats Results/04-trim/$trimmed
	printf "===\n"
    done | more
fi

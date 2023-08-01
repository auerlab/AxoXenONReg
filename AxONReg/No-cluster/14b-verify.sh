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

cmd="head Results/14-fasda-fc-hisat2/*.tsv"
printf "Running '$cmd'...\n"
pause
$cmd | more

cmd="head Results/14-fasda-fc-hisat2/*.txt"
printf "Running '$cmd'...\n"
pause
$cmd | more

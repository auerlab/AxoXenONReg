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

cmd='more Logs/13-fasda-abundance-hisat2/*.out'
printf "Running $cmd...\n"
pause
$cmd

cmd='more Logs/13-fasda-abundance-hisat2/*.err'
printf "Running $cmd...\n"
pause
$cmd

cmd='head Results/13-fasda-abundance-hisat2/*-abundance.tsv'
printf "Running $cmd...\n"
pause
$cmd | more

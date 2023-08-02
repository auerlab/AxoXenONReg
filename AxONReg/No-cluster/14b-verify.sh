#!/bin/sh -e

##########################################################################
#   Example verification script.  Edit to add additional checks or
#   tune the existing checks to your personal taste.
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

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
printf "Running '$cmd'...\n"
pause
$cmd | more


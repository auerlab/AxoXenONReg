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

cmd="webbrowser Results/06-multiqc-trimmed/multiqc_report.html"
printf "Running '$cmd'...\n"
pause
$cmd

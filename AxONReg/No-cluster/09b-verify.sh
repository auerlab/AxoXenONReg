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

cmd='more Logs/09-kallisto-quant/*.out'
printf "Running $cmd...\n"
pause
$cmd

cmd='more Logs/09-kallisto-quant/*.err'
printf "Running $cmd...\n"
pause
$cmd

cmd='more Results/09-kallisto-quant/*/run_info.json'
printf "Running $cmd...\n"
pause
$cmd

cmd='head Results/09-kallisto-quant/*/abundance.tsv'
printf "Running $cmd...\n"
pause
$cmd | more

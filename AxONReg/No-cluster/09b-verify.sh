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

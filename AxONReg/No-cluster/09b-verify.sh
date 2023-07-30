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

more Logs/09-kallisto-quant/*.out
more Logs/09-kallisto-quant/*.err
head Results/09-kallisto-quant/*/abundance.tsv | more

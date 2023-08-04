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

cat Logs/12-hisat2-align/*.err | fgrep 'aligned concordantly exactly' 
pause

cat Logs/12-hisat2-align/*.err | fgrep 'aligned concordantly >'
pause

cat Logs/12-hisat2-align/*.err | fgrep 'overall alignment'
pause

cmd='ls -lh Results/12-hisat2-align/*.bam'
printf "Running $cmd...\n"
pause

$cmd | more
cmd='more Logs/12-hisat2-align/*.err'
printf "Running $cmd...\n"
pause
$cmd

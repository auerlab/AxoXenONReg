#!/bin/sh -e

##########################################################################
#   Description:
#       Report whether each stage has been submitted/run, so that
#       the user knows where they left off after being distracted.
##########################################################################

for script in [012][0-9]-*.s*; do
    dir=Logs/${script%.*}
    if [ -d $dir ]; then
	printf "%-40s: " $script
	if [ -z "$(ls $dir)" ]; then
	    printf "Not submitted.\n"
	else
	    printf "Submitted\n"
	fi
    fi
done

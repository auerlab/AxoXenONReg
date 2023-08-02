#!/bin/sh -e

##########################################################################
#   Example verification script.  Edit to add additional checks or
#   tune the existing checks to your personal taste.
##########################################################################

more Logs/10-fasda-kallisto/*.out
more Logs/10-fasda-kallisto/*.err
more Results/10-fasda-kallisto/*.tsv
more Results/10-fasda-kallisto/*.txt

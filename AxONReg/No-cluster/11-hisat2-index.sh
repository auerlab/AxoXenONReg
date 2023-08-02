#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw and trimmed data for comparison
#
#   Dependencies:
#       Requires directory structure.  Run after *-organize.sh.
#
#   History:
#   Date        Name        Modification
#   2023-06     Jason Bacon Begin
##########################################################################

Sh/11-hisat2-index.sh \
    2>&1 | tee Logs/11-hisat2-index/hisat2-index.out

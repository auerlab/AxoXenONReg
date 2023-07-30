#!/bin/sh -e

#   Dependencies:
#       Requires trimmed reads.

# This script is single threaded, so there is no difference between
# running on a cluster or a workstation.

Sh/07-reference.sh 2>&1 | tee Logs/07-reference/07-reference.out

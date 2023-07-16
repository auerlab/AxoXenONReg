#!/bin/sh -e

build=$(Reference/genome-build.sh)
release=$(Reference/genome-release.sh)

echo Xenopus_tropicalis.UCB_Xtro_$build.$release.chr.gtf


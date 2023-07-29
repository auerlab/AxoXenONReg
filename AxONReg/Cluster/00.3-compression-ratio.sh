#!/bin/sh -e

gz_total=0
xz_total=0
ls Raw/*/*.fq.gz    # Terminate script if gz files have been removed
for file in Raw/*/*.fq.gz; do
    gz_size=$(du -sm $file | cut -f 1)
    xz_size=$(du -sm ${file%.gz}.xz | cut -f 1)
    printf "$file:\t$gz_size\t$xz_size\t"
    gz_total=$(($gz_total + $gz_size))
    xz_total=$(($xz_total + $xz_size))
    printf "$xz_size / $gz_size\n" | bc -l
done
printf "Gzip total:\t$(($gz_total / 1000)) GB\n"
printf "Xz total:  \t$(($xz_total / 1000)) GB\n"
printf "Space saved:\t$((($gz_total - $xz_total) / 1000)) GB\n"

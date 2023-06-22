#!/bin/sh -e

if [ $# != 1 ]; then
    printf "Usage: ./00.2-gz2xz.sh file.gz\n"
    exit 1
fi

my_file_gz=$1
my_file_xz=${my_file_gz%.gz}.xz

if [ ! -e $my_file_xz ]; then
    printf "Converting $my_file_gz to $my_file_xz...\n"
    #time gunzip -c $my_file_gz | xz > $my_file_xz
fi

printf "Checking $my_file_gz...\n"
gunzip -c $my_file_gz | md5 > $my_file_gz-ck.txt
printf "Checking $my_file_xz...\n"
xzcat $my_file_xz | md5 > $my_file_xz-ck.txt

# diff returns success when files are not different
if diff $my_file_gz-ck.txt $my_file_xz-ck.txt; then
    cat << EOM2
    
Files are the same.  If should be safe to remove this copy of
$my_file_gz, assuming you have another copy of the file readily
accessible as a backup.

EOM2
else
    printf "Files are different.  $my_file_xz should not be trusted.\n"
    printf "Inspect it to see what the problem is and remove it.\n"
fi


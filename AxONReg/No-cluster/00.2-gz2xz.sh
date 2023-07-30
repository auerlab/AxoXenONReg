#!/bin/sh -e

cat << EOM > Sh/00.2-gz2xz.sh
#!/bin/sh -e

if [ \$# != 1 ]; then
    printf "Usage: $0 file.gz\n"
    exit 1
fi

my_gz=\$1
my_xz=\${my_gz%.gz}.xz

if [ ! -e \$my_xz ]; then
    printf "Converting \$my_gz to \$my_xz...\n"
    time gunzip -c \$my_gz | xz > \$my_xz
else
    printf "\$my_xz already exists.\n"
fi

printf "Checking \$my_gz...\n"
gunzip -c \$my_gz | md5 > \$my_gz-ck.txt
printf "Checking \$my_xz...\n"
xzcat \$my_xz | md5 > \$my_xz-ck.txt

# diff returns success when files are not different
if diff \$my_gz-ck.txt \$my_xz-ck.txt; then
    printf "Files are the same.  Removing \$my_gz.\n"
    rm \$my_gz
else
    printf "Files are different.  \$my_xz should not be trusted.\n"
    printf "Inspect it to see what the problem is and remove it.\n"
fi

EOM

chmod 755 Sh/00.2-gz2xz.sh
hw_threads=$(../Common/get-hw-threads.sh)
jobs=$hw_threads    # gunzip takes almost nothing

# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Raw/*/*.fq.gz | xargs -n 1 -P $jobs Sh/00.2-gz2xz.sh

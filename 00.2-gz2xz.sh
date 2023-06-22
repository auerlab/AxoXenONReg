#!/bin/sh -e

if which sbatch; then
    printf "Generating sbatch file...\n"

    files=$(ls Raw/*/*.fq.gz)
    count=$(echo $files | awk '{ print NF }')

    cat << EOM > SLURM/00.2-gz2xz.sbatch
#!/bin/sh -e

#SBATCH --array=1-$count
#SBATCH --output=Logs/00-gz2xz/slurm-%A_%a.out
#SBATCH --error=Logs/00-gz2xz/slurm-%A_%a.err

: \${SLURM_ARRAY_TASK_ID:=1}

files=\$(ls Raw/*/*.fq.gz)
my_file_gz=\$(echo \$files | awk -v index=\$SLURM_ARRAY_TASK_ID '{ print \$index }')
my_file_xz=\${my_file_gz%.gz}.xz

printf "Converting \$my_file_gz...\n"
if [ ! -e \$my_file_xz ]; then
    time gunzip -c \$my_file_gz | xz > \$my_file_xz
fi

printf "Checking \$my_file_gz...\n"
gunzip -c \$my_file_gz | md5 > \$my_file_gz-ck.txt
printf "Checking \$my_file_xz...\n"
xzcat \$my_file_xz | md5 > \$my_file_xz-ck.txt

# diff returns success when files are not different
if diff \$my_file_gz-ck.txt \$my_file_xz-ck.txt; then
    cat << EOM2
    
Files are the same.  If should be safe to remove this copy of
\$my_file_gz, assuming you have another copy of the file readily
accessible as a backup.

EOM2
else
    printf "Files are different.  \$my_file_xz should not be trusted.\n"
    printf "Inspect it to see what the problem is and remove it.\n"
fi
EOM

    chmod 755 SLURM/00.2-gz2xz.sbatch
    mkdir -p Logs/00.2-gz2xz
    # sbatch SLURM/00.2-gz2xz.sbatch

else

    printf "No SLURM.  Using xargs...\n"
    cat << EOM > Xargs/00.2-gz2xz.sh
#!/bin/sh -e

if [ \$# != 1 ]; then
    printf "Usage: $0 file.gz\n"
    exit 1
fi

my_file_gz=\$1
my_file_xz=\${my_file_gz%.gz}.xz

if [ ! -e \$my_file_xz ]; then
    printf "Converting \$my_file_gz to \$my_file_xz...\n"
    #time gunzip -c \$my_file_gz | xz > \$my_file_xz
fi

printf "Checking \$my_file_gz...\n"
gunzip -c \$my_file_gz | md5 > \$my_file_gz-ck.txt
printf "Checking \$my_file_xz...\n"
xzcat \$my_file_xz | md5 > \$my_file_xz-ck.txt

# diff returns success when files are not different
if diff \$my_file_gz-ck.txt \$my_file_xz-ck.txt; then
    cat << EOM2
    
Files are the same.  If should be safe to remove this copy of
\$my_file_gz, assuming you have another copy of the file readily
accessible as a backup.

EOM2
else
    printf "Files are different.  \$my_file_xz should not be trusted.\n"
    printf "Inspect it to see what the problem is and remove it.\n"
fi

EOM
    #ape Xargs/00.2-gz2xz.sh
    #exit
    chmod 755 Xargs/00.2-gz2xz.sh
    hw_threads=$(./get-hw-threads.sh)
    jobs=$hw_threads    # gunzip takes almost nothing
    
    # Tried GNU parallel and ran into bugs.  Xargs just works.
    ls Raw/*/*.fq.gz | xargs -n 1 -P $jobs Xargs/00.2-gz2xz.sh
fi

#!/bin/sh -e

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
my_gz=\$(echo \$files | awk -v index=\$SLURM_ARRAY_TASK_ID '{ print \$index }')
my_xz=\${my_gz%.gz}.xz

printf "Converting \$my_gz...\n"
if [ ! -e \$my_xz ]; then
    time gunzip -c \$my_gz | xz > \$my_xz
fi

printf "Checking \$my_gz...\n"
gunzip -c \$my_gz | md5 > \$my_gz-ck.txt
printf "Checking \$my_xz...\n"
xzcat \$my_xz | md5 > \$my_xz-ck.txt

# diff returns success when files are not different
if diff \$my_gz-ck.txt \$my_xz-ck.txt; then
    cat << EOM2
    
Files are the same.  If should be safe to remove this copy of
\$my_gz, assuming you have another copy of the file readily
accessible as a backup.

EOM2
else
    printf "Files are different.  \$my_xz should not be trusted.\n"
    printf "Inspect it to see what the problem is and remove it.\n"
fi
EOM

chmod 755 SLURM/00.2-gz2xz.sbatch
mkdir -p Logs/00.2-gz2xz
sbatch SLURM/00.2-gz2xz.sbatch

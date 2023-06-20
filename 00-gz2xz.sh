#!/bin/sh -e

files=$(ls Raw/*/*.fq.gz)
count=$(echo $files | awk '{ print NF }')

cat << EOM > gz2xz.sbatch
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
    gunzip -c \$my_file_gz | xz > \$my_file_xz
fi

printf "Checking \$my_file_gz...\n"
gunzip -c \$my_file_gz | wc > \$my_file_gz-wc.txt
printf "Checking \$my_file_xz...\n"
xzcat \$my_file_xz | wc > \$my_file_xz-wc.txt

# diff returns success when files are not different
if diff \$my_file_gz-wc.txt \$my_file_xz-wc.txt; then
    printf "Files are the same.\n"
    # rm \$my_file_gz
else
    printf "Files are different.\n"
    # rm \$my_file_xz
fi
EOM

chmod 755 gz2xz.sbatch
mkdir -p Logs/00-gz2xz
cat gz2xz.sbatch


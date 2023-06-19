#!/bin/sh -e

files=$(ls Raw/*/*.fq.gz)
count=$(echo $files | awk '{ print NF }')

cat << EOM > gz2xz.sbatch
#!/bin/sh -e

#SBATCH --array=1-$count
#SBATCH --output=Logs/00-gz2xz/slurm-%A_%a.out
#SBATCH --error=Logs/00-gz2xz/slurm-%A_%a.err

files=\$(ls Raw/*/*.fq.gz)
my_file=\$(echo \$files | awk -v index=\$SLURM_ARRAY_TASK_ID '{ print \$index }')
set -x
gzcat \$my_file | xz > \${my_file%.gz}.xz
EOM

mkdir -p Logs/00-gz2xz
cat gz2xz.sbatch


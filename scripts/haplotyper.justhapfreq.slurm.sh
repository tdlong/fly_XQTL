#!/bin/bash
#SBATCH --job-name=haplotyper
#SBATCH -A tdlong_lab        ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --cpus-per-task=1 

module load R/3.6.2
file=$1
folder=$2
SNPtable=$3
foundernames=$4

pool=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f 1` 
OutName=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f 2`

echo "Rscript scripts/haplotyper.justhapfreq.R $pool $OutName $folder $SNPtable $foundernames"
Rscript --verbose scripts/haplotyper.justhapfreq.R $pool $OutName $folder $SNPtable $foundernames


#!/bin/bash
#SBATCH --job-name=power     ## Name of the job.
#SBATCH -A ecoevo283         ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --array=1-1000        ## normal wc -l $file is helpful ... but only the first 17 chromosomes are not junkola 
#SBATCH --cpus-per-task=1    ## number of cores the job needs, can the programs you run make used of multiple cores?

module load R/3.6.2
file="helperfile/seeds.txt"
rseed=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1` 
echo "Rscript scripts/power2.R $SLUM_ARRAY_TASK_ID $rseed"
Rscript --verbose scripts/power2.R $SLURM_ARRAY_TASK_ID $rseed


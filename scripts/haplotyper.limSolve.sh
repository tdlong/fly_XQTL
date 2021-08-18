#!/bin/bash
#$ -N lhap
#$ -q adl
#$ -pe openmp 4 
#$ -R y

module load R/3.4.1
file=$1
folder=$2
SNPtable=$3
foundernames=$4

Cpool=`head -n $SGE_TASK_ID $file | tail -n 1 | cut -f 1` 
Tpool=`head -n $SGE_TASK_ID $file | tail -n 1 | cut -f 2` 
OutName=`head -n $SGE_TASK_ID $file | tail -n 1 | cut -f 3`

echo "Rscript scripts/haplotyper.limSolve.R $Cpool $Tpool $OutName $folder $SNPtable $foundernames"
Rscript --verbose scripts/haplotyper.limSolve.R $Cpool $Tpool $OutName $folder $SNPtable $foundernames


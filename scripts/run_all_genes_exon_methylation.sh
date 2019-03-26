#!/bin/bash --login
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --job-name all_genes_exon_methylation
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
#Set tmp directories
export TMPDIR=$PBS_O_WORKDIR
export TMP=$PBS_O_WORKDIR
export TEMP=$PBS_O_WORKDIR

#variables
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)

#get total weighted mC
echo "Get gene exon methylation data for $sample"
cd combined
python ../../../scripts/all_genes_exon_methylation.py $sample


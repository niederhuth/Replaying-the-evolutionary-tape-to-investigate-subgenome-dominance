#!/bin/bash --login
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --job-name total_mC
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/mC/bin:$PATH"

#Set tmp directories
export TMPDIR=$PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/mC/bin:$PATH"

export TMP=$PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/mC/bin:$PATH"

export TEMP=$PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/mC/bin:$PATH"


#variables
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)

#get total weighted mC
echo "Get total methylation of $sample"
cd combined
python ../../../scripts/total_mC.py $sample


#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120GB
#SBATCH --job-name 100kb_windows
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
echo "Get 100kb windows data for $sample"
cd combined
python ../../../scripts/100kb_windows.py $sample


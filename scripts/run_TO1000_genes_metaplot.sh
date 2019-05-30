#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --job-name TO1000_genes_metaplot
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
echo "Get gene metaplot data for $sample"
cd combined
python ../../../scripts/TO1000_genes_metaplot.py $sample


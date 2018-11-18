#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60GB
#SBATCH --job-name total_mC
#SBATCH --output=%x-%j.SLURMout

export TMPDIR=$PBS_O_WORKDIR
export TMP=$PBS_O_WORKDIR
export TEMP=$PBS_O_WORKDIR

cd $PBS_O_WORKDIR
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)

echo "MethylC Analysis of $sample"
python ../../scripts/total_mC.py $sample


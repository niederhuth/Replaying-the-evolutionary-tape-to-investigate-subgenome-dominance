#!/bin/bash -login
#PBS -l walltime=3:59:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=20gb
#PBS -N methylC_analysis

export TMPDIR=$PBS_O_WORKDIR
export TMP=$PBS_O_WORKDIR
export TEMP=$PBS_O_WORKDIR

cd $PBS_O_WORKDIR
module load BEDTools/2.24.0
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)
mkdir results

echo "MethylC Analysis of $sample"
python ../../scripts/methylC_analysis.py $sample

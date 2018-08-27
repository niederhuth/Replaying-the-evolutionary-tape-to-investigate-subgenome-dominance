#!/bin/bash -login
#PBS -l walltime=3:59:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=50gb
#PBS -N chromosome_map

export TMPDIR=$PBS_O_WORKDIR
export TMP=$PBS_O_WORKDIR
export TEMP=$PBS_O_WORKDIR

cd $PBS_O_WORKDIR
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)

echo "MethylC Analysis of $sample"
python ../../scripts/gene_methylation.py $sample

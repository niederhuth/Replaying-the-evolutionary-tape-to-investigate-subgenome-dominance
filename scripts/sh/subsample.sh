#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120GB
#SBATCH --job-name subsample
#SBATCH --output=%x-%j.SLURMout

cd $PBS_O_WORKDIR

seqtk sample -2 -s100 BS-To_H357WCCXY_L5_1.fastq 194909329 > To_1.fastq
seqtk sample -2 -s100 BS-To_H357WCCXY_L5_2.fastq 194909329 > To_2.fastq


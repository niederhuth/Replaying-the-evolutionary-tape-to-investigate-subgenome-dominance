#!/bin/bash --login
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25GB
#SBATCH --job-name setup
#SBATCH --output=%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/mC/bin:$PATH"

i=$(pwd | sed s/^.*\\///)

#Prep genome
echo "Setting up $i"
samtools faidx $i.fa
cut -f1,2 $i.fa.fai > $i.genome
methylpy build-reference --input-files $i.fa \
	--output-prefix $i --aligner bowtie2

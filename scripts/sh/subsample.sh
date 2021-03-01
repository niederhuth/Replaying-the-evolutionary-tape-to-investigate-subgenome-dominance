#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120GB
#SBATCH --job-name subsample
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/Bnapus-polyploidy/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/Bnapus-polyploidy/lib:$LD_LIBRARY_PATH"

seqtk sample -2 -s100 ../TO1000/fastq/_1.fastq 194909329 > fastq/To_1.fastq
seqtk sample -2 -s100 ../TO1000/fastq/_2.fastq 194909329 > fastq/To_2.fastq


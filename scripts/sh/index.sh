#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=25GB
#SBATCH --job-name setup
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/Bnapus-polyploidy/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/Bnapus-polyploidy/lib:$LD_LIBRARY_PATH"


samtools faidx combined.fa
cut -f1,2 combined.fa.fai > combined.genome

methylpy build-reference \
	--input-files combined.fa \
	--output-prefix combined \
	--num-procs 4 \
	--aligner bowtie2

echo "Done"
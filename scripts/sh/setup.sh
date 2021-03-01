#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25GB
#SBATCH --job-name setup
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/Bnapus-polyploidy/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/Bnapus-polyploidy/lib:$LD_LIBRARY_PATH"

sample=$(pwd | sed s/^.*\\///)

mkdir combined combined/results fastq

#Download SRA
cd fastq
echo "Downloading SRA files"
sra_list=$(grep ${sample} ../../../misc/samples.csv | cut -d ',' -f 5)
for i in $sra_list
do
	echo "Downloading $i"
	prefetch --max-size 100000000 $i
	mv $i/*sra ./
	rmdir $i
	echo "Converting $i"
	fastq-dump --gzip --split-3 "$i".sra
	rm "$i".sra
done

echo "Done"
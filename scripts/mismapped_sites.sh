#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60GB
#SBATCH --job-name find_mismapped_sites
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/mC/bin:$PATH"


#Get IMB218 sites mapped to Bo genome
cd IMB218/combined
zcat allc_IMB218.tsv.gz | awk '/^Bo/' > mismapped
#Adjust for data size
$size=$(zcat allc_IMB218.tsv.gz | awk '{sum+=$6} END {print sum}')
awk -v OFS="\t" -v x=$size 'print $1,$2,$3,$4,$5/x,$6/x,0' mismapped > size_adjusted


#Get TO1000 sites mapped to Br genome
cd ../../TO1000/combined
zcat allc_TO1000.tsv.gz | awk '/^Br/' > mismapped
#Adjust for data size
$size=$(zcat allc_TO1000.tsv.gz | awk '{sum+=$6} END {print sum}')
awk -v OFS="\t" -v x=$size 'print $1,$2,$3,$4,$5/x,$6/x,0' mismapped > size_adjusted

#Combined IMB218 & TO1000
cd ../../
header="##fileformat=VCFv4.1"
echo $header | cat - IMB218/combined/size_adjusted TO1000/combined/size_adjusted > misaligned_sites


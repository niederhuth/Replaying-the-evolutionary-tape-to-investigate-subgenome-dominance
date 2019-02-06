#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60GB
#SBATCH --job-name correct_sites
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR

#List Variables
sample=$(pwd | sed s/^.*\\///)
header="##fileformat=VCFv4.1"

#unzip allc file
cd combined2
zcat allc_"$sample".tsv.gz > tmp 
echo $header | cat - tmp > tmp2

#count mapped sites
z=$(awk '{sum+=$6} END {print sum}' tmp)

#adjust by number of sites and round up
awk -v z=$z '{printf "%s\t%s\t%s\t%s\t%.0f\t%.0f\t%s\n", $1,$2,$3,$4,($5*z)+0.5,($6*z)+0.5,$7}' ../../misaligned_sites > tmp3

echo $header | cat - tmp3 > misaligned_sites

#bedtools intersect -v -a tmp2 -b misaligned_sites > no_match
bedtools intersect -wa -wb -a tmp2 -b misaligned_sites > match
awk -v OFS="\t" '{print $1,$2,$3,$4,$5-$12,$6-$13,$7}' match | awk -v OFS="\t" '{if ($6 > 0) print $0}' | awk -v OFS="\t" '{if ($5 > $6) print $1,$2,$3,$4,$6,$6,$7; else print $0}' | awk -v OFS="\t" '{if ($5 < 0) print $1,$2,$3,$4,0,$6,$7; else print $0}' > corrected_tmp
cat no_match corrected_tmp | grep -v \# | sort -k1,2h > corrected_allc.tsv

#remove tmp files
rm tmp tmp2 tmp3 misaligned_sites corrected_tmp #no_match match



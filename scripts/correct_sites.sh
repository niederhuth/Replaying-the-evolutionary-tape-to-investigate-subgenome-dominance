#!/bin/bash --login
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --job-name correct_sites
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR

#List Variables
sample=$(pwd | sed s/^.*\\///)

#unzip allc file
cd combined2
zcat allc_"$sample".tsv.gz > tmp 

#count mapped sites
z=$(awk '{sum+=$6} END {print sum}' tmp)

#adjust by number of sites and round up
awk -v z=$z '{printf "%s\t %s\t %s\t %s\t %.0f\t %.0f\t %s\n", $1,$2,$3,$4,($5*z)+0.5,($6*z)+0.5,$7}' ../../misaligned_sites > misaligned_sites

#correct sites
while read line
do
a=$(echo $line | tr ' ' '\t' | cut -f1)
b=$(echo $line | tr ' ' '\t' | cut -f2)
c=$(echo $line | tr ' ' '\t' | cut -f3)
d=$(echo $line | tr ' ' '\t' | cut -f4)
e=$(awk -v a=$a -v b=$b -v c=$c -v d=$d '{if ($1==a && $2==b && $3==c && $4==d) print $0}' misaligned_sites)
f=$(echo $e | tr ' ' '\t' | cut -f5)
g=$(echo $e | tr ' ' '\t' | cut -f6)
echo $line | awk -v f=$f -v g=$g -v OFS="\t" '{print $1,$2,$3,$4,$5-f,$6-g,$7}' >> corrected_sites
done < tmp

#remove tmp files
rm tmp misaligned_sites


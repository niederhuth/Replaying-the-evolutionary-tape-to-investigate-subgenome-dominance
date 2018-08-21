#!/bin/bash -login
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=6gb
#PBS -N setup

cd $PBS_O_WORKDIR
i=$(pwd | sed s/^.*\\///)
module load bowtie2/2.3.1
module load SAMTools/1.5

#Prep genome
echo "Setting up $i"
samtools faidx $i.fa
cut -f1,2 $i.fa.fai > $i.genome
methylpy build-reference --input-files $i.fa \
	--output-prefix $i --bowtie2 True

#!/bin/bash --login
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name methylpy
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR

sample=$(pwd | sed s/^.*\\///)

if ls fastq/*_2.fastq >/dev/null 2>&1
then
	methylpy paired-end-pipeline \
	--read1-files fastq/*_1.fastq \
	--read2-files fastq/*_2.fastq \
	--sample $sample \
	--forward-ref ../ref/combined/combined_f.mmi \
	--reverse-ref ../ref/combined/combined_r.mmi \
	--ref-fasta ../ref/combined/combined.fa \
	--libraries "libA" \
	--path-to-output "" \
	--pbat False \
	--check-dependency False \
	--num-procs 20 \
	--sort-mem 5G \
	--num-upstream-bases 0 \
	--num-downstream-bases 2 \
	--generate-allc-file True \
	--generate-mpileup-file True \
	--compress-output True \
	--bgzip False \
	--path-to-bgzip "" \
	--path-to-tabix "" \
	--trim-reads True \
	--path-to-cutadapt "" \
	--path-to-aligner "" \
	--aligner "minimap2" \
	--aligner-options "-ax sr --secondary=no --heap-sort no" \
	--merge-by-max-mapq True \
	--remove-clonal True \
	--path-to-picard /mnt/home/niederhu/anaconda3/share/picard-2.18.16-0 \
	--keep-clonal-stats True \
	--java-options "" \
	--path-to-samtools "" \
	--adapter-seq-read1 AGATCGGAAGAGCACACGTCTGAAC \
	--adapter-seq-read2 AGATCGGAAGAGCGTCGTGTAGGGA \
	--remove-chr-prefix False \
	--add-snp-info False \
	--unmethylated-control "37_Plastid" \
	--binom-test True \
	--sig-cutoff .01 \
	--min-mapq 30 \
	--min-cov 3 \
	--max-adapter-removal 1 \
	--overlap-length 3 \
	--error-rate 0.1 \
	--min-qual-score 10 \
	--min-read-len 30 \
	--min-base-quality 1 \
	--keep-temp-files False
fi

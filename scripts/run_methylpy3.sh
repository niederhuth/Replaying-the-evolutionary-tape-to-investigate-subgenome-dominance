#!/bin/bash --login
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB
#SBATCH --job-name R500_methylpy
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
mkdir R500

#List Variables
sample=$(pwd | sed s/^.*\\///)
f_ref="../../ref/R500/R500_f"
r_ref="../../ref/R500/R500_r"
fasta="../../ref/R500/R500.fa"
read1="../fastq/*_1.fastq"
read2="../fastq/*_2.fastq"
unmethylated_control="37_Plastid"
adaptor1="AGATCGGAAGAGCACACGTCTGAAC"
adaptor2="AGATCGGAAGAGCGTCGTGTAGGGA"
aligner="bowtie2"
aligner_options="--very-sensitive -X 1000"
picard="/mnt/home/niederhu/anaconda3/share/picard-2.18.16-0"

#Unzip fastq files
echo "Decompressing fastq files"
cd fastq
for i in *fastq.gz
do
	gunzip $i
done 
cd ../R500

#Run Methylpy
echo "Running methylpy"
methylpy paired-end-pipeline \
	--read1-files $read1 \
	--read2-files $read2 \
	--sample $sample \
	--forward-ref $f_ref \
	--reverse-ref $r_ref \
	--ref-fasta $fasta \
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
	--aligner $aligner \
	--aligner-options "$aligner_options" \
	--merge-by-max-mapq True \
	--remove-clonal True \
	--path-to-picard $picard \
	--keep-clonal-stats True \
	--java-options "" \
	--path-to-samtools "" \
	--adapter-seq-read1 $adaptor1 \
	--adapter-seq-read2 $adaptor2 \
	--remove-chr-prefix False \
	--add-snp-info False \
	--unmethylated-control $unmethylated_control \
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
	--keep-temp-files True

#rm *mpileup_output.tsv *_reads_no_clonal*.bam* *_libA.metric

#Compress fastq files
#echo "Compressing fastqs"
#cd ../fastq
#for i in *fastq
#do
#	gzip $i
#done

echo "Done"

#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --job-name methylpy

cd $PBS_O_WORKDIR

sample=$(pwd | sed s/^.*\\///)

cd fastq
for i in *fq.gz
do
        gunzip $i
        name=$(echo $i | sed s/\.gz//)
        new=$(echo $name | sed s/fq$/fastq/)
        mv $name $new
done

for i in *fastq.gz
do
	gunzip $i
done 
cd ../

if ls fastq/*_2.fastq >/dev/null 2>&1
then
	echo "Data is paired-end"
	echo "Running methylpy"
	methylpy paired-end-pipeline \
	--read1-files fastq/*_1.fastq \
	--read2-files fastq/*_2.fastq \
	--sample $sample \
	--forward-ref ../ref/combined/combined_f \
	--reverse-ref ../ref/combined/combined_r \
	--ref-fasta ../ref/combined/combined.fa \
	--libraries "libA" \
	--path-to-output "" \
	--pbat False \
	--check-dependency False \
	--num-procs 10 \
	--sort-mem 5 \
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
	--aligner "bowtie2" \
	--aligner-options "" \
	--merge-by-max-mapq False \
	--remove-clonal True \
	--path-to-picard /mnt/home/niederhu/anaconda3/share/picard-2.18.11-0/ \
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
else
	echo "Data is single-end"
	echo "Running methylpy"
	methylpy single-end-pipeline \
	--read-files fastq/*fastq \
	--sample $sample \
        --forward-ref ../ref/combined/combined_f \
        --reverse-ref ../ref/combined/combined_r \
        --ref-fasta ../ref/combined/combined.fa \
        --libraries "libA" \
        --path-to-output "" \
        --pbat False \
        --check-dependency False \
        --num-procs 10 \
        --sort-mem 5 \
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
        --aligner "bowtie2" \
        --aligner-options "" \
        --merge-by-max-mapq True \
        --remove-clonal True \
        --path-to-picard /mnt/home/niederhu/anaconda3/share/picard-2.18.11-0/ \
        --keep-clonal-stats True \
        --java-options "" \
        --path-to-samtools "" \
	--adapter-seq AGATCGGAAGAGCACACGTCTG \
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

#rm *mpileup_output.tsv *_reads_no_clonal*.bam* *_libA.metric

echo "Compressing fastqs"
cd fastq
for i in *fastq
do
	gzip $i
done

echo "Done"

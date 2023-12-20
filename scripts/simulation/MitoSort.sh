#!/bin/bash

##PBS configure
#PBS -N MitoSort
#PBS -j oe
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=8

## Set variables
runtime_output_dir=~/MitoSort/simulation/MitoSort/runtime/
input_dir=~/MitoSort/simulation/bam_data/simulate_depth/
output_dir=~/MitoSort/simulation/MitoSort/
nreads=(10000 15000 20000 25000 30000 35000 40000 45000 50000)

## Activate conda environment
source /public/home/shipy3/miniconda3/bin/activate
conda activate MitoSort

for i in "${nreads[@]}"
do
	# make output directory
	mkdir -p ${output_dir}${i}
	mkdir -p ${output_dir}${i}/MitoSort
	mkdir -p ${output_dir}${i}/MitoSort/BAM

	# Subset reads mapped to chrM
	samtools view -b ~/MitoSort/simulation/souporcell/${i}/100_${i}_1_6_simulated.sorted.bam chrM -o ${output_dir}${i}/MitoSort/BAM/possorted_chrM_realign.bam

	# index BAM file
	samtools index ${output_dir}${i}/MitoSort/BAM/possorted_chrM_realign.bam

	# process barcode file
	cut -f 1 ${input_dir}100_${i}_1_6_simulated.bam_barcodes.txt | cut -d "_" -f 2 > ${output_dir}${i}/MitoSort/100_${i}_1_6_simulated.bam_barcodes.txt

	# generate-snp-matrix
	START=$(date +%s.%N)
	python ~/MitoSort/bin/MitoSort_pipeline.py generate-snp-matrix -b ${output_dir}${i}/MitoSort/BAM/possorted_chrM_realign.bam -f ~/MitoSort/simulation/bam_data/simulate_depth/MitoSort_test/tools_test/50_15000_8/reheader_genome.fa -c ${output_dir}${i}/MitoSort/100_${i}_1_6_simulated.bam_barcodes.txt -m ~/MitoSort/data/hg19_chrM.bed --varscan_path /public/home/chenbzh5/Tools/VarScan.v2.3.7.jar -o ${output_dir}${i}
	END=$(date +%s.%N)
	Duration=$(echo "$END - $START" | bc)
	formatted_duration=$(date -u -d @$Duration +%H:%M:%S)
	echo "generate-snp-matrix : $formatted_duration" >> ${runtime_output_dir}${i}_runtime.txt

	# demultiplex
	START=$(date +%s.%N)
	python ~/MitoSort/bin/MitoSort_pipeline.py demultiplex -o ${output_dir}${i} -k 6
	END=$(date +%s.%N)
	Duration=$(echo "$END - $START" | bc)
	formatted_duration=$(date -u -d @$Duration +%H:%M:%S)
	echo "demultiplex : $formatted_duration" >> ${runtime_output_dir}${i}_runtime.txt
done













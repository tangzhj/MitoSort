#!/bin/bash

##PBS configure
#PBS -N Vireo
#PBS -j oe
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=8

## Set variables
runtime_output_dir=~/MitoSort/simulation/Vireo/runtime/
input_dir=~/MitoSort/simulation/bam_data/simulate_depth/
output_dir=~/MitoSort/simulation/Vireo/
nreads=(10000 15000 20000 25000 30000 35000 40000 45000 50000)

## Activate conda environment
source /public/home/shipy3/miniconda3/bin/activate
conda activate python3


for i in "${nreads[@]}"
do
	# make output directory
	mkdir -p ${output_dir}${i}

	# remove duplicated barcodes
	sort ~/MitoSort/simulation/souporcell/${i}/100_${i}_1_6_simulated.bam_barcodes.txt | uniq > ${output_dir}${i}/100_${i}_1_6_simulated.bam_barcodes.txt

	# genotyping
	START=$(date +%s.%N)
	cellsnp-lite -s /md01/shipy3/tzj/SMILE/simulation/souporcell/${i}/100_${i}_1_6_simulated.sorted.bam -b ${output_dir}${i}/100_${i}_1_6_simulated.bam_barcodes.txt -O ${output_dir}${i} -p 8 --minMAF 0.1 --minCOUNT 20 --gzip --UMItag None
	END=$(date +%s.%N)
	Duration=$(echo "$END - $START" | bc)
	formatted_duration=$(date -u -d @$Duration +%H:%M:%S)
	echo "genotyping : $formatted_duration" >> ${runtime_output_dir}${i}_runtime.txt

	# Demultiplexing for donors
	START=$(date +%s.%N)
	vireo -c ${output_dir}${i} -N 6 -o ${output_dir}${i}
	END=$(date +%s.%N)
	Duration=$(echo "$END - $START" | bc)
	formatted_duration=$(date -u -d @$Duration +%H:%M:%S)
	echo "demultiplexing : $formatted_duration" >> ${runtime_output_dir}${i}_runtime.txt
done










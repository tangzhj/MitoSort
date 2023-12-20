#!/bin/bash

##PBS configure
#PBS -N Freemuxlet
#PBS -j oe
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=8

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
set -x  ## Log everything

## Set variables
runtime_output_dir=~/MitoSort/simulation/Freemuxlet/runtime/
input_dir=~/MitoSort/simulation/bam_data/simulate_depth/
output_dir=~/MitoSort/simulation/Freemuxlet/
nreads=(10000 15000 20000 25000 30000 35000 40000 45000 50000)


for i in "${nreads[@]}"
do
	# make output directory
	mkdir -p ${output_dir}${i}

	# Popscle Pileup: identify the number of reads from each allele at each of the common SNP location
	START=$(date +%s.%N)
	echo "`date +%Y/%m/%d_%H:%M:%S` Start Popscle Pileup"
	singularity exec /md01/shipy3/software/Demuxafy.sif popscle dsc-pileup --sam ~/MitoSort/simulation/souporcell/${i}/100_${i}_1_6_simulated.sorted.bam --tag-group CB --vcf /md01/shipy3/ori/DB/Demuxafy/reheader_keep_main_chr_GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf --group-list ~/MitoSort/simulation/souporcell/${i}/100_${i}_1_6_simulated.bam_barcodes.txt --out ${output_dir}${i}
	echo "`date +%Y/%m/%d_%H:%M:%S` Popscle Pileup completed"
	END=$(date +%s.%N)
	Duration=$(echo "$END - $START" | bc)
	formatted_duration=$(date -u -d @$Duration +%H:%M:%S)
	echo "Popscle Pileup: $formatted_duration" >> ${runtime_output_dir}${i}_runtime.txt

	# demultiplex
	START=$(date +%s.%N)
	echo "`date +%Y/%m/%d_%H:%M:%S` Start demultiplexing"
	singularity exec /md01/shipy3/software/Demuxafy.sif popscle freemuxlet --plp ${output_dir}${i} --out ${output_dir}${i} --group-list  ~/MitoSort/simulation/souporcell/${i}/100_${i}_1_6_simulated.bam_barcodes.txt --nsample 6
	echo "`date +%Y/%m/%d_%H:%M:%S` Demultiplexing completed"
	END=$(date +%s.%N)
	Duration=$(echo "$END - $START" | bc)
	formatted_duration=$(date -u -d @$Duration +%H:%M:%S)
	echo "demultiplexing: $formatted_duration" >> ${runtime_output_dir}${i}_runtime.txt
done
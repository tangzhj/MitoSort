#!/bin/bash

##PBS configure
#PBS -N soupercell
#PBS -j oe
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=8

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
set -x  ## Log everything

## Set variables
runtime_output_dir=~/MitoSort/simulation/souporcell/runtime/
input_dir=/~/MitoSort/simulation/bam_data/simulate_depth/
output_dir=~/MitoSort/simulation/souporcell/
nreads=(10000 15000 20000 25000 30000 35000 40000 45000 50000)

for i in "${nreads[@]}"
do
	# make output directory
	mkdir -p ${output_dir}${i}

	# sort BAM files
	samtools sort -o ${output_dir}${i}/100_${i}_1_6_simulated.sorted.bam ${input_dir}100_${i}_1_6_simulated.bam

	# index BAM file
	samtools index ${output_dir}${i}/100_${i}_1_6_simulated.sorted.bam

	# process barcode file
	cut -f 2 ${input_dir}100_${i}_1_6_simulated.bam_barcodes.txt | cut -d "_" -f 2 | sed "1d" > ${output_dir}${i}/100_${i}_1_6_simulated.bam_barcodes.txt

	# run soupercell
	START=$(date +%s.%N)
	echo "`date +%Y/%m/%d_%H:%M:%S` Run start"
	singularity exec --bind ${input_dir},${output_dir}${i},/md01/chenbzh5/bottleneck2/simulation/bam_data/simulate_depth/MitoSort_test/tools_test/50_15000_8/ ~/software/souporcell_latest.sif souporcell_pipeline.py -i ${output_dir}${i}/100_${i}_1_6_simulated.sorted.bam -b ${output_dir}${i}/100_${i}_1_6_simulated.bam_barcodes.txt -f ~/MitoSort/simulation/bam_data/simulate_depth/MitoSort_test/tools_test/50_15000_8/reheader_genome.fa -t 8 -o ${output_dir}${i} -k 6 --no_umi True --min_alt 2 --min_ref 2
	echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
	END=$(date +%s.%N)
	Duration=$(echo "$END - $START" | bc)
	formatted_duration=$(date -u -d @$Duration +%H:%M:%S)
	echo "${i}\t$formatted_duration" >> ${runtime_output_dir}soupercell_runtime.txt
done
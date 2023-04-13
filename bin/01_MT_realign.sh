#!/bin/bash


## Parse input parameters
while getopts "b:g:t:o:" opt; do
  case $opt in
    b) bam_file="$OPTARG"
    ;;
    g) genome_fasta="$OPTARG"
    ;;
    t) gatk_path="$OPTARG"
    ;;
    o) output_dir="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done


## Make output directory
clean_output_dir=`dirname $output_dir`
BAM_output_dir=$clean_output_dir/MitoSort/BAM
mkdir -p $BAM_output_dir

cd $BAM_output_dir

## Sort BAM file by position if it is not sorted
if `samtools view -H $bam_file | grep -q "SO:coordinate"`;then
  echo "[Log] The input BAM file has been sorted by position"
else
  echo "[`date +"%Y-%m-%d %H:%M:%S"`] Sort BAM file by position"
  samtools sort -o possorted.bam $bam_file
  ## Index sorted BAM file
  echo "[`date +"%Y-%m-%d %H:%M:%S"`] Index sorted BAM file"
  samtools index possorted.bam
fi

## Subset reads mapped to chrM
if `samtools view -H $bam_file | grep -q "SO:coordinate"`;then
  echo "[`date +"%Y-%m-%d %H:%M:%S"`] Subset reads mapped to chrM"
  samtools view -b $bam_file chrM > possorted_chrM.bam
else
  echo "[`date +"%Y-%m-%d %H:%M:%S"`] Subset reads mapped to chrM"
  samtools view -b possorted.bam chrM > possorted_chrM.bam
fi

## Index BAM file that containing reads mapped to chrM
echo "[`date +"%Y-%m-%d %H:%M:%S"`] Index BAM file that containing reads mapped to chrM"
samtools index possorted_chrM.bam

## Create target intervals list using RealignerTargetCreator
echo "[`date +"%Y-%m-%d %H:%M:%S"`] Create target intervals list using GATK RealignerTargetCreator"
java -Xmx8g -jar $gatk_path \
        -R $genome_fasta \
        -T RealignerTargetCreator  \
        -nt 10 \
        -I possorted_chrM.bam \
        -o chrM_realignertargetcreator.intervals


## Realign reads using IndelRealigner
echo "[`date +"%Y-%m-%d %H:%M:%S"`] Realign reads using GATK IndelRealigner"
java -Xmx100g -jar $gatk_path \
        -R $genome_fasta \
        -T IndelRealigner \
        -filterNoBases \
        -maxReads 10000000 \
        -I possorted_chrM.bam \
        -targetIntervals chrM_realignertargetcreator.intervals \
        -o possorted_chrM_realign.bam


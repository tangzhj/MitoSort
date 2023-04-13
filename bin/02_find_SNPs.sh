#!/bin/bash

## Parse input parameters
while getopts "g:m:v:o:" opt; do
  case $opt in
    g) genome_fasta="$OPTARG"
    ;;
    m) chrM_length="$OPTARG"
    ;;
    v) VarScan_path="$OPTARG"
    ;;
    o) output_dir="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

## Set output directory
clean_output_dir=`dirname $output_dir`
BAM_output_dir=$clean_output_dir/MitoSort/BAM


## Generate text pileup output for chrM BAM file, with each line representing the pileup of reads at a single genomic position
echo "[`date +"%Y-%m-%d %H:%M:%S"`] Pileup reads at chrM region"
samtools mpileup -l $chrM_length -q 30 -Q 30 -f $genome_fasta -x $BAM_output_dir/possorted_chrM_realign.bam > $BAM_output_dir/possorted_chrM_realign.mpileup


## Call SNPs by VarScan pileup2snp 
echo "[`date +"%Y-%m-%d %H:%M:%S"`] Call SNPs by VarScan pileup2snp"
java -Xmx4g -jar $VarScan_path pileup2snp \
  $BAM_output_dir/possorted_chrM_realign.mpileup \
  --min-var-freq 0.01 \
  --min-reads2 2 >$BAM_output_dir/possorted_chrM_realign.snv
#!/bin/python

import os
import subprocess 
import datetime
import argparse

## Set paramter
parser = argparse.ArgumentParser()
parser.add_argument('--genome_fasta','-f', type=str,required=True, help='FASTA file of reference genome')
parser.add_argument('--chrm_length','-m', type=str,required=True, help='a BED file containing chrM region')
parser.add_argument('--varscan_path',type=str,required=True, help='Path to VarScan')
parser.add_argument('--output_dir','-o',type=str,required=True, help='Output parent directory')
args = parser.parse_args()

def record_time():
    today_date = datetime.date.today()
    shown_date = today_date.strftime("%Y/%m/%d/")
    now = datetime.datetime.now()
    shown_time = now.strftime("%H:%M:%S")
    value = shown_date+" "+shown_time
    return value


BAM_output_dir = os.path.join(args.output_dir, "MitoSort", "BAM")
os.chdir(BAM_output_dir)

# subset part of reads for faster computation
print("["+record_time()+"] Subset reads for faster computation")
subset_cmd = ["samtools", "view", "-s", "0.2", "-b",os.path.join(BAM_output_dir, "possorted_chrM_realign.bam"),">",os.path.join(BAM_output_dir, "possorted_chrM_realign_0.2.bam")]
subprocess.call(" ".join(subset_cmd), shell=True)


# Generate text pileup output for chrM BAM file, with each line representing the pileup of reads at a single genomic position
print("["+record_time()+"] Pileup reads at chrM region")
mpileup_cmd = ["samtools", "mpileup", "-l", args.chrm_length, "-q", "30", "-Q", "30", "-f", args.genome_fasta, "-x", os.path.join(BAM_output_dir, "possorted_chrM_realign_0.2.bam"), "-o", os.path.join(BAM_output_dir, "possorted_chrM_realign_0.2.mpileup")]
subprocess.call(" ".join(mpileup_cmd), shell=True)

# Call SNPs by VarScan pileup2snp 
print("["+record_time()+"] Call SNPs by VarScan pileup2snp")
pileup2snp_cmd = ["java", "-Xmx4g", "-jar", args.varscan_path, "pileup2snp", os.path.join(BAM_output_dir, "possorted_chrM_realign_0.2.mpileup"), "--min-var-freq", "0.01", "--min-reads2", "2", ">", os.path.join(BAM_output_dir, "possorted_chrM_realign.snv")]
subprocess.call(" ".join(pileup2snp_cmd), shell=True)

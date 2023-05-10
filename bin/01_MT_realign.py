#!/bin/python

import subprocess 
import os
import datetime
import argparse

## Set paramter
parser = argparse.ArgumentParser()
parser.add_argument('--bam_file','-b', type=str,required=True, help='BAM file output by Cellranger')
parser.add_argument('--genome_fasta','-f', type=str,required=True, help='FASTA file of reference genome')
parser.add_argument('--gatk_path',type=str,required=True, help='Path to GenomeAnalysisTK(GATK)')
parser.add_argument('--output_dir','-o', type=str,required=True, help='Output parent directory')
parser.add_argument('--data_type',type=str,required=False, default="ATAC",help='Type of data. Default is ATAC. Set to RNA if your data is scRNA')
parser.add_argument('--picard_path',type=str,required=False,help='picard.jar is required when the type of data is RNA')
args = parser.parse_args()

def record_time():
    today_date = datetime.date.today()
    shown_date = today_date.strftime("%Y/%m/%d/")
    now = datetime.datetime.now()
    shown_time = now.strftime("%H:%M:%S")
    value = shown_date+" "+shown_time
    return value

# make parent directory
MitoSort_output_dir = os.path.join(args.output_dir, "MitoSort")
if not os.path.exists(MitoSort_output_dir):
    os.mkdir(MitoSort_output_dir)


# make output directory of BAM files
BAM_output_dir = os.path.join(args.output_dir, "MitoSort", "BAM")
if not os.path.exists(BAM_output_dir):
    os.mkdir(BAM_output_dir)
os.chdir(BAM_output_dir)

# Sort BAM file by position if it is not sorted
if "SO:coordinate" in subprocess.check_output(["samtools", "view", "-H", args.bam_file]).decode():
    print("[Log] The input BAM file has been sorted by position")
else:
    print("["+record_time()+"] Sort BAM file by position")
    subprocess.call(["samtools", "sort", "-o", "possorted.bam", args.bam_file], check=True)
    print("["+record_time()+"] Index sorted BAM file")
    subprocess.call(["samtools", "index", "possorted.bam"])

# Subset reads mapped to chrM
if "SO:coordinate" in subprocess.check_output(["samtools", "view", "-H", args.bam_file]).decode():
    print("["+record_time()+"] Subset reads mapped to chrM")
    #subprocess.call(["samtools", "view", "-b", args.bam_file, "chrM", "-o", "possorted_chrM.bam"])
else:
    print("["+record_time()+"] Subset reads mapped to chrM")
    subprocess.call(["samtools", "view", "-b", "possorted.bam", "chrM", "-o", "possorted_chrM.bam"])

# Index BAM file that containing reads mapped to chrM
print("["+record_time()+"] Index BAM file that containing reads mapped to chrM")
#subprocess.call(["samtools", "index", "possorted_chrM.bam"])

if args.data_type == "ATAC":
    print("["+record_time()+"] Create target intervals list using GATK RealignerTargetCreator")
    subprocess.call(["java", "-Xmx8g", "-jar", args.gatk_path, "-R", args.genome_fasta, "-T", "RealignerTargetCreator", "-nt", "10", "-I", "possorted_chrM.bam", "-o", "chrM_realignertargetcreator.intervals"])

    # Realign reads using IndelRealigner
    print("["+record_time()+"] Realign reads using GATK IndelRealigner")
    subprocess.call(["java", "-Xmx100g", "-jar", args.gatk_path, "-R", args.genome_fasta, "-T", "IndelRealigner", "-filterNoBases", "-maxReads","10000000","-I","possorted_chrM.bam","-targetIntervals","chrM_realignertargetcreator.intervals","-o","possorted_chrM_realign.bam"])

if args.data_type == "RNA":
    print("["+record_time()+"] Splits reads that contain Ns in their cigar string")
    #subprocess.call(["/data/R04/chenbzh5/ori/Tools/gatk-4.2.5.0/gatk","SplitNCigarReads","-R", args.genome_fasta, "-I", "possorted_chrM.bam", "-O", "possorted_chrM_splitN.bam"])

    print("["+record_time()+"] Replace read groups in a BAM file")
    #subprocess.call(["java", "-Xmx4g", "-jar", args.picard_path, "AddOrReplaceReadGroups","I=possorted_chrM_splitN.bam","O=possorted_chrM_splitN_addgp.bam","RGID=4","RGLB=lib1","RGPL=illumina","RGPU=unit1","RGSM=20"])

    print("["+record_time()+"] Index BAM file")
    subprocess.call(["samtools", "index", "possorted_chrM_splitN_addgp.bam"])

    print("["+record_time()+"] Create target intervals list using GATK RealignerTargetCreator")
    subprocess.call(["java", "-Xmx8g", "-jar", args.gatk_path, "-R", args.genome_fasta, "-T", "RealignerTargetCreator", "--filter_reads_with_N_cigar","-nt", "32", "-I", "possorted_chrM_splitN_addgp.bam", "-o", "chrM_realignertargetcreator.intervals"])

    print("["+record_time()+"] Realign reads using GATK IndelRealigner")
    subprocess.call(["java", "-Xmx100g", "-jar", args.gatk_path, "-R", args.genome_fasta, "-T", "IndelRealigner","-maxReads","10000000","-I","possorted_chrM_splitN_addgp.bam","-targetIntervals","chrM_realignertargetcreator.intervals","-o","possorted_chrM_realign.bam"])
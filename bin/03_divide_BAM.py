#!/bin/python

import pysam as ps
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial
import os
import datetime
import argparse

## Set paramter
parser = argparse.ArgumentParser()
parser.add_argument('--output_dir','-o', type=str,required=True, help='Output file path')
args = parser.parse_args()

output_dir=args.output_dir


bam_file=os.path.join(output_dir,"MitoSort/BAM/possorted_chrM_realign.bam")
temp_dir=os.path.join(output_dir,"MitoSort/temp")

## Make directory 
if not os.path.exists(temp_dir):
	os.mkdir(temp_dir)


## Defind function to divide BAM file based on MT positions
def fetch_region(regions):
	unsplit_file=bam_file
	start = regions
	end = regions + 2100
	out_dir=temp_dir
	outname = "_".join((str(start),str(end)))
	filename = out_dir + "/" +outname+ ".bam"
	samfile = ps.AlignmentFile( unsplit_file, "rb")
	split_file = ps.AlignmentFile( filename, "wb", template=samfile)
	for read in samfile.fetch('chrM',start,end):
		split_file.write(read)
	split_file.close()
	os.system("samtools index "+filename)


if __name__=='__main__':
	today_date = datetime.date.today()
	shown_date = today_date.strftime("%Y/%m/%d/")
	now = datetime.datetime.now()
	shown_time = now.strftime("%H:%M:%S")
	print("["+ shown_date+" "+shown_time+"] Generate temporary BAM files for faster computation")
	region_lit = [item*2100 for item in range(8)]
	pool = Pool(8)
	pool.map(fetch_region,region_lit)
	pool.close()
	pool.join()





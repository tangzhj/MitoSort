#!/bin/python

import os
import re
import datetime
import pandas as pd
import argparse

## Set paramters
parser = argparse.ArgumentParser()
parser.add_argument('--cell_barcode','-c', type=str,required=True, help='A CSV file containing barcode metrics')
parser.add_argument('--output_dir','-o', type=str,required=True, help='Output file path')
args = parser.parse_args()


cell_barcode_file=args.cell_barcode
output_dir=args.output_dir

## Make directory 
barcode_dir=os.path.join(output_dir,"MitoSort/barcode")
if not os.path.exists(barcode_dir):
	os.mkdir(barcode_dir)


today_date = datetime.date.today()
shown_date = today_date.strftime("%Y/%m/%d/")
now = datetime.datetime.now()
shown_time = now.strftime("%H:%M:%S")
print("["+ shown_date+" "+shown_time+"] Collecting cell barcodes")

barcode_dic = {}
mito_counts_dic = {}
genome_counts_dic = {}

bars = pd.read_csv(cell_barcode_file,header = 0)
if 'singlecell' in cell_barcode_file:
	for j in bars.index:
		if bars.loc[j,'is__cell_barcode'] != 1:
			continue
		ap_name =  bars.loc[j,'barcode']
		mito_counts_dic[ap_name] = bars.loc[j,'mitochondrial']
		genome_counts_dic[ap_name] = bars.loc[j,'passed_filters']

		if ap_name not in barcode_dic.keys():
			barcode_dic[ap_name] = 0
		else:
			barcode_dic[ap_name] += 1


if 'metrics' in cell_barcode_file:
	for j in bars.index:
		if bars.iloc[j,3] != 1:
			continue
		ap_name =  bars.loc[j,'barcode']
		mito_counts_dic[ap_name] = bars.loc[j,'atac_mitochondrial_reads'] #mito reads
		genome_counts_dic[ap_name] = bars.loc[j,'atac_fragments'] #genome reads

		if ap_name not in barcode_dic.keys():
			barcode_dic[ap_name] = 0
		else:
			barcode_dic[ap_name] += 1


output_file=barcode_dir + "/" + "barcode_result.txt"
f = open(output_file,"w")
f.write("barcode"+"\t"+"mito_reads"+"\t"+"genome_reads"+"\n")
for z in  barcode_dic.keys():
	final = str(z)+"\t"+str(mito_counts_dic[z])+"\t"+str(genome_counts_dic[z])+"\n"
	if  barcode_dic[z]==0:
		f.write(final)
	if  barcode_dic[z] > 0:
		print(final)
f.close()

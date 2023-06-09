
# -*- coding: utf-8 -*-
import time
import pysam
import random
from random import sample
from multiprocessing import Pool
import multiprocessing as mp
import os
import sys
# Define function to simulate reads for a single cell barcode

# Open a file to write the selected barcodes
def write_passed_cell(want_reads):

    with open("{}_{}_{}_{}_simulated.bam".format(num_cells_per_sample, want_reads,repeat1,sample_num)+"_barcodes.txt", "w") as outfile:
        outfile.write("sample"+"\t"+"barcode"+"\t"+"matching_reads"+"\t"+"want_reads"+"\t"+"mt_reads"+"\n")
        for j in dict_cell.keys():
            wanted = dict_cell[j]
            for i in range(need_cell_num):
                sample, barcodes,num_matching_reads ,want_reads,mt_reads = wanted[i]
                h = str(sample)+"_"+str(barcodes)
                #outfile.write("Sample: {}\n".format(sample))
                outfile.write("{}\t{}\t{}\t{}\t{}\n".format(sample, barcodes,num_matching_reads ,want_reads,mt_reads))  

    outfile.close()


def simulate_reads_for_cell_barcode(args,args2):

    if args2 == "silglet":
        # Open the merged BAM file
        sample,cell_barcode = args.split("_")
        sample_path1 = "../splitbam_raw/" + cell_barcode + ".bam"

        try:

            if os.path.isfile(sample_path1):
                # File exists, do something with it
                #os.system("samtools index "+sample_path1)
                cell_bam = pysam.AlignmentFile(sample_path1, "rb")
            else:
                print("######no bam:"+sample_path1)
                return 'a'
        except Exception as e:
            # Handle any other exceptions that may occur
            print("Error:", str(e))
        

        # Search for reads in the merged BAM file with the matching cell barcode

        matching_reads = []
        for read in cell_bam.fetch(until_eof=True):

            if read.has_tag("CB") and read.is_paired  :
            #if :
                matching_reads.append(read)
        cell_bam.close()




        if len(matching_reads) < want_reads :
            print('+++++++++++++++++++++++++++++++ no enough matching reads')
            return 'a'
        else:




            # Select a random set of matching reads
            selected_reads = random.sample(matching_reads, int(want_reads))


            mt_reads = 0
            for read_index, read in enumerate(selected_reads):
                if read.reference_name=="chrM":
                    mt_reads = mt_reads + 1
            if mt_reads< need_mt_reads:
                print('+++++++++++++++++++++++++++++++ no enough mitochondrial reads')
                return 'a'
            print("mt:",mt_reads)
            a = (sample,cell_barcode,len(matching_reads),want_reads,mt_reads)
        
            print(a)

            
            get_a_count.append(a)

            if sample not in dict_cell:
                dict_cell[sample] = []
            dict_cell[sample] = dict_cell[sample] + [a]

            if len(dict_cell[args.split("_")[0]]) > num_cells_per_sample:
                print('#######################')
                return 'enough'


        # Modify the read name to include the cell barcode and return the read

        for read_index, read in enumerate(selected_reads):
            #read.query_name = read.query_name + "_cell" + cell_barcode + "_read" + str(read_index)
            
            simulated_bam.write(read)

        return a


get_sample_count = []
get_a_count = []
dict_cell = {}
global a
(num_cells_per_sample1,want_reads1,repeat1,sample_num1) = (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
global want_reads
global need_cell_num
global repeat
global sample_num
global doublet_num
global need_mt_reads
num_cells_per_sample = int(num_cells_per_sample1)
need_cell_num = int(num_cells_per_sample)
want_reads = int(want_reads1)
repeat = int(repeat1)
sample_num = int(sample_num1)
need_mt_reads = 500
merged_bam = pysam.AlignmentFile("../zUMIs/CB_sorted_each_donor_sample_200_filtered.tagged.Aligned.out.bam", "rb")
# Create a new BAM file and write the simulated reads to it
#simulated_bam = pysam.AlignmentFile("{}_{}_simulated.bam".format(num_cells_per_sample, want_reads), "wb", template=merged_bam)





if __name__ == "__main__":


   



    start_time = time.time()









    #sample_names = sample(['CD34','BMMC','15#16','CRC','CCL1','sample1','sample7','lib2'],sample_num)
    sample_names = sample(['donor6','donor7','donor8'],sample_num)

    

    # Open the barcode list file and read the barcodes and read counts for each cell
    sample_barcodes = {}
    with open("../Smart-seq3xpress/PBMCs_run6_QCpass_each_donor_sample_200_barcodes_with_donor_annotation.txt", "r") as barcode_file:
        for line in barcode_file:
            if "barcode" in line:
                continue
            fields = line.strip().split("\t")
            sample = fields[1]
            if sample not in sample_names: # with random sample
                continue

            cell_barcode = fields[0]
            if sample not in sample_barcodes :
                sample_barcodes[sample] = []
            sample_barcodes[sample].append(cell_barcode)


    simulated_bam = pysam.AlignmentFile("{}_{}_{}_{}_simulated.bam".format(num_cells_per_sample,want_reads,repeat1,sample_num1), "wb", template=merged_bam)


    for sample, barcodes in sample_barcodes.items():
           
            # Select up to 1200 cell barcodes for this sample
            #args = []  
            selected_barcodes = random.sample(barcodes, min(num_cells_per_sample+100, len(barcodes)))

            for barcode in selected_barcodes:
                args = sample+'_'+str(barcode)
                # Simulate reads for the selected cell barcodes using multiple processes
                results = simulate_reads_for_cell_barcode(args,"silglet")


    for j in dict_cell.keys():
        wanted = dict_cell[j]
        print('################################\n')
        print(wanted)

    print(dict_cell.keys())

    write_passed_cell(want_reads)
    simulated_bam.close()
    end_time = time.time()
    running_time = end_time - start_time
    print("Running time: {:.2f} seconds".format(running_time))


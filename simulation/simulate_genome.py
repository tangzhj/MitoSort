# -*- coding: utf-8 -*-
import time
import pysam
import random
from random import sample
from multiprocessing import Pool
#from pathos.multiprocessing import ProcessingPoll as Pool
import multiprocessing as mp
import os
import sys

################
#This is used to simulate genome bam file.
#(num_cells_per_sample1,doublet_num1,want_reads1,repeat1,sample_num1)
#python simulate_genome.py 300 30 20000 1 8
'''
sh lib_run_snv.sh 300 20000 1 8
python generate_split_bam.py 300_20000_1_8
python standard_pool_test_split.py ./splitbam_raw/ ./bamfile/300_20000_1_8_simulated.sorted.snv ./barcode/300_20000_1_8_simulated.bam_barcodes.txt ./300_20000_1_8/
python standard_Demultiplex.py 300_20000_1_8 8
'''
################

# Open a file to write the selected barcodes
def write_passed_cell(want_reads):

    with open("{}_{}_{}_{}_simulated.bam".format(num_cells_per_sample, want_reads,repeat1,sample_num)+"_barcodes.txt", "w") as outfile:
        outfile.write("barcode"+"\t"+"CB_barcode"+"\t"+"matching_reads"+"\t"+"old_read"+"\t"+"need_reads"+"\n")
        for j in dict_cell.keys():
            wanted = dict_cell[j]
            for i in range(need_cell_num):
                sample, barcodes,num_matching_reads ,read_count,want_reads = wanted[i]
                h = str(sample)+"_"+str(barcodes)
                #outfile.write("Sample: {}\n".format(sample))
                outfile.write("{}_{}\t{}\t{}\t{}\t{}\n".format(sample, barcodes,CB_barcode_dic[h],num_matching_reads ,read_count,want_reads))  

    outfile.close()



# Generate {sample：barcode in splitbam of genome }
def generate_dic():
    samplename_barcode_of_genome = {} 
    #genome_spiltbam = ['CD34_genome_bam', 'BMMC_genome_bam', '15#16_genome_bam', 'CRC_genome_bam', 'CCL1_genome_bam', 'sample1_genome_bam', 'sample7_genome_bam', 'lib2_genome_bam']
    genome_spiltbam = ['CD34_genome_bam', 'BMMC_genome_bam', 'CRC_genome_bam', 'CCL1_genome_bam', 'sample1_genome_bam','lib2_genome_bam']
    for i in genome_spiltbam:
        k = i.split("_")[0]
        for j in  os.listdir("../simulate_reads/"+i):
            if k not in samplename_barcode_of_genome:
                print(k)
                samplename_barcode_of_genome[k] = []
            else:
                samplename_barcode_of_genome[k].append(j.replace(".bam",""))
    return(samplename_barcode_of_genome)



# Define function to simulate reads for a single cell barcode
def simulate_reads_for_cell_barcode(args,args2):

    if args2 == "silglet":
        # Open the merged BAM file
        sample,cell_barcode , read_count = args.split("_")
        sample_path1 = sample_dir[sample] + cell_barcode + ".bam"
        sample_path2 = sample_dir[sample] + cell_barcode.replace("-1","") + ".bam"
        try:

            if os.path.isfile(sample_path1):
                # File exists, do something with it
                #os.system("samtools index "+sample_path1)
                cell_bam = pysam.AlignmentFile(sample_path1, "rb")
            elif os.path.isfile(sample_path2):
                #os.system("samtools index "+sample_path2)
                cell_bam = pysam.AlignmentFile(sample_path2, "rb")
            else:
                print("######no bam:"+sample_path1)
                return 'a'
        except Exception as e:
            # Handle any other exceptions that may occur
            print("Error:", str(e))
        

        # Search for reads in the merged BAM file with the matching cell barcode

        matching_reads = []
        for read in cell_bam.fetch(until_eof=True):

            if read.has_tag("CR") and read.is_paired  :
            #if :
                matching_reads.append(read)
                #if args.split("_")[0] == "sample1" or sample == "sample7" or sample == "lib2":
                #    read.reference_id = 24
                 #   read.reference_name = "chrM"
        cell_bam.close()



        a = (sample,cell_barcode,len(matching_reads),read_count,want_reads)
        
        print(a)
        if len(matching_reads) < want_reads :
            print('+++++++++++++++++++++++++++++++ no enough matching reads')
            return a
        else:

            # Select a random set of matching reads
            selected_reads = random.sample(matching_reads, min(int(read_count), int(want_reads)))


            mt_reads = 0
            for read_index, read in enumerate(selected_reads):
                if read.reference_name=="chrM":
                    mt_reads = mt_reads + 1
            if mt_reads< need_mt_reads:
                print('+++++++++++++++++++++++++++++++ no enough matching reads')
                return a
            print("mt:",mt_reads)

            
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

    #-------------------------------------------#doublet
    else:
        print("come to doublet !!!!!!")
        # Open the merged BAM file
        sample,cell_barcode , read_count = args.split("_")
        sample_path1 = sample_dir[sample] + cell_barcode + ".bam"
        sample_path2 = sample_dir[sample] + cell_barcode.replace("-1","") + ".bam"
        try:

            if os.path.isfile(sample_path1):
                # File exists, do something with it
                #os.system("samtools index "+sample_path1)
                cell_bam = pysam.AlignmentFile(sample_path1, "rb")
            elif os.path.isfile(sample_path2):
                #os.system("samtools index "+sample_path2)
                cell_bam = pysam.AlignmentFile(sample_path2, "rb")
            else:
                print("######no bam:"+sample_path1)
                return 'a'
        except Exception as e:
            # Handle any other exceptions that may occur
            print("Error:", str(e))

        # Search for reads in the merged BAM file with the matching cell barcode

        matching_reads1 = []
        mt_reads1 = 0
        for read in cell_bam.fetch(until_eof=True):
            if read.reference_name=="chrM":
                mt_reads1 = mt_reads1 + 1
            if read.has_tag("CR") and read.is_paired  :
            #if :
                matching_reads1.append(read)
                #if args.split("_")[0] == "sample1" or sample == "sample7" or sample == "lib2":
                #    read.reference_id = 24
                #    read.reference_name = "chrM"
        cell_bam.close()


        #-------------------------------------------
        #second cell
        sample,cell_barcode , read_count = args2.split("_")
        sample_path1 = sample_dir[sample] + cell_barcode + ".bam"
        sample_path2 = sample_dir[sample] + cell_barcode.replace("-1","") + ".bam"
        try:

            if os.path.isfile(sample_path1):
                # File exists, do something with it
                #os.system("samtools index "+sample_path1)
                cell_bam2 = pysam.AlignmentFile(sample_path1, "rb")
            elif os.path.isfile(sample_path2):
                #os.system("samtools index "+sample_path2)
                cell_bam2 = pysam.AlignmentFile(sample_path2, "rb")
            else:
                print("######no bam:"+sample_path1)
                return 'a'
        except Exception as e:
            # Handle any other exceptions that may occur
            print("Error:", str(e))

         # Search for reads in the merged BAM file with the matching cell barcode

        matching_reads2 = []
        
        for read in cell_bam2.fetch(until_eof=True):
            
            if read.has_tag("CR") and read.is_paired:
            #if :
                matching_reads2.append(read)
                #if args2.split("_")[0] == "sample1" or sample == "sample7" or sample == "lib2":
                #    read.reference_id = 24
                #    read.reference_name = "chrM"

        cell_bam2.close()

        #-------------------------------------------


        new_read_count = args1.split("_")[2]+"+"+args2.split("_")[2]+"+"+args2.split("_")[0]+"_"+args2.split("_")[1]
        new_matching_reads = str(len(matching_reads1)) +"+"+str(len(matching_reads2))
        a = (args.split("_")[0],args.split("_")[1],new_matching_reads,new_read_count,want_reads)


        print(a)
        if len(matching_reads1) < want_reads or len(matching_reads2) < want_reads :
            print('+++++++++++++++++++++++++++++++ no enough matching reads')
            return a


        else:
        
            # Select a random set of two matching reads
            selected_reads1 = random.sample(matching_reads1,  int(want_reads))
            selected_reads2 = random.sample(matching_reads2,  int(want_reads))

            #check mt reads not too low
            mt_reads1 = 0
            for read_index, read in enumerate(selected_reads1):
                if read.reference_name=="chrM":
                    mt_reads1 = mt_reads1 + 1
            mt_reads2 = 0
            for read_index, read in enumerate(selected_reads2):
                #print(read.reference_name)
                if read.reference_name=="chrM":
                    mt_reads2 = mt_reads2 + 1
            print("mt1:g1:mt2:g2:{}\t{}\t{}\t{}".format(mt_reads1,len(selected_reads1),mt_reads2,len(selected_reads2)))
            d1 = mt_reads1/float(len(selected_reads1))
            d2 = mt_reads2/float(len(selected_reads2))
            print("d1d2",d1,d2)
            if  (d1/d2)<0.2 or (d1/d2) > 5.0:
                print('+++++++++++++++++++++++++++++++ no enough matching reads')
                return a




            
            if args.split("_")[0] not in dict_cell:
                dict_cell[args.split("_")[0]] = []
            dict_cell[args.split("_")[0]] = dict_cell[args.split("_")[0]] + [a]

            if len(dict_cell[args.split("_")[0]]) > num_cells_per_sample:
                print('#######################enough')
                return 'enough'

            # Modify the read name to include the cell barcode and return the read

            for read_index, read in enumerate(selected_reads1):
                #read.query_name = read.query_name + "_cell" + cell_barcode + "_read" + str(read_index)

                
                simulated_bam.write(read)


            for read_index, read in enumerate(selected_reads2):
                read.set_tag("CR",args.split("_")[1])
                read.query_name = read.query_name + "_cell" + cell_barcode + "_read" + str(read_index)
                
                simulated_bam.write(read)

            return a


##

get_sample_count = []
get_a_count = []
dict_cell = {}
global a
(num_cells_per_sample1,doublet_num1,want_reads1,repeat1,sample_num1) = (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
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
doublet_num = int(doublet_num1)
need_mt_reads = 500
merged_bam = pysam.AlignmentFile("../merged_10people.bam", "rb")
# Create a new BAM file and write the simulated reads to it
#simulated_bam = pysam.AlignmentFile("{}_{}_simulated.bam".format(num_cells_per_sample, want_reads), "wb", template=merged_bam)



if __name__ == "__main__":

    start_time = time.time()

    sample_dir = {
    'CD34':'../CD34_genome_bam/',
    'BMMC':'../BMMC_genome_bam/',
    '15#16':'../15#16_genome_bam/',
    'CRC':'../CRC_genome_bam/',
    'CCL1':'../CCL1_genome_bam/',
    'sample1':'../sample1_genome_bam/',
    'sample7':'../sample7_genome_bam/',
    'lib2':'../lib2_genome_bam/'
    }    


    #sample_names = sample(['CD34','BMMC','15#16','CRC','CCL1','sample1','sample7','lib2'],sample_num)
    sample_names = sample(['CD34','BMMC','CRC','CCL1','sample1','lib2'],sample_num)

    samplename_barcode_of_genome = generate_dic()
    

    # Open the barcode list file and read the barcodes and read counts for each cell
    sample_barcodes = {}
    CB_barcode_dic = {}
    with open("../barcode_result2.txt", "r") as barcode_file:
        for line in barcode_file:
            if "barcode" in line:
                continue
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            sample = fields[0].split("_")[0]
            if sample not in sample_names: # with random sample
                continue

            cell_barcode = fields[0].split("_")[1].replace("-1","")
            CB_barcode = fields[1].split("_")[1]
            read_count = int(fields[2])
            genome_read_count = int(fields[3])
            read_count = genome_read_count + read_count
            dumplicate = float(fields[4])
            #print(cell_barcode,samplename_barcode_of_genome[sample])

            if cell_barcode not in samplename_barcode_of_genome[sample]:#only passed filtered genome bam of cell
                continue

            if read_count*3 > want_reads and dumplicate < 0.6 :
                if sample not in sample_dir.keys():
                    continue
                if sample not in sample_barcodes :
                    sample_barcodes[sample] = []
                sample_barcodes[sample].append((cell_barcode, read_count))


            if fields[0] not in CB_barcode_dic.keys():
                CB_barcode_dic[fields[0].replace("-1","")] = fields[1]

    # Define the number of reads per cell
    num_reads_per_cell = want_reads
    # count passed cell num
    #print(CB_barcode_dic.keys())
    for i in range(len(sample_names)):
        print(sample_names[i],len(sample_barcodes[sample_names[i]]),"Open the barcode list file")




    # Loop through each sample and simulate reads for the selected barcodes

    #pool = Pool(processes=5)

    simulated_bam = pysam.AlignmentFile("{}_{}_{}_{}_simulated.bam".format(num_cells_per_sample,want_reads,repeat1,sample_num1), "wb", template=merged_bam)

    if 1:
        for sample, barcodes in sample_barcodes.items():
           
            # Select up to 1200 cell barcodes for this sample
            #args = []  
            selected_barcodes2 = random.sample(barcodes, min(num_cells_per_sample+500, len(barcodes)))
            selected_barcodes = list(set(selected_barcodes2))

            ## random generate doublet barcode for selected_barcodes

            selet_to_be_doublet = random.sample(selected_barcodes,int(doublet_num))


            for barcode in selected_barcodes:
                #generate doublet
                if  barcode in selet_to_be_doublet:
                    #get pooled cell from other sample
                    #a = ['CD34','BMMC','15#16','CRC','CCL1','sample1','sample7','lib2']
                    a = ['CD34','BMMC','CRC','CCL1','sample1','lib2']
                    a.remove(sample)

                    doublet_sample = random.sample(a,1)[0]

                    doublet_barcode = random.sample(samplename_barcode_of_genome[doublet_sample],1)[0]
                    if doublet_barcode not in selected_barcodes:    

                        args2 = doublet_sample+'_'+str(doublet_barcode)+"_"+str(random.randint(5000,5000)) # this is doublet barcode
                        


                        args1 = sample+'_'+str(barcode[0])+"_"+str(barcode[1]) # this is singlet barcode
                        results = simulate_reads_for_cell_barcode(args1,args2)
                    else:
                        continue


                #generate singlet
                else:
                    print(sample)
                    args = sample+'_'+str(barcode[0])+"_"+str(barcode[1])
                    # Simulate reads for the selected cell barcodes using multiple processes
                    results = simulate_reads_for_cell_barcode(args,"silglet")

        

    #pool.close()
    #pool.join()
    
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

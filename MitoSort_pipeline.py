#!/data/R04/chenbzh5/miniconda2/envs/MitoSort/bin/python

import os
import pysam as ps
import pandas as pd
import numpy as np
import datetime



def record_time():
    today_date = datetime.date.today()
    shown_date = today_date.strftime("%Y/%m/%d/")
    now = datetime.datetime.now()
    shown_time = now.strftime("%H:%M:%S")
    value = shown_date+" "+shown_time
    return value

def fetch_region(bam_file, temp_dir,regions):
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

def fetch_region_wrapper(args):
    return fetch_region(*args)

def site_to_file(Position):##caculate site from which file
    region = [item*2100 for item in range(8)]
    diff = [abs(item*2100-Position) for item in range(8)]
    min_region = min([abs(item*2100-Position) for item in range(8)])
    seletect_region = region[diff.index(min_region)]
    if Position >= seletect_region : 
        return(seletect_region,seletect_region+2100)
    if Position < seletect_region :
        return(seletect_region-2100,seletect_region)

def fetch_site(temp_dir,barcodes,tag,site):
    Position = int(site.split("-")[0])##site position
    start_p,end_p = site_to_file(Position)
    outname = "_".join((str(start_p),str(end_p)))
    filename = temp_dir + "/" +outname+ ".bam"
    in_sam = ps.AlignmentFile(filename, 'rb')
    Ref = site.split("-")[1]
    VarAllele = site.split("-")[2]
    ref_base_calls_mtx = pd.DataFrame(0, index=[site], columns=barcodes, dtype=np.int16)
    alt_base_calls_mtx = pd.DataFrame(0, index=[site], columns=barcodes, dtype=np.int16)
    count = 0
    Position = Position-1
    test_dic1 = {}
    test_dic2= {}
    for i in barcodes:
        test_dic1[i] = 0
        test_dic2[i] = 0
    for read in in_sam.fetch('chrM', Position, Position+1):
        if read.has_tag(tag):
            if (Position) in read.get_reference_positions():
                base_index =  read.get_reference_positions().index(Position)
            # if the read aligned positions cover the SNV position
                barcode = read.get_tag(tag)
                if barcode in barcodes:
                
                    count += 1
                    base = read.query_sequence[ read.get_aligned_pairs(True)[base_index][0]]
                    if base == Ref:
                        test_dic1[barcode] += 1
                    elif base == VarAllele:
                        test_dic2[barcode] += 1
    ref_base_calls_mtx  = pd.DataFrame([test_dic1],index=[site])
    alt_base_calls_mtx  = pd.DataFrame([test_dic2],index=[site])
    in_sam.close()
    return(ref_base_calls_mtx,alt_base_calls_mtx)

def fetch_site_wrapper(args):
    return fetch_site(*args)


import click
@click.group()
def cli():
    pass

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--bam_file', '-b', required=True, type=click.Path(exists=True),help='BAM file output by Cellranger')
@click.option('--genome_fasta', '-f', required=True, type=click.Path(exists=True),help='FASTA file of reference genome')
@click.option('--gatk_path',required=True, type=click.Path(exists=True),help='Path to GenomeAnalysisTK(GATK)')
@click.option('--output_dir', '-o', required=True, type=click.Path(),
              help='Output parent directory')
@click.option('--data_type',required=False, default="ATAC",help='Type of data. Default is ATAC. Set to RNA if your data is scRNA')
def MT_realign(bam_file,genome_fasta,gatk_path,output_dir,data_type):
    """
    Realign mitochondrial reads using GATK
    """
    if not bam_file:
        click.echo('Error: BAM file is required. Please specify it by `-b <bam_file>`')
        sys.exit(1)

    if not genome_fasta:
        click.echo('Error: FASTA file of reference genome is required. Please specify it by `-f <genome_fasta>`')
        sys.exit(1)

    if not gatk_path:
        click.echo('Error: Path to GenomeAnalysisTK(GATK) is required. Please specify it by `--gatk_path <gatk_path>`')
        sys.exit(1)

    if not output_dir:
        click.echo('Error: Output directory is required. Please specify it by `-o <output_dir>`')
        sys.exit(1)

    click.echo('BAM file: {}'.format(bam_file))
    click.echo('FASTA file of reference genome: {}'.format(genome_fasta))
    click.echo('Path to GenomeAnalysisTK(GATK): {}'.format(gatk_path))
    click.echo('Output parent directory: {}'.format(output_dir))
    click.echo('Data type: {}'.format(data_type))

    print("checking modules")
    import subprocess 
    import os
    import datetime
    import time
    print("imports done")

    start = time.time()
    # make parent directory
    MitoSort_output_dir = os.path.join(output_dir, "MitoSort")
    if not os.path.exists(MitoSort_output_dir):
        os.mkdir(MitoSort_output_dir)


    # make output directory of BAM files
    BAM_output_dir = os.path.join(output_dir, "MitoSort", "BAM")
    if not os.path.exists(BAM_output_dir):
        os.mkdir(BAM_output_dir)
    os.chdir(BAM_output_dir)

    # Sort BAM file by position if it is not sorted
    if "SO:coordinate" in subprocess.check_output(["samtools", "view", "-H", bam_file]).decode():
        print("[Log] The input BAM file has been sorted by position")
    else:
        print("["+record_time()+"] Sort BAM file by position")
        subprocess.call(["samtools", "sort", "-o", "possorted.bam", bam_file], check=True)
        print("["+record_time()+"] Index sorted BAM file")
        subprocess.call(["samtools", "index", "possorted.bam"])

    # Subset reads mapped to chrM
    if "SO:coordinate" in subprocess.check_output(["samtools", "view", "-H", bam_file]).decode():
        print("["+record_time()+"] Subset reads mapped to chrM")
        subprocess.call(["samtools", "view", "-b", bam_file, "chrM", "-o", "possorted_chrM.bam"])
    else:
        print("["+record_time()+"] Subset reads mapped to chrM")
        subprocess.call(["samtools", "view", "-b", "possorted.bam", "chrM", "-o", "possorted_chrM.bam"])

    # Index BAM file that containing reads mapped to chrM
    print("["+record_time()+"] Index BAM file that containing reads mapped to chrM")
    subprocess.call(["samtools", "index", "possorted_chrM.bam"])
    '''
    if data_type == "ATAC":
        print("["+record_time()+"] Create target intervals list using GATK RealignerTargetCreator")
        subprocess.call(["java", "-Xmx8g", "-jar", gatk_path, "-R", genome_fasta, "-T", "RealignerTargetCreator", "-nt", "10", "-I", "possorted_chrM.bam", "-o", "chrM_realignertargetcreator.intervals"])

        # Realign reads using IndelRealigner
        print("["+record_time()+"] Realign reads using GATK IndelRealigner")
        subprocess.call(["java", "-Xmx100g", "-jar", gatk_path, "-R", genome_fasta, "-T", "IndelRealigner", "-filterNoBases", "-maxReads","10000000","-I","possorted_chrM.bam","-targetIntervals","chrM_realignertargetcreator.intervals","-o","possorted_chrM_realign.bam"])
    '''

    end = time.time()
    elapsed = end-start
    print("Total time:",str(time.strftime("%Hh%Mm%Ss", time.gmtime(elapsed))))


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--bam_file', '-b', required=True, type=click.Path(exists=True),help='BAM file containing realigned MT reads `possorted_chrM_realign.bam`')
@click.option('--genome_fasta', '-f', required=True, type=click.Path(exists=True),help='FASTA file of reference genome')
@click.option('--cell_barcode', '-c', required=True, type=click.Path(exists=True),help="A CSV file containing barcode metrics output by Cellranger (singlecell.csv/per_barcode_metrics.csv). If it isn't generated by Cellranger, please input a file cotaining cell barcodes with header named `barcode`")
@click.option('--chrm_length', '-m', required=True, type=click.Path(exists=True),help='a BED file containing chrM region')
@click.option('--varscan_path', required=True, type=click.Path(exists=True),help='Path to VarScan')
@click.option('--output_dir', '-o', required=True, type=click.Path(),
              help='Output parent directory')
@click.option("--cell_tag", required = False, default = "CR", help = "Default is CR. Set if your cell barcode tag is not CR")
def Generate_SNP_matrix(bam_file,genome_fasta,chrm_length,varscan_path,cell_barcode,output_dir,cell_tag):
    """
    Generate SNP matrices
    """
    if not bam_file:
        click.echo('Error: BAM file containing realigned MT reads is required. Please specify it by `-b <bam_file>`')
        sys.exit(1)

    if not genome_fasta:
        click.echo('Error: FASTA file of reference genome is required. Please specify it by `-f <genome_fasta>`')
        sys.exit(1)

    if not chrm_length:
        click.echo('Error: a BED file containing chrM region is required. Please specify it by `-m <chrm_length>`')
        sys.exit(1)

    if not varscan_path:
        click.echo('Error: Path to VarScan is required. Please specify it by `--varscan_path <varscan_path>`')
        sys.exit(1)

    if not cell_barcode:
        click.echo('Error: A CSV file containing barcode metrics is required. Please specify it by `-c <cell_barcode>`')
        sys.exit(1)

    if not output_dir:
        click.echo('Error: Output directory is required. It should be the same directory when running `MT_realign`command. Please specify it by `-o <output_dir>`')
        sys.exit(1)

    print("checking modules")
    import subprocess 
    import os
    import re
    import pysam as ps
    import pandas as pd
    import numpy as np
    from multiprocessing import Pool
    from functools import partial
    from collections import defaultdict
    import sys
    import datetime
    import functools
    import dill
    import time 
    print("imports done")

    start = time.time()
    BAM_output_dir = os.path.join(output_dir, "MitoSort", "BAM")
    os.chdir(BAM_output_dir)

    # subset part of reads for faster computation
    print("["+record_time()+"] Subset reads for faster computation")
    subset_cmd = ["samtools", "view", "-s", "0.2", "-b",bam_file,">",os.path.join(BAM_output_dir, "possorted_chrM_realign_0.2.bam")]
    subprocess.call(" ".join(subset_cmd), shell=True)

    # Generate text pileup output for chrM BAM file, with each line representing the pileup of reads at a single genomic position
    print("["+record_time()+"] Pileup reads at chrM region")
    mpileup_cmd = ["samtools", "mpileup", "-l", chrm_length, "-q", "30", "-Q", "30", "-f", genome_fasta, "-x", os.path.join(BAM_output_dir, "possorted_chrM_realign_0.2.bam"), "-o", os.path.join(BAM_output_dir, "possorted_chrM_realign_0.2.mpileup")]
    subprocess.call(" ".join(mpileup_cmd), shell=True)

    # Call SNPs by VarScan pileup2snp 
    print("["+record_time()+"] Call SNPs by VarScan pileup2snp")
    pileup2snp_cmd = ["java", "-Xmx4g", "-jar", varscan_path, "pileup2snp", os.path.join(BAM_output_dir, "possorted_chrM_realign_0.2.mpileup"), "--min-var-freq", "0.01", "--min-reads2", "2", ">", os.path.join(BAM_output_dir, "possorted_chrM_realign.snv")]
    subprocess.call(" ".join(pileup2snp_cmd), shell=True)

    # make temp directory 
    temp_dir=os.path.join(output_dir,"MitoSort/temp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    # divide BAM file based on MT positions for faster computation
    print("["+record_time()+"] Generate temporary BAM files for faster computation")
    region_lit = [(bam_file, temp_dir, item*2100) for item in range(8)]
    pool = Pool(8)
    results = pool.imap(fetch_region_wrapper, region_lit)
    pool.close()
    pool.join()

    cell_barcode_file=cell_barcode
    # make barcode directory 
    barcode_dir=os.path.join(output_dir,"MitoSort/barcode")
    if not os.path.exists(barcode_dir):
        os.mkdir(barcode_dir)

    # collect cell barcodes
    print("["+record_time()+"] Collecting cell barcodes")
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
    elif 'metrics' in cell_barcode_file:
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
    else:
        output_file=barcode_dir + "/" + "barcode_result.txt"
        cp_cmd = ["cp", cell_barcode_file, output_file]
        subprocess.call(" ".join(cp_cmd), shell=True)

    # Make SNP directory 
    SNP_matrix_dir=os.path.join(output_dir,"MitoSort/SNP_matrix")
    if not os.path.exists(SNP_matrix_dir):
        os.mkdir(SNP_matrix_dir)

    # generate SNP matrix
    print("["+record_time()+"] Generating SNP matrices")
    bam_dir=os.path.join(output_dir,"MitoSort/temp/")
    temp_dir=os.path.join(output_dir,"MitoSort/temp")
    snv_file=os.path.join(output_dir,"MitoSort/BAM/possorted_chrM_realign.snv")
    barcode_file=os.path.join(output_dir,"MitoSort/barcode/barcode_result.txt")
    filtered_vcf = pd.read_csv(snv_file,sep="\t") ##germline snp
    row_lit = []
    blacklit = [302,309,312,313,316,514,515,523,524,3106,3107,3109,3110,16181]


    f = pd.read_csv(barcode_file,sep="\t") ##cellbarcodes
    barcodes = {}
    if cell_tag == "CR":


        for position, row in f.iterrows():
            i = str(row['barcode']).replace("-1","")##cell barcode in multi is CR and not '-1' in cell barcodes.
            barcodes[i.strip()] = 0  
    else:
        for position, row in f.iterrows():
            i = str(row['barcode'])
            barcodes[i.strip()] = 0
    barcodes_lit = list(barcodes)

    for position, row in filtered_vcf.iterrows():
        site_infor = str(row['Position']),row['Ref'],row['VarAllele']
        fre = float(row['VarFreq'].split("%")[0])
        if fre > 1.0 and fre <99.0 and row['Position'] not in blacklit:
            row_lit.append((temp_dir,barcodes,cell_tag,"-".join(site_infor)))

    #tag=cell_tag # default is CB
    bam = bam_dir ##8 bam dir
    result_frequency = os.path.join(SNP_matrix_dir,"frequency.csv")
    result_alt = os.path.join(SNP_matrix_dir,"alt.csv")
    result_ref = os.path.join(SNP_matrix_dir,"ref.csv")

    pool = Pool(30)
    result_list = pool.map(fetch_site_wrapper,row_lit)
    pool.close()
    pool.join()
    ref_matrix = pd.DataFrame(columns=barcodes, dtype=np.int16)
    alt_matrix = pd.DataFrame(columns=barcodes, dtype=np.int16)
    for j in result_list:
        ref_matrix = pd.concat([ref_matrix,j[0]],axis=0,sort=True)
        alt_matrix = pd.concat([alt_matrix,j[1]],axis=0,sort=True)

    frequency_matrix = alt_matrix/(ref_matrix+alt_matrix)
    frequency_matrix.to_csv(result_frequency)
    ref_matrix.to_csv(result_ref)
    alt_matrix.to_csv(result_alt)
    end = time.time()
    elapsed = end-start
    print("Total time:",str(time.strftime("%Hh%Mm%Ss", time.gmtime(elapsed))))


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--output_dir', '-o', required=True, type=click.Path(),
              help='Output parent directory')
@click.option('--clusters', '-k', required=True,help='number of pooled individuals')
@click.option("--p1_cutoff", required = False, default = "0.9", help = "maximum cutoff of p1 for doublet identification. Default to be 0.9")
@click.option("--p2_cutoff", required = False, default = "0.1", help = "minimun cutoff of p2 for doublet identification. Default to be 0.1")
@click.option("--method", required = False, default = "full", help = "if not 'full' it would be 'direct',which direct caculate p_value. ")
@click.option("--depth_cutoff", required = False, default = "1", help = " the depth percell lower than depth_percell get 'unassign' ")
def demultiplex(output_dir,clusters,p1_cutoff,p2_cutoff,method,depth_cutoff):
    """
    Identifying cross-genotype doublets and demultiplexing samples
    """
    if not output_dir:
        click.echo('Error: Output directory is required. It should be the same directory when running `Generate_SNP_matrix` command. Please specify it by `-o <output_dir>`')
        sys.exit(1)

    if not clusters:
        click.echo('Error: the number of pooled individuals is required. Please specify it by `-k <clusters>`. If the number of pooled individuals is unkown, you can first run the step of demultiplex with random k. Based on the silhouette score, you can rerun the the step of demultiplex with optimal k' )
        sys.exit(1)

    print("checking modules")
    import os
    import numpy as np
    from numpy import pi
    import matplotlib
    matplotlib.use('Agg')
    from sklearn.cluster import KMeans
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from sklearn import datasets
    from scipy.cluster.hierarchy import linkage, dendrogram
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    from scipy import stats
    from sklearn.mixture import GaussianMixture as GMM
    from sklearn.decomposition import PCA
    import random
    from numpy import random
    import time
    from random import sample
    from sklearn.metrics import confusion_matrix   
    from sklearn.metrics import accuracy_score, average_precision_score,precision_score,f1_score,recall_score
    from sklearn.metrics import silhouette_score
    import sys
    from sklearn.neighbors import KNeighborsClassifier
    from jinja2 import Template
    import argparse
    import datetime
    import time
    sys.setrecursionlimit(100000)
    print("imports done")
    start = time.time()
    demultiplex_dir=os.path.join(output_dir,"MitoSort/Demultiplex_output")
    if not os.path.exists(demultiplex_dir):
        os.mkdir(demultiplex_dir)
    output_path = demultiplex_dir

    SNP_matrix_dir=os.path.join(output_dir,"MitoSort/SNP_matrix")
    path1 = os.path.join(SNP_matrix_dir,"alt.csv")
    path2 = os.path.join(SNP_matrix_dir,"ref.csv")
    path3 = os.path.join(SNP_matrix_dir,"frequency.csv")
    barcode_dir=os.path.join(output_dir,"MitoSort/barcode")
    path4 = os.path.join(barcode_dir,"barcode_result.txt")

    sample_num = int(clusters)
    confident_germline_ratio = 0.01 #0.5
    D_x = float(p1_cutoff)
    D_y = float(p2_cutoff)

    def read_file(path1,path2,path3,path4):

        #read csv
        alt_data = pd.read_csv(path1)
        alt_data.head()

        ref_data = pd.read_csv(path2,header=0)
        ref_data.head()

        data3 = pd.read_csv(path3,sep = ",",header = 0)
        data3.head()

        alt_data = alt_data.drop_duplicates(keep='first')
        ref_data = ref_data.drop_duplicates(keep='first')
        data3 = data3.drop_duplicates(keep='first')

        matrix_ref = ref_data.T
        matrix_ref.columns = matrix_ref.iloc[0,:]
        matrix_ref = matrix_ref.iloc[1:,0:]
        matrix_ref = pd.DataFrame(matrix_ref,dtype=np.int64)
        matrix_ref = matrix_ref[matrix_ref>0]
        matrix_ref = matrix_ref.where(matrix_ref.notnull(), 0)
        matrix_ref.dtypes

        matrix_alt = alt_data.T
        matrix_alt.columns = matrix_alt.iloc[0,:]
        matrix_alt = matrix_alt.iloc[1:,0:]
        matrix_alt = pd.DataFrame(matrix_alt,dtype=np.int64)
        matrix_alt = matrix_alt[matrix_alt>0]
        matrix_alt = matrix_alt.where(matrix_alt.notnull(), 0)
        matrix_alt.dtypes

        matrix_fre = matrix_alt/(matrix_ref+matrix_alt)
        matrix_fre = matrix_fre.where(matrix_fre.notnull(), 0)

        matrix_depth = (matrix_ref+matrix_alt)
        matrix_depth = matrix_depth.where(matrix_depth.notnull(), 0)
        matrix_depth = matrix_depth[matrix_fre<=1]
        matrix_depth = matrix_depth.where(matrix_depth.notnull(), 0)

        fragment_data = pd.read_csv(path4,sep = "\t",header = 0,index_col='barcode')
        fragment_data.head()

        return(matrix_alt,matrix_ref,matrix_fre,matrix_depth,fragment_data)

    def run_pip(matrix_alt,matrix_ref,matrix_fre,matrix_depth,sample_num,confident_germline_ratio,D_x,D_y):
        print("Filter variants!")
        matrix_fre_pre = matrix_fre #copy
        matrix_fre  = pd.DataFrame(matrix_fre)

        barcode_lit =  matrix_fre.index

        ##methods2
        a = np.sum(matrix_fre.loc[barcode_lit,:]>0.95,axis=0)/np.sum(matrix_fre.loc[barcode_lit,:]>0.05,axis=0)##>0.05
        matrix_alt = matrix_alt.where(matrix_alt.notnull(), 0)
        seletec_germ = []
        for i in a.index:
            #print(i,a[i])
            if a[i] > confident_germline_ratio:
                seletec_germ.append(i)
        
        #generate final matrix.
        #need to filter low depth
        def filtered_depth_again(matrix_depth):
            depth_percell = np.mean(matrix_depth.loc[:,:],axis = 1)

            low_barcode = []
            for i in range(len(depth_percell)):
                if depth_percell[i] < 5:
                    low_barcode.append(depth_percell.index[i])

            #can filter barcode
            barcode_lit = matrix_depth.index
          
            barcode_lit =  list(set(barcode_lit) -set(low_barcode)) #barcode_lit 
            print(len(matrix_fre.index),len(barcode_lit),len(low_barcode))
            return(low_barcode,barcode_lit)
        ##
        (low_barcode,barcode_lit2) = filtered_depth_again(matrix_depth)

        To_GMM_barcode = barcode_lit
        To_GMM_barcode2 = barcode_lit2

        To_GMM_matrix = matrix_fre.loc[barcode_lit,seletec_germ]

        high_depth_matrix = matrix_fre.loc[barcode_lit2,seletec_germ]


        #final_matrix to GMM or kmean.

        X = np.array(high_depth_matrix) #filtered doublet & low_depth

        ##define cluster by kmeans.
        sample_num = sample_num
        print('sample_num:',sample_num)


        sil = []
        kmax = 12

        # dissimilarity would not be defined for a single cluster, thus, minimum number of clusters should be 2
        for k in range(2, kmax+1):
            kmeans = KMeans(n_clusters = k).fit(X)
            labels = kmeans.labels_
            sil.append(silhouette_score(X, labels, metric = 'euclidean'))
        fig,ax = plt.subplots(figsize = (8,6))
        ax.set_facecolor("white")
        ax.grid(ls="")
         # Convert bottom-left and top-right to display coordinates
        x_0, y_0 = ax.transAxes.transform((0, 0))
        x_1, y_1 = ax.transAxes.transform((1, 1))

        # Convert back to Axes coordinates
        x_0, y_0 = ax.transAxes.inverted().transform((x_0, y_0))
        x_1, y_1 = ax.transAxes.inverted().transform((x_1, y_1))
        width = 1
        rect = plt.Rectangle(
            (x_0, y_0), x_1-x_0, y_1-y_0,
            color="black",
            transform=ax.transAxes,
            zorder=-1,
            lw=2*width+0.02,
            fill=None,
        )
        ax.patches.append(rect)
        sil_max_index = sil.index(max(sil))
        if sil[sil_max_index] - sil[sil_max_index-1] < 0.001:
            
            #print("sil:",sil[sil_max_index],sil[sil_max_index-1])
            sil_max_index = sil_max_index-1
            
        plt.axvline(x = sil_max_index+2, ymin = 0.05, ymax = 1.0, color = 'r',ls='--',label = "Best k="+str(sil_max_index+2))
        
        plt.legend(loc = "best",fontsize=15)
        plt.plot(range(2,kmax+1,1), sil, '-p', color='gold')
        plt.xticks(range(1,kmax+1,1))
        plt.xlabel('number of clusters, k')
        plt.ylabel('silhouette_score')
        plt.savefig("{}/sihouette_score.pdf".format(output_path),dpi = 400,bbox_inches = 'tight') # save figure 
        plt.show()
        
        if len(barcode_lit)>2000:
            ab = sns.clustermap(To_GMM_matrix.iloc[:2000,:], metric="euclidean")
            ab.ax_heatmap.set_xlabel("Site")
            ab.savefig("{}/raw_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') #save figure
        else:
            ab = sns.clustermap(To_GMM_matrix, metric="euclidean")
            ab.ax_heatmap.set_xlabel("Site")
            ab.savefig("{}/raw_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') #save figure


        kmeans = KMeans(n_clusters = sample_num,max_iter=1000,tol=1e-10,algorithm='full',n_init=100)
        kmeans.fit(X)
        labels = kmeans.labels_
        kmeans.labels_


        ##get specific germline
        print("Get specific germline!")
        createVar = locals()
        myVar1 = [] # classified samples
        for i in range(sample_num):
            createVar['cell_'+ str(i)] = []
            myVar1.append('cell_'+ str(i)) 
        #print(len(labels))
        #print(len(To_GMM_barcode))
        for i in range(len(labels)):
            for j in range(sample_num):
                if labels[i]==j:
                    locals()['cell_'+ str(j)].append(To_GMM_barcode2[i])



        myVar2 = [] # specific variants for each sample
        for i in range(sample_num):
            createVar['cell_'+ str(i)+"_site"] = []
            myVar2.append('cell_'+ str(i)+"_site")

        for k in range(sample_num):
            #print(locals()['cell_'+ str(j)][1:20])
            jkl = np.sum(high_depth_matrix.loc[locals()["cell_"+str(k)],:][high_depth_matrix.loc[locals()["cell_"+str(k)],:]>=0.99],axis = 0)/(len(locals()["cell_"+str(k)]))
            for i in range(len(high_depth_matrix.columns)):
                if jkl[i]>0.1:
                    locals()["cell_"+str(k)+"_site"].append(high_depth_matrix.columns[i])



        site_dic = []
        for k in range(sample_num):
            site_dic = site_dic + locals()["cell_"+str(k)+"_site"]


        dicts = {}
        for key in site_dic:
            dicts[key] = dicts.get(key, 0) + 1
        double_lit = []
        for key in dicts.keys():
            if dicts[key] >1:
                double_lit.append(key)

        useful_site=list(set(dicts)-(set(double_lit)))
        print(len(useful_site))

        #get clean_matrix
        clean_matrix = To_GMM_matrix.reindex(columns = useful_site)

        if len(barcode_lit)>2000:
            ad = sns.clustermap(clean_matrix.iloc[:2000,:], metric="euclidean")
            ad.ax_heatmap.set_xlabel("Site")
            ad.savefig("{}/clean_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') 
        else:
            ad = sns.clustermap(clean_matrix, metric="euclidean")
            ad.ax_heatmap.set_xlabel("Site")
            ad.savefig("{}/clean_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') 

        myVar3 = [] # specific variants for each sample
        for i in range(sample_num):
            createVar['single_cell_'+ str(i)] = []
            myVar3.append('single_cell_'+ str(i)) 
            locals()['single_cell_'+ str(i)] = set(locals()["cell_"+str(i)+"_site"]).difference(set(double_lit))
            print(locals()['single_cell_'+ str(i)])
            #print( locals()['cell_'+ str(i)][1:5],locals()['single_cell_'+ str(i)])


        myVar4 = [] # martix for each sample
        for i in range(sample_num):
            createVar['cell_'+ str(i)+"_clean_matrix"] = []  #cell_0_clean_matrix 
            myVar4.append('cell_'+ str(i)+"_clean_matrix")

        for j in range(sample_num):
            locals()['cell_'+ str(j)+"_clean_matrix"] = clean_matrix.loc[:,locals()['single_cell_'+ str(j)]]

        ##caculate p-value in first fitted cluster and second fitted
        print("caculate p-value!")
        tag = 0
        var_lit = myVar3
        cell_lit = myVar1
        for j in low_barcode:
            locals()['cell_'+ str(sample_num-1)].append(j)
        #def caculate_p_each_k(cell_lit,var_lit,matrix_alt,matrix_ref,):
        if 1:
            cell_pvalue_dic = {}
            cell_pvalue_dic1 = {}
            unassign_cell  = []
            for c in cell_lit:
                cells = locals()[str(c)]
                for a_cell in cells:
                    if a_cell in clean_matrix.index:
                        cell_pvalue_dic[a_cell] = []
                        cell_pvalue_dic1[a_cell] = []
                    else:
                        continue
                    tag = tag + 1
                    #print(tag)
                    if tag > 1000000:
                        break

                    k_alt_lit = []
                    k_ref_lit = []
                    for l in var_lit:
                        if c.split("_")[1] == l.split("_")[2]:
                            


                            ppsum = 1
                            #for a_site in locals()[str(l)]:
                            #    k_alt = np.mean(matrix_alt.loc[locals()[str(c)],a_site])
                            #    k_ref = np.mean(matrix_ref.loc[locals()[str(c)],a_site])
                            #    k_alt_lit.append(k_alt)
                            #    k_ref_lit.append(k_ref)
                            #k_alt = np.mean(k_alt_lit)  
                            #k_ref = np.mean(k_ref_lit) 

                            #................................................. 
                            ll = list(locals()[str(l)])

                            #site_fine = np.sum(matrix_fre.loc[a_cell,ll]>0.05)
                            #if site_fine/float(len(ll)) < 0.01:
                                #unassign_cell.append(a_cell)
                            #................................................. 
                            depth_percell = np.mean(matrix_depth.loc[a_cell,:])
                            if depth_percell<float(depth_cutoff):
                                unassign_cell.append(a_cell)


                            c_alt = np.sum(matrix_alt.loc[a_cell,ll])
                            c_ref = np.sum(matrix_ref.loc[a_cell,ll])
                            if c_alt > c_ref:

                                #print(ll)
                                k_alt = np.mean(np.sum(matrix_depth.loc[locals()[str(c)],ll],axis = 1))
                                #k_ref = np.mean(np.sum(matrix_ref.loc[locals()[str(c)],locals()[str(l)]],axis = 1))
                                k_ref = 0
                            else:
                                k_alt = 0
                                k_ref = 10


                            alpha_post = k_alt+c_alt
                            beta_post = k_ref+c_ref
                            #print("############",c,l,"####################")
                            #print(k_alt,k_ref,c_alt,c_ref)
                            post_mean = alpha_post/(alpha_post + beta_post)
                            post_mode = (alpha_post-1)/(alpha_post+beta_post-2)
                                                            
                            if post_mode>1:
                                post_mode = 1
                            if post_mode<0:
                                post_mode = 0
                                #print(post_mean)
                            #ppsum = ppsum * abs(post_mean)
                            cell_pvalue_dic[a_cell].append(post_mean)
                            cell_pvalue_dic1[a_cell].append(post_mode)
                        else:
                            ppsum = 1
                            ll = list(locals()[str(l)])
                            c_alt = np.sum(matrix_alt.loc[a_cell,ll])
                            c_ref = np.sum(matrix_ref.loc[a_cell,ll])
                            if c_alt > c_ref:
                                k_alt = 10
                                k_ref = 0
                            #for a_site in locals()[str(l)]:
                            #    k_alt = np.mean(matrix_alt.loc[locals()[str(c)],a_site])
                            #    k_ref = np.mean(matrix_ref.loc[locals()[str(c)],a_site])
                            #    k_alt_lit.append(k_alt)
                            #    k_ref_lit.append(k_ref)
                            #k_alt = np.mean(k_alt_lit)  
                            #k_ref = np.mean(k_ref_lit)  
                            #k_alt = np.mean(np.sum(matrix_alt.loc[locals()[str(c)],locals()[str(l)]],axis = 1))
                            else:
                                k_alt = 0
                                k_ref = np.mean(np.sum(matrix_depth.loc[locals()[str(c)],ll],axis = 1))


                            alpha_post = k_alt+c_alt
                            beta_post = k_ref+c_ref
                            #print(k_alt,k_ref,c_alt,c_ref)
                            '''
                            if (alpha_post + beta_post) == 0:
                                #print(a_cell,ll)
                                post_mean = 0
                            if (alpha_post+beta_post-2) <0:
                                post_mode = 0
                            else:
                            '''
                            post_mean = alpha_post/(alpha_post + beta_post)
                            post_mode = (alpha_post-1)/(alpha_post+beta_post-2)
                                
                            if post_mode>1:
                                post_mode = 1
                            if post_mode<0:
                                post_mode = 0
                                
                                #print(post_mean)
                            #ppsum = ppsum * abs(post_mean)
                            cell_pvalue_dic[a_cell].append(post_mean)
                            cell_pvalue_dic1[a_cell].append(post_mode)
            #return(cell_pvalue_dic)
        #test_dic = caculate_p_each_k(myVar1,myVar2,matrix_alt,matrix_ref)    

        
        lit = []
        raw_lit1 = []
        raw_lit2 = []
        bar_lit = [] 
        
        wrong_cell = []
        
        lit1 = []
        lit2 = []
        lit_doublet1 = []
        lit_doublet2 = []
        lit_con = []
        index_lit1 = []
        index_lit2 = []
        psort_dic = {}
        min_sample = 0
        cell_pvalue_dic = cell_pvalue_dic1
        bar = []
        for i in cell_pvalue_dic.keys():

            #k = sorted(cell_pvalue_dic[i],reverse=True)
            indices = sorted(range(len(cell_pvalue_dic[i])),key=lambda index: cell_pvalue_dic[i][index],reverse=True)
            k = [cell_pvalue_dic[i][index] for index in indices]
            #print(k)
            raw_lit1.append(k[0])
            raw_lit2.append(k[1])
            
            bar_lit.append(i)
            psort_dic[i] = k
            '''
            if k[0]>0.99 and (k[0]+k[1])<1.1:
                min_sample = min_sample + 1
            if k[0]>0.95 and k[1]>0.05:
                wrong_cell.append(i)


            if i not in ture_doublet:
                bar.append(i)

                lit1.append(k[0])
                lit2.append(k[1])
            else:
                index_lit1.append(indices[0])
                index_lit2.append(indices[1])

                lit_doublet1.append(k[0])
                lit_doublet2.append(k[1]) 

            '''
        min_sample_ratio = min_sample/len(cell_pvalue_dic.keys())

        ## define doublet
        ##use knn

        data_P = list(zip(raw_lit1, raw_lit2))
        df3 = pd.DataFrame(np.array([raw_lit1,raw_lit2]).T,columns=['a', 'b'],index=bar_lit)
        #print(data_P)
        
        known_class = []
        new_lit1 = []
        new_lit2 = []
        no_label_lit1 = []
        no_label_lit2 = []
        known_bar = []
        prediction_bar = []
        for i in range(len(raw_lit1)):
            if raw_lit1[i] > 0.99 and raw_lit2[i]<0.01: ## give confident cell "S" label
                new_lit1.append(raw_lit1[i])
                new_lit2.append(raw_lit2[i])
                known_class.append("S")
                known_bar.append(bar_lit[i])

            elif raw_lit1[i] < D_x and raw_lit2[i] > D_y:## give doublet cell "D" label
                new_lit1.append(raw_lit1[i])
                new_lit2.append(raw_lit2[i])
                known_class.append("D")
                known_bar.append(bar_lit[i])
            else: 
                no_label_lit1.append(raw_lit1[i])
                no_label_lit2.append(raw_lit2[i])
                prediction_bar.append(bar_lit[i])
        label_data_P = list(zip(new_lit1, new_lit2))
        nolabel_data_P = list(zip(no_label_lit1, no_label_lit2))

        knn = KNeighborsClassifier(n_neighbors=5) ## n = 5
        knn.fit(label_data_P, known_class)    
        if nolabel_data_P:
            prediction = knn.predict(nolabel_data_P)
        else:
            prediction = []
        #print(prediction)
        all_class = list(known_class)+list(prediction)
        all_bar = list(known_bar)+list(prediction_bar)

        
        # Generate scatter plot for training data

        #colors = list(map(lambda x: '#3b4cc0' if x == "S" else '#b40426', list(known_class)+list(prediction)))
        kp = list(known_class)+list(prediction)
        lp = label_data_P+nolabel_data_P
        s_lit = []
        d_lit = []
        for i in range(len(kp)):
            if kp[i] == 'S':
                s_lit.append(lp[i])
            if kp[i] == 'D':
                d_lit.append(lp[i])

        #print("Singlet:{}\nDoublet:{}".format,len(s_lit),len(d_lit))
        plt.figure(figsize=(10,10))
        plt.scatter(np.array(s_lit[1:500])[:,0], np.array(s_lit[1:500])[:,1], c="#3b4cc0", marker="o", picker=True)
        if len(d_lit) > 0:
            plt.scatter(np.array(d_lit)[:,0], np.array(d_lit)[:,1], c="#b40426", marker="o", picker=True)

        plt.xlim(0,1.1,0.01)
        plt.title("first and second fitted p_value for each cell")
        plt.xlabel('P_value_1')
        plt.ylabel('P_value_2')
        plt.legend(["Singlet", "Doublet"], loc='upper left', markerscale=2, fontsize="large")
        plt.axvline(x = D_x, ymin = 0.0, ymax = 1.1, color = 'r',ls='--')
        plt.axhline(y = D_y, xmin = 0.0, xmax = 1.1, color = 'r',ls='--')

        plt.savefig("{}/p_value.pdf".format(output_path),dpi = 400,bbox_inches = 'tight') 
        plt.show()
        all_p = (label_data_P+nolabel_data_P)
        
        ## find maximum cutoff for p1
        singlet_blue1= []
        doublet_red1 = []

        c = 0
        clean_barcode = []
        for i in range(len(all_class)):
            if all_class[i] == "D" and  df3.loc[all_bar[i],'a'] < 0.99:
                clean_barcode.append(all_bar[i])
                c += 1
                #if clean_matrix.index[i] in true_doublet:
                    #print(clean_matrix.index[i],litss_lit[i])
        
        if len(clean_barcode) > 10:
            aa = sns.clustermap(clean_matrix.loc[clean_barcode,:], metric="euclidean")
            aa.ax_heatmap.set_xlabel("Site")
            aa.savefig("{}/Doublet_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') 

        ## final cluster
        print("Demultiplex singlet!")
        matrix_fre = matrix_fre_pre

        final_barcode = list(set(set(barcode_lit)-set(clean_barcode)))
        #final_barcode = list(low_barcode)
        final_barcode_sorted = sorted(final_barcode)
        matrix_final = matrix_fre.loc[set(final_barcode),useful_site].reindex(index = final_barcode_sorted)

        #sns.clustermap(matrix_final.loc[low_barcode,useful_site], metric="euclidean",figsize=(45, 45),yticklabels=False,col_cluster=True,row_cluster=True)

        if len(final_barcode)>2000:
            ac = sns.clustermap(clean_matrix.loc[final_barcode[:2000],:], metric="euclidean")
            ac.ax_heatmap.set_xlabel("Site")
            ac.savefig("{}/Singlet_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') 
        else:
            ac = sns.clustermap(clean_matrix.loc[final_barcode,:], metric="euclidean")
            ac.ax_heatmap.set_xlabel("Site")
            ac.savefig("{}/Singlet_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') 

        X = matrix_final
        print("matrix_final:",matrix_final.shape)

        ##define cluster by kmeans.
        kmeans = KMeans(n_clusters = sample_num,max_iter=500,tol=1e-10,algorithm='full',n_init=1000)
        kmeans.fit(X)
        labels = kmeans.labels_
        kmeans.labels_
        
        final_labels = labels
       

        finalVar1 = [] 
        for i in range(sample_num):
            createVar['end_cell_'+ str(i)] = []
            finalVar1.append('end_cell_'+ str(i)) 

        for i in range(len(final_labels)):
            for j in range(sample_num):
                if labels[i]==j:
                    locals()['end_cell_'+ str(j)].append(final_barcode_sorted[i])

        finalVar2 = [] 
        for i in range(sample_num):
            createVar['end_cell_'+ str(i)+"_site"] = []
            finalVar2.append('end_cell_'+ str(i)+"_site") 

        for k in range(sample_num):
            jkl = np.sum(matrix_final.loc[locals()["end_cell_"+str(k)],:][matrix_final.loc[locals()["end_cell_"+str(k)],:]>=0.99],axis = 0)/(len(locals()["end_cell_"+str(k)]))
            for i in range(len(matrix_final.columns)):
                if jkl[i]>0.1:
                    locals()["end_cell_"+str(k)+"_site"].append(matrix_final.columns[i])


        ## get mutipul p valuse index from two sample specific germline set: 1.myVar3  2.finalVar2
        def calculate_similarity(list1, list2):
            similarity_count = len(set(list1) & set(list2))
            similarity_ratio = similarity_count / len(list1)
            return similarity_ratio

        new_index_of_germline_set = []

        for i in range(len(myVar3)):
            for l in range(len(finalVar2)):
                similarity_ratio = calculate_similarity(myVar3, finalVar2)
                if sorted(locals()[str(myVar3[i])]) == sorted(locals()[str(finalVar2[l])]) or similarity_ratio > 0.6:
                    new_index_of_germline_set.append(l)


         ## change sample of cell and need to recheck the cluster wiht p value?
        final_dict = {}
        for i in all_bar:

            if i in unassign_cell:
                final_dict[i] = "Unassign"
            if i in clean_barcode and i not in unassign_cell:
                final_dict[i] = "Doublet"

            if i in final_barcode_sorted and i not in unassign_cell:
                p_value = np.array(cell_pvalue_dic[i])
                all_p_value_list = list(p_value[new_index_of_germline_set])##get p_value and reorder by new index.
                max_p_index = all_p_value_list.index(max(all_p_value_list))

                #sample_of_cell = final_labels[final_barcode_sorted.index(i)]
                #changed_cell_sample = new_index_of_germline_set[sample_of_cell]
                final_dict[i] = "Sample"+str(max_p_index)


            


        specific_germline = {}
        for k in range(sample_num):
            k2 = new_index_of_germline_set[k]
            specific_germline["Sample"+str(k2)] = locals()["end_cell_"+str(k2)+"_site"]
            
        return(cell_pvalue_dic,new_index_of_germline_set,specific_germline,all_p,all_class,all_bar,clean_matrix,useful_site,var_lit,wrong_cell,clean_barcode,matrix_fre,final_dict,matrix_final)



    ##--------------##
    ##run_pip_simple##
    ##              ##
    ##              ##
    ##              ##
    ##--------------##
    def run_pip_simple(matrix_alt,matrix_ref,matrix_fre,matrix_depth,sample_num,confident_germline_ratio,D_x,D_y):
        print("Filter variants!")
        matrix_fre_pre = matrix_fre #copy
        matrix_fre  = pd.DataFrame(matrix_fre)

        barcode_lit =  matrix_fre.index

        ##methods2
        a = np.sum(matrix_fre.loc[barcode_lit,:]>0.95,axis=0)/np.sum(matrix_fre.loc[barcode_lit,:]>0.05,axis=0)##>0.05
        b = np.sum(matrix_fre.loc[barcode_lit,:]>0.95,axis=0)/np.sum(matrix_fre.loc[barcode_lit,:]>=0.0,axis=0)##>0.05
        matrix_alt = matrix_alt.where(matrix_alt.notnull(), 0)
        seletec_germ = []
        for i in a.index:
            #print(i,a[i])
            
            if a[i] > confident_germline_ratio and b[i] < (2.0/sample_num):
                seletec_germ.append(i)
        #generate final matrix.

        To_GMM_barcode = barcode_lit

        To_GMM_matrix = matrix_fre.loc[barcode_lit,seletec_germ]


        #final_matrix to GMM or kmean.

        X = np.array(To_GMM_matrix) #filtered doublet & low_depth

        ##define cluster by kmeans.
        sample_num = sample_num
        print('sample_num:',sample_num)

        sil = []
        kmax = 12

        # dissimilarity would not be defined for a single cluster, thus, minimum number of clusters should be 2
        for k in range(2, kmax+1):
            kmeans = KMeans(n_clusters = k).fit(X)
            labels = kmeans.labels_
            sil.append(silhouette_score(X, labels, metric = 'euclidean'))
        fig,ax = plt.subplots(figsize = (8,6))
        ax.set_facecolor("white")
        ax.grid(ls="")
         # Convert bottom-left and top-right to display coordinates
        x_0, y_0 = ax.transAxes.transform((0, 0))
        x_1, y_1 = ax.transAxes.transform((1, 1))

        # Convert back to Axes coordinates
        x_0, y_0 = ax.transAxes.inverted().transform((x_0, y_0))
        x_1, y_1 = ax.transAxes.inverted().transform((x_1, y_1))
        width = 1
        rect = plt.Rectangle(
            (x_0, y_0), x_1-x_0, y_1-y_0,
            color="black",
            transform=ax.transAxes,
            zorder=-1,
            lw=2*width+0.02,
            fill=None,
        )
        ax.patches.append(rect)
        sil_max_index = sil.index(max(sil))
        if sil[sil_max_index] - sil[sil_max_index-1] < 0.0001:
            
            #print("sil:",sil[sil_max_index],sil[sil_max_index-1])
            sil_max_index = sil_max_index-1
            
        plt.axvline(x = sil_max_index+2, ymin = 0.05, ymax = 1.0, color = 'r',ls='--',label = "Best k="+str(sil_max_index+2))
        
        plt.legend(loc = "best",fontsize=15)
        plt.plot(range(2,kmax+1,1), sil, '-p', color='gold')
        plt.xticks(range(1,kmax+1,1))
        plt.xlabel('number of clusters, k')
        plt.ylabel('silhouette_score')
        plt.savefig("{}/sihouette_score.pdf".format(output_path),dpi = 400,bbox_inches = 'tight') # save figure 
        plt.show()




        
        if len(barcode_lit)>2000:
            ab = sns.clustermap(To_GMM_matrix.loc[barcode_lit[1:2000],:], metric="euclidean")
            ab.ax_heatmap.set_xlabel("Site")
            ab.savefig("{}/raw_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') #save figure
        else:
            ab = sns.clustermap(To_GMM_matrix.loc[:,:], metric="euclidean")
            ab.ax_heatmap.set_xlabel("Site")
            ab.savefig("{}/raw_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') #save figure


        if len(barcode_lit)>2000:
            ad = sns.clustermap(To_GMM_matrix.loc[barcode_lit[1:2000],seletec_germ], metric="euclidean")
            ad.ax_heatmap.set_xlabel("Site")
            ad.savefig("{}/clean_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') 
        else:
            ad = sns.clustermap(To_GMM_matrix.loc[:,seletec_germ], metric="euclidean")
            ad.ax_heatmap.set_xlabel("Site")
            ad.savefig("{}/clean_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') 




        kmeans = KMeans(n_clusters = sample_num,max_iter=1000,tol=1e-10,algorithm='full',n_init=100)
        kmeans.fit(X.T)
        labels = kmeans.labels_
        kmeans.labels_



        createVar = locals()
        myVar3 = [] # specific variants for each sample
        for i in range(sample_num):
            createVar['single_cell_'+ str(i)] = []
            myVar3.append('single_cell_'+ str(i)) 
            site = []
            for j in range(len(seletec_germ)):
                if labels[j] == i:
                    site.append(seletec_germ[j]) 
            
            locals()['single_cell_'+ str(i)] = site
            #print( locals()['cell_'+ str(i)][1:5],locals()['single_cell_'+ str(i)])




        ##caculate p-value in first fitted cluster and second fitted
        print("caculate p-value!")
        tag = 0
        var_lit = myVar3
        cell_lit = barcode_lit
        #def caculate_p_each_k(cell_lit,var_lit,matrix_alt,matrix_ref,):
        if 1:
            cell_pvalue_dic = {}
            cell_pvalue_dic1 = {}
            unassign_cell  = []
            for c in cell_lit:
                a_cell = c
                cell_pvalue_dic[a_cell] = []
                cell_pvalue_dic1[a_cell] = []
                if 1:

                    k_alt_lit = []
                    k_ref_lit = []
                    for l in var_lit:
                        if 1:
                            ppsum = 1
                            #for a_site in locals()[str(l)]:
                            #    k_alt = np.mean(matrix_alt.loc[locals()[str(c)],a_site])
                            #    k_ref = np.mean(matrix_ref.loc[locals()[str(c)],a_site])
                            #    k_alt_lit.append(k_alt)
                            #    k_ref_lit.append(k_ref)
                            #k_alt = np.mean(k_alt_lit)  
                            #k_ref = np.mean(k_ref_lit)  
                            ll = list(locals()[str(l)])


                            c_alt = np.sum(matrix_alt.loc[a_cell,ll])
                            c_ref = np.sum(matrix_ref.loc[a_cell,ll])
                            if c_alt > c_ref:
                            #print(ll)
                                k_alt = 10
                            #k_ref = np.mean(np.sum(matrix_ref.loc[locals()[str(c)],locals()[str(l)]],axis = 1))
                                k_ref = 0
                            else:
                                k_alt = 0
                            #k_ref = np.mean(np.sum(matrix_ref.loc[locals()[str(c)],locals()[str(l)]],axis = 1))
                                k_ref = 10


                            alpha_post = k_alt+c_alt
                            beta_post = k_ref+c_ref

                            post_mean = alpha_post/(alpha_post + beta_post)
                            post_mode = (alpha_post-1)/(alpha_post+beta_post-2)
                            if post_mode>1:
                                post_mode = 1
                            if post_mode<0:
                                post_mode = 0

                                #print(post_mean)
                            #ppsum = ppsum * abs(post_mean)
                            cell_pvalue_dic[a_cell].append(post_mean)
                            cell_pvalue_dic1[a_cell].append(post_mode)

        #print(cell_pvalue_dic)
        lit = []
        raw_lit1 = []
        raw_lit2 = []
        bar_lit = [] 

        wrong_cell = []

        lit1 = []
        lit2 = []
        lit_doublet1 = []
        lit_doublet2 = []
        lit_con = []
        index_lit1 = []
        index_lit2 = []
        psort_dic = {}
        min_sample = 0
        cell_pvalue_dic = cell_pvalue_dic1
        bar = []
        for i in cell_pvalue_dic.keys():

            #k = sorted(cell_pvalue_dic[i],reverse=True)
            indices = sorted(range(len(cell_pvalue_dic[i])),key=lambda index: cell_pvalue_dic[i][index],reverse=True)
            k = [cell_pvalue_dic[i][index] for index in indices]
            #print(k)
            raw_lit1.append(k[0])
            raw_lit2.append(k[1])

            bar_lit.append(i)
            psort_dic[i] = k
            '''
            if k[0]>0.99 and (k[0]+k[1])<1.1:
                min_sample = min_sample + 1
            if k[0]>0.95 and k[1]>0.05:
                wrong_cell.append(i)


            if i not in ture_doublet:
                bar.append(i)

                lit1.append(k[0])
                lit2.append(k[1])
            else:
                index_lit1.append(indices[0])
                index_lit2.append(indices[1])

                lit_doublet1.append(k[0])
                lit_doublet2.append(k[1]) 

            '''
        min_sample_ratio = min_sample/len(cell_pvalue_dic.keys())

        ## define doublet
        ##use knn

        data_P = list(zip(raw_lit1, raw_lit2))
        df3 = pd.DataFrame(np.array([raw_lit1,raw_lit2]).T,columns=['a', 'b'],index=bar_lit)
        #print(data_P)

        known_class = []
        new_lit1 = []
        new_lit2 = []
        no_label_lit1 = []
        no_label_lit2 = []
        known_bar = []
        prediction_bar = []
        for i in range(len(raw_lit1)):
            if raw_lit1[i] > 0.99 and raw_lit2[i]<0.01: ## give confident cell "S" label
                new_lit1.append(raw_lit1[i])
                new_lit2.append(raw_lit2[i])
                known_class.append("S")
                known_bar.append(bar_lit[i])

            elif (raw_lit1[i] < D_x and raw_lit2[i] > D_y) :## give doublet cell "D" label
                new_lit1.append(raw_lit1[i])
                new_lit2.append(raw_lit2[i])
                known_class.append("D")
                known_bar.append(bar_lit[i])
            else: 
                no_label_lit1.append(raw_lit1[i])
                no_label_lit2.append(raw_lit2[i])
                prediction_bar.append(bar_lit[i])
        label_data_P = list(zip(new_lit1, new_lit2))
        nolabel_data_P = list(zip(no_label_lit1, no_label_lit2))

        knn = KNeighborsClassifier(n_neighbors=5) ## n = 5
        knn.fit(label_data_P, known_class)    
        if nolabel_data_P:
            prediction = knn.predict(nolabel_data_P)
        else:
            prediction = []
        #print(prediction)
        all_class = list(known_class)+list(prediction)
        all_bar = list(known_bar)+list(prediction_bar)


        # Generate scatter plot for training data

        #colors = list(map(lambda x: '#3b4cc0' if x == "S" else '#b40426', list(known_class)+list(prediction)))
        kp = list(known_class)+list(prediction)
        lp = label_data_P+nolabel_data_P
        s_lit = []
        d_lit = []
        for i in range(len(kp)):
            if kp[i] == 'S':
                s_lit.append(lp[i])
            if kp[i] == 'D':
                d_lit.append(lp[i])

        #print("Singlet:{}\nDoublet:{}".format,len(s_lit),len(d_lit))
        plt.figure(figsize=(10,10))
        plt.scatter(np.array(s_lit[1:5000])[:,0], np.array(s_lit[1:5000])[:,1], c="#3b4cc0", marker="o", picker=True)
        if len(d_lit) > 0:
            plt.scatter(np.array(d_lit)[:,0], np.array(d_lit)[:,1], c="#b40426", marker="o", picker=True)

        plt.xlim(0,1.1,0.01)
        plt.title("first and second fitted p_value for each cell")
        plt.xlabel('P_value_1')
        plt.ylabel('P_value_2')
        plt.legend(["Singlet", "Doublet"], loc='upper left', markerscale=2, fontsize="large")
        plt.axvline(x = D_x, ymin = 0.0, ymax = 1.1, color = 'r',ls='--')
        plt.axhline(y = D_y, xmin = 0.0, xmax = 1.1, color = 'r',ls='--')

        plt.savefig("{}/p_value.pdf".format(output_path),dpi = 400,bbox_inches = 'tight') 
        plt.show()
        all_p = (label_data_P+nolabel_data_P)

        ## find maximum cutoff for p1
        singlet_blue1= []
        doublet_red1 = []

        c = 0
        clean_barcode = []
        for i in range(len(all_class)):
            if all_class[i] == "D" and  df3.loc[all_bar[i],'a'] < 0.99:
                clean_barcode.append(all_bar[i])
                c += 1
                #if clean_matrix.index[i] in true_doublet:
                    #print(clean_matrix.index[i],litss_lit[i])
        final_barcode = list(set(set(all_bar)-set(clean_barcode)))

        if len(clean_barcode) > 10:
            aa = sns.clustermap(To_GMM_matrix.loc[clean_barcode,:], metric="euclidean")
            aa.ax_heatmap.set_xlabel("Site")
            aa.savefig("{}/Doublet_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') 

        if len(final_barcode)>2000:
            ac = sns.clustermap(To_GMM_matrix.loc[final_barcode[:2000],:], metric="euclidean")
            ac.ax_heatmap.set_xlabel("Site")
            ac.savefig("{}/Singlet_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') 
        else:
            ac = sns.clustermap(To_GMM_matrix.loc[final_barcode,:], metric="euclidean")
            ac.ax_heatmap.set_xlabel("Site")
            ac.savefig("{}/Singlet_heatmap.png".format(output_path),dpi = 400,bbox_inches = 'tight') 

        ## change sample of cell and need to recheck the cluster wiht p value?
        final_dict = {}
        for i in all_bar:

            if i in clean_barcode:
                final_dict[i] = "Doublet"

            if i in barcode_lit and i not in clean_barcode:
                p_value = np.array(cell_pvalue_dic[i])
                all_p_value_list = list(p_value)##get p_value and reorder by new index.
                max_p_index = all_p_value_list.index(max(all_p_value_list))

                #sample_of_cell = final_labels[final_barcode_sorted.index(i)]
                #changed_cell_sample = new_index_of_germline_set[sample_of_cell]
                final_dict[i] = "Sample"+str(max_p_index)


            


        specific_germline = {}
        for k in range(sample_num):
            specific_germline["Sample"+str(k)] = locals()["single_cell_"+str(k)]

        return(cell_pvalue_dic,specific_germline,all_p,all_class,all_bar,clean_barcode,matrix_fre,final_dict)









    def generate_html(output_path):
        # Output directory where the PNG and PDF files are located
        output_dir = output_path

        # List all PNG and PDF files in the output directory
        Discription_png = ["raw_heatmap.png: Initial germline muations.", "clean_heatmap.png: Only sample specific germline muations.", "Singlet_heatmap.png: Only Singlet cells.","Doublet_heatmap.png: Only Doublet cells."]
        Discription_pdf = ['sihouette_score.pdf: The optimal classification',"p_value.pdf: p-value for each cells. The cells in the red box in the upper left section are labeled as the original doublet."]
        png_files = ['raw_heatmap.png', 'clean_heatmap.png', 'Singlet_heatmap.png','Doublet_heatmap.png']
        pdf_files = ['sihouette_score.pdf', 'p_value.pdf']

        # Read the data from the text file
        with open(output_dir+'/result_pvalue.txt', 'r') as f:
            data = [line.strip().split('\t') for line in f]

        # Define a template for the HTML file



        html_template = '''
        <html>
        <head>
        <title>My Figures and Data</title>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.25/css/jquery.dataTables.min.css"/>

        <style>
            .image-container {
                display: inline-block;
                width: 45%;
                margin-right: 2%;
                margin-bottom: 10px;
            }
            .image-container:last-child {
                margin-right: 0;
            }
            .image-title {
                font-size: 14px;
                font-weight: bold;
                margin-bottom: 5px;
            }
        </style>

        </head>
        <body>
        <h1>Output figures</h1>

        {% for i in range(1) %}

        <div class="image-container">
                <div class="image-title">{{ Discription_png[i] }}</div>
                <img src="{{ png_files[i] }}" width="100%">
            </div>

        <div class="image-container">
                <div class="image-title">{{ Discription_png[i+1] }}</div>
                <img src="{{ png_files[i+1] }}" width="100%">
            </div>


        <div class="image-container">
                <div class="image-title">{{ Discription_png[i+2] }}</div>
                <img src="{{ png_files[i+2] }}" width="100%">
            </div>

        <div class="image-container">
                <div class="image-title">{{ Discription_png[i+3] }}</div>
                <img src="{{ png_files[i+3] }}" width="100%">
            </div>
            


        {% endfor %}

        {% for i in range(1) %}
            <h2> {{ Discription_pdf[i]}}</h2>
            <iframe src="{{ pdf_files[i] }}" width="90%" height="600px"></iframe>
            <h2>{{Discription_pdf[i+1]}}</h2>
            <iframe src="{{ pdf_files[i+1] }}" width="90%" height="600px"></iframe>
        {% endfor %}


        <h2>p-value for each cell</h2>
        <table id="my-table">
        <thead>
        <tr>
        {% for header in data[0] %}
            <th>{{ header }}</th>
        {% endfor %}
        </tr>
        </thead>
        <tbody>
        {% for row in data[1:] %}
            <tr>
            {% for cell in row %}
                <td>{{ cell }}</td>
            {% endfor %}
            </tr>
        {% endfor %}
        </tbody>
        </table>
        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script src="https://cdn.datatables.net/1.10.25/js/jquery.dataTables.min.js"></script>
        <script>
        $(document).ready(function() {
            $('#my-table').DataTable();
        } );
        </script>





        </body>
        </html>
        ''' 

        # Create a Jinja2 template object
        template = Template(html_template)

        # Render the template with the file lists and data
        html_content = template.render(png_files=png_files, pdf_files=pdf_files, data=data,Discription_png=Discription_png,Discription_pdf=Discription_pdf)

        # Write the HTML file to disk
        with open(output_dir+'/my_figures_and_data.html', 'w') as f:
            f.write(html_content)
        print("Generate html!")



######################-------------------------------##########################
    print("["+record_time()+"] Start demultiplexing")
    (matrix_alt,matrix_ref,matrix_fre,matrix_depth,fragment_data) = read_file(path1,path2,path3,path4)

    if method=="full":
        (cell_pvalue_dic,new_index_of_germline_set,specific_germline,all_p,all_class,all_bar,clean_matrix,useful_site,var_lit,wrong_cell,clean_barcode,matrix_fre,final_dict,matrix_final) = run_pip(matrix_alt,matrix_ref,matrix_fre,matrix_depth,sample_num,confident_germline_ratio,D_x,D_y)

        def rechange_p_value_order(cell):
            p_value = np.array(cell_pvalue_dic[cell])
            all_p_value_list = list(p_value[new_index_of_germline_set])
            all_p_value_list
            return(all_p_value_list)






        write_first_line = "Barcode\tDemultiplex\tP_value_1\tP_value_2"
        for i in range(sample_num):
            write_first_line = write_first_line+"\t"+"Sample"+str(i)+"_Pr"
        write_first_line = write_first_line + "\n"  

        with open("{}/result_pvalue.txt".format(demultiplex_dir),"w") as ff:
            ff.write(write_first_line)
            for i in range(len(all_bar)):

                if all_p[i][0] > 1:
                    p1 = 1.0
                else:
                    p1 = all_p[i][0]
                if all_p[i][1] > 1:
                    p2 = 1.0
                else:
                    p2 = all_p[i][1]


                all_p_value_list = rechange_p_value_order(all_bar[i])
                ff.write("{}\t{}\t{}\t{}".format(all_bar[i],final_dict[all_bar[i]],p1,p2))#,all_class[i]
                for i in range(sample_num):
                    sp = all_p_value_list[i]
                    if sp > 1:
                        sp = 1.0
                    ff.write("\t{}".format(sp))
                ff.write("\n")
        ff.close()



        with open("{}/specific_germline.txt".format(demultiplex_dir),"w") as ffff:
            ffff.write("Specific germline\tSample\n")
            for i in specific_germline.keys():
                for j in specific_germline[i]:
                    ffff.write("{}\t{}\n".format(j,i))
        ffff.close()

        generate_html(demultiplex_dir)
        print("["+record_time()+"] Demultiplexing done!")
        end = time.time()
        elapsed = end-start
        print("Total time:",str(time.strftime("%Hh%Mm%Ss", time.gmtime(elapsed))))





######################-------------------------------##########################
    else:#run_pip_simple
        (cell_pvalue_dic,specific_germline,all_p,all_class,all_bar,clean_barcode,matrix_fre,final_dict) = run_pip_simple(matrix_alt,matrix_ref,matrix_fre,matrix_depth,sample_num,confident_germline_ratio,D_x,D_y)

        write_first_line = "Barcode\tDemultiplex\tP_value_1\tP_value_2"
        for i in range(sample_num):
            write_first_line = write_first_line+"\t"+"Sample"+str(i)+"_Pr"
        write_first_line = write_first_line + "\n"  

        with open("{}/result_pvalue.txt".format(demultiplex_dir),"w") as ff:
            ff.write(write_first_line)
            for i in range(len(all_bar)):

                if all_p[i][0] > 1:
                    p1 = 1.0
                else:
                    p1 = all_p[i][0]
                if all_p[i][1] > 1:
                    p2 = 1.0
                else:
                    p2 = all_p[i][1]


                all_p_value_list = cell_pvalue_dic[all_bar[i]]
                ff.write("{}\t{}\t{}\t{}".format(all_bar[i],final_dict[all_bar[i]],p1,p2))#,all_class[i]
                for i in range(sample_num):
                    sp = all_p_value_list[i]
                    if sp > 1:
                        sp = 1.0
                    if sp < 0:
                        sp = 0.0
                    ff.write("\t{}".format(sp))
                ff.write("\n")
        ff.close()

        with open("{}/specific_germline.txt".format(demultiplex_dir),"w") as ffff:
            ffff.write("Specific germline\tSample\n")
            for i in specific_germline.keys():
                for j in specific_germline[i]:
                    ffff.write("{}\t{}\n".format(j,i))
        ffff.close()


        generate_html(demultiplex_dir)
        print("["+record_time()+"] Demultiplexing done!")
        end = time.time()
        elapsed = end-start
        print("Total time:",str(time.strftime("%Hh%Mm%Ss", time.gmtime(elapsed))))



if __name__ == '__main__':
    cli()

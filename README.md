# MitoSort
Preprint manuscript of this method available at XXX
## Overview 
- MitoSort is an efficient computational method to demultiplex samples from different individuals and detect cross-genotype doublets using endogenous mtDNA germline variants. 
- MitoSort is comprised of 6 steps with the first 2 using external tools and other using in-house script. Users can run each step one by one:
	1. Realign MT sequences (using GATK, 01_MT_realign.py)
	2. Find SNP (using Varscan2, 02_find_SNPs.py)
	3. Divide bam (03_divide_BAM.py)
	4. Retain cell barcode (04_retain_cell_barcodes.py)
	5. Generate SNP matrix(05_generate_SNP_matrix.py)
	6. Demultiplex(06_demultiplex.py)
- For the convenience of users, we also encapsulate MitoSort pipeline into three subcommands including mt-realign, generate-snp-matrix and demultiplex. Users can input the name of subcommand for further help information : 
```shell
Usage: MitoSort_pipeline.py [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  demultiplex          Identifying cross-genotype doublets and demultiplexing samples
  generate-snp-matrix  Generate SNP matrices
  mt-realign           Realign mitochondrial reads using GATK
```

## Installation
- Users can create a conda environment to install all the required python packages of MitoSort
```shell
# clone MitoSort repository 
git clone https://github.com/tangzhj/MitoSort.git

# create a conda environment and install all the required python packages
conda env create -f /path/to/MitoSort/MitoSort_env.yaml
conda activate MitoSort
```
- In addition to required python packages. MitoSort also requires [GATK](https://github.com/broadgsa/gatk/releases) for MT realignment and [VarScan2](https://github.com/Jeltje/varscan2/releases) for variant calling. Users should install them and specify the path of tool when running MitoSort. 


## Usage 
### Realign mitochondrial reads using GATK
- First, mitochondrial reads are realigned by GATK by calling `mt-realign` command. The following are options for `mt-realign` command. It takes a BAM file as input and output a BAM file containing realigned MT reads (possorted_chrM_realign.bam).
```shell
Usage: MitoSort_pipeline.py mt-realign [OPTIONS]

  Realign mitochondrial reads using GATK

Options:
  -b, --bam_file PATH      BAM file output by Cellranger  [required]
  -f, --genome_fasta PATH  FASTA file of reference genome  [required]
  --gatk_path PATH         Path to GenomeAnalysisTK(GATK)  [required]
  -o, --output_dir PATH    Output parent directory  [required]
  --data_type TEXT         Type of data. Default is ATAC. Set to RNA if your
                           data is scRNA

  -h, --help               Show this message and exit.
```
- A typical `mt-realign` command looks like
```shell
python MitoSort_pipeline.py -b /path/to/possorted_bam.bam -f /path/to/reference.fasta --gatk_path /path/to/GenomeAnalysisTK_3.5-0.jar -o /path/to/output_dir
```

### Generate SNP matrices
- Given a list of barcodes, SNP matrices containg alternative counts, reference counts and allele frequency are generated for these barcodes by calling `generate-snp-matrix` command. The following are options for `generate-snp-matrix` command. Note that the `output_dir` option should be the same as that of `mt-realign` command.
```shell
Usage: MitoSort_pipeline.py generate-snp-matrix [OPTIONS]

  Generate SNP matrices

Options:
  -b, --bam_file PATH      BAM file containing realigned MT reads
                           `possorted_chrM_realign.bam`  [required]

  -f, --genome_fasta PATH  FASTA file of reference genome  [required]
  -c, --cell_barcode PATH  A CSV file containing barcode metrics output by
                           Cellranger
                           (singlecell.csv/per_barcode_metrics.csv). If it
                           isn't generated by Cellranger, please input a file
                           cotaining cell barcodes with header named `barcode`
                           [required]

  -m, --chrm_length PATH   a BED file containing chrM region  [required]
  --varscan_path PATH      Path to VarScan  [required]
  -o, --output_dir PATH    Output parent directory  [required]
  --cell_tag TEXT          Default is CR. Set if your cell barcode tag is not
                           CR

  -h, --help               Show this message and exit.
```
- A typical `generate-snp-matrix` command looks like
```shell
python MitoSort_pipeline.py generate-snp-matrix -b /path/to/possorted_chrM_realign.bam -f /path/to/reference.fasta -c /path/to/singlecell.csv -m /path/to/MitoSort/data/hg38_chrM.bed --varscan_path /path/to/VarScan.v2.3.7.jar -o /path/to/output_dir
```

### Doublet identification and sample demultiplexing
- The core of MitoSort is to identify cross-genotype doublets and demultiplex samples from different individuals. It can be achieved by calling `demultiplex` command. Users can tune the parameters and run it multiple times to achieve best performance on specific datasets. The following are options for `demultiplex` command. Note that the `output_dir` option should be the same as that of `mt-realign` command.
```shell
Usage: MitoSort_pipeline.py demultiplex [OPTIONS]

  Identifying cross-genotype doublets and demultiplexing samples

Options:
  -o, --output_dir PATH  Output parent directory  [required]
  -k, --clusters TEXT    number of pooled individuals  [required]
  --p1_cutoff TEXT       maximum cutoff of p1 for doublet identification.
                         Default to be 0.9

  --p2_cutoff TEXT       minimun cutoff of p2 for doublet identification.
                         Default to be 0.1

  -h, --help             Show this message and exit.
```
- A typical `demultiplex` command looks like
```shell
python MitoSort_pipeline.py demultiplex -o /path/to/output_dir -k number_of_pooled_individuals
```
- The important output files are 
	1. specific_gemeline.txt (specific mtDNA germline variants for each individual)
	2. The result_pvalue.txt  
	- The following is example of result_pvalue.txt, which includes the cell barcode, singlet/doublet status, assigned sample, p1 and p2:
```
barcode             S_or_D    Demultiplex     P_value_1              P_value_2
TGCGTAAAGTCGGGAT-1  S         Sample2         1.0                    0.0011289436655
CCTAAAGCACCATTCC-1  S         Sample2         0.99454523553          0.00941613048339
AGCCAGCTCTGGCCAG-1  S         Sample3         0.991007704872         0.00158759496457
TTCAACTTCGTTGTTT-2  S         Sample0         0.998110251201         0.0017108811068
CATTCATAGTACGCGA-4  S         Sample0         0.998505613353         0.00139568480579
TGCACCTGTATTCGCA-1  D         Doublet         0.722414353581         0.253402082068
GCGATCGAGTACAGAT-1  S         Sample0         0.997696759297         0.00915089076023
CAACCAACAGTTACAC-4  S         Sample0         0.997044832437         0.00447870195437
TACGGATGTGTTTCTT-1  S         Sample1         1.0                    0.00186859425923
GCAGCCATCAGCCGGT-1  S         Sample1         0.996103622578         0.00234697252993
TGTACGAGTATCCTTT-1  S         Sample1         0.993533729465         0.00477584068789

```
- We also generate a HTML file cotaining output figures
- Raw allele frequency matrix
![Image text](https://github.com/tangzhj/MitoSort/blob/main/example_output_figures/raw_heatmap.png)

- Clean allele frequency matrix for all cell barcodes (removing common variants)
![Image text](https://github.com/tangzhj/MitoSort/blob/main/example_output_figures/clean_heatmap.png)

- Scatter plot of p1 and p2 for doubelt identification. Users can modify the cutoff of p1 and p2 based on the plot
![Image text](https://github.com/tangzhj/MitoSort/blob/main/example_output_figures/p_value.png)

- Clean allele frequency matrix for singlets
![Image text](https://github.com/tangzhj/MitoSort/blob/main/example_output_figures/Singlet_heatmap.png)

- Clean allele frequency matrix for doublets
![Image text](https://github.com/tangzhj/MitoSort/blob/main/example_output_figures/Doublet_heatmap.png)

- https://github.com/tangzhj/MitoSort/blob/main/example_output_figures/Doublet_heatmap.png


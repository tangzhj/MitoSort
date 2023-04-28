# MitoSort
Preprint manuscript of this method available at XXX
## Overview 
- MitoSort is an efficient computational method to demultiplex samples from different individuals and detect cross-genotype doublets using endogenous mtDNA germline variants. 
MitoSort is comprised of 6:
1. Realign ([]())
2. Find SNP ([Varscan2]())
3. Divide bam ([generate_csv]())
4. Retain cell barcode ()
5. Generate SNP matrix()
6. Demultiplex()

## Installation
### Easy Installation (Linux) (recommended) 
Wait for setting.
### Hard Installation 
Wait for setting.
### Input and Output
The output files are:
1. specific_gemeline.txt
2. The result_pvalue.txt 
result_pvalue.txt looks like:
```
barcode	S_or_D	Demultiplex	P_value_1	P_value_2
TGCGTAAAGTCGGGAT-1	S	Sample2	1.0 		0.0011289436655
CCTAAAGCACCATTCC-1	S	Sample2	0.99454523553	0.00941613048339
AGCCAGCTCTGGCCAG-1	S	Sample3	0.991007704872	0.00158759496457
TTCAACTTCGTTGTTT-2	S	Sample0	0.998110251201	0.0017108811068
CATTCATAGTACGCGA-4	S	Sample0	0.998505613353	0.00139568480579
TGCACCTGTATTCGCA-1	D	Doublet	0.722414353581	0.253402082068
GCGATCGAGTACAGAT-1	S	Sample0	0.997696759297	0.00915089076023
CAACCAACAGTTACAC-4	S	Sample0	0.997044832437	0.00447870195437
TACGGATGTGTTTCTT-1	S	Sample1	1.0 		0.00186859425923
GCAGCCATCAGCCGGT-1	S	Sample1	0.996103622578	0.00234697252993
TGTACGAGTATCCTTT-1	S	Sample1	0.993533729465	0.00477584068789

```



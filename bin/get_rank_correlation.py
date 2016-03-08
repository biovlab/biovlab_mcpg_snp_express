#!/bin/python

import sys

# parameter check
if len(sys.argv) < 5 :
	print "[USAGE : python get_rank_correlation.py [gene list file] [expression level file list] [methylation level file list] [TSS_info_file]"
	sys.exit(0)

# parameters
gene_list_fd = open(sys.argv[1], "r")
gene_info_fd = open(sys.argv[5], "r")


#variales
gene_list=[]
# exp array
# [filename][chr][starat][end][count]

# sample file list assume ';'seperated
expression_sample_list = sys.argv[2].split(';')
methylation_sample_list = sys.argv[3].split(';')
sample_num = len(expression_smple_list)

# dictionaries
refseq2chr={}
refseq2strand={}
refseq2TSS={}
refseq2TSE={}
refseq2geneSymbol_list={}

# extract gene list
for line in gene_list_fd.readslines():
	gene_list.append(line[:-1])

# extract gene expression file
for exp_file in expression_sample_list:
	exp_fd = open(exp_file, "r")
	for line in exp_fd.readlines():

	close(exp_file)
		

# extract methylation level file

# extract TSS file
# $1:refseqID, $2:chr, $3:strand, $4:TSS, $5:TSE, $11:Gene symbol(possiblly multiple gene symbol)

for line in gene_info_fd.readliens():
	tokens = line[:-1].split('\t')

	# extract data
	refseqID = tokens[0]
	chr_num = tokens[1]
	strand = tokens[2]
	TSS = tokens[3]
	TSE = tokens[4]
	geneSymbol_list = tokens[10].split(',')

	# store to dictionaries
	refseq2chr[refseqID] = chr_num
	refseq2strand[refseqID] = strand
	refseq2TSS[refseqID] = TSS
	refseq2TSE[refseqID] = TSE
	refseq2geneSymbol_list[refseqID] = geneSymbol_list


# For each gene, compute correlation of each bins
for refseq_id in refseq2chr.keys():
	for (i=0; i<sample_num; i++):


# close files
close(gene_list_fd)
close(gene_info_fd)







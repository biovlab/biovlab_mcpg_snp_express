import sys
import Bio.Cluster

# parameter check
if len(sys.argv) < 8 :
	print "[USAGE : python bin/merge_snp.py [snp sample_file_list] [gene expression sample file list] [correlation threshold] [merge out_file] [stat out file] [out matrix file for theshold] [out matrix file for PC NC ratio]"
        sys.exit(0)


# variables
snp_file_list	 		= sys.argv[1].split(';')
exp_file_list 			= sys.argv[2].split(';')
corr_threshold			= float(sys.argv[3])
merge_out_file_fd		= open(sys.argv[4], "w")
stat_out_file_fd		= open(sys.argv[5], "w")
matrix_out_file_threshold_fd = open(sys.argv[6], "w")
matrix_out_file_PC_NC_fd = open(sys.argv[7], "w")
sample_num 			= len(snp_file_list)
data				= {}				# double dictionary data = { sample : value }
exp_data			= {}
common_position_count		= {}
header				= ""

# stat variables
# snp
num_snp_over_threshold		= {}
num_snp_having_positive_corr 	= {}
num_snp_having_negative_corr 	= {}
num_snps_within_ranges	= {}
num_snps_having_corr	 	= {}
# gene 
num_gene_having_snp_over_threshold = 0




# result from above
#chr    bin_start       bin_end ref/alt snp_strand      info    range_start  range_end   strand  gene_symbol	range_kind	refseq
#chr10  108709  108710  a/aAAA  +       INDEL;IS=3,0.600000;DP=5;QS=0.000000,1.000000,0.000000,0.000000;VDB=2.063840e-02;AF1=1;AC1=2;DP4=0,0,3,0;MQ=20;FQ=-43.5 0       345178  -       TUBB8   SNP_TSS_TSE_flanking_range+-250kb       NM_177987

# find common position of snp

# extract snp files
for file_num in range(len(snp_file_list)):
	temp_file = snp_file_list[file_num]
	in_fd = open(temp_file, "r")

	# per file
	for line in in_fd:
		tokens 		= line[:-1].split('\t')

		if (len(tokens) < 11): 
			print "[Error] Input snp sample has less then 11 column! Please check the input file : " + snp_file_list[file_num]
			print "[Error] The input line is : " + line[:-1]
			sys.exit(0)

		if line.startswith('#') : 
			header = line[:-1]
			continue

		chr_num 	= tokens[0]
		start 		= str(tokens[1])
		end		= str(tokens[2])
		ref_alt		= tokens[3]
		af1		= str(tokens[5].split("AF1")[1].split('=')[1].split(';')[0])
		snp_range_start	= str(tokens[6])
		snp_range_end	= str(tokens[7])
		strand		= tokens[8]
		refseq		= tokens[11]
		gene_symbol	= tokens[9]
		range_kind	= tokens[10]

		# store data with double dictionary
		# data={}
		# data[fixed_key][sample/count]={count : value, sample : ref/alt"_"af1}, so with samekey, data has multiple(sample number of) ref/alt"_"af1 values
		fixed_key 	= chr_num+"~"+ start +"~"+ snp_range_start + "~" + snp_range_end + "~" + refseq + "~" + gene_symbol

		# count
		try:
			temp = data[range_kind]
		except KeyError:
			data[range_kind]={}

		try:	
			data[range_kind][fixed_key]['count'] += + 1
		except KeyError:
			data[range_kind][fixed_key]={}
			data[range_kind][fixed_key]['count'] = 1

		# data
		data[range_kind][fixed_key][file_num] = ref_alt + "~" + af1

	in_fd.close()

# extract gene expression file

# gene expression file (.cell file) : # 100730_s_4_export.txt.CEL.exp
# 200000_s_at     PRPF8   9.416984390952
# 200000_s_at     PRPF8   8.416984390952
# use average

for file_num in range(len(exp_file_list)):
	temp_file = exp_file_list[file_num]
	in_exp_fd = open(temp_file, "r")

	# per file
	for line in in_exp_fd:
		tokens = line[:-1].split('\t')

		if (len(tokens) < 2): 
			print "[Error] Input exp sample has less then 3 column! Please check the input file : " + exp_file_list[file_num]
			print "[Error] The input line is : " + line[:-1]
			sys.exit(0)

		gene_symbol	= tokens[0]
		exp_level	= float(tokens[1])

		try:
			exp_data[gene_symbol,file_num]['value'] = (exp_data[gene_symbol,file_num]['value'] * exp_data[gene_symbol,file_num]['count'] + exp_level) / (exp_data[gene_symbol,file_num]['count'] + 1)
			exp_data[gene_symbol,file_num]['count'] = exp_data[gene_symbol,file_num]['count'] + 1
		except KeyError:
			exp_data[gene_symbol,file_num] ={}
			exp_data[gene_symbol,file_num]['count'] = 1
			exp_data[gene_symbol,file_num]['value'] = exp_level

	in_exp_fd.close()

# find snps with count over threshold (that is, common snps), and merge snps

threshold=1

for range_kind in data.keys():
	for key in data[range_kind].keys():

		# print only snps over threshold
		if (data[range_kind][key]['count'] > threshold):

			tokens 		= key.split('~')

			chr_num		= tokens[0]
			start		= tokens[1]
			end		= str(int(start) + 1)
			snp_range_start	= tokens[2]
			snp_range_end	= tokens[3]
			refseq		= tokens[4]
			gene_symbol	= tokens[5]

			# print as original format
			#chr10  108709  108710  a/aAAA  +       INDEL;IS=3,0.600000;DP=5;QS=0.000000,1.000000,0.000000,0.000000;VDB=2.063840e-02;AF1=1;AC1=2;DP4=0,0,3,0;MQ=20;FQ=-43.5 0       345178  -       TUBB8   SNP_TSS_TSE_flanking_range+-250kb       NM_177987
			merge_out_file_fd.write(chr_num + "\t" + start + "\t" + end + "\t")

			# increase snp count
			try:
				num_snps_within_ranges[range_kind] += 1
			except KeyError:
				num_snps_within_ranges[range_kind] = 1

			# print ref/alt values per sample
			for file_num in range(len(snp_file_list)):

				try:
					tokens2 	= data[range_kind][key][file_num].split('~')
					ref_alt		= tokens2[0] 
				except KeyError:
					ref_alt		= "-/-"
					
				if file_num == 0:
					merge_out_file_fd.write(ref_alt)
				else:
					merge_out_file_fd.write("\t" + ref_alt)

			merge_out_file_fd.write("\t")

			# print af1 values per sample
			rank_correlation_x=[]
			for file_num in range(len(snp_file_list)):
				try:
					tokens2 	= data[range_kind][key][file_num].split('~')
					af1		= tokens2[1] 
				except KeyError:
					af1		= 0

				if file_num == 0:
					merge_out_file_fd.write(str(af1))
				else:
					merge_out_file_fd.write("\t" + str(af1))

				# store af1 value for rank correlation
				rank_correlation_x.append(af1)

			merge_out_file_fd.write("\t" + snp_range_start + "\t" + snp_range_end + "\t" + strand + "\t" + gene_symbol + "\t" + range_kind + "\t" + refseq)

			# print out exp level per sample
			rank_correlation_y=[]
			for file_num in range(len(exp_file_list)):
				try:
					exp_level=str(exp_data[gene_symbol,file_num]['value'])
				except KeyError:
					exp_level=0

				merge_out_file_fd.write("\t" + str(exp_level))

				# store exp_level for rank correlation
				rank_correlation_y.append(exp_level)


			# remove '-' element on both
			remove_index=[]
			for i in range(len(rank_correlation_x)):
				if rank_correlation_x[i] == '-':
					remove_index.append(i)
				if rank_correlation_y[i] == '-':
					remove_index.append(i)
			# remove duplicate
			remove_index=list(set(remove_index))

			# sort index to remove the list from backward
			remove_index.sort(reverse=True)

			# remove list by index
			#print rank_correlation_x
			#print rank_correlation_y
			#print remove_index
			for index in remove_index:
				del rank_correlation_x[index]
				del rank_correlation_y[index]

			#print "here"
			#print len(rank_correlation_x)
			#print len(snp_file_list)*0.3

			# check if the list is empty 
			if (len(rank_correlation_x) > 0) and (len(rank_correlation_x) > (len(snp_file_list)*0.3)) :
				# compute spearman'rho 
				rank_corr=1 - Bio.Cluster.distancematrix((rank_correlation_x,rank_correlation_y), dist="s")[1][0]
				merge_out_file_fd.write("\t" + str(rank_corr) + "\n")

				try:
					num_snps_having_corr[range_kind] += 1
				except KeyError:
					num_snps_having_corr[range_kind] = 1

				if rank_corr > corr_threshold or rank_corr < -corr_threshold : 
					try:
						num_snp_over_threshold[range_kind] += 1
					except KeyError:
						num_snp_over_threshold[range_kind] = 1

					if rank_corr > 0 :
						try:
							num_snp_having_positive_corr[range_kind] += 1
						except KeyError:
							num_snp_having_positive_corr[range_kind] = 1

					if rank_corr < 0 :
						try:
							num_snp_having_negative_corr[range_kind] += 1
						except KeyError:
							num_snp_having_negative_corr[range_kind] = 1

			else:
				merge_out_file_fd.write("\t-\n")

matrix_out_file_threshold_fd.write("stat\tvariable\tvalue\n")
matrix_out_file_PC_NC_fd.write("stat\tvariable\tvalue\n")

# print out statistics
for range_kind in num_snps_within_ranges.keys():
	stat_out_file_fd.write("<< Range : " + range_kind + " >>\n")
	try:
		stat_out_file_fd.write("(1). Number of SNPs within " + range_kind + " range :" + str(num_snps_within_ranges[range_kind]) + "\n")
	except KeyError:
		stat_out_file_fd.write("(1). Number of SNPs within " + range_kind + " range : 0\n")
	try:
		stat_out_file_fd.write("(2). Number of SNPs having Spermans's rho with gene expression among (1) :" + str(num_snps_having_corr[range_kind]) + "\n")
	except KeyError:
		stat_out_file_fd.write("(2). Number of SNPs having Spermans's rho with gene expression among (1) : 0\n")
	try:
		stat_out_file_fd.write("(3). Number of SNPs having Spermans's rho(rank correlation) over +- " + str(corr_threshold) + " among (2) :" + str(num_snp_over_threshold[range_kind]) + "\n")
	except KeyError:
		stat_out_file_fd.write("(3). Number of SNPs having Spermans's rho(rank correlation) over +- " + str(corr_threshold) + " among (2) : 0\n")
	try:
		stat_out_file_fd.write("(4). % of SNP over rank correlation threshold +-" + str(corr_threshold) + " (1)/(3) : " + str(float(num_snp_over_threshold[range_kind])/float(num_snps_having_corr[range_kind]) * 100) + "%\n")
	except KeyError:
		stat_out_file_fd.write("(4). % of SNP over rank correlation threshold +-" + str(corr_threshold) + " (1)/(3) : 0%\n")
	try:
		stat_out_file_fd.write("(5). % of SNP having positive correlation among (4) : " + str((float(num_snp_having_positive_corr[range_kind])/float(num_snps_having_corr[range_kind])*100)) + "%\n")
	except KeyError:
		stat_out_file_fd.write("(5). % of SNP having positive correlation among (4) : 0%\n")
	try:
		stat_out_file_fd.write("(6). % of SNP having negative correlation among (4) : " + str((float(num_snp_having_negative_corr[range_kind])/float(num_snps_having_corr[range_kind]))*100) + "%\n\n")
	except KeyError:
		stat_out_file_fd.write("(6). % of SNP having negative correlation among (4) : 0%\n\n")

	try:
		temp_value='%.2f'%(float(num_snp_over_threshold[range_kind])/float(num_snps_having_corr[range_kind]) * 100)
	except KeyError:
		temp_value='%.2f'%(float(0.00))
	temp_value2=100-float(temp_value)
	matrix_out_file_threshold_fd.write("Corr_>_threshold\t"+range_kind+"\t" + str(temp_value)+"\n")
	matrix_out_file_threshold_fd.write("Corr_<_threshold\t"+range_kind+"\t" + str(temp_value2)+"\n")

	try:
		temp_value='%.2f'%(float(num_snp_having_positive_corr[range_kind])/float(num_snp_over_threshold[range_kind])*100)
	except KeyError:
		temp_value='%.2f'%(float(0.00))
	temp_value2=100-float(temp_value)
	matrix_out_file_PC_NC_fd.write("Negative_correlation\t"+range_kind+"\t" + str(temp_value2) + "\n")
	matrix_out_file_PC_NC_fd.write("Positive_correlation\t"+range_kind+"\t" + str(temp_value) + "\n")
# gene
#print "Number of genes having SNP which has correlation over threshold : " + str(num_gene_having_snp_over_threshold)

#print exp_data
		

					

stat_out_file_fd.close()
matrix_out_file_threshold_fd.close()
matrix_out_file_PC_NC_fd.close()
		
		



		






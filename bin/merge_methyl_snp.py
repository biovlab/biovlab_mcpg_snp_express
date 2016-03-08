import sys
import Bio.Cluster

# parameter check
if len(sys.argv) < 7 :
	print "[USAGE : python bin/merge_methyl_snp.py [sample_file_list] [gene expression sample file list] [correlation threshold] [merge out_file] [stat out file] [methyl/snp]"
        sys.exit(0)


# variables
file_list	 		= sys.argv[1].split(';')
exp_file_list 			= sys.argv[2].split(';')
corr_threshold			= float(sys.argv[3])
merge_out_file_fd		= open(sys.argv[4], "w")
stat_out_file_fd		= open(sys.argv[5], "w")
sample_num 			= len(file_list)
data				= {}				# double dictionary data = { sample : value }
exp_data			= {}
common_position_count		= {}
header				= ""
methyl_snp=sys.argv[6]

# stat variables
# snp 
num_count_in_ranges		= {}
num_having_corr_in_ranges	= {}
num_over_threshold_in_ranges	= {}
num_having_positive_corr_in_ranges	= {}
num_having_negative_corr_in_ranges	= {}
num_snps_having_corr	 	= 0
snp_count_within_TSS_range	= 0
num_snp_over_threshold		= 0
num_snp_having_positive_corr 	= 0
num_snp_having_negative_corr 	= 0
num_snps_within_tss_range	= 0
num_snps_having_corr	 	= 0
# gene 
num_gene_having_snp_over_threshold = 0




# result from above
# methyl file1
#chr    bin_start       bin_end methyl_count    methyl_level(rms)       range_start     range_end       strand  refseq  gene_symbol     range_kind
#chr1    10401   10500   29      0.334728794189335       10373   11873   +       NR_046018,      DDX11L1 TSS_range

# methyl file2
#chr    bin_start       bin_end methyl_count    methyl_level(rms)       range_start     range_end       strand  refseq  gene_symbol     range_kind
#chr1    10401   10500   29      0.334728794189335       10373   11873   +       NR_046018,      DDX11L1 TSS_range

# snp
#chr    bin_start       bin_end ref/alt snp_strand      info    range_start   range_end   strand  refseq  gene_symbol	range_kind
#chr1    14929   14930   A/G     +       DP=5;QS=0.000000,1.000000,0.000000,0.000000;VDB=1.165196e-02;AF1=1;AC1=2;DP4=0,0,2,2;MQ=20;FQ=-39       0       320008  +       NM_001005484,   OR4F5	TSS_range
#chr1    14929   14930   A/G     +       DP=5;QS=0.000000,1.000000,0.000000,0.000000;VDB=1.165196e-02;AF1=1;AC1=2;DP4=0,0,2,2;MQ=20;FQ=-39       0       264409  +       NR_046018,      DDX11L1	genebody

# find common position of snp

# extract files
for file_num in range(len(file_list)):
	temp_file = file_list[file_num]
	in_fd = open(temp_file, "r")

	# per file
	for line in in_fd:
		tokens 		= line[:-1].split('\t')

		if (methyl_snp == "snp"):
			if (len(tokens) < 11): 
				print "[Error] Input snp sample has less then 11 column! Please check the input file : " + file_list[file_num]
				print "[Error] The input line is : " + line[:-1]
				sys.exit(0)
		else:
			if (len(tokens) < 10): 
				print "[Error] Input sample has less then 10 column! Please check the input file : " + file_list[file_num]
				print "[Error] The input line is : " + line[:-1]
				sys.exit(0)
			

		if line.startswith('#') : 
			header = line[:-1]
			continue

		chr_num 	= tokens[0]
		start 		= str(tokens[1])
		end		= str(tokens[2])

		if (methyl_snp == "snp"):
			ref_alt		= tokens[3]
			af1_rms		= str(tokens[5].split("AF1")[1].split('=')[1].split(';')[0])
			range_start	= str(tokens[7])
			range_end	= str(tokens[8])
			strand		= tokens[9]
			refseq		= tokens[10]
			gene_symbol	= tokens[11]
			range_kind	= tokens[12]
		else:
			ref_alt		= "-/-"
			af1_rms		= tokens[4]
			range_start	= str(tokens[5])
			range_end	= str(tokens[6])
			strand		= tokens[7]
			refseq		= tokens[8]
			gene_symbol	= tokens[9]
			range_kind	= tokens[10]


		# store data with double dictionary
		# data={}
		# data[fixed_key][sample/count]={count : value, sample : ref/alt"_"af1_rms}, so with samekey, data has multiple(sample number of) ref/alt"_"af1_rms values
		fixed_key 	= chr_num+"~"+ start +"~"+ range_start + "~" + range_end + "~" + refseq + "~" + gene_symbol + "~" + range_kind

		# count
		try:
			data[fixed_key]['count'] = data[fixed_key]['count'] + 1
		except KeyError:
			data[fixed_key]={}
			data[fixed_key]['count'] = 1

		# data
		data[fixed_key][file_num] = ref_alt + "~" + af1_rms

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

		# get average exp value for multiple exp values of single gene
		try:
			exp_data[gene_symbol,file_num]['value'] = (exp_data[gene_symbol,file_num]['value'] * exp_data[gene_symbol,file_num]['count'] + exp_level) / (exp_data[gene_symbol,file_num]['count'] + 1)
			exp_data[gene_symbol,file_num]['count'] = exp_data[gene_symbol,file_num]['count'] + 1
		except KeyError:
			exp_data[gene_symbol,file_num] ={}
			exp_data[gene_symbol,file_num]['count'] = 1
			exp_data[gene_symbol,file_num]['value'] = exp_level

	in_exp_fd.close()


		


# find over threshold (that is, common snps), and merge snps

threshold=1

for key in data.keys():

	# print only snps over threshold
	if (data[key]['count'] > threshold):

		tokens 		= key.split('~')

		chr_num		= tokens[0]
		start		= tokens[1]
		end		= str(int(start) + 1)
		range_start	= tokens[2]
		range_end	= tokens[3]
		refseq		= tokens[4]
		gene_symbol	= tokens[5]
		range_kind	= tokens[6]

		# print as original format
		merge_out_file_fd.write(chr_num + "\t" + start + "\t" + end + "\t")

		# increase snp count
		num_snps_within_tss_range = num_snps_within_tss_range + 1
		num_bin_within_tss_range = num_bins_within_tss_range + 1

		try:
			num_count_in_ranges[range_kind] = num_count_in_ranges[range_kind] + 1
		except KeyError:
			num_count_in_ranges[range_kind] = 1



		# print ref/alt values per sample
		if (methyl_snp == "snp"):
			for file_num in range(len(file_list)):
				try:
					tokens2 	= data[key][file_num].split('~')
					ref_alt		= tokens2[0] 
				except KeyError:
					ref_alt		= "-/-"
					
				if file_num == 0:
					merge_out_file_fd.write(ref_alt)
				else:
					merge_out_file_fd.write("\t" + ref_alt)

			merge_out_file_fd.write("\t")

		# print af1_rms values per sample
		rank_correlation_x=[]
		for file_num in range(len(file_list)):
			try:
				tokens2 	= data[key][file_num].split('~')
				af1_rms		= tokens2[1] 
			except KeyError:
				af1_rms		= '-'

			if file_num == 0:
				merge_out_file_fd.write(af1_rms)
			else:
				merge_out_file_fd.write("\t" + af1_rms)

			# store af1_rms value for rank correlation
			rank_correlation_x.append(af1_rms)

		merge_out_file_fd.write("\t" + range_start + "\t" + range_end + "\t" + strand + "\t" + refseq + "\t" + gene_symbol + "\t" + range_kind + "\t")


		# print out exp level per sample
		rank_correlation_y=[]
		for file_num in range(len(exp_file_list)):
			try:
				exp_level=str(exp_data[gene_symbol,file_num]['value'])
			except KeyError:
				exp_level="-"

			if (file_num == 0):
				merge_out_file_fd.write(exp_level)
			else:
				merge_out_file_fd.write("\t" + exp_level)

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
		if (len(rank_correlation_x) > 0) and (len(rank_correlation_x) > (len(file_list)*0.3)) :
			# compute spearman'rho 
			rank_corr=1 - Bio.Cluster.distancematrix((rank_correlation_x,rank_correlation_y), dist="s")[1][0]
			merge_out_file_fd.write("\t" + str(rank_corr) + "\n")

			num_snps_having_corr = num_snps_having_corr + 1

			try:
				num_having_corr_in_ranges[range_kind] = num_having_corr_in_ranges[range_kind] + 1
			except KeyError:
				num_having_corr_in_ranges[range_kind] = 1


			if rank_corr > corr_threshold or rank_corr < -corr_threshold : 
				num_snp_over_threshold = num_snp_over_threshold + 1
				try:
					num_over_threshold_in_ranges[range_kind] = num_over_threshold_in_ranges[range_kind] + 1
				except KeyError:
					num_over_threshold_in_ranges[range_kind] = 1

				if rank_corr > 0 :
					num_snp_having_positive_corr = num_snp_having_positive_corr + 1
					try:
						num_having_positive_corr_in_ranges[range_kind] = num_having_positive_corr_in_ranges[range_kind]  + 1
					except KeyError:
						num_having_positive_corr_in_ranges[range_kind] = 1

				if rank_corr < 0 :
					num_snp_having_negative_corr = num_snp_having_negative_corr + 1

					try:
						num_having_negative_corr_in_ranges[range_kind] = num_having_negative_corr_in_ranges[range_kind]  + 1
					except KeyError:
						num_having_negative_corr_in_ranges[range_kind] = 1
				#num_gene_having_snp_over_threshold++
				


		else:
			merge_out_file_fd.write("\t-\n")

# print out statistics
if (methyl_snp=="snp"):

	# snp
	stat_out_file_fd.write("(1). Number of SNPs within TSS +- range :" + str(num_snps_within_tss_range) + "\n")
	stat_out_file_fd.write("(2). Number of SNPs having Spermans's rho with gene expression among (1) :" + str(num_snps_having_corr) + "\n")
	stat_out_file_fd.write("(3). Number of SNPs having Spermans's rho(rank correlation) over +- " + str(corr_threshold) + " among (2) :" + str(num_snp_over_threshold) + "\n")
	stat_out_file_fd.write("(4). % of SNP over rank correlation threshold +-" + str(corr_threshold) + " (1)/(3) : " + str(float(num_snp_over_threshold)/float(num_snps_having_corr) * 100) + "%\n")
	stat_out_file_fd.write("(5). % of SNP having positive correlation among (4) : " + str((float(num_snp_having_positive_corr)/float(num_snps_having_corr)*100)) + "%\n")
	stat_out_file_fd.write("(6). % of SNP having negative correlation among (4) : " + str((float(num_snp_having_negative_corr)/float(num_snps_having_corr))*100) + "%\n")

else:
	# methyl
	for range_kind in num_count_in_ranges.keys():
		stat_out_file_fd.write("(1). Number of bin within " + range_kind + " range :" + str(num_count_in_ranges[range_kind]) + "\n")
		stat_out_file_fd.write("(2). Number of bin having Spermans's rho with gene expression among (1) :" + str(num_having_corr_in_ranges[range_kind]) + "\n")
		stat_out_file_fd.write("(3). Number of bins having Spermans's rho(rank correlation) over +- " + str(corr_threshold) + " among (2) :" + str(num_over_threshold_in_ranges[range_kind]) + "\n")
		stat_out_file_fd.write("(4). % of bins over rank correlation threshold +-" + str(corr_threshold) + " (1)/(3) : " + str(float(num_having_corr_in_ranges[range_kind])/float(num_having_corr_in_ranges[range_kind]) * 100) + "%\n")
		stat_out_file_fd.write("(5). % of bins having positive correlation among (4) : " + str((float(num_having_positive_corr_in_ranges[range_kind])/float(num_having_corr_in_ranges[range_kind])*100)) + "%\n")
		stat_out_file_fd.write("(6). % of SNP having negative correlation among (4) : " + str((float(num_snp_having_negative_corr)/float(num_snps_having_corr))*100) + "%\n")

# gene
#print "Number of genes having SNP which has correlation over threshold : " + str(num_gene_having_snp_over_threshold)

#print exp_data
		

					

		
		



		






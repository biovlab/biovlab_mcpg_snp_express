import sys
import Bio.Cluster

# parameter check
if len(sys.argv) < 7 :
	print "[USAGE : python bin/merge_methyl.py [rms methyl sample_file_list] [gene expression sample file list] [correlation threshold] [merge out_file] [stat out file] [gene_symbol]"
        sys.exit(0)


# variables
in_fd_list			= []
methyl_file_list	 		= sys.argv[1].split(';')
exp_file_list 			= sys.argv[2].split(';')
corr_threshold			= float(sys.argv[3])
merge_out_file_fd		= open(sys.argv[4], "w")
stat_out_file_fd		= open(sys.argv[5], "w")
sample_num 			= len(methyl_file_list)
data				= {}				# double dictionary data = { sample : value }
exp_data			= {}
common_position_count		= {}
header				= ""
gene_symbol_fd			= open(sys.argv[6], "r")
gene_symbol_list		= []

# stat variables
# methyl
methyl_count_within_TSS_range	= 0
num_methyl_over_threshold		= 0
num_methyl_having_positive_corr 	= 0
num_methyl_having_negative_corr 	= 0
num_methyls_within_tss_range	= 0
num_methyls_having_corr	 	= 0
# gene 
num_gene_having_methyl_over_threshold = 0

# extract gene symbol 
for gene_symbol in gene_symbol_fd.readlines():
	gene_symbol_list.append(gene_symbol)

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

		#print gene_symbol
		#print exp_level
		try:
			exp_data[gene_symbol,file_num]['value'] = (exp_data[gene_symbol,file_num]['value'] * exp_data[gene_symbol,file_num]['count'] + exp_level) / (exp_data[gene_symbol,file_num]['count'] + 1)
			exp_data[gene_symbol,file_num]['count'] += 1
		except KeyError:
			exp_data[gene_symbol,file_num] ={}
			exp_data[gene_symbol,file_num]['count'] = 1
			exp_data[gene_symbol,file_num]['value'] = exp_level

		#if (gene_symbol == "TP73"):
			#print file_num
			#print exp_data[gene_symbol,file_num]['value']


	in_exp_fd.close()


		


# result from above
# methyl file1
#chr    bin_start       bin_end methyl_count    methyl_level(rms)       range_start     range_end       strand  refseq  gene_symbol     range_kind
#chr1    10401   10500   29      0.334728794189335       10373   11873   +       NR_046018,      DDX11L1 TSS_range

# methyl file2
#chr    bin_start       bin_end methyl_count    methyl_level(rms)       range_start     range_end       strand  refseq  gene_symbol     range_kind
#chr1    10401   10500   29      0.334728794189335       10373   11873   +       NR_046018,      DDX11L1 TSS_range

# find common position of methyl

# since all file has same coordinate read each line from files and compute correlation to reduce memory usage for large data
# extract methyl files

# open all files


for file_num in range(len(methyl_file_list)):
	temp_file = methyl_file_list[file_num]
	in_fd_list.append(open(temp_file, "r"))

count=-1
while 1:
	# init
	count +=1
	rank_correlation_x=[]
	rank_correlation_y=[]

	for file_num in range(len(in_fd_list)):

				
		in_fd=in_fd_list[file_num]

		line = in_fd.readline()

		# exit eof
		if not line: break

		rms 		= line[:-1]

		if (rms == ""): 
			print "[Error] The input line is : " + line[:-1]
			sys.exit(0)

		if line.startswith('#') : 
			header = line[:-1]
			continue

		gene_symbol	= gene_symbol_list[count]

		rank_correlation_x.append(rms)

		# store exp_level for rank correlation
		try:
			rank_correlation_y.append(str(exp_data[gene_symbol,file_num]['value']))
		except KeyError:
			if (gene_symbol=="TP73"):
				print "bad"
				print file_num
				print exp_data[gene_symbol,file_num]['value']
			rank_correlation_y.append("-")

	if len(rank_correlation_x) ==0: continue

	for temp_rms in rank_correlation_x:
		merge_out_file_fd.write(temp_rms + "\t")
	for temp_exp in rank_correlation_y:
		merge_out_file_fd.write(temp_exp + "\t")
		

	remove_index=[]
	for i in range(len(rank_correlation_x)):
		if rank_correlation_x[i] == '-' or rank_correlation_x[i] == "NA":
			remove_index.append(i)
		if rank_correlation_y[i] == '-':
			remove_index.append(i)
	# remove duplicate
	remove_index=list(set(remove_index))

	# sort index to remove the list from backward
	remove_index.sort(reverse=True)

	# remove list by index
	for index in remove_index:
		del rank_correlation_x[index]
		del rank_correlation_y[index]

#	print rank_correlation_x
#	print rank_correlation_y
	if (len(rank_correlation_x) > 0) and (len(rank_correlation_x) > (len(methyl_file_list)*0.3)) :
		# compute spearman'rho 
		rank_corr=1 - Bio.Cluster.distancematrix((rank_correlation_x,rank_correlation_y), dist="s")[1][0]
		merge_out_file_fd.write(str(rank_corr) + "\n")

		num_methyls_having_corr += 1

		if rank_corr > corr_threshold or rank_corr < -corr_threshold : 
			num_methyl_over_threshold += 1

			if rank_corr > 0 :
				num_methyl_having_positive_corr += 1
			elif rank_corr < 0 :
				num_methyl_having_negative_corr += 1
			#num_gene_having_methyl_over_threshold++

	else:
		merge_out_file_fd.write("\t-\n")



	# exit eof
	if not line: break

for temp_fd in in_fd_list:
	temp_fd.close()

sys.exit(1)
# find methyls with count over threshold (that is, common methyls), and merge methyls

threshold=1

for key in data.keys():

	# print only methyls over threshold
	if (data[key]['count'] > threshold):

		tokens 		= key.split('~')

		chr_num		= tokens[0]
		start		= tokens[1]
		end		= str(int(start) + 1)
		methyl_range_start	= tokens[2]
		methyl_range_end	= tokens[3]
		refseq		= tokens[4]
		gene_symbol	= tokens[5]

		# print as original format
		merge_out_file_fd.write(chr_num + "\t" + start + "\t" + end + "\t")

		# increase methyl count
		num_methyls_within_tss_rang = num_methyls_within_tss_range + 1


		# print ref/alt values per sample
		for file_num in range(len(methyl_file_list)):

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

		# print af1 values per sample
		rank_correlation_x=[]
		for file_num in range(len(methyl_file_list)):
			try:
				tokens2 	= data[key][file_num].split('~')
				af1		= tokens2[1] 
			except KeyError:
				af1		= '-'

			if file_num == 0:
				merge_out_file_fd.write(af1)
			else:
				merge_out_file_fd.write("\t" + af1)

			# store af1 value for rank correlation
			rank_correlation_x.append(af1)

		merge_out_file_fd.write("\t" + methyl_range_start + "\t" + methyl_range_end + "\t" + strand + "\t" + refseq + "\t" + gene_symbol + "\t")


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
		#print len(methyl_file_list)*0.3

		# check if the list is empty 
		if (len(rank_correlation_x) > 0) and (len(rank_correlation_x) > (len(methyl_file_list)*0.3)) :
			# compute spearman'rho 
			rank_corr=1 - Bio.Cluster.distancematrix((rank_correlation_x,rank_correlation_y), dist="s")[1][0]
			merge_out_file_fd.write("\t" + str(rank_corr) + "\n")

			num_methyls_having_corr = num_methyls_having_corr + 1


			if rank_corr > corr_threshold or rank_corr < -corr_threshold : 
				num_methyl_over_threshold = num_methyl_over_threshold + 1

				if rank_corr > 0 :
					num_methyl_having_positive_corr = num_methyl_having_positive_corr + 1
				if rank_corr < 0 :
					num_methyl_having_negative_corr = num_methyl_having_negative_corr + 1
				#num_gene_having_methyl_over_threshold++
				


		else:
			merge_out_file_fd.write("\t-\n")

# print out statistics
# methyl
stat_out_file_fd.write("(1). Number of METHYLs within TSS +- range :" + str(num_methyls_within_tss_range) + "\n")
stat_out_file_fd.write("(2). Number of METHYLs having Spermans's rho with gene expression among (1) :" + str(num_methyls_having_corr) + "\n")
stat_out_file_fd.write("(3). Number of METHYLs having Spermans's rho(rank correlation) over +- " + str(corr_threshold) + " among (2) :" + str(num_methyl_over_threshold) + "\n")
stat_out_file_fd.write("(4). % of METHYL over rank correlation threshold +-" + str(corr_threshold) + " (1)/(3) : " + str(float(num_methyl_over_threshold)/float(num_methyls_having_corr) * 100) + "%\n")
stat_out_file_fd.write("(5). % of METHYL having positive correlation among (4) : " + str((float(num_methyl_having_positive_corr)/float(num_methyls_having_corr)*100)) + "%\n")
stat_out_file_fd.write("(6). % of METHYL having negative correlation among (4) : " + str((float(num_methyl_having_negative_corr)/float(num_methyls_having_corr))*100) + "%\n")



import sys

# parameter check
if len(sys.argv) < 7 :
	print "[USAGE : python bin/profile_ME_GE.py [merged ME_GE_COR file] [out file] [sample num] [correlation threshold] [out matrix file for theshold] [out matrix file for PC NC ratio]"
        sys.exit(0)

# variables
range_count={}
range_having_corr_count={}
range_corr_over_threshold_count={}
range_corr_overed_pos_count={}
range_corr_overed_neg_count={}
corr=""

corr_threshold = float(sys.argv[4])
stat_out_file_fd = open(sys.argv[2], "w")
matrix_out_file_threshold_fd = open(sys.argv[5], "w")
matrix_out_file_PC_NC_fd = open(sys.argv[6], "w")
sample_num=int(sys.argv[3])
#range_kind_index=sample_num+8-1
#corr_index=(sample_num*2)+10-1
range_kind_index=sample_num+9-1
corr_index=(sample_num*2)+11-1

in_fd = open(sys.argv[1], "r")

for line in in_fd.readlines():
	# skip header
	if line[0] == '#': continue

	tokens = line[:-1].split('\t')

	try:
		range_kind=tokens[range_kind_index]
		corr=tokens[corr_index]

	except IndexError:
		print tokens
		sys.exit(0)
	#print corr

	try:
		range_count[range_kind] +=1
	except KeyError:
		range_count[range_kind] =1

	if corr !='-' and corr !="NA" and corr !="nan":
		try:
			range_having_corr_count[range_kind] +=1
		except KeyError:
			range_having_corr_count[range_kind] =1

		if float(corr) > corr_threshold or float(corr) < -corr_threshold:
			try:
				range_corr_over_threshold_count[range_kind] +=1
			except KeyError:
				range_corr_over_threshold_count[range_kind] =1

			if float(corr) > 0:
				try:
					range_corr_overed_pos_count[range_kind] +=1
				except KeyError:
					range_corr_overed_pos_count[range_kind] =1

			elif float(corr) < 0:
				try:
					range_corr_overed_neg_count[range_kind] +=1
				except KeyError:
					range_corr_overed_neg_count[range_kind] =1

matrix_out_file_threshold_fd.write("stat\tvariable\tvalue\n")
matrix_out_file_PC_NC_fd.write("stat\tvariable\tvalue\n")
		
for range_kind in range_count.keys():

	stat_out_file_fd.write("(1). Number of bins within " + range_kind + " range :" + str(range_count[range_kind]) + "\n")
	try:
		stat_out_file_fd.write("(2). Number of bins having Correlation with gene expression among (1) :" + str(range_having_corr_count[range_kind]) + "\n")
	except KeyError:
		stat_out_file_fd.write("(2). Number of bins having Correlation with gene expression among (1) : 0\n")
	try:
		stat_out_file_fd.write("(3). Number of bins having Correlation over +- " + str(corr_threshold) + " among (2) :" + str(range_corr_over_threshold_count[range_kind]) + "\n")
	except KeyError:
		stat_out_file_fd.write("(3). Number of bins having Correlation over +- " + str(corr_threshold) + " among (2) : 0\n")
	try:
		stat_out_file_fd.write("(4). % of bins over Correlation threshold +-" + str(corr_threshold) + " (3)/(2) : " + str(float(range_corr_over_threshold_count[range_kind])/float(range_having_corr_count[range_kind]) * 100) + "%\n")
	except KeyError:
		stat_out_file_fd.write("(3). Number of bins having Correlation over +- " + str(corr_threshold) + " among (2) : 0\n")
	try:
		stat_out_file_fd.write("(5). % of bins having positive correlation among (4) : " + str((float(range_corr_overed_pos_count[range_kind])/float(range_corr_over_threshold_count[range_kind])*100)) + "%\n")
	except KeyError:
		stat_out_file_fd.write("(5). % of bins having positive correlation among (4) : 0%\n") 
	try:
		stat_out_file_fd.write("(6). % of bins having negative correlation among (4) : " + str((float(range_corr_overed_neg_count[range_kind])/float(range_corr_over_threshold_count[range_kind]))*100) + "%\n")
	except KeyError:
		stat_out_file_fd.write("(5). % of bins having positive correlation among (4) : 0%\n") 

	stat_out_file_fd.write("\n")
	try:
		temp_value='%.2f'%(float(range_corr_over_threshold_count[range_kind])/float(range_having_corr_count[range_kind]) * 100)
	except KeyError:
		temp_value='%.2f'%(0.00)
	temp_value2=100-float(temp_value)
	matrix_out_file_threshold_fd.write("Corr_>_threshold\t"+range_kind+"\t" + str(temp_value)+"\n")
	matrix_out_file_threshold_fd.write("Corr_<_threshold\t"+range_kind+"\t" + str(temp_value2)+"\n")

	try:
		temp_value='%.2f'%(float(range_corr_overed_pos_count[range_kind])/float(range_corr_over_threshold_count[range_kind])*100)
	except KeyError:
		temp_value='%.2f'%(0.00)
	temp_value2=100-float(temp_value)
	matrix_out_file_PC_NC_fd.write("Negative_correlation\t"+range_kind+"\t" + str(temp_value2) + "\n")
	matrix_out_file_PC_NC_fd.write("Positive_correlation\t"+range_kind+"\t" + str(temp_value) + "\n")

stat_out_file_fd.close()
matrix_out_file_threshold_fd.close()
matrix_out_file_PC_NC_fd.close()

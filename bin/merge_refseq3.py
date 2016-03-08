#!bin/python

import sys

in_fd = open(sys.argv[1], "r")
refseq_column=int(sys.argv[2])

refseq_column_index=refseq_column-1


# assume refseq column 9

temp={}

# sample 
# methyl
# chr1    10401   10500   15      0.450757106523661       10373   11873   +       NR_046018       DDX11L1	TSS_range
# chr1    10501   10600   13      0.337175188241805       10373   11873   +       NR_046018       DDX11L1	TSS_range
# chr1    10601   10700   8       0.0679522028884321      10373   11873   +       NR_046018       DDX11L1	TSS_range

# snp
# chr10   95447   95448   G/A     +       DP=2;QS=0.000000,1.000000,0.000000,0.000000;VDB=5.960000e-02;AF1=1;AC1=2;DP4=0,0,2,0;MQ=20;FQ=-33       chr10   0       545729  +       NM_212479       ZMYND11
# chr10   95447   95448   G/A     +       DP=2;QS=0.000000,1.000000,0.000000,0.000000;VDB=5.960000e-02;AF1=1;AC1=2;DP4=0,0,2,0;MQ=20;FQ=-33       chr10   0       550577  +       NM_006624       ZMYND11

# ASSUME data is sorted. Possible to reduce memory usage
# extract 
previous_position_key=""
for line in in_fd.readlines():
	tokens = line[:-1].split('\t')

	# header	
	if tokens[0].startswith('#'):
		print line[:-1]
		continue

	# remove n/a data
	if tokens[refseq_column_index + 1] == "n/a":
		continue

	# column num check
	try:
		refseq = tokens[refseq_column_index]
	except IndexError:
		print "[INFO] IndexError"
		print tokens
		sys.exit(0)

	# init variable
	position_key=""

	# check within same chr~start~end key since TSS TSE is not sorted
	for i in range(0, 3):
		position_key = position_key + tokens[i] + "~"


	# create key
	data_key=""
	gene_symbol=tokens[refseq_column_index + 1]
	range_kind=tokens[refseq_column_index + 2]

	for i in range(0, refseq_column_index):
		data_key = data_key + tokens[i] + "~"
	data_key = data_key + gene_symbol +"~"+ range_kind

	# for first data
	if previous_position_key=="":
		previous_position_key=position_key
		temp[data_key]=[]
		temp[data_key].append(refseq)
		continue

	# different position key -> merge & print within stored data
	if previous_position_key != position_key:
		
		for key in temp.keys():
			tokens2 = key.split('~')

			# print data keys
			for i in range(0,refseq_column_index):
				sys.stdout.write(tokens2[i] + "\t")

			# print merged refseq
			for key_index in range(len(temp[key])):
				refseq2 = temp[key][key_index]
				if key_index==0:
					sys.stdout.write(refseq2)
				else:
					sys.stdout.write("," + refseq2)


			gene_symbol=tokens2[refseq_column_index]	# since in the key there is no refseq, so refseq index point gene symbol
			range_kind=tokens2[refseq_column_index+1]	# since in the key there is no refseq, so refseq index point gene symbol

			sys.stdout.write("\t"+gene_symbol+"\t")
			print range_kind

		# init
		temp={}

	# store data
	try:
		if refseq not in temp[data_key]:
			temp[data_key].append(refseq)
	except KeyError:
		temp[data_key]=[]
		temp[data_key].append(refseq)

	previous_position_key = position_key 

# for last 
for key in temp.keys():
	tokens2 = key.split('~')

	for i in range(0,refseq_column_index):
		sys.stdout.write(tokens2[i] + "\t")

	for key_index in range(len(temp[key])):
	#for refseq in temp[key]:
		refseq2 = temp[key][key_index]
		if key_index==0:
			sys.stdout.write(refseq2)
		else:
			sys.stdout.write("," + refseq2)



	#	sys.stdout.write(refseq+",")

	gene_symbol=tokens2[refseq_column_index]	# since in the key there is no refseq, so refseq index point gene symbol
	range_kind=tokens2[refseq_column_index+1]	# since in the key there is no refseq, so refseq index point gene symbol

	sys.stdout.write("\t"+gene_symbol+"\t")
	print range_kind

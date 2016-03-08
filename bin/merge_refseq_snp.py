#!bin/python

import sys

in_fd = open(sys.argv[1], "r")
refseq_column=int(sys.argv[2])

refseq_column_index=refseq_column-1


# assume refseq column 9

temp={}

# sample 
# chr1    10401   10500   15      0.450757106523661       10373   11873   +       NR_046018       DDX11L1
# chr1    10501   10600   13      0.337175188241805       10373   11873   +       NR_046018       DDX11L1
# chr1    10601   10700   8       0.0679522028884321      10373   11873   +       NR_046018       DDX11L1

# snp example
# chr10   95352   95353   G/A     +       DP=2;QS=0.000000,1.000000,0.000000,0.000000;VDB=5.960000e-02;AF1=1;AC1=2;DP4=0,0,0,2;MQ=20;FQ=-33       chr10   0       550577  +       NM_001202465    ZMYND11


# extract 
for line in in_fd.readlines():
	tokens = line[:-1].split('\t')
	if tokens[0].startswith('#'):
		header=line[:-1]
		continue

	try:
		refseq = tokens[refseq_column_index]
	except IndexError:
		print "[INFO] IndexError"
		print tokens
		sys.exit(0)

	temp_key=""

	
	for i in range(0, refseq_column_index):
		temp_key = temp_key + tokens[i] + "~"

	temp_key = temp_key + tokens[refseq_column]


	try:
		if refseq not in temp[temp_key]:
			temp[temp_key].append(refseq)
	except KeyError:
		temp[temp_key]=[]
		temp[temp_key].append(refseq)




# print merged refseq

print header

for key in temp.keys():
	tokens = key.split('~')

	for i in range(0,refseq_column_index):
		sys.stdout.write(tokens[i] + "\t")

	
	for refseq in temp[key]:
		sys.stdout.write(refseq+",")
	sys.stdout.write("\t")
	

	print tokens[refseq_column_index]	









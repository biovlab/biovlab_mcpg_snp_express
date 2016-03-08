#!bin/python

import sys

in_fd = open(sys.argv[1], "r")
refseq_column=int(sys.argv[2])

refseq_column_index=refseq_column-1


# assume refseq column 9

#temp={}
temp=[]

# to avoid sorting at the end
temp_key_list=[]


# sample 
# methyl
# chr1    10401   10500   15      0.450757106523661       10373   11873   +       NR_046018       DDX11L1
# chr1    10501   10600   13      0.337175188241805       10373   11873   +       NR_046018       DDX11L1
# chr1    10601   10700   8       0.0679522028884321      10373   11873   +       NR_046018       DDX11L1

# snp
# chr10   95447   95448   G/A     +       DP=2;QS=0.000000,1.000000,0.000000,0.000000;VDB=5.960000e-02;AF1=1;AC1=2;DP4=0,0,2,0;MQ=20;FQ=-33       chr10   0       545729  +       NM_212479       ZMYND11
# chr10   95447   95448   G/A     +       DP=2;QS=0.000000,1.000000,0.000000,0.000000;VDB=5.960000e-02;AF1=1;AC1=2;DP4=0,0,2,0;MQ=20;FQ=-33       chr10   0       550577  +       NM_006624       ZMYND11

# ASSUME data is sorted. Possible to reduce memory usage
# extract 
previous_key=""
previous_refseq=""
for line in in_fd.readlines():
	tokens = line[:-1].split('\t')
	if tokens[0].startswith('#'):
		header=line[:-1]
		print header
		continue

	if tokens[refseq_column_index + 1] == "n/a":
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

	temp_key = temp_key + ''.join(tokens[refseq_column].split(','))

	if previous_key=="":
		previous_key=temp_key
	#	previous_refseq=refseq
		temp.append(refseq)
		continue


	#print temp_key
	# check need to be merge
	if previous_key != temp_key:
		
		#if refseq not in temp:
		#	temp.append(refseq)

	#else:
		# different key, so no need to be merge -> print!

		tokens2 = previous_key.split('~')

		for i in range(0,refseq_column_index):
			sys.stdout.write(tokens2[i] + "\t")

		
		for refseq2 in temp:
			sys.stdout.write(refseq2+",")
		sys.stdout.write("\t")
		

		print tokens2[refseq_column_index]	

		temp=[]
# python ../../bin/merge_refseq2.py 100730_s_8.fq.bin_TSS_range 9 | grep -v "n/a" > t
	if refseq not in temp:
		temp.append(refseq)
	previous_key = temp_key

# for last 
tokens2 = previous_key.split('~')
for i in range(0,refseq_column_index):
	sys.stdout.write(tokens2[i] + "\t")
for refseq2 in temp:
	sys.stdout.write(refseq2+",")
sys.stdout.write("\t")
print tokens2[refseq_column_index]	




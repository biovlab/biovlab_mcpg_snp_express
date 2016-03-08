import sys
import Bio.Cluster

# parameter check
if len(sys.argv) < 3 :
	print "[USAGE : python bin/get_GDSNP.py [merged snp file] [class num] [class list]"
        sys.exit(0)

# variables
merged_snp_file_fd		= open(sys.argv[1], "r")
sample_class_num		= int(sys.argv[2])
sample_class_list		= sys.argv[3].split(';')
sample_num			= len(sample_class_list)
sample_class_set 	= list(set(sample_class_list))		

af_value_start_index 	= 3 + sample_num 
af_value_end_index 	= 3 + sample_num + sample_num

# threshold
GDSNP_exist_pct_in_others=0.1
GDSNP_exist_pct_in_the_subtype=0.3

# init variables
# chr19	53849632	53849633	-/-	-/-	c/cC	-/-	c/cC	-/-	-/-	-/-	-	-	1	-	1	-	-	-	53837001	53858122	+	ZNF845	genebody	NM_138374-	-	-	-	-	-	-	-	-

# get total count for each class 
class_total_count={}
class_total_count_mul_GDSNP_exist_pct_in_the_subtype={}
class_total_count_mul_GDSNP_exist_pct_in_others={}

for sample_class in sample_class_list:
	try:
		class_total_count[sample_class] +=1
	except KeyError:
		class_total_count[sample_class] = 1

	class_total_count_mul_GDSNP_exist_pct_in_the_subtype[sample_class] = class_total_count[sample_class] * GDSNP_exist_pct_in_the_subtype
	class_total_count_mul_GDSNP_exist_pct_in_others[sample_class] = class_total_count[sample_class] * GDSNP_exist_pct_in_others


# find GDSNP: should follow rules 1. $GDSNP_exist_pct_in_the_subtype 2. $GDSNP_exist_pct_in_others
for line in merged_snp_file_fd.readlines():
	tokens 			= line[:-1].split('\t')

	# get af1 value list
	af1_list 		= tokens[af_value_start_index:af_value_end_index]

	line_class_exist_count={}

	# init 
	for sample_class in sample_class_set:
		line_class_exist_count[sample_class]=0

	# example : - - 1 - 1 - - -
	for af1_index in range(len(af1_list)):

		# example : a a a b b b c c
		current_index_class = sample_class_list[af1_index]

		if af1_list[af1_index] !='-' and af1_list[af1_index] != '0':
			line_class_exist_count[current_index_class] += 1

	#check GDSNP 
	# example : a b c
	is_GDSNP=0
	for sample_class_first in sample_class_set:

		# check $GDSNP_exist_pct_in_the_subtype
		if (line_class_exist_count[sample_class_first] >= class_total_count_mul_GDSNP_exist_pct_in_the_subtype[sample_class_first]):
			is_GDSNP=1
			for sample_class_second in sample_class_set:
				if sample_class_first == sample_class_second : continue

				# check $GDSNP_exist_pct_in_others
				if (line_class_exist_count[sample_class_second] >= class_total_count_mul_GDSNP_exist_pct_in_others[sample_class_second]):
					is_GDSNP = 0
					break

			if (is_GDSNP == 1):
				print line[:-1]
				break

import sys
import scipy
from scipy.stats import pearsonr
from scipy.stats import spearmanr

first_value_start_col=int(sys.argv[2])
second_value_start_col=int(sys.argv[3])
sample_num=int(sys.argv[4])
pearson_spearman=int(sys.argv[5])

first_value_end_col=first_value_start_col + sample_num
second_value_end_col=second_value_start_col + sample_num

for line in sys.stdin:
	tokens = line[:-1].split('\t')

	try:
		x=scipy.array(map(float,tokens[first_value_start_col:first_value_end_col]))
		y=scipy.array(map(float,tokens[second_value_start_col:second_value_end_col]))
	except ValueError:
		print line[:-1] + "\tNA"
		continue


	#print x
	#print y
	try:
		if pearson_spearman == 0:
			corr,pval=pearsonr(x, y)
		else:
			corr,pval=spearmanr(x, y)
			
	except ValueError:
		print "[INFO] Error in values"
		print "[INFO] methylation"
		print x
		print "[INFO] gene expression"
		print y
		sys.exit(0)
	print line[:-1] + "\t" + str(corr)


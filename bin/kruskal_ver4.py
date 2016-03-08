import sys
import scipy.stats as stats
import csv
from multiprocessing import Pool
import statsmodels.stats.multitest as smm
from functools import partial
import math
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

def split_list(val_l, sample_num_l): 
	index = 0
	result_l = []

	for i in sample_num_l:
		result_l.append(val_l[index:index+i])
		index+=i
	
	return result_l

def calc_kruskal(x, sample_num_l, alpha):
	tmp_input_l = split_list(x[1:],sample_num_l) #ignore id column

	try:
#		h,p = stats.kruskal(*tmp_input_l) #run kruskal-wallist test
		h,p = stats.f_oneway(*tmp_input_l)
	except ValueError:
		return x+[float('nan')]
	
#	if math.isnan(p) :
#		return x+['1.00','0']	

	result = [`p`]

	return x+result
		


#stats.kruskal(*a)

# u,p=stats.mannwhitneyu([1,2,3,4],[5,6,7,9])

#input1 : input file
#input2 : class num (ex. 13,7,10)
#input3 : p-value
#input4 : cpu_num
#input5 : output file
#input6 : correction methods ('bonferroni', 'fdr', 'holm' and so on. see R stats.p.adjust)

with open(sys.argv[1], 'r') as f:
	reader = csv.reader(f, delimiter="\t")
	
	header_line = next(f)

	header_line = header_line.rstrip() + '\tp_value\tp.adj\n'

	sample_num_l = map(int,sys.argv[2].split(","))
	
	alpha = float(sys.argv[3])

	partial_kruskal = partial(calc_kruskal, sample_num_l=sample_num_l, alpha=alpha)
	
	pool = Pool(processes=int(sys.argv[4]))

	result = pool.map(partial_kruskal,[row for row in reader])

	p_val_list=[]

	for elem in result:
		p_val_list += [float(elem[-1])]
	
	stats = importr('stats')

	p_adjust = stats.p_adjust(FloatVector(p_val_list), method = sys.argv[6])

#	rej, pval_corr = smm.multipletests(p_val_list, alpha=alpha, method=sys.argv[6])[:2]

	for index in range(len(result)):
			result[index] = result[index] + [`p_adjust[index]`]
	
	with open(sys.argv[5], 'w') as f_out:
		f_out.write(header_line)
		f_out.writelines('\t'.join(i) + '\n' for i in result)
	
#	with open(sys.argv[5], 'r') as correc:
#		correc_reader = csv.reader(correc, delimiter="\t")
		
#		correc_header_line = next(correc)
#		correc_header_line = correc_header_line.rstrip() + '\tp.adj'
#
#		p_val_list=[]

#		for row in correc_reader:
#			p_val_list += [row[-2]]
		
#		w = csv.writer(f_out, dialect='excel-tab')
#		w.writerow(result)

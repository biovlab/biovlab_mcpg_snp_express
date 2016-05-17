#!/usr/local/bin/python
import os
import json
import sys

def main(base_dir):
	#target_dir = base_dir + "/integrated"

		
	#file_tf_tg = target_dir + "/Cyto_panel_TF_TG.txt"
	#file_tg_exp = target_dir + "/Cyto_panel_TG_exp.txt"
	#file_tg_info = target_dir + "/Cyto_panel_TG_Info.txt"
	#file_tg_pm = target_dir + "/Cyto_panel_TG_promoter_methyl.txt"

	file_tf_tg = sys.argv[1]
	file_tg_exp = sys.argv[2]
	file_tg_info = sys.argv[3]
	file_tg_pm = sys.argv[4]
	file_json = sys.argv[5]

	#file_output = base_dir + "/final_result/methylation/cyto_extended.txt"
	file_output = sys.argv[6]
	# open json
	with open(file_output, "w") as f_out:
		with open(file_json, "r") as f_json:
			i_json = json.load(f_json)
			## tf-tg
			i_tf_tg = {} # tf -> array(tg) / tg -> array(tf)
			with open(file_tf_tg, "r") as f_tf_tg:
				lines = f_tf_tg.readlines()
				# skip first line
				for line in lines[1:]:
					(tf, tg) = line.split()
					item = { 'name': tg }
					if (tf in i_tf_tg):
						i_tf_tg[tf].append( item )
					else: 
						i_tf_tg[tf] = [ item ]

					item = { 'name': tf }
					if (tg in i_tf_tg):
						i_tf_tg[tg].append( item )
					else: 
						i_tf_tg[tg] = [ item ]

			## tg_exp
			i_tg_exp = {} # gene -> (class_pair, exp_value)
			with open(file_tg_exp, "r") as f_tg_exp:
				lines = f_tg_exp.readlines()
				# split headers, except first item
				headers = lines[0].split()[1:]
				# from second line
				for line in lines[1:]:
					entries = line.split()
					gene = entries[0]
					i_tg_exp[gene] = {}
					for idx, exp in enumerate( entries[1:] ):
						i_tg_exp[gene][ headers[idx] ] = exp

			## tg_info
			i_tg_info = {} # gene/promoter -> (promoter/gene: corr)
			with open(file_tg_info, "r") as f_tg_info:
				lines = f_tg_info.readlines();
				for line in lines[1:]:
					(gene, promoter, corr) = line.split()
					i_tg_info[gene] = { promoter: corr }
					i_tg_info[promoter] = { gene: corr }

			## tg_pm
			classes = []
			i_tg_pm = {} # promoter -> ( class: value )
			with open(file_tg_pm, "r") as f_tg_pm:
				lines = f_tg_pm.readlines();
				# split headers, except first two
				classes = headers = lines[0].split()[1:]
				# from second line
				for line in lines[1:]:
					entries = line.split()
					promoter = entries[0]
					i_tg_pm[promoter] = {}
					for idx, meth in enumerate( entries[1:] ):
						i_tg_pm[promoter][ headers[idx] ] = meth

			res = {}
			key_set = set( list(i_tf_tg.keys()) + list(i_tg_exp.keys()) + list(i_tg_pm.keys()) )
			for key in key_set:
				res[key] = {
					'target': i_tf_tg.get(key),
					'exp': i_tg_exp.get(key),
					'corr': i_tg_info.get(key),
					'meth': i_tg_pm.get(key)
				}
			
			i_json['info'] = res
			i_json['classes'] = classes
			f_out.write( json.dumps(i_json) )




			



if __name__ == '__main__':
	main(sys.argv[1])

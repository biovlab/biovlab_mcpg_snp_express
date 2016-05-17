#!/usr/local/bin/python
import os
import json
import sys

def check_exist_add(info_dict, temp_dict, data_type, class_type):
	print data_type
	if data_type not in info_dict:
		info_dict[data_type]=[]
		temp_dict_local={}
		info_dict[data_type].append(temp_dict_local)

	if class_type not in info_dict[data_type][0]:
		info_dict[data_type][0][class_type]=[]

	info_dict[data_type][0][class_type].append(temp_dict)

	print info_dict
	return info_dict


# get argument as
	# data_type;class_type;print_type;datas...,
	# print type
		# normal : name, grid, type
		# tag : name, grid, type, tag
		# link : name, grid, type, link
		# dummy : grid

def add_data(info_dict, data_string):
	data_list=data_string.split(',')
	for data in data_list:
		tokens=data.split(';')

		# variables
		data_type=tokens[0]
		class_type=tokens[1]
		print_type=tokens[2]

		temp_dict={}

		# set each print types

		# dummy
		if print_type=="dummy":
			temp_dict["grid"]=tokens[3]

			
			info_dict=check_exist_add(info_dict, temp_dict, data_type, class_type)
			
		#	if data_type not in info_dict:
		#		info_dict[data_type]=[]
		#		info_dict[data_type].append({})
		#	else
		#		if class_type not in info_dict[data_type[0]]
		#			info_dict[data_type][0][class_type]=[]
		#			info_dict[data_type][0][class_type].append(temp_dict)
		#		else					
		#			info_dict[data_type][0][class_type].append(temp_dict)
			continue;

		# normal + others
		temp_dict["name"]=tokens[3]
		temp_dict["grid"]=tokens[4]
		temp_dict["type"]=tokens[5]

		# others
		if print_type!="normal":
			temp_dict[print_type]=tokens[6]
		
		print info_dict
		print temp_dict
		print data_type
		print class_type
		info_dict=check_exist_add(info_dict, temp_dict, data_type, class_type)

	#	if data_type not in info_dict:
	#		info_dict[data_type]=[]
	#		info_dict[data_type].append({})
	#	else
	#		if class_type not in info_dict[data_type[0]]
	#			info_dict[data_type][0][class_type]=[]
	#			info_dict[data_type][0][class_type].append(temp_dict)
	#		else					
	#			info_dict[data_type][0][class_type].append(temp_dict)

		return info_dict

def main():

	# variables	
	file_json = sys.argv[1]
	#data_type=["gene_expression", "methylation", "mutation"]
	class_kind_list=sys.argv[2].split(';')
	data_string=sys.argv[3] # data_type;class_type;print_type;datas...,
	final_json_info={}

	print "[INFO] input clases : " + sys.argv[2]

	final_json_info["result"]=[]

	# init data type dict
	temp_dict={}
	#for type in data_type:
	#	temp_dict[type]=[]

	# init clas pairs
	#temp_sub_dict={}
	#temp_sub_dict["ALL_class"]=[]
	#for i in range(0, len(class_kind_list)):
	#	for j in range(i+1, len(class_kind_list)):
	#		pair_prefix=class_kind_list[i]+"_vs_"+class_kind_list[j]
	#		temp_sub_dict[pair_prefix]=[]

	# add cass pairs to data type dict
	#for type in data_type:
	#	temp_dict[type].append(temp_sub_dict)

	temp_dict=add_data(temp_dict, data_string)	

	# add data type dict to final json
	final_json_info["result"].append(temp_dict)
	print json.dumps(final_json_info, sort_keys=True, indent=4, separators=(',', ': ')) 

	sys.exit()
	temp={}
	with open(file_json, "w") as f_out:
		with open(file_hall_mark_gene_set, "r") as f_hall:


			for line in f_hall:
				tokens=line.strip().split('\t')

				hall_id=tokens[0]
				hall_site=tokens[1]
				hall_genes=tokens[2:]

				temp[hall_id]=hall_genes

			final_json_info["hallmark"]=temp

			f_out.write( json.dumps(final_json_info, sort_keys=True, indent=4, separators=(',', ': ')) )

if __name__ == '__main__':
	main()

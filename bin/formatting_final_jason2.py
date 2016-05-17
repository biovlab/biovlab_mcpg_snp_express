#!/usr/local/bin/python
import os
import json
import sys

def check_exist_add(info_dict, temp_dict, data_type, class_type):

	# check data exist
	if data_type not in info_dict:
		info_dict[data_type]=[]
		temp_dict_local={}
		info_dict[data_type].append(temp_dict_local)

	# check class exist
	if class_type not in info_dict[data_type][0]:
		info_dict[data_type][0][class_type]=[]

	# add to info_dict
	info_dict[data_type][0][class_type].append(temp_dict)

	return info_dict


# get argument as
	# data_type;class_type;print_type;datas...,
	# print type
		# normal : name, grid, type
		# tag : name, grid, type, tag
		# link : name, grid, type, link
		# dummy : grid, type

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
			temp_dict["grid"]=int(tokens[3])
			temp_dict["type"]=tokens[4]
			info_dict=check_exist_add(info_dict, temp_dict, data_type, class_type)
			continue;

		# normal + others
		temp_dict["name"]=tokens[3]
		temp_dict["grid"]=int(tokens[4])
		temp_dict["type"]=tokens[5]

		# others
		if print_type!="normal":
			temp_dict[print_type]=tokens[6]
		
		info_dict=check_exist_add(info_dict, temp_dict, data_type, class_type)

	return info_dict

def main():
	# variables	
	file_json = sys.argv[1]
	data_string=sys.argv[2] # data_type;class_type;print_type;datas...,
	final_json_info={}

	print data_string
	#print "[INFO] input clases : " + sys.argv[2]

	final_json_info["result"]=[]

	# init data type dict
	temp_dict={}
	temp_dict=add_data(temp_dict, data_string)	

	# add data type dict to final json
	final_json_info["result"].append(temp_dict)

	with open(file_json, "w") as f_out:
			f_out.write( json.dumps(final_json_info, sort_keys=True, indent=4, separators=(',', ': ')) )

	print json.dumps(final_json_info, sort_keys=True, indent=4, separators=(',', ': ')) 


if __name__ == '__main__':
	main()

#!/usr/local/bin/python
import os
import json

def main():
	file_onco_tsg1="/data/project/MSG/lib/driver_genes_by_CNV.csv"
	file_onco_tsg2="/data/project/MSG/lib/driver_genes_by_subtle_mutations.csv"
	file_json = "/data/project/MSG/lib/onco_tsg.json"

	onco_tsg_info={}
	temp={}
	with open(file_json, "w") as f_out:
		with open(file_onco_tsg1, "r") as f_onco_tsg1:
			with open(file_onco_tsg2, "r") as f_onco_tsg2:

				header1 = f_onco_tsg1.readline();
				header1 = f_onco_tsg2.readline();

				for line in f_onco_tsg1:
					tokens=line.strip().split(',')

					onco_tsg_id=tokens[0]
					onco_tsg_class=tokens[2]
					onco_tsg_core_pathway=tokens[3]
					onco_tsg_process=tokens[4]
	
					try:
						temp[onco_tsg_class].append(onco_tsg_id)
					except KeyError:
						temp[onco_tsg_class]=[]
						temp[onco_tsg_class].append(onco_tsg_id)

				for line in f_onco_tsg2:
					tokens=line.strip().split(',')

					onco_tsg_id=tokens[0]
					onco_tsg_class=tokens[4]
					onco_tsg_core_pathway=tokens[5]
					onco_tsg_process=tokens[6]

					try:
						temp[onco_tsg_class].append(onco_tsg_id)
					except KeyError:
						temp[onco_tsg_class]=[]
						temp[onco_tsg_class].append(onco_tsg_id)

		f_out.write( json.dumps(temp) )

if __name__ == '__main__':
	main()

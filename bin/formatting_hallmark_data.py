#!/usr/local/bin/python
import os
import json
import sys

def main():
	file_hall_mark_gene_set="/data/project/MSG/lib/hallmark_genes_msigdb.txt"
	file_json = "/data/project/MSG/lib/hallmark_genes_msigdb.json"

	hall_info={}
	temp={}
	with open(file_json, "w") as f_out:
		with open(file_hall_mark_gene_set, "r") as f_hall:

			for line in f_hall:
				tokens=line.strip().split('\t')

				hall_id=tokens[0]
				hall_site=tokens[1]
				hall_genes=tokens[2:]

				temp[hall_id]=hall_genes

			hall_info["hallmark"]=temp

			f_out.write( json.dumps(hall_info, sort_keys=True, indent=4, separators=(',', ': ')) )

if __name__ == '__main__':
	main()

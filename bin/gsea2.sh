#!/bin/sh

# create DAVAID links

# Referenced by http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html
# Usage
# http://david.abcc.ncifcrf.gov/api.jsp?type=xxxxx&ids=XXXXX,XXXXX,XXXXXX,&tool=xxxx&annot=xxxxx,xxxxxx,xxxxx,

# type  =  one of DAVID recognized gene types
# annot  = a list of desired annotation  categories separated by ","
# ids  = a list of user's gene IDs separated by ","
# tool  = one of DAVID tool names

# examples
# 1. functional annotation summary page
# http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=2919,6347,6348,6364&tool=summary

# 2. functional annotation chart
# http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=2919,6347,6348,6364&tool=chartReport&annot=GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE

#########################################################################################
# NOTE : DAVID HAVE LIMITATION UP TO 400 GENES & HTTP URL HAS 2048 LENGTH LIMIT !!!
#########################################################################################


# variables
URL_length_limit=2048
gene_list_file=$1
result_dir=$2
id_type=$3
out_prefix=$4

type="?type=$id_type"
http_prefix="http://david.abcc.ncifcrf.gov/api.jsp"
tool="&tool=summary"
annot="GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE"

# fixed query
fixed_query=$http_prefix$type"&ids="$tool$annot
# fixed query length & remaining length
((remain_length=URL_length_limit-${#fixed_query}))
#echo $remain_length


# read gene list to array
IFS=";" input_gene_list=($(awk '{if (substr($1, 0, 1)=="#" || $1=="") next; printf "%s", $1";"}' $gene_list_file))
#echo $temp
#IFS=";" read -a input_gene_list <<< "$temp"

# create gene_id_list
for (( i=0; i<${#input_gene_list[@]}; i++ ));  
do
	if [ $i -eq 0 ]; then
		temp_ids=${input_gene_list[$i]}
	else
		temp_ids=$temp_ids,${input_gene_list[$i]}
	fi

	#echo $temp_ids
	#echo ${#temp_ids}

	# total query length check
	#temp_query=$http_prefix$type"&ids="$temp_ids$tool$annot

	#if [ ${#temp_query} -qt $URL_length_limit ]; then
	if [ ${#temp_ids} -gt $remain_length ]; then
		break
	else
		temp_ids2=$temp_ids
	fi
done

query=$http_prefix$type"&ids="$temp_ids2$tool
echo "Functional_annotation_summary="$query
echo $query > $result_dir/$out_prefix"Functional_annotation_summary.txt"

tool="&tool=chartReport"
annot="&annot=GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE" 
query=$http_prefix$type"&ids="$temp_ids2$tool$annot
echo "Functional_annotation_chart="$query
echo $query > $result_dir/$out_prefix"Functional_annotation_chart.txt"

tool="&tool=annotationReport"
annot="&annot=GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE"
query=$http_prefix$type"&ids="$temp_ids2$tool$annot
echo "Functional_annotation_table="$query
echo $query > $result_dir/$out_prefix"Functional_annotation_table.txt"

tool="&tool=term2term"
annot="&annot=GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE"
query=$http_prefix$type"&ids="$temp_ids2$tool$annot
echo "Functional_annotation_clustering="$query
echo $query > $result_dir/$out_prefix"Functional_annotation_clustering.txt"

#!/bin/bash

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
temp_ids=""
for (( i=0; i<${#input_gene_list[@]}; i++ ));  
do
	if [ $i -eq 0 ]; then
		temp_ids=${input_gene_list[$i]}
	else
		temp_ids=$temp_ids,${input_gene_list[$i]}
	fi

	if [ $i -ge 400 ]; then 
		echo "[INFO] The input gene list over 400 genes. Only top 400 are used to generate DAVID link. Please directly use DAVID site to input more than 3000 genes"
		break
	fi

	
#	echo $temp_ids 
#	echo $i
#	echo ${#temp_ids}
#	echo $remain_length

	# total query length check
	#temp_query=$http_prefix$type"&ids="$temp_ids$tool$annot

	#if [ ${#temp_query} -qt $URL_length_limit ]; then
	if [ ${#temp_ids} -gt $remain_length ]; then
		break
	fi
done
echo "[INFO] temp_ids:$temp_ids"

query=$http_prefix$type"&ids="$temp_ids$tool
david_html_file=$result_dir/$out_prefix".DAVID_summary.html"
david_txt_file=$result_dir/$out_prefix".DAVID_summary.txt"
echo "Functional_annotation_summary="$query
echo $query > $david_txt_file

# create david_html_file
echo '<html><head>' > $david_html_file
echo '<meta http-equiv="refresh" content="0;url='$query'" />' >> $david_html_file
echo '</head></html>' >> $david_html_file



tool="&tool=chartReport"
annot="&annot=GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE" 
query=$http_prefix$type"&ids="$temp_ids$tool$annot
david_html_file=$result_dir/$out_prefix".DAVID_chart.html"
david_txt_file=$result_dir/$out_prefix".DAVID_chart.txt"
echo "Functional_annotation_chart="$query
echo $query > $david_txt_file

# create david_html_file
echo '<html><head>' > $david_html_file
echo '<meta http-equiv="refresh" content="0;url='$query'" />' >> $david_html_file
echo '</head></html>' >> $david_html_file



tool="&tool=annotationReport"
annot="&annot=GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE"
query=$http_prefix$type"&ids="$temp_ids$tool$annot
david_html_file=$result_dir/$out_prefix".DAVID_table.html"
david_txt_file=$result_dir/$out_prefix".DAVID_table.txt"
echo "Functional_annotation_table="$query
echo $query > $david_txt_file

# create david_html_file
echo '<html><head>' > $david_html_file
echo '<meta http-equiv="refresh" content="0;url='$query'" />' >> $david_html_file
echo '</head></html>' >> $david_html_file



tool="&tool=term2term"
annot="&annot=GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE"
query=$http_prefix$type"&ids="$temp_ids$tool$annot
david_html_file=$result_dir/$out_prefix".DAVID.clustering.html"
david_txt_file=$result_dir/$out_prefix".DAVID_clustering.txt"
echo "Functional_annotation_clustering="$query
echo $query > $david_txt_file

# create david_html_file
echo '<html><head>' > $david_html_file
echo '<meta http-equiv="refresh" content="0;url='$query'" />' >> $david_html_file
echo '</head></html>' >> $david_html_file

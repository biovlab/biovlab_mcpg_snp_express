#!/bin/bash

delimeter=$6
IFS=$delimeter read -ra file_list_with_path <<< "$1"
IFS=$delimeter read -ra type_list <<< "$2"
IFS=$delimeter read -ra track_name_list <<< "$3"

output_dir=$4
genome_version=$5
# NOTE : this url should contains '/' at the end, so no more additional '/' required when use
base_url=$7
out_file_prefix=$8


file_list=()

# get file name only list and create symbolic links
for ((i=0;i<${#file_list_with_path[@]};i++)); do
	file_with_path=${file_list_with_path[$i]}

	# get file name list
	temp_file=`basename \$file_with_path`
	file_list+=($temp_file)
	extension="${temp_file##*.}"

	# create symbolic links
	ln -s $file_with_path $output_dir/$temp_file

	# indexed file? create symbolic links for index
	if [ ${type_list[$i]} == "vcfTabix" ]; then
		ln -s $file_with_path".tbi" $output_dir/$temp_file".tbi"
	fi
done



# variables
hub_name="BioVLAB_mCpG_SNP_EXPRESS"
short_label="BioVLAB_results"
long_labbel="BioVLAB_mCpG_SNP_EXPRESS_results"
hub_file_name=$out_file_prefix".UCSC_hub.txt"
genome_file=$out_file_prefix"_genomes.txt"
description_url="biovlab_description_url.html"
email_addr="snu.biovlab@gmail.com"
trackDb=$out_file_prefix"_trackDb.txt"
visibility="dense"

UCSC_hub_link="http://genome.ucsc.edu/cgi-bin/hgTracks?db=$genome_version&hubUrl=$base_url$hub_file_name"


# create hub.txt
{
echo "hub $hub_name"
echo "shortLabel $short_label"
echo "longLabel $long_labbel"
echo "genomesFile $genome_file"
echo "email $email_addr"
echo "descriptionUrl $description_url"
} > $output_dir/$hub_file_name

# crete genomes.txt
{
echo "genome  $genome_version"
echo "trackDb $trackDb"
} > $output_dir/$genome_file

# create trackDb file
{
for ((i=0;i<${#file_list[@]};i++)); do
	echo "track ${track_name_list[$i]}"
	echo "bigDataUrl ${file_list[$i]}"
	echo "shortLabel ${track_name_list[$i]}"
	echo "longLabel ${track_name_list[$i]}"
	echo "visibility $visibility"
	echo "type ${type_list[$i]}"
	echo ""
done
}>$output_dir/$trackDb

# create URL
echo $UCSC_hub_link




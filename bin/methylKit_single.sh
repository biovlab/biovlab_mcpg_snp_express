#!/bin/bash
source `dirname $0`/../env.sh

# parameters
sample_paths=$1
sample_ids=$2
class_ids=$3
result_dir=$4
window=$5

################################################
# Script for methylKit analysis
################################################

##run methylKit for single file statstics
#statstics will be saved in result_dir/file_name + ".stat.txt"
#base methylation level will be saved in result_dir/file_name + ".bed"

Rscript ./single_stat.r -s $sample_paths -i $sample_ids -c $class_ids -r $result_dir -w $window

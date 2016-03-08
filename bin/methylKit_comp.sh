#!/bin/bash
source `dirname $0`/../env.sh

# parameters
sample_paths1=$1
sample_ids1=$2
class_ids1=$3

sample_paths2=$4
sample_ids2=$5
class_ids2=$6

result_dir=$7
ref_path=$8
cpgi_ref_path=$9

window=${10}
normal=${11}


################################################
# Script for methylKit analysis
################################################

sample_paths=$sample_paths1";"$sample_paths2
sample_ids=$sample_ids1";"$sample_ids2
class_ids=$class_ids1";"$class_ids2


#run methylKit for sample comparision
Rscript ./comp_stat.r -s $sample_paths -i $sample_ids -c $class_ids -r $result_dir --reference $ref_path --cpgi_reference $cpgi_ref_path -w $window -n $normal


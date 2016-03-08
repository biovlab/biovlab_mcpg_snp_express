#!/bin/bash   
input_matrix=$1			# matrix format : (id	value1	value2	...)		
class_list=$2			# NOTE :  THIS VALUE SHOULD BE DELIMINATED BY ','. example : "a,a,a,b,b,b,c,c"
class_num=$3			# number of kind of class. for above class list the value whould be 3
num_process=$4			

cat $input_matrix |  parallel --no-notice --pipe -j$num_process -L1000000 -k python /data/project/mcpg/bin/get_mi_score_heejoon2.py - $class_list $class_num

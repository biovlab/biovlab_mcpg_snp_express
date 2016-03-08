#!/bin/bash
source `dirname $0`/../env.sh

# parameters
input_fq1=$1
input_fq2=$2

single_pair=$3
result_dir=$4
adaptor=$5

	
# variables
input_fq1_filename_only=`basename \$input_fq1`
input_fq2_filename_only=`basename \$input_fq2`
file_extention=$(echo $input_fq1_filename_only |awk -F . '{if (NF>1) {print $NF}}')
filename1_wo_ext=`basename $input_fq1_filename_only "."$file_extention`
filename2_wo_ext=`basename $input_fq2_filename_only "."$file_extention`

# temp testing
#ref_dir="/data/project/mcpg/ref_data/human/ucsc/masked_for_targeted_bs_seq/"


################################################
# Script for BS-seq data analysis
################################################


# adaptor clipping & align with bismark
echo "[INFO] Adaptor trimming and aligning"
echo "[INFO] Adaptor: $adaptor"
if [ "$adaptor" == "NA" ]; then
	echo "[INFO] Adaptor is not specified. Use standard adaptor"
       adaptor="AGATCGGAAGAGC"	# universal adaptor
fi

if [ "$single_pair" == "pair" ]; then
	echo "[INFO] This is paired end input"
	trim_galore --paired -a $adaptor -o $result_dir $input_fq1 $input_fq2	

	# NOTE : THIS REFERENCE GENOME IS MASKED FOR ICBP TARGETED BS SEQ!!! NEED TO USE DIFFERENT ONE FOR GENERAL USEAGE
	bismark --bam -o $result_dir $BISMARK_REF_DIR_HUMAN  -1 $result_dir/$filename1_wo_ext"_val_1.fq" -2 $result_dir/$filename2_wo_ext"_val_2.fq"

	# sort sam by samtools
	#samtools view -@ $NUM_CPUS -ubS $result_dir/$filename1_wo_ext"_val_1.fq_bismark_pe.sam" | samtools sort -@ $NUM_CPUS -m 3G  - $result_dir/$filename1_wo_ext"_val_1.fq_bismark_pe.sorted"
	samtools sort -@ $NUM_CPUS -m 3G $result_dir/$filename1_wo_ext"_val_1.fq_bismark_bt2_pe.bam" $result_dir/$filename1_wo_ext"_bismark_bt2_pe.sorted"
	samtools view -@ $NUM_CPUS $result_dir/$filename1_wo_ext"_bismark_bt2_pe.sorted.bam" > $result_dir/$filename1_wo_ext"_bismark_bt2_pe.sorted.sam"

else
	echo "[INFO] This is single end input"
	trim_galore -a $adaptor -o $result_dir $input_fq1
	bismark -o $result_dir $BISMARK_REF_DIR_HUMAN $result_dir/$filename1_wo_ext"_trimmed.fq" 

	# sort sam by samtools
	samtools sort -@ $NUM_CPUS -m 3G $result_dir/$filename1_wo_ext"_trimmed.fq_bismark_bt2_pe.bam" $result_dir/$filename1_wo_ext"_bismark_bt2_pe.sorted"
	samtools view -@ $NUM_CPUS $result_dir/$filename1_wo_ext"_bismark_bt2_pe.sorted.bam" > $result_dir/$filename1_wo_ext"_bismark_bt2_pe.sorted.sam" 
	
fi 

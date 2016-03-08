#!/bin/bash
source `dirname $0`/../../env.sh

data_dir=""

if [ "$is_array" -eq "1" ]; then
	data_dir=`dirname \${cel_list_temp[0]}`
else
	data_dir=`dirname \${sample_fq1_list_temp[0]}`
fi

cel_list=()
sample_fq1_list=()
sample_fq2_list=()

if [ "$is_array" -eq "1" ]; then
	for ((i=0;i<${#cel_list_temp[@]};i++)); do
		file_extension1=$(echo `basename \${cel_list_temp[\$i]}` | awk -F . '{if(NF>1){print $NF}}' | tr '[:upper:]' '[:lower:]')
	
		if [ $file_extension1 == "cel" ]; then
			cel_list+=(${cel_list_temp[$i]})
		elif [ $file_extension1 == "gz" ]; then
			filename1=$data_dir"/"`basename ${cel_list_temp[$i]} ".gz"`
			
			cel_list+=($filename1)
			
			if [ "$unzip" -eq "1" ]; then
				gunzip -d -c ${cel_list_temp[$i]} > $filename1
			fi
		else 
			echo "Invalid file format : "$file_extension1" (Allowed : CEL, gz)"
			exit
		fi
	done

else
	for ((i=0;i<${#sample_fq1_list_temp[@]};i++)); do
		file_extension1=$(echo `basename \${sample_fq1_list_temp[\$i]}` | awk -F . '{if(NF>1){print $NF}}' | tr '[:upper:]' '[:lower:]')
	
		if [ $file_extension1 == "fq" ]; then
			sample_fq1_list+=(${sample_fq1_list_temp[$i]})
		elif [ $file_extension1 == "gz" ]; then
			filename1=$data_dir"/"`basename ${sample_fq1_list_temp[$i]} ".gz"`
		
			sample_fq1_list+=($filename1)
			
			if [ "$unzip" -eq "1" ];then
				gunzip -d -c ${sample_fq1_list_temp[$i]} > $filename1
			fi
		else
			echo "Invalid file format : "$file_extension1" (Allowed : fq, gz)"
			exit
		fi

		if [ $single_pair == "pair" ]; then
			file_extension2=$(echo `basename \${sample_fq2_list_temp[\$i]}` | awk -F . '{if(NF>1){print $NF}}' | tr '[:upper:]' '[:lower:]')
			if [ $file_extension2 == "fq" ]; then
				sample_fq2_list+=(${sample_fq2_list_temp[$i]})
			elif [ $file_extension2 == "gz" ]; then
				filename2=$data_dir"/"`basename ${sample_fq2_list_temp[$i]} ".gz"`
		
				sample_fq2_list+=($filename2)
			
				if [ "$unzip" -eq "1" ]; then
					gunzip -d -c ${sample_fq2_list_temp[$i]} > $filename2
				fi
			else
				echo "Invalid file format : "$file_extension2" (Allowed : fq, gz)"
				exit
			fi
		else
			sample_fq2_list+=("NA")
		fi
	done
fi




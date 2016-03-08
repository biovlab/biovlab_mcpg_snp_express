#!/bin/bash
source `dirname $0`/../../env.sh

####################################################
# parse all experiment parameters
####################################################
declare -A exp_info

# MySQL query to extract all parameters by uid 
uid=$1
mysql_db_name=$MYSQL_DB_NAME
mysql_user_id=$MSQL_USER_ID
mysql_passwd=$MYSQL_PASSWD
query="$MYSQL_QUERY'$uid';"
mysql_hostname=$MYSQL_HOSTNAME
mysql_host_port=$MYSQL_HOST_PORT

echo "[INFO] PGA_UID:"$uid

# variables
field_name_list=()

ge_replicate_exist=-1
me_replicate_exist=-1
mu_replicate_exist=-1

ge_is_array=-1
me_is_array=-1
mu_is_array=-1

ge_single_pair=""
me_single_pair=""
mu_single_pair=""

ge_array_list=() # for Microarray
me_array_list=() # for Microarray
ge_array_list=() # for Microarray

ge_sample_fq1_list=() # for RNA-seq pair-1
ge_sample_fq2_list=() # for RNA-seq pair-2
me_sample_fq1_list=() # for RNA-seq pair-1
me_sample_fq2_list=() # for RNA-seq pair-2
mu_sample_fq1_list=() # for RNA-seq pair-1
mu_sample_fq2_list=() # for RNA-seq pair-2

ge_type_list=()
me_type_list=()
mu_type_list=()

ge_type_kind=()
me_type_kind=()
mu_type_kind=()

me_sample_format=-1

# read query result by line
tuple_num=-1 # init iter
while read line; do 
	# split single line data
	IFS=$'\t' read -ra fields <<< "$line"	

	echo "$line"
	# parse field names
	if [ ${#field_name_list[@]} -eq 0 ]; then
		field_name_list=("${fields[@]}"); continue
	fi

	((tuple_num++))
	# for each field assign value
	for (( field_index=0; field_index<${#field_name_list[@]}; field_index++ )); do
		# assign value with [tuple_num, field name] as key
		exp_info[$tuple_num,${field_name_list[$field_index]}]=${fields[$field_index]}	

	done
done < <(mysql --host=$mysql_hostname --port=$mysql_host_port --user="$mysql_user_id" --password="$mysql_passwd" --column-names=TRUE "$mysql_db_name"  --execute="$query")
total_num_tuple=$tuple_num	# this started with 0

# above will extract following tables. NOTE : sample_format is added. below is old version
#######################################################################################################################################################
#+------------+--------------+-------------+---------+------------+-----------------------------------------------------------------+-----------+-------------+------------------+----+------------------------------------------------+---------------+
#| class_name | patient_name | sample_type | is_pair | replicates | file_path                                                       | file_name | pair_number | replicate_number | id | uid                                            | analysis_type |
#+------------+--------------+-------------+---------+------------+-----------------------------------------------------------------+-----------+-------------+------------------+----+------------------------------------------------+---------------+
#| F          | S1           | Expression  |       1 |          1 | /usr/local/apache2/htdocs/airavata-php-gateway/app/storage/temp | 2.txt     |           1 |                1 |  3 | expdbtest_216c621d_4c1d_46e0_b1c0_72c1ad92a40d |             0 |
#| F          | S1           | Expression  |       1 |          1 | /usr/local/apache2/htdocs/airavata-php-gateway/app/storage/temp | 5.txt     |           2 |                1 |  3 | expdbtest_216c621d_4c1d_46e0_b1c0_72c1ad92a40d |             0 |
#| F          | S1           | Methylation |       0 |          1 | /usr/local/apache2/htdocs/airavata-php-gateway/app/storage/temp | 1.txt     |           1 |                1 |  3 | expdbtest_216c621d_4c1d_46e0_b1c0_72c1ad92a40d |             0 |
#| F          | S1           | Mutation    |       0 |          1 | /usr/local/apache2/htdocs/airavata-php-gateway/app/storage/temp | 3.txt     |           1 |                1 |  3 | expdbtest_216c621d_4c1d_46e0_b1c0_72c1ad92a40d |             0 |
#+------------+--------------+-------------+---------+------------+-----------------------------------------------------------------+-----------+-------------+------------------+----+------------------------------------------------+---------------+
#######################################################################################################################################################
# NOTE : for the field name, please refer env.sh

####################################################
# assign parsed value to each list and variables
####################################################

#######################
# setup parameters 
#######################

# NOTE : WE ASSUMES ALL THE CLASS, REPLICATES, SINGLE/PAIR-end information comes in order
for (( i=0; i<=$total_num_tuple; i++ )); do
	work_file=${exp_info[$i,$FILE_PATH]}"/"${exp_info[$i,$FILE_NAME]}
	sample_type=${exp_info[$i,$SAMPLE_TYPE]}
	sample_format=${exp_info[$i,$SAMPLE_FORMAT]}
	class_name=${exp_info[$i,$CLASS_NAME]}
	pair_num=${exp_info[$i,$PAIR_NUMBER]}
	ge_type=${exp_info[$i,$EXPRESSION_TYPE]}
	me_type=${exp_info[$i,$METHYLATION_TYPE]}
	mu_type=${exp_info[$i,$MUTATION_TYPE]}
	replicates=${exp_info[$i, $REPLICATES]}

	####################################################
	# Gene expression
	####################################################
	if [ $sample_type == "Expression" ]; then

		if [ $sample_format == "Microarray" ]; then
			ge_is_array=1; ge_single_pair="single"; ge_array_list+=($work_file)			# Microarray
			ge_type_list+=($class_name)
		else
			ge_is_array=0																									# RNA-seq

			# 0:exclude, 1:single-end, 2:paired_end
			if [ "$ge_type" == "1" ]; then																		# single-end	
				ge_single_pair="single"; ge_sample_fq1_list+=($work_file)
				ge_type_list+=($class_name)
			else																														# paired-end
				ge_single_pair="pair"

				if [ "$pair_num" == "1" ]; then
					ge_sample_fq1_list+=($work_file) 																# for pair-1
					ge_type_list+=($class_name)
				else
					ge_sample_fq2_list+=($work_file) 																# for pair-2
				fi
			fi
		fi

		# replicate	
		if [ "$ge_replicate_exist" -eq -1 ] && [ "$replicates" != "1" ]; then
			ge_replicate_exist=1
		fi

	####################################################
	# Methylation
	####################################################
	elif [ $sample_type == "Methylation" ]; then	

		if [ "$sample_format" == "MBD-Seq / MeDIP-Seq" ]; then						# MBD/MeDIP	
			me_sample_format=0

			# 0:exclude, 1:single-end, 2:paired_end
			if [ $me_type == "1" ]; then																		# single-end	
				me_single_pair="single"; me_sample_fq1_list+=($work_file)
				me_sample_fq2_list+=("NA")
				me_type_list+=($class_name)
			else																														# paired-end
				me_single_pair="pair"

				if [ "$pair_num" == "1" ]; then
					me_sample_fq1_list+=($work_file) 																# for pair-1
					me_type_list+=($class_name)
				else
					me_sample_fq2_list+=($work_file) 																# for pair-2
				fi
			fi

	 	elif [ "$sample_format" == "Infinium 27k" ]; then									# Infinium 27/450k
			me_sample_format=1; me_single_pair="single"; 
			me_array_list+=($work_file)
		  me_type_list+=($class_name)

		elif [ "$sample_format" == "Infinium 450k" ]; then									# Infinium 27/450k
			# TODO : set for Infinium
			#me_is_array=1; 
			me_sample_format=2
			me_single_pair="single"; me_array_list+=($work_file)
		  me_type_list+=($class_name)
		else																															# BS-seq
			me_sample_format=3

			# 0:exclude, 1:single-end, 2:paired_end
			if [ "$me_type" == "1" ]; then																		# single-end	
				me_single_pair="single"; me_sample_fq1_list+=($work_file)
				me_sample_fq2_list+=("NA")
				me_type_list+=($class_name)
			else																														# paired-end
				me_single_pair="pair"

				if [ "$pair_num" == "1" ]; then
					me_sample_fq1_list+=($work_file) 																# for pair-1
					me_type_list+=($class_name)
				else
					me_sample_fq2_list+=($work_file) 																# for pair-2
				fi
			fi
		fi

		# replicate	
		if [ "$me_replicate_exist" -eq -1 ] && [ "$replicates" != "1" ]; then
			me_replicate_exist=1
		fi

	####################################################
	# Mutation
	####################################################
	elif [ $sample_type == "Mutation" ]; then		

		if [ $sample_format == "DNA-Seq" ] || [ $sample_format == "RNA-Seq" ]; then	# DNA-seq/RNA-seq
			if [ $sample_format == "DNA-Seq" ]; then
				mu_sample_format=0
			else
				mu_sample_format=1
			fi
			mu_is_array=0 

			# 0:exclude, 1:single-end, 2:paired_end
			if [ "$mu_type" == "1" ]; then																		# single-end	
				mu_single_pair="single"; mu_sample_fq1_list+=($work_file)
				mu_type_list+=($class_name)
			else																														# paired-end
				mu_single_pair="pair"

				if [ "$pair_num" == "1" ]; then
					mu_sample_fq1_list+=($work_file) 																# for pair-1
					mu_type_list+=($class_name)
				else
					mu_sample_fq2_list+=($work_file) 																# for pair-2
				fi
			fi

		elif [ $sample_format == "Humanhap 550" ]; then									# Humanhap 550
			# TODO : set for Infinium
			mu_sample_format=2
			mu_is_array=1; mu_single_pair="single"; mu_array_list+=($work_file)
			mu_type_list+=($class_name)
		fi

		# replicate
		if [ "$mu_replicate_exist" -eq -1 ] && [ "$replicates" != "1" ]; then
			mu_replicate_exist=1
		fi
	fi

done

# NOTE : THIS ASSUMES ALL DATA IS EXIST. THIS MAY NEED TO BE PROCESSED IN EACH SCRIPT.
ge_type_kind=($(echo ${ge_type_list[@]} | tr ' ' '\n' | sort | uniq | tr '\n' ' '))
me_type_kind=($(echo ${me_type_list[@]} | tr ' ' '\n' | sort | uniq | tr '\n' ' '))
mu_type_kind=($(echo ${mu_type_list[@]} | tr ' ' '\n' | sort | uniq | tr '\n' ' '))

# for testing
echo ""
echo "[INFO] Gene expression info"
echo "[INFO] Type kind : ${ge_type_kind[@]}"
echo "[INFO] Type list : ${ge_type_list[@]}"
echo "[INFO] Sample fq1 list: ${ge_sample_fq1_list[@]}"
echo "[INFO] Sample fq2 list: ${ge_sample_fq2_list[@]}"
echo "[INFO] Array list: ${ge_array_list[@]}"
echo "[INFO] Single or Pair: $ge_single_pair"
echo "[INFO] Is array?: $ge_is_array"

echo ""
echo "[INFO] Methylation info"
echo "[INFO] Type kind : ${me_type_kind[@]}"
echo "[INFO] Type list : ${me_type_list[@]}"
echo "[INFO] Sample fq1 list: ${me_sample_fq1_list[@]}"
echo "[INFO] Sample fq2 list: ${me_sample_fq2_list[@]}"
echo "[INFO] Array list: ${me_array_list[@]}"
echo "[INFO] Single or Pair: $me_single_pair"
echo "[INFO] Sample format: $me_sample_format"

echo ""
echo "[INFO] Mutation info"
echo "[INFO] Type kind : ${mu_type_kind[@]}"
echo "[INFO] Type list : ${mu_type_list[@]}"
echo "[INFO] Sample fq1 list: ${mu_sample_fq1_list[@]}"
echo "[INFO] Sample fq2 list: ${mu_sample_fq2_list[@]}"
echo "[INFO] Array list: ${mu_array_list[@]}"
echo "[INFO] Single or Pair: $mu_single_pair"
echo "[INFO] Is array?: $mu_is_array"




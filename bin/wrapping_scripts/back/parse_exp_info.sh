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

echo "[INFO] UID:"$uid

# variables
field_name_list=()

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
done < <(mysql --user="$mysql_user_id" --password="$mysql_passwd" --column-names=TRUE "$mysql_db_name"  --execute="$query")
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



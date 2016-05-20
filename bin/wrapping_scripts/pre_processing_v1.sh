#!/bin/bash

# TODO : check if this is local or cloud


uid=$1
output=$2	# result_root
resource=$3

# check local or cloud
if [ $resource == "localhost" ]; then
	echo "[INFO] This is in bhi2 local server"
else
	echo "[INFO] This is in cloud"
	cp `dirname $0`/../../../env.sh `dirname $0`/../../
fi
	
source `dirname $0`/../../env.sh
source $WORK_DIR/bin/wrapping_scripts/parse_exp_info.sh

for (( i=0; i<=$total_num_tuple; i++ )); do

	# Check if the input is url based or direct input (local only)
	if [ ${exp_info[$i,$URL]} == "NULL" ] && [ $resource == "localhost" ]; then
		echo "[INFO] User file name is : ${exp_info[$i,$FILE_NAME]}"
		echo "[INFO] User input is direct uploaded"
	else # URL based (local or cloud)
		echo "[INFO] User input is url based"
		echo "[INFO] file with url : $work_file_with_url"

		work_file_with_url=${exp_info[$i,$URL]}
		file_name_only=`basename $work_file_with_url`

		# get user data
		echo "[INFO] url based User data downloading"

		wget -c -N  $work_file_with_url -O $USER_DATA_UPLOAD_DIR/$file_name_only
		#lftp -c 'set xfer:clobber true;pget -n 10 -O /data/project/MSG/tmp/ http://147.46.15.115/temp/icbp_link/100730_s_4.fq'
		#command="set xfer:clobber true;pget -n 10 -O $USER_DATA_UPLOAD_DIR $work_file_with_url"
		#lftp -c "$command"
		#lftp -c 'set xfer:clobber true;pget -n 5 http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz'
	fi
done

<<'COMMENT'
else
#######################################
# cloud only
#######################################
	# NOTE : env.sh in the cloud is different with local one. It is pre-deployed, so need to be copy after update all the scripts.

	# NOTE : System assumes All inputs for cloud are url based

	echo "[INFO] This is in cloud"
	cp `dirname $0`/../../../env.sh `dirname $0`/../../
	source `dirname $0`/../../env.sh

	# parse experiment information
	source $WORK_DIR/bin/wrapping_scripts/parse_exp_info.sh

											# replaced with github
											# copy all latest scripts to amazon
											# maybe use bbcp for faster copying
											# bbcp -P 10 -f -T 'ssh -i "$AMAZON_SSH_PRIVATE_KEY" -x -a %I -l %U %H bbcp' \ $WORK_DIR/bin $WORK_DIR/env.sh airavata@$amazon_ip:$AMAZON_WORK_DIR


											#latest_scripts_file=$WORK_DIR/latest_scripts.tar.gz

											#tar cfz $latest_scripts_file -C $WORK_DIR bin 
											#scp -i $AMAZON_SSH_PRIVATE_KEY $latest_scripts_file airavata@$amazon_ip:$AMAZON_WORK_DIR
											#ssh -i $AMAZON_SSH_PRIVATE_KEY airavata@$amazon_ip tar xfz $AMAZON_WORK_DIR/latest_scripts.tar.gz

											#scp -r -i $AMAZON_SSH_PRIVATE_KEY $WORK_DIR/bin $WORK_DIR/env.sh airavata@$amazon_ip:$AMAZON_WORK_DIR

	# download URL based inpts
	# NOTE: Cloud assuems only accepting URL based input
	for (( i=0; i<=$total_num_tuple; i++ )); do
		# FIXME
		if [ ${exp_info[$i,$URL]} != "NULL" ]; then
			work_file_url=${exp_info[$i,$FILE_PATH]}"/"${exp_info[$i,$FILE_NAME]}
			work_file=${exp_info[$i,$FILE_PATH]}"/"${exp_info[$i,$FILE_NAME]}

			# excute download inputs command on cloud
			# maybe use lftp for fast file transfer?
			# lftp -c 'set xfer:clobber true;pget -n 5 http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz'
			echo "[INFO] User input is url based"

			work_file_with_url=${exp_info[$i,$URL]}
			echo "[INFO] file with url : $work_file_with_url"
			file_name_only=`basename $work_file_with_url`

			# get user data
			echo "[INFO] url based User data download"
			wget -c -N  $work_file_with_url -O $USER_DATA_UPLOAD_DIR/$file_name_only
		fi
		#wget -c -N $work_file_url -O $work_file
	done

	# set cloud specific parameters in env.sh

		#ssh -i $AMAZON_SSH_PRIVATE_KEY airavata@$amazon_ip

fi
COMMENT

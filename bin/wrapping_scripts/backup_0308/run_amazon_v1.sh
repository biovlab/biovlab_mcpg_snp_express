#!/bin/bash
source `dirname $0`/../../env.sh

# parse experiment information
source `dirname $0`/parse_exp_info.sh

#amazon_ip=52.21.9.63
amazon_ip=$1

# copy all latest scripts to amazon
# maybe use bbcp for faster copying
# bbcp -P 10 -f -T 'ssh -i "$AMAZON_SSH_PRIVATE_KEY" -x -a %I -l %U %H bbcp' \ $WORK_DIR/bin $WORK_DIR/env.sh airavata@$amazon_ip:$AMAZON_WORK_DIR


latest_scripts_file=$WORK_DIR/latest_scripts.tar.gz

tar cfz $latest_scripts_file -C $WORK_DIR bin 
scp -i $AMAZON_SSH_PRIVATE_KEY $latest_scripts_file airavata@$amazon_ip:$AMAZON_WORK_DIR
ssh -i $AMAZON_SSH_PRIVATE_KEY airavata@$amazon_ip tar xfz $AMAZON_WORK_DIR/latest_scripts.tar.gz

#scp -r -i $AMAZON_SSH_PRIVATE_KEY $WORK_DIR/bin $WORK_DIR/env.sh airavata@$amazon_ip:$AMAZON_WORK_DIR

# download URL based inpts
for (( i=0; i<=$total_num_tuple; i++ )); do
	# FIXME
	work_file_url=${exp_info[$i,$FILE_PATH]}"/"${exp_info[$i,$FILE_NAME]}
	work_file=${exp_info[$i,$FILE_PATH]}"/"${exp_info[$i,$FILE_NAME]}

	# excute download inputs command on cloud
	# maybe use lftp for fast file transfer?
	# lftp -c 'set xfer:clobber true;pget -n 5 http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz'
	ssh -i $AMAZON_SSH_PRIVATE_KEY airavata@$amazon_ip wget $work_file_url -o $work_file
done

# set cloud specific parameters in env.sh

	#ssh -i $AMAZON_SSH_PRIVATE_KEY airavata@$amazon_ip



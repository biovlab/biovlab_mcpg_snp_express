#!/bin/bash
source `dirname $0`/../../env.sh

# parse experiment information
source `dirname $0`/parse_exp_info.sh


NUM_CPUS=$me_cpu

# variables
bin_dir="$WORK_DIR/bin"
<<'COMMENT'
result_dir="$WORK_DIR/result"
COMMENT
#result_dir="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/small_test_data/targeted_bs_breast_cancer_80genes/small_test/"

#functions
function my_join { local IFS="$1"; shift; echo "$*"; }

# directories
bs_result_dir=$result_dir"/methyl/bs"
methylkit_from_bs_result_dir=$result_dir"/methyl/bs/methylkit"

mbd_result_dir=$result_dir"/methyl/mbd"
medips_from_mbd_result_dir=$result_dir"/methyl/mbd/medips"
final_result_dir=$final_result_root_dir/methylation
mkdir -p $final_result_dir
####################################################################################
# 0. Parse input samples
####################################################################################

# MBD analysis

# check BS or MBD or ....
me_type_kind=($(echo ${ge_type_list[@]} | tr ' ' '\n' | sort | uniq | tr '\n' ' '))

# assign parsed values to each variable
cel_list_temp=(${me_array_list[@]})
type_kind=(${me_type_kind[@]})
type_list=(${me_type_list[@]})
sample_fq1_list_temp=(${me_sample_fq1_list[@]})
sample_fq2_list_temp=(${me_sample_fq2_list[@]})
sample_format=$me_sample_format	# 0: MBD-Seq / MeDIP-Seq, 1: "Infinium 27k", 2:"Infinium 450k", 3:"BS-seq"
pair_single_flag=$me_single_pair

# preprocess input / unzip
unzip=1
source `dirname $0`/preprocess_me_input.sh


bs_r1_list=(${sample_fq1_list[@]})
bs_r2_list=(${sample_fq2_list[@]})
BS_TYPE_KIND=(${me_type_kind[@]})
BS_TYPE_LIST=(${me_type_list[@]})

mbd_r1_list=(${sample_fq1_list[@]})
mbd_r2_list=(${sample_fq2_list[@]})
MBD_TYPE_KIND=(${me_type_kind[@]})
MBD_TYPE_LIST=(${me_type_list[@]})


####################################################################################
# 1. DNA methylation
####################################################################################


#TODO : option setting (thanks to SJ)
####################################################################################

# temp sample list
#test_dir="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/small_test_data/targeted_bs_breast_cancer_80genes/"

# bs or mbd
#sample_format="mbd"

# bs inputs
#bs_r1_list=($test_dir"AU565_ATCACG_L001_R1_001.fastq.head100000.fq" $test_dir"BT549_CGATGT_L001_R1_001.fastq.head100000.fq")
#bs_r2_list=($test_dir"AU565_ATCACG_L001_R2_001.fastq.head100000.fq" $test_dir"BT549_CGATGT_L001_R2_001.fastq.head100000.fq")

# mbd inputs
#mbd_r1_list=($test_dir"../icbp/BrCa-02_head10000000.fastq" $test_dir"../icbp/BrCa-03_head10000000.fastq")
#mbd_r2_list=("NA" "NA")

# adapter seq
#adapter_seq="NA"

# pair? single?
#pair_single_flag="single"

# bin size
window_size=100

# class kind list
#GROUP_TYPE_KIND=("class1" "class2")

# class file list
#GROUP_TYPE_LIST=("class1" "class2")


#BS_TYPE_KIND=("${GROUP_TYPE_KIND[@]}")
#BS_TYPE_LIST=("${GROUP_TYPE_LIST[@]}")

#MBD_TYPE_KIND=("${GROUP_TYPE_KIND[@]}")
#MBD_TYPE_LIST=("${GROUP_TYPE_LIST[@]}")

pval_threshold=$P_VALUE_CUT

####################################################################################

if [ "$sample_format" -eq "3" ]; then # BS

	# 1.1 BS-seq

	# create result directory
	mkdir -p $bs_result_dir
	
	count=0
	#trimming & mapping
	for (( i=0; i<${#bs_r1_list[@]}; i++ )); do
		bash $bin_dir/bs_ver2.sh ${bs_r1_list[$i]} ${bs_r2_list[$i]} $pair_single_flag $bs_result_dir $NUM_CPUS

<<'COMMENT'
&
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait
COMMENT
	done

	# run methylkit
	mkdir -p $methylkit_from_bs_result_dir

	# run parallel using parallel? Use background and wait process!!
	count=0
	for (( i=0; i<${#bs_r1_list[@]}; i++ )); do
		# variables
		input_fq1_filename_only=`basename \${bs_r1_list[$i]}`
		file_extention=$(echo $input_fq1_filename_only |awk -F . '{if (NF>1) {print $NF}}')
		filename1_wo_ext=`basename $input_fq1_filename_only "."$file_extention`
		temp_input=$bs_result_dir/$filename1_wo_ext"_bismark_bt2_pe.sorted.sam"

		# run METHYLKIT
		$R_DIR/Rscript $bin_dir/single_stat.r -s $temp_input -i $input_fq1_filename_only -c $i -r $methylkit_from_bs_result_dir -w $window_size &

		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait

	echo "[INFO] methylkit single sample profiling part done" 

	#TODO : MAYBE SUBTYPE(GROUP) COMPARIONS NEEDED?

	# methylkit pair comparison part 
	count=0

	for (( k=0; k<${#BS_TYPE_KIND[@]}; k++)); do
		work_file=""
		for (( l=$k+1; l<${#BS_TYPE_KIND[@]}; l++ )); do
			echo "[INFO] Currently processing ${BS_TYPE_KIND[$k]} and ${BS_TYPE_KIND[$l]}"

			work_file=${BS_TYPE_KIND[$k]}"_vs_"${BS_TYPE_KIND[$l]}
		
			# init file list
			first_sample_class_list=()
			first_sample_class_name=()
			second_sample_class_list=()
			second_sample_class_name=()

			for (( i=0; i<${#bs_r1_list[@]}; i++ )); do

				if [ ${BS_TYPE_LIST[$i]} == ${BS_TYPE_KIND[$k]} ]; then
					input_fq1_filename_only=`basename \${bs_r1_list[$i]}`

					first_sample_class_list+=($methylkit_from_bs_result_dir"/"$input_fq1_filename_only"_CpG.txt")
					first_sample_class_name+=($input_fq1_filename_only)

				elif [ ${BS_TYPE_LIST[$i]} == ${BS_TYPE_KIND[$l]} ]; then
					input_fq1_filename_only=`basename \${bs_r1_list[$i]}`

					second_sample_class_list+=($methylkit_from_bs_result_dir"/"$input_fq1_filename_only"_CpG.txt")
					second_sample_class_name+=($input_fq1_filename_only)
				fi
			done
			
			class_ids="1"

			for(( i=0; i<${#first_sample_class_name[@]}-1; i++)); do
				class_ids+=",1"
			done
			
			for(( i=0; i<${#second_sample_class_name[@]}; i++)); do
				class_ids+=",0"
			done
			$R_DIR/Rscript $bin_dir/comp_stat_ver2.r -s `my_join "," ${first_sample_class_list[@]}`","`my_join "," ${second_sample_class_list[@]}` -i `my_join "," ${first_sample_class_name[@]}`","`my_join "," ${second_sample_class_name[@]}` -c $class_ids -r $methylkit_from_bs_result_dir -w $window_size -n FALSE --out_prefix $work_file --pval $P_VALUE_CUT &
			let count+=1
			[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
		done
	done; wait

	count=0
	for (( k=0; k<${#BS_TYPE_KIND[@]}; k++)); do
		work_file=""
		for (( l=$k+1; l<${#BS_TYPE_KIND[@]}; l++ )); do
			work_file=${BS_TYPE_KIND[$k]}"_vs_"${BS_TYPE_KIND[$l]}

			cut -f1-3 $methylkit_from_bs_result_dir"/"$work_file".DMR_list.txt" > $methylkit_from_bs_result_dir"/"$work_file".DMR.bed" &
			let count+=1
			[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
		done
	done; wait
		
	echo "[INFO] methylkit pair comparison part done"

	# convert beds(methylation levels) to bigWig file to visualize on UCSC
	echo "[INFO] Convert bed to bigWig to visualize methylation levels on UCSC"
	ls $methylkit_from_bs_result_dir/*.sam.bed | xargs -I {} bash -c 'tail -n +2 $1 | sort -k1,1 -k2,2n - > $1.sorted' -- {}
	ls $methylkit_from_bs_result_dir/*.sorted | xargs -I {} $bin_dir/bedGraphToBigWig {} $REF_HUMAN_CHR_SIZE {}".bw"
	
	# make met_level files
	for (( i=0; i<${#bs_r1_list[@]}; i++ )); do
	{
		# variables
		input_fq1_filename_only=`basename \${bs_r1_list[$i]}`
		file_extention=$(echo $input_fq1_filename_only |awk -F . '{if (NF>1) {print $NF}}')
		filename1_wo_ext=`basename $input_fq1_filename_only "."$file_extention`
		
		MR_input=$methylkit_from_bs_result_dir/$filename1_wo_ext"_bismark_bt2_pe.sorted.sam_MR.bed"
		MR_coverage_input=$methylkit_from_bs_result_dir/$filename1_wo_ext"_bismark_bt2_pe.sorted.sam_MR_coverage.bed"

		methylation_level_file=$methylkit_from_bs_result_dir/$input_fq1_filename_only".met_level"

		echo -e "chr\tbin_start\tbin_end\tmethyl_count\tmethyl_level" > $methylation_level_file

		paste <(tail -n+2 $MR_coverage_input) <(tail -n+2 $MR_input) | cut -f1-4,8 >> $methylation_level_file
<<'COMMENT_TEMP'
		awk '
 				 BEGIN{FS=OFS="\t"}
  			 FILENAME == ARGV[1] {
   				 my_index=$1"/"$2"/"$3;
    			 arr[my_index]=$4;
   				 next
  			 }
  			 FILENAME == ARGV[2] {
    			 my_index=$1"/"$2"/"$3;
    			 arr[my_index]=arr[my_index]"\t"$4;
    			 next
  			 }
  			 FILENAME == ARGV[3] {
    			 my_index=$1"/"$2"/"$3;
    			 if(my_index in arr){
        			 print $1,$2,$3,arr[my_index];
    			 }
    			 else{
        			 print $1,$2,$3,0,0;
    			 }
  			 }
			' $MR_coverage_input $MR_input $REF_HUMAN_GENOME_100_BIN_BED >> $methylation_level_file
COMMENT_TEMP
	} &
	done; wait

	# copy all final results
#	cp $methylkit_from_bs_result_dir/*.png $methylkit_from_bs_result_dir/*.DMR.bed $methylkit_from_bs_result_dir/*.DMR_list.txt $final_result_dir

	# move to final results
	mv $methylkit_from_bs_result_dir/*CoverageStats.png $methylkit_from_bs_result_dir/*.MethylStats.png $methylkit_from_bs_result_dir/*.corr.png $final_result_dir	
	

	# for testing
	# small_set_dir="/data/project/mcpg/test_data/small_test_data"
	#sh bin/bs.sh $samll_set_dir"/targeted_bs_breast_cancer_80genes/ZR7530_AAAGCA_L001_R1_001.fastq.head100000" $samll_set_dir"/targeted_bs_breast_cancer_80genes/ZR7530_AAAGCA_L001_R2_001.fastq.head100000" pair $result_dir NA
	#sh bin/bs.sh $small_set_dir"/targeted_bs_breast_cancer_80genes/ZR751_GAAACC_L001_R1_001.fastq.head100000" $small_set_dir"/targeted_bs_breast_cancer_80genes/ZR751_GAAACC_L001_R2_001.fastq.head100000" pair $result_dir NA

NUM_CPUS=$(( $SYS_NUM_CPUS / 3 ))
elif [ "$sample_format" -eq "0" ]; then
	# 1.2 MBD-seq
	# create result directory
	mkdir -p $mbd_result_dir
	echo "[INFO] Start MBDseq analysis"

	echo "[INFO] Trimming in parallel"
	count=0
	for (( i=0; i<${#mbd_r1_list[@]}; i++ )); do
		if [ "$pair_single_flag" == "single" ]; then
			trim_galore -o $mbd_result_dir ${mbd_r1_list[$i]} &
		else
			trim_galore --paired -o $mbd_result_dir ${mbd_r1_list[$i]} ${mbd_r2_list[$i]} &
		fi
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done;wait

	echo "[INFO] Alignning by STAR"
	for (( i=0; i<${#mbd_r1_list[@]}; i++ )); do
			# variables
			input_fq1_filename_only=`basename \${mbd_r1_list[$i]}`
			file_extention=$(echo $input_fq1_filename_only |awk -F . '{if (NF>1) {print $NF}}')
			filename1_wo_ext=`basename $input_fq1_filename_only "."$file_extention`

		if [ "$pair_single_flag" == "single" ]; then
			cd $mbd_result_dir
			STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $mbd_result_dir/$filename1_wo_ext"_trimmed.fq" --runThreadN $SYS_NUM_CPUS --alignIntronMax 1 --outFileNamePrefix $mbd_result_dir/$input_fq1_filename_only"." --alignEndsType EndToEnd --seedSearchStartLmax 20 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 # --outFilterMatchNmin 15 # --alignEndsType EndToEnd --seedSearchStartLmax 30 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 options are for increasing sensitivity for less than 50bp reads. MBDseq ICBP 36bp length 
			cd $WORK_DIR

		else # for paired-end
			#variables
			input_fq2_filename_only=`basename \${mbd_r2_list[$i]}`
			filename2_wo_ext=`basename $input_fq2_filename_only "."$file_extention`

			cd $mbd_result_dir
			STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $mbd_result_dir/$filename1_wo_ext"_val_1.fq" $mbd_result_dir/$filename2_wo_ext"_val_2.fq" --runThreadN $SYS_NUM_CPUS --alignIntronMax 1 --outFileNamePrefix $mbd_result_dir/$input_fq1_filename_only"." --alignEndsType EndToEnd --seedSearchStartLmax 20 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 
			cd $WORK_DIR
		fi
	done


	echo "[INFO] Convert SAM to Sorted BAM"
	count=0
	for (( i=0; i<${#mbd_r1_list[@]}; i++ )); do
		# covert sam to sorted bam
		samtools view -ubS -@ $SYS_NUM_CPUS $mbd_result_dir/$input_fq1_filename_only".Aligned.out.sam" | samtools sort -@ $SYS_NUM_CPUS -m 3G - $mbd_result_dir/$input_fq1_filename_only".sorted"
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done;wait

	echo "[INFO] Index sorted BAM"
	count=0
	for (( i=0; i<${#mbd_r1_list[@]}; i++ )); do
		# covert sam to sorted bam
		samtools index $mbd_result_dir/$input_fq1_filename_only".sorted.bam" &
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done;wait


# for testing
	#sh bin/mbd.sh $fastq_sample1 NA single $result_dir 
	#sh bin/mbd.sh $fastq_sample2 NA single $result_dir

	#TODO : Create input list with merging biological replicate by ';'
	#TODO : Create sample num list with merging biological replicate by ';'

	# create result directory
	mkdir -p $medips_from_mbd_result_dir

	echo "[INFO] Start MEDIPS based on MBDseq analysis"

	# run MEDIPS to extract DMRs, ****SUB TYPE COMPARISON****
	for (( k=0; k<${#MBD_TYPE_KIND[@]}; k++ )); do 
		work_file=""
		for (( l=$k+1; l<${#MBD_TYPE_KIND[@]}; l++ )); do
			echo "[INFO] Currently processing ${MBD_TYPE_KIND[$k]} and ${MBD_TYPE_KIND[$l]}"

			work_file=${MBD_TYPE_KIND[$k]}"_vs_"${MBD_TYPE_KIND[$l]}

			# init list
			first_sample_class_name=()
			second_sample_class_name=()

			# TODO : create strand specific BAM file since MEDIPS do not care about strand, divide bam to strand specific bam and run MEDIPS and merge with strand information. However, the medips dmr extract process will be considered as well.
			# create sample list per subtype
			for (( i=0; i<${#mbd_r1_list[@]}; i++ )); do
				if [ ${MBD_TYPE_LIST[$i]} == ${MBD_TYPE_KIND[$k]} ]; then
					first_sample_class_name+=($mbd_result_dir"/"`basename \${mbd_r1_list[$i]}`".sorted.bam")
				elif [ ${MBD_TYPE_LIST[$i]} == ${MBD_TYPE_KIND[$l]} ]; then
					second_sample_class_name+=($mbd_result_dir"/"`basename \${mbd_r1_list[$i]}`".sorted.bam")
				fi 
			done

			echo "[INFO] Input1 sample lists are : ${first_sample_class_name[@]}"
			echo "[INFO] Input2 sample lists are : ${second_sample_class_name[@]}"

			# run medips between 2 classes with replicate information
			# NOTE : For each process, it takes 10GB, so check memory availability before assign number of cpu.

			# TODO : Parallelize MEDIPS.meth() function to speed up
#			$R_DIR/Rscript $bin_dir/medips_pair_parallel_v2.R -P $NUM_CPUS -a `my_join ";" ${first_sample_class_name[@]}` -b `my_join ";" ${second_sample_class_name[@]}` -r $medips_from_mbd_result_dir -p ${MBD_TYPE_KIND[$k]}"_vs_"${MBD_TYPE_KIND[$l]} -w $window_size --pvalue $pval_threshold
		
			$R_DIR/Rscript $bin_dir/medips_pair_v2.R -a `my_join ";" ${first_sample_class_name[@]}` -b `my_join ";" ${second_sample_class_name[@]}` -r $medips_from_mbd_result_dir -p ${MBD_TYPE_KIND[$k]}"_vs_"${MBD_TYPE_KIND[$l]} --pvalue $P_VALUE_CUT 
			
			# DMR file format setting
			
			tail -n+2 $medips_from_mbd_result_dir"/"${MBD_TYPE_KIND[$k]}"_vs_"${MBD_TYPE_KIND[$l]}".mr.edgeR.s.txt" | cut -f1-3 > $medips_from_mbd_result_dir"/"${MBD_TYPE_KIND[$k]}"_vs_"${MBD_TYPE_KIND[$l]}".DMR.bed"		
	
			# create per sample methyl level file (based on rms)
			# get rms start end colums
			rms_start_column=$(( 5 + ( (${#first_sample_class_name[@]} + ${#second_sample_class_name[@]}) * 2 ) + ${#first_sample_class_name[@]} ))
			rms_end_column=$(( $rms_start_column + ${#second_sample_class_name[@]} -1 ))

			if [ $k -eq 0 ]; then
				echo "[INFO] Extract methylation level per sample"
				# NOTE : Since MEDIPS put X if the file name start with number, remove X from header
				filename_string=`my_join ";" ${second_sample_class_name[@]}`
				#awk -v $filename_string -v result_dir=$medips_from_mbd_result_dir -v first_class_num=${#first_sample_class_name[@]} -v second_class_num=${#second_sample_class_name[@]} -v rms_start_column=$rms_start_column -v rms_end_column=$rms_end_column 'NR==1{for (i=rms_start_column;i<=rms_end_column; i++){split($i, tokens, ".sorted.bam.rms");filename[i]=substr(tokens[1],1)".met_level";};next} {for (i=rms_start_column;i<=rms_end_column; i++){print $1"\t"$2"\t"$3"\t"$(i-((first_class_num + second_class_num)*2))"\t"$i> result_dir"/"filename[i]}}' $medips_from_mbd_result_dir"/"$work_file".mr.edgeR.txt"
				awk -v filename_string=$filename_string -v result_dir=$medips_from_mbd_result_dir -v first_class_num=${#first_sample_class_name[@]} -v second_class_num=${#second_sample_class_name[@]} -v rms_start_column=$rms_start_column -v rms_end_column=$rms_end_column 'NR==1{split(filename_string, filename_list, ";");temp=rms_end_column-rms_start_column+1;for (i=1 ;i<=temp; i++){split(filename_list[i], tokens, ".sorted.bam");filename[i]=tokens[1]".met_level";};next}{OFS="\t";j=1;for (i=rms_start_column;i<=rms_end_column; i++){print $1,$2,$3,$(i-((first_class_num + second_class_num)*2)),$i> result_dir"/"filename[j]; j=j+1}}' $medips_from_mbd_result_dir"/"$work_file".mr.edgeR.txt"

<<'COMMENT'
				awk -v filename_string=$filename_string -v result_dir="." -v first_class_num=2 -v second_class_num=2 -v rms_start_column=15 -v rms_end_column=16 'NR==1{split(filename_string, filename_list, ";");temp=rms_end_column-rms_start_column+1;for (i=1 ;i<=temp; i++){split(filename_list[i], tokens, ".sorted.bam");filename[i]=tokens[1]".met_level";print filename[i];};next}'

				# convert wig to bigWig to visualize on UCSC genome browser
				echo "[INFO] create bigwig for visualization on UCSC genome browser"
				wig_file=$medips_from_mbd_result_dir"/"$work_file".file2.wig"
				tail -n +2 $wig_file > $medips_from_mbd_result_dir"/"${MBD_TYPE_KIND[$l]}".notrack.wig"
				wigToBigWig $medips_from_mbd_result_dir"/"${MBD_TYPE_KIND[$l]}".notrack.wig" $REF_HUMAN_CHR_SIZE $medips_from_mbd_result_dir"/"${MBD_TYPE_KIND[$l]}".bw"
COMMENT
			fi


		done

		# for first class
		rms_start_column=$(( 5 + ( (${#first_sample_class_name[@]} + ${#second_sample_class_name[@]}) * 2 ) ))
		rms_end_column=$(( $rms_start_column + ${#first_sample_class_name[@]} -1 ))

		if [ $k -eq 0 ]; then
			echo "[INFO] Extract methylation level per sample"
			filename_string=`my_join ";" ${first_sample_class_name[@]}`
			awk -v filename_string=$filename_string -v result_dir=$medips_from_mbd_result_dir -v first_class_num=${#first_sample_class_name[@]} -v second_class_num=${#second_sample_class_name[@]} -v rms_start_column=$rms_start_column -v rms_end_column=$rms_end_column 'NR==1{split(filename_string, filename_list, ";");temp=rms_end_column-rms_start_column+1;for (i=1 ;i<=temp; i++){split(filename_list[i], tokens, ".sorted.bam");filename[i]=tokens[1]".met_level";};next}{OFS="\t";j=1;for (i=rms_start_column;i<=rms_end_column; i++){print $1,$2,$3,$(i-((first_class_num + second_class_num)*2)),$i> result_dir"/"filename[j]; j=j+1}}' $medips_from_mbd_result_dir"/"$work_file".mr.edgeR.txt"



<<'COMMENT'
			awk -v result_dir=$medips_from_mbd_result_dir -v first_class_num=${#first_sample_class_name[@]} -v second_class_num=${#second_sample_class_name[@]} -v rms_start_column=$rms_start_column -v rms_end_column=$rms_end_column 'NR==1{for (i=rms_start_column;i<=rms_end_column; i++){split($i, tokens, ".sorted.bam.rms");filename[i]=substr(tokens[1],1)".met_level";};next} {for (i=rms_start_column;i<=rms_end_column; i++){print $1"\t"$2"\t"$3"\t"$(i-((first_class_num + second_class_num)*2))"\t"$i> result_dir"/"filename[i]}}' $medips_from_mbd_result_dir"/"$work_file".mr.edgeR.txt" 

			# convert wig to bigWig to visualize on UCSC genome browser
			echo "[INFO] create bigwig for visualization on UCSC genome browser"
			wig_file=$medips_from_mbd_result_dir"/"$work_file".file1.wig"
			tail -n +2 $wig_file > $medips_from_mbd_result_dir"/"${MBD_TYPE_KIND[$k]}".notrack.wig"
			$bin_dir/wigToBigWig $medips_from_mbd_result_dir"/"${MBD_TYPE_KIND[$k]}".notrack.wig" $REF_HUMAN_CHR_SIZE $medips_from_mbd_result_dir"/"${MBD_TYPE_KIND[$k]}".bw"
COMMENT
		fi
	done

else
	>&2 echo "[INFO] Not Specited Option. Only bs or mbd are allowed"
	exit
fi

<<'COMMENT'
	count=0
	# DO NOT USE FOR SUBTYPE PAIR COMPARISON THIS IS SINGLE FILE PAIR COMPARISON
	# run MEDIPS to extract DMRs, ****ALL PAIR-WISE**** 
	for (( i=0; i<${#mbd_list[@]}; i++ )); do
		for (( j=$i+1; j<${#mbd_list[@]}; j++ )); do
			sample_num_list=$i";"$j
			temp_filename_only1=`basename \${mbd_list[$i]}`
			temp_filename_only2=`basename \${mbd_list[$j]}`

			#TODO : Create input list with merging biological replicate by ';'

			$R_DIR/Rscript $bin_dir/medips_pair.R -a $mbd_result_dir"/"$temp_filename_only1".sorted.bam" -b $mbd_result_dir"/"$temp_filename_only2".sorted.bam" -r $medips_from_mbd_result_dir -p $temp_filename_only1"_vs_"$temp_filename_only2 &

			let count+=1
			#[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
			# NOTE : each process takes 30GB memory,so reisrict the # of parallel process to 3 
			[[ $((count%3)) -eq 0 ]] && wait

			# create per sample methyl level file (based on rms)
			if [ $i -eq 0 ]; then
				cut -f1,2,3,6,10 $mbd_result_dir/$temp_filename_only1"_vs_"$temp_filename_only2".mr.edgeR.txt" > $temp_filename_only2".met_level"
			fi
		done

		# create per sample methyl level file (based on rms)
		if [ $i -eq 0]; then
			cut -f1,2,3,5,9 $mbd_result_dir/$temp_filename_only1"_vs_"$temp_filename_only2".mr.edgeR.txt" > $temp_filename_only1".met_level"
		fi

		# convert wig to bigWig to visualize on UCSC genome browser
		temp_filename_only1=`basename \${mbd_list[$i]}`
		temp_filename_only2=`basename \${mbd_list[$i+1]}`

		wig_file=$temp_filename_only1"_vs_"$temp_filename_only2".file1.wig"
		tail -n +2 $wig_file > $temp_filename_only1".notrack.wig"
		wigToBigWig $temp_filename_only1".notrack.wig" $REF_HUMAN_CHR_SIZE $temp_filename_only1".bw"

		last_i=$(( ${#mbd_list[@]} - 1 ))

		if [ $last_i -eq i ]; then
			wig_file=$temp_filename_only1"_vs_"$temp_filename_only2".file2.wig"
			tail -n +2 $wig_file > $temp_filename_only2".notrack.wig"
			wigToBigWig $temp_filename_only2".notrack.wig" $REF_HUMAN_CHR_SIZE $temp_filename_only2".bw"
		fi
	done

	wait

COMMENT
	#TODO : Create sample num list with merging biological replicate by ';'
	#sample_num_list="1;1;2;2;"

	# for testing
		#$R_DIR/Rscript /data/project/mcpg/bin/medips_pair.R -a $sample1_list -b $sample2_list -r $medips_from_mbd_result_dir -p "test"

# copy to final result
cp $methylkit_from_bs_result_dir/*.DMR.bed $final_result_dir

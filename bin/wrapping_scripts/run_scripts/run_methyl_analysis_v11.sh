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

###################################################################################

#######################################
# BS-seq
#######################################
if [ "$sample_format" -eq "3" ]; then # BS

	# 1.1 BS-seq

	# create result directory
	mkdir -p $bs_result_dir
	
	count=0
	#trimming & mapping
	for (( i=0; i<${#bs_r1_list[@]}; i++ )); do
		bash $bin_dir/bs_ver2.sh ${bs_r1_list[$i]} ${bs_r2_list[$i]} $pair_single_flag $bs_result_dir $NUM_CPUS
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
			first_sample_class_list=(); first_sample_class_name=()
			second_sample_class_list=(); second_sample_class_name=()

			for (( i=0; i<${#bs_r1_list[@]}; i++ )); do

				input_fq1_filename_only=`basename \${bs_r1_list[$i]}`
				input_fq1_filename_with_cpg=$methylkit_from_bs_result_dir"/"$input_fq1_filename_only"_CpG.txt"
				
				if [ ${BS_TYPE_LIST[$i]} == ${BS_TYPE_KIND[$k]} ]; then
					first_sample_class_list+=($input_fq1_filename_with_cpg)
					first_sample_class_name+=($input_fq1_filename_only)
				elif [ ${BS_TYPE_LIST[$i]} == ${BS_TYPE_KIND[$l]} ]; then
					second_sample_class_list+=($input_fq1_filename_with_cpg)
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
	} &
	done; wait

	# move to final results
	mv $methylkit_from_bs_result_dir/*CoverageStats.png $methylkit_from_bs_result_dir/*.MethylStats.png $methylkit_from_bs_result_dir/*.corr.png $final_result_dir	
	

	# for testing
	# small_set_dir="/data/project/mcpg/test_data/small_test_data"
	#sh bin/bs.sh $samll_set_dir"/targeted_bs_breast_cancer_80genes/ZR7530_AAAGCA_L001_R1_001.fastq.head100000" $samll_set_dir"/targeted_bs_breast_cancer_80genes/ZR7530_AAAGCA_L001_R2_001.fastq.head100000" pair $result_dir NA
	#sh bin/bs.sh $small_set_dir"/targeted_bs_breast_cancer_80genes/ZR751_GAAACC_L001_R1_001.fastq.head100000" $small_set_dir"/targeted_bs_breast_cancer_80genes/ZR751_GAAACC_L001_R2_001.fastq.head100000" pair $result_dir NA

	# copy to final result
	cp $methylkit_from_bs_result_dir/*.DMR.bed $final_result_dir


#######################################
# MBD-seq
#######################################
#NUM_CPUS=$(( $SYS_NUM_CPUS / 3 ))
elif [ "$sample_format" -eq "0" ]; then
	NUM_CPUS=$(( $SYS_NUM_CPUS / 3 ))
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

		cd $mbd_result_dir
		if [ "$pair_single_flag" == "single" ]; then
			STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $mbd_result_dir/$filename1_wo_ext"_trimmed.fq" --runThreadN $SYS_NUM_CPUS --alignIntronMax 1 --outFileNamePrefix $mbd_result_dir/$input_fq1_filename_only"." --alignEndsType EndToEnd --seedSearchStartLmax 20 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 # --outFilterMatchNmin 15 # --alignEndsType EndToEnd --seedSearchStartLmax 30 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 options are for increasing sensitivity for less than 50bp reads. MBDseq ICBP 36bp length 

		else # for paired-end
			#variables
			input_fq2_filename_only=`basename \${mbd_r2_list[$i]}`
			filename2_wo_ext=`basename $input_fq2_filename_only "."$file_extention`

			STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $mbd_result_dir/$filename1_wo_ext"_val_1.fq" $mbd_result_dir/$filename2_wo_ext"_val_2.fq" --runThreadN $SYS_NUM_CPUS --alignIntronMax 1 --outFileNamePrefix $mbd_result_dir/$input_fq1_filename_only"." --alignEndsType EndToEnd --seedSearchStartLmax 20 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 
		fi
		cd $WORK_DIR
	done



	echo "[INFO] Convert SAM to Sorted BAM"
	count=0
	for (( i=0; i<${#mbd_r1_list[@]}; i++ )); do
		input_fq1_filename_only=`basename \${mbd_r1_list[$i]}`
		# covert sam to sorted bam
		samtools view -ubS -@ $SYS_NUM_CPUS $mbd_result_dir/$input_fq1_filename_only".Aligned.out.sam" | samtools sort -@ $SYS_NUM_CPUS -m 10G - $mbd_result_dir/$input_fq1_filename_only".sorted"
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done;wait

	echo "[INFO] Index sorted BAM"
	count=0
	for (( i=0; i<${#mbd_r1_list[@]}; i++ )); do
		# covert sam to sorted bam
		input_fq1_filename_only=`basename \${mbd_r1_list[$i]}`
		samtools index $mbd_result_dir/$input_fq1_filename_only".sorted.bam" &
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done;wait


# for testing
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
			work_file_with_dir=$medips_from_mbd_result_dir"/"${MBD_TYPE_KIND[$k]}"_vs_"${MBD_TYPE_KIND[$l]}

			# init list
			first_sample_class_name=()
			second_sample_class_name=()
			first_sample_class_filename_only=()
			second_sample_class_filename_only=()

			# TODO : create strand specific BAM file since MEDIPS do not care about strand, divide bam to strand specific bam and run MEDIPS and merge with strand information. However, the medips dmr extract process will be considered as well.
			# create sample list per subtype
			for (( i=0; i<${#mbd_r1_list[@]}; i++ )); do
				if [ ${MBD_TYPE_LIST[$i]} == ${MBD_TYPE_KIND[$k]} ]; then
					first_sample_class_name+=($mbd_result_dir"/"`basename \${mbd_r1_list[$i]}`".sorted.bam")
					first_sample_class_filename_only+=(`basename \${mbd_r1_list[$i]}`)
				elif [ ${MBD_TYPE_LIST[$i]} == ${MBD_TYPE_KIND[$l]} ]; then
					second_sample_class_name+=($mbd_result_dir"/"`basename \${mbd_r1_list[$i]}`".sorted.bam")
					second_sample_class_filename_only+=(`basename \${mbd_r1_list[$i]}`)
				fi 
			done

			echo "[INFO] Input1 sample lists are : ${first_sample_class_name[@]}"
			echo "[INFO] Input2 sample lists are : ${second_sample_class_name[@]}"

			# run medips between 2 classes with replicate information
			# NOTE : For each process, it takes 10GB, so check memory availability before assign number of cpu.

			# TODO : Parallelize MEDIPS.meth() function to speed up
			$R_DIR/Rscript $bin_dir/medips_pair_v2.R -a `my_join ";" ${first_sample_class_name[@]}` -b `my_join ";" ${second_sample_class_name[@]}` -r $medips_from_mbd_result_dir -p $work_file --pvalue $P_VALUE_CUT 
		
			# Manual sorting by lexico order
			medips_result_file=$work_file_with_dir".mr.edgeR.txt"
			medips_sig_result_file=$work_file_with_dir".mr.edgeR.s.txt"
			sorted_medips_result=$work_file_with_dir".mr.edgeR.sorted.txt"

			head -n1 $medips_result_file > $sorted_medips_result

			# generate sort list
			sort_string=()
			for chr in "1" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2" "20" "21" "22" "3" "4" "5" "6" "7" "8" "9" "X" "Y"; do
				sort_string+=($medips_from_mbd_result_dir"/M_chr"$chr)				
			done

			# manual sort. NOTE: Each chr is required to be sorted within the chr!!!
			awk -v output_dir=$medips_from_mbd_result_dir '{print >> output_dir"/M_"$1}' $medips_result_file
			cat `my_join " " ${sort_string[@]}` >> $sorted_medips_result

			# remove temp files
			rm -rf `my_join " " ${sort_string[@]}`

			# DMR file format setting (This is edgeR.s.txt!!! not edgeR.txt)
			tail -n+2 $medips_sig_result_file | cut -f1-3 > $work_file_with_dir".DMR.bed"		
	
			# create per sample methyl level file (based on rms)
			# get rms start end colums
			rms_start_column=$(( 5 + ( (${#first_sample_class_name[@]} + ${#second_sample_class_name[@]}) * 2 ) + ${#first_sample_class_name[@]} ))
			rms_end_column=$(( $rms_start_column + ${#second_sample_class_name[@]} -1 ))

			if [ $k -eq 0 ]; then
				echo "[INFO] Extract methylation level per sample"
				# NOTE : Since MEDIPS put X if the file name start with number, remove X from header
				filename_string=`my_join ";" ${second_sample_class_filename_only[@]}`
				awk -v f_str=$filename_string -v result_dir=$medips_from_mbd_result_dir -v first_class_num=${#first_sample_class_name[@]} -v second_class_num=${#second_sample_class_name[@]} -v rms_start_column=$rms_start_column -v rms_end_column=$rms_end_column 'NR==1{split(f_str, f_list, ";");print "chr\tbin_start\tbin_end\tmethyl_count\tmethyl_level";next;}{OFS="\t";j=1;for (i=rms_start_column;i<=rms_end_column; i++){print $1,$2,$3,$(i-((first_class_num + second_class_num)*2)),$i> result_dir"/"f_list[j]".met_level"; j=j+1}}' $sorted_medips_result
			fi
		done

		# for first class
		rms_start_column=$(( 5 + ( (${#first_sample_class_name[@]} + ${#second_sample_class_name[@]}) * 2 ) ))
		rms_end_column=$(( $rms_start_column + ${#first_sample_class_name[@]} -1 ))

		if [ $k -eq 0 ]; then
			echo "[INFO] Extract methylation level per sample"
			filename_string=`my_join ";" ${first_sample_class_filename_only[@]}`
			awk -v f_str=$filename_string -v result_dir=$medips_from_mbd_result_dir -v first_class_num=${#first_sample_class_name[@]} -v second_class_num=${#second_sample_class_name[@]} -v rms_start_column=$rms_start_column -v rms_end_column=$rms_end_column 'NR==1{split(f_str, f_list, ";");print "chr\tbin_start\tbin_end\tmethyl_count\tmethyl_level";next;}{OFS="\t";j=1;for (i=rms_start_column;i<=rms_end_column; i++){print $1,$2,$3,$(i-((first_class_num + second_class_num)*2)),$i> result_dir"/"f_list[j]".met_level"; j=j+1}}' $sorted_medips_result
		fi
	done

	# copy to final result
	cp $medips_from_mbd_result_dir/*.DMR.bed $final_result_dir

else
	echo "[INFO] Not Specited Option. Only bs or mbd are allowed"
	exit
fi


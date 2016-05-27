#!/bin/bash
source `dirname $0`/../../env.sh

# parse experiment information
source `dirname $0`/parse_exp_info.sh


NUM_CPUS=$SYS_NUM_CPUS

# directories
bin_dir="$WORK_DIR/bin"
#<<'COMMENT'
#result_dir="$WORK_DIR/result"
#COMMENT
#profile_result_dir="$WORK_DIR/profile"
integ_result_dir=$result_dir"/integrated"
mbd_result_dir=$result_dir"/methyl/mbd"
bs_result_dir=$result_dir"/methyl/bs"
methylkit_from_bs_result_dir=$result_dir"/methyl/bs/methylkit"
medips_from_mbd_result_dir=$result_dir"/methyl/mbd/medips"
cel_result_dir=$result_dir"/gene_exp/cel"
rnaseq_result_dir=$result_dir"/gene_exp/rna_seq/deg/"
mu_dnaseq_result_dir=$result_dir"/mu/dna_seq/"
mu_rnaseq_result_dir=$result_dir"/mu/rna_seq/"
runDir_2pass=$rnaseq_result_dir"/2pass"
viz_dir=$result_dir"/visualization"
final_result_dir=$final_result_root_dir/integrated
final_result_dir_ge=$final_result_root_dir/gene_expression
final_result_dir_me=$final_result_root_dir/methylation
final_result_dir_mu=$final_result_root_dir/mutation
mkdir -p $final_result_dir
mkdir -p $final_result_dir_ge
mkdir -p $final_result_dir_me
mkdir -p $final_result_dir_mu

# directory based variables
ME_all_merged=$integ_result_dir/ME_MGD		# prepare initiall prefix for all merged file
ME_GE_all_merged=$integ_result_dir/ME_GE_MGD
ME_GE_COR_all_merged=$integ_result_dir/ME_GE_COR_MGD
ME_GE_COR_all_merged_THRESH=$integ_result_dir/ME_GE_COR_MGD_THRESH.txt
#REF_HUMAN_GENE_RNAGE_INFO=$integ_result_dir/REF_HUMAN_GENE_RNAGE_INFO.txt


# variables
window_size=100 
TSS_range_for_methyl=2000
TSS_range_for_snp=250000
gene_exp_fold_chage=2
corr_threshold=0.3
kw_threshold=$P_VALUE_CUT
no_replicate=0

mbdseq_bs_corr=1

unzip=0

# functions
function my_join { local IFS="$1"; shift; echo "$*"; }

ge_list=()
me_list=()
mu_list=()

###################
#######################
# Gene Expreession Data Parse
ge_type_kind=($(echo ${ge_type_list[@]} | tr ' ' '\n' | sort | uniq | tr '\n' ' '))

# assign parsed values to each variable
cel_list_temp=(${ge_array_list[@]})
type_kind=(${ge_type_kind[@]})
type_list=(${ge_type_list[@]})
sample_fq1_list_temp=(${ge_sample_fq1_list[@]})
sample_fq2_list_temp=(${ge_sample_fq2_list[@]})
is_array=$ge_is_array
single_pair=$ge_single_pair
replicate_exist=$ge_replicate_exist

source `dirname $0`/preprocess_ge_input.sh

# GE Data type
ge_data_type=$is_array

# GE Data list
if [ "$is_array" -eq "1" ]; then
	ge_list=(${cel_list[@]})
else
	ge_list=(${sample_fq1_list[@]})
fi

# GE result dir
ge_result_dir=""
if [ "$ge_data_type" -eq "1" ]; then
	ge_result_dir=$cel_result_dir
else
	ge_result_dir=$rnaseq_result_dir
fi



#######################
# Methylation Data Parse
########################
# check BS or MBD or ....
me_type_kind=($(echo ${ge_type_list[@]} | tr ' ' '\n' | sort | uniq | tr '\n' ' '))

# assign parsed values to each variable
cel_list_temp=(${me_array_list[@]})
sample_fq1_list_temp=(${me_sample_fq1_list[@]})
sample_fq2_list_temp=(${me_sample_fq2_list[@]})
sample_format=$me_sample_format	# 0: MBD-Seq / MeDIP-Seq, 1: "Infinium 27k", 2:"Infinium 450k", 3:"BS-seq"
pair_single_flag=$me_single_pair

source `dirname $0`/preprocess_me_input.sh


# ME Data type
me_data_type=$sample_format

# ME Data list
if [ "$sample_format" -eq "1" ] || [ "$sample_format" -eq "2" ]; then
	me_list=(${cel_list[@]})
else
	me_list=(${sample_fq1_list[@]})
fi

# ME result dir
me_result_dir=""
if [ "$me_data_type" -eq "0" ]; then
	me_result_dir=$medips_from_mbd_result_dir
elif [ "$me_data_type" -eq "3" ]; then
	me_result_dir=$methylkit_from_bs_result_dir
else
	me_result_dir=""
fi


#######################
# Mutation Data Parse
########################
# type kind
mu_type_kind=($(echo ${mu_type_list[@]} | tr ' ' '\n' | sort | uniq | tr '\n' ' '))


# assign parsed values to each variable
cel_list_temp=(${mu_array_list[@]})
sample_fq1_list_temp=(${mu_sample_fq1_list[@]})
sample_fq2_list_temp=(${mu_sample_fq2_list[@]})
is_array=$mu_is_array
single_pair=$mu_single_pair
data_type=$mu_sample_format					# 0:DNAseq, 1:RNAseq, 2:HumanHap

source `dirname $0`/preprocess_mu_input.sh
# MU Data type
mu_data_type=$data_type

# MU Data list
if [ "$mu_data_type" -eq "0" ]; then # 0 : DNA-seq
	mu_list=(${sample_fq1_list[@]})
elif [ "$mu_data_type" -eq "1" ]; then # 1 : RNA-seq
	mu_list=(${sample_fq1_list[@]})
else # 2 : HumanMap 550
	mu_list=(${cel_list[@]})
fi



# MU result dir
mu_result_dir=""

if [ "$mu_data_type" -eq "2" ]; then # Human Map
	mu_result_dir=""
elif [ "$mu_data_type" -eq "1" ]; then # RNA-seq
	mu_result_dir=$runDir_2pass
else #DNA-seq
	mu_result_dir=$mu_dnaseq_result_dir
fi

##################


####################################################################################
# 0. set each of data type 
####################################################################################
# NOTE : THIS SHOULD BE SAME FOR ALL ME, GE, MU
type_kind=(${ge_type_kind[@]})
type_list=(${ge_type_list[@]})
type_len=($(echo ${ge_type_list[@]} | tr ' ' '\n' | sort | uniq -c | awk '{print $1}'| tr '\n' ' '))

# methyl
#me_data_type=1			# 1:MBD, 2:BS, 3:infinium

# gene express
#ge_type=1
# mutation
#mu_type=2					# 1:RNAseq, 2:DNAseq

# input class
#type_kind=("lu" "baa" "bab")
#type_list=("lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "baa" "baa" "baa" "baa" "baa" "baa" "baa" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab")
#TYPE_LEN=(13 7 10)

# multiple correciton
multiple_correction="fdr_bh" #'b', 'fdr_bh', 'hs' and so on. see statsmodels.stats.multitest.multipletests

# pearson or spearman
pear_sp_cor=0 # 0 : pearson / 1: spearman


##############################################################


# TODO
# check all the data type follows same format 

####################################################################################
# 4. Profiling : DEG, DMuG, DMeG association
####################################################################################
	# create result directory
	mkdir -p $integ_result_dir

#################################
# 4.1 Within Subtype
#################################
# 4.1.1 ME vs GE : Correlation, % over threshold, % of NC and PC : 1. TSS +- range, 2. genebody, 3. CpG island, shore, shelf, intergenic, intragenic
	
	# TODO : subtype iteration (for package)	

	# mbd methyl level file : methyl/mbd/medips/100730_s_1.fq_vs_100812_s_1.fq.mr.edgeR.txt
	# 200037_s_at     CBX3    8.49603579205156

	# merge all met_level file	
	paste_list=()
#	for (( k=0; k<${#type_kind[@]}; k++ )); do 
	header_list=()
#		paste_sub_list=()
		for (( i=0; i<${#me_list[@]}; i++ )); do
			temp_filename_only=`basename \${me_list[$i]}`
			met_level_file=$me_result_dir/$temp_filename_only".met_level"
			met_level_file_value_only=$integ_result_dir//$temp_filename_only".only.met_level"

#			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				cut -f1,2,3,5 $met_level_file | tail -n+2 > $met_level_file_value_only 

				paste_list+=($met_level_file_value_only)
				header_list+=(`basename \${me_list[$i]}`"-methyl")
#				paste_sub_list+=($met_level_file_value_only)
#			fi
		done
		# union all the region and merge values for subtype with filling 0 for empty
#		bedtools unionbedg -i `my_join " " ${paste_sub_list[@]}` | bedtools intersect -wa -wb -a - -b $REF_HUMAN_GENE_RANGE_INFO_GROUP_BY_REFSEQ > $ME_all_merged"_"${type_kind[$k]}
#	done

	# union all the region and merge values with filling 0 for empty
#	header_string=$(IFS=$'\t'; echo "${header_list[*]}")
	header_string=`my_join "/" ${header_list[@]}`

	echo -e "#chr\tbin_start\tbin_end\t$header_string\trange_chr\trange_start\trange_end\tstrand\tgene_symbol\trange_kind\trefseq" | tr '/' '\t' > $ME_all_merged".txt"
	bedtools unionbedg -i `my_join " " ${paste_list[@]}` | bedtools intersect -wa -wb -a - -b $REF_HUMAN_GENE_RANGE_INFO_GROUP_BY_REFSEQ >> $ME_all_merged".txt"

	sample_pos_list=()
	temp_sum=4
	
	sample_pos_list+=($temp_sum)

	for ((i=0;i<${#type_len[@]}; i++)); do
		let temp_sum+=${type_len[$i]}
		sample_pos_list+=($temp_sum)
	done

<<'COMMENT'
	for ((k=0;k<${#type_kind[@]};k++)); do
		for ((l=$k+1;l<${#type_kind[@]};l++)); do
			cut -f1-3,${sample_pos_list[$k]}-$((${sample_pos_list[$k]}+${type_len[$k]}-1)),${sample_pos_list[$l]}-$((${sample_pos_list[$l]}+${type_len[$l]}-1)),${sample_pos_list[${#type_kind[@]}]}- > 
COMMENT

	for ((k=0;k<${#type_kind[@]};k++)); do	
		 cut -f1-3,${sample_pos_list[$k]}-$((${sample_pos_list[$k]}+${type_len[$k]}-1)),${sample_pos_list[${#type_kind[@]}]}-  $ME_all_merged".txt" > $ME_all_merged"_"${type_kind[$k]}
	done

<<'COMMENT'

	#################################################################
	# extract only bins in range & merge refseq and remove reseq with no gene symbol
	#################################################################
	count=0
	for (( i=0; i<${#me_list[@]}; i++ )); do
		temp_filename_only=`basename \${me_list[$i]}`
		met_level_file=$me_result_dir/$temp_filename_only".met_level"

		echo -e "#chr\tbin_start\tbin_end\tmethyl_count\tmethyl_level\trange_start\trange_end\tstrand\tgene_symbol\trange_kind\trefseq" > $integ_result_dir/$temp_filename_only".RNGs.MGDrefS"
		awk 'NR==1{next;}{print}' $met_level_file | bedtools intersect -wa -wb -a - -b $REF_HUMAN_GENE_RANGE_INFO  | grep -v "n/a"  | bedtools groupby -g 1,2,3,4,5,7,8,9,11,12 -c 10 -o collapse >> $integ_result_dir/$temp_filename_only".RNGs.MGDrefS" &
		let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait
	# result from above
	# chr	bin_start	bin_end	count	rms	TSS+-range/TSS	TSS+-range/TSS	strand	GeneSymbol	range_kind	refseq
	# chr1    10301   10400   31      0      10373   15909   +	DDX11L1,DDX11L2	TSS_range	NR_046018,NR_046019  
	#################################################################
	# merge all methylation result to one file
	#################################################################

	## prefix
	cut -f1-3 $integ_result_dir"/"`basename \${me_list[0]}`".RNGs.MGDrefS" > $integ_result_dir/temp_prefix.txt &
	cut -f6- $integ_result_dir"/"`basename \${me_list[0]}`".RNGs.MGDrefS" > $integ_result_dir/temp_pos.txt &
	wait


	paste_all_list=()
	for (( k=0; k<${#type_kind[@]}; k++ )); do 

		## rms per file
		paste_list=(); count=0
		for (( i=0; i<${#me_list[@]}; i++ )); do
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				work_file=$integ_result_dir/`basename \${me_list[$i]}`.methyl

				awk -v header=`basename \${me_list[$i]}`"-methyl" 'NR==1{print header;next;}{print $5}' $integ_result_dir"/"`basename \${me_list[$i]}`".RNGs.MGDrefS" > $work_file &
				paste_list+=($work_file); paste_all_list+=($work_file)

				let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
			fi
		done; wait

		## paste all
 		paste $integ_result_dir/temp_prefix.txt ${paste_list[@]} $integ_result_dir/temp_pos.txt > $ME_all_merged"_"${type_kind[$k]}
	done

	# for all merged file
 	paste $integ_result_dir/temp_prefix.txt ${paste_all_list[@]} $integ_result_dir/temp_pos.txt > $ME_all_merged".txt"

COMMENT

	#################################################################
	# merge gene expression to merged methyl file
	#################################################################

	## merge each sample's avg GE to pre_merge_GE_all file
	paste_all_list=()
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		
		paste_list=();count=0
		for (( i=0; i<${#ge_list[@]}; i=i+1 )); do
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				work_file=$integ_result_dir/`basename \${ge_list[$i]}`".AVG_GE"

				echo `basename \${ge_list[$i]}`"-gene_exp" > $work_file
			
				temp_exp_input=""

				if [ "$ge_data_type" -eq "1" ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
					temp_exp_input=$ge_result_dir"/"`basename \${ge_list[\$i]}`".exp.geneSymbol_avg_exp.txt"
				else #Li1_head_100000_1.htseq.geneSymbol
					file_extension=$(echo `basename \${ge_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
					temp_exp_input=$ge_result_dir"/"`basename \${ge_list[0]} "."$file_extension`".htseq.geneSymbol"
				fi
				
				# merge exp to based on methyl merged file gene order
				awk -v input_file=$temp_exp_input -f $bin_dir/append_avg_GE.awk $ME_all_merged"_"${type_kind[$k]} >> $work_file &
				paste_list+=($work_file);paste_all_list+=($work_file)

				let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
			fi
		done; wait

		## merge ME GE
		paste $ME_all_merged"_"${type_kind[$k]} ${paste_list[@]} > $ME_GE_all_merged"_"${type_kind[$k]}
	done

	# for all merged file
	paste $ME_all_merged".txt" ${paste_all_list[@]} > $ME_GE_all_merged".txt"

#################################################################
	# compute correlation between each bin and gene expression
	#################################################################

	# FOR PACKAGE, THIS IS NOT USEFUL, so commented
<<'COMMENT'
	echo "[INFO] compute correlation"
	for (( k=0; k<${#type_kind[@]}; k++ )); do 

		work_file=$ME_GE_all_merged"_"${type_kind[$k]}
		work_file2=$ME_GE_COR_all_merged"_"${type_kind[$k]}

		echo "[INFO] pearson/spearman correlation for subtype "${type_kind[$k]}
		echo "[INFO] subtype length : "${type_len[$k]}

		head -n1 $work_file | xargs -I {} echo -e "{}\tcorrelation" > $work_file2

		# Pearson
		# awk 'NR==1{next;}{print $0}' $work_file | parallel --no-notice --pipe -j"$NUM_CPUS" -L1000000 -k python $bin_dir/pearson.py - 3 $((9 + ${type_len[$k]} )) ${type_len[$k]} >> $work_file2

		# Spearman : Score is bad, so not use currently
		#Rscript p_spearman.r -i test.txt -n 3 --s1 2 --s2 8 -p 1 -o ./test_result.txt 

		# compute % of over correlation threshold
		#echo "[INFO] Compute % of correlation over threshold"
		#python $bin_dir/profile_ME_GE.py $work_file2 $integ_result_dir"/ME_GE_stat_"${type_kind[$k]}".txt" ${type_len[$k]} $corr_threshold $integ_result_dir/ME_GE_MAT_THRESH_${type_kind[$k]}.txt $integ_result_dir/ME_GE_MAT_PNC_${type_kind[$k]}.txt
	done
COMMENT

# 4.1.2 SNP vs GE : TSS, TSE +- range, Correlation(Spearman's rho), % over threshold, % of NC and PC
# TODO: Change as ME vs GE to speed up
	# extract only snps in range
	count=0
	
	snp_file_extension=$(echo `basename \${mu_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')

	paste_list=()
	header_list=()
	header_list2=()

	for (( i=0; i<${#mu_list[@]}; i++ )); do
		temp_filename_only=`basename \${mu_list[\$i]}`
    temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`
    snp_level_file=""
    if [ "$mu_data_type" -eq "0" ] || [ "$mu_data_type" -eq "1" ]; then
      snp_level_file=$mu_result_dir/$temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf.bed"
    else
      snp_level_file=""
    fi

	  snp_level_af_file=$integ_result_dir"/"$temp_filename_only".snp"
		paste_list+=($snp_level_af_file)
		header_list+=(`basename \${mu_list[$i]}`"-snp")
		header_list2+=(`basename \${mu_list[$i]}`"-AF")

 		
		awk 'BEGIN{FS=OFS="\t"}{num=split($6,value_arr,";");num2=split(value_arr[2],af_arr,"="); af_val=af_arr[2]; print $1,$2,$3,$4"|"af_val}' $snp_level_file > $snp_level_af_file &
		let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait

	done; wait


	header_string=`my_join "/" ${header_list[@]}`
	header_string2=`my_join "/" ${header_list2[@]}`

	echo -e "#chr\tbin_start\tbin_end\t$header_string\t$header_string2\trange_chr\trange_start\trange_end\tstrand\tgene_symbol\trange_kind\trefseq" | tr '/' '\t' > $integ_result_dir/SNP_MGD.txt
	
	bedtools unionbedg -i `my_join " " ${paste_list[@]}` -filler "-" | awk -v sample_num=${#type_list[@]} 'BEGIN{FS=OFS="\t"}{for(i=4;i<(4+sample_num);i++){if( $(i) == "-"){chr_arr[i]="-/-";af_arr[i]=0}else{num=split($(i),arr,"|");chr_arr[i]=arr[1];af_arr[i]=arr[2]}}; printf "%s\t%s\t%s", $1,$2,$3; for(i=4;i<(4+sample_num);i++){printf "\t%s", chr_arr[i]}; for(i=4;i<(4+sample_num);i++){printf "\t%s", af_arr[i]};printf "\n"}' |bedtools intersect -wa -wb -a - -b $REF_HUMAN_GENE_RANGE_INFO_GROUP_BY_REFSEQ >> $integ_result_dir/SNP_MGD.txt



 # for class(subtype) file
  for ((k=0;k<${#type_kind[@]};k++)); do
     cut -f1-3,${sample_pos_list[$k]}-$((${sample_pos_list[$k]}+${type_len[$k]}-1)),$((${sample_pos_list[$k]}+${#type_list[@]}))-$((${sample_pos_list[$k]}+${type_len[$k]}-1 + ${#type_list[@]})),$((${sample_pos_list[${#type_kind[@]}]}+${#type_list[@]}))- $integ_result_dir/SNP_MGD.txt > $integ_result_dir/SNP_MGD"_"${type_kind[$k]}
	done 


	# merge each sample's avg GE to SNP merged file

	paste_all_list=()

	for (( k=0; k<${#type_kind[@]}; k++ )); do
	
		paste_list=();count=0
		for ((i=0; i<${#ge_list[@]}; i++)); do
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				work_file=$integ_result_dir/`basename \${ge_list[$i]}`".AVG_GE_for_SNP"
		
				temp_exp_input=""

				echo `basename \${ge_list[$i]}`"-gene_exp" > $work_file
		
    		if [ "$ge_data_type" -eq "1" ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
      		 temp_exp_input=$ge_result_dir"/"`basename \${ge_list[\$i]}`".exp.geneSymbol_avg_exp.txt"
    		else #Li1_head_100000_1.htseq.geneSymbol
      		 file_extension=$(echo `basename \${ge_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
      		 temp_exp_input=$ge_result_dir"/"`basename \${ge_list[0]} "."$file_extension`".htseq.geneSymbol"
    		fi

				# merge exp to based on snp merged file gene order
				awk -v input_file=$temp_exp_input -f $bin_dir/append_avg_GE.awk $integ_result_dir/SNP_MGD"_"${type_kind[$k]} >> $work_file &

				paste_list+=($work_file);paste_all_list+=($work_file)

				let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
			fi
		done; wait

		# merge for class file

		paste $integ_result_dir/SNP_MGD"_"${type_kind[$k]} ${paste_list[@]} > $integ_result_dir/SNP_GE_MGD_${type_kind[$k]}.txt
	done

	# for all merged file
	paste $integ_result_dir/SNP_MGD.txt ${paste_all_list[@]} > $integ_result_dir/SNP_GE_MGD.txt


<<'COMMENT'
	for (( i=0; i<${#mu_list[@]}; i++ )); do
		temp_filename_only=`basename \${mu_list[\$i]}`
		temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`
		snp_level_file=""
		if [ "$mu_data_type" -eq "0" ] || [ "$mu_data_type" -eq "1" ]; then
			snp_level_file=$mu_result_dir/$temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf.bed"
		else
			snp_level_file=""
		fi
			
		echo -e "#chr\tbin_start\tbin_end\tref/alt\tsnp_strand\tinfo\trange_start\trange_end\tstrand\tgene_symbol\trange_kind\trefseq" > $integ_result_dir/$temp_filename_only".SNP_RNGs.MGDrefS"
		bedtools intersect -wa -wb -a $snp_level_file -b $REF_HUMAN_GENE_RANGE_INFO | grep -v "n/a" | bedtools groupby -g 1,2,3,4,5,6,8,9,10,12,13 -c 11 -o collapse >> $integ_result_dir/$temp_filename_only".SNP_RNGs.MGDrefS" &

		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait

	# merge refseq and remove reseq with no gene symbol
	echo "[INFO] Merging same subtype SNP data"
	count=0; SNP_MGD_LIST_all=(); GE_MGD_LIST_all=()
	for (( k=0; k<${#type_kind[@]}; k++ )); do 

		# get subtype snp, GE file list
		SNP_MGD_LIST=(); GE_MGD_LIST=();
		for (( i=0; i<${#mu_list[@]}; i++ )); do
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				SNP_MGD_LIST+=($integ_result_dir"/"`basename \${mu_list[$i]}`".SNP_RNGs.MGDrefS")
				SNP_MGD_LIST_all+=($integ_result_dir"/"`basename \${mu_list[$i]}`".SNP_RNGs.MGDrefS")
		
				temp_exp_input=""

				if [ $ge_data_type -eq 1 ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
					temp_exp_input=$ge_result_dir"/"`basename \${ge_list[\$i]}`".exp.geneSymbol_avg_exp.txt"
				else #Li1_head_100000_1.htseq.geneSymbol
					file_extension=$(echo `basename \${ge_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
					temp_exp_input=$ge_result_dir"/"`basename \${ge_list[0]} "."$file_extension`".htseq.geneSymbol"
				fi

				GE_MGD_LIST+=($temp_exp_input)
				GE_MGD_LIST_all+=($temp_exp_input)
			fi
		done

		# result from above: #chr10	108709	108710	a/aAAA	+	INDEL;IS=3,0.600000;DP=5;QS=0.000000,1.000000,0.000000,0.000000;VDB=2.063840e-02;AF1=1;AC1=2;DP4=0,0,3,0;MQ=20;FQ=-43.5	0	345178	-	TUBB8	SNP_TSS_TSE_flanking_range+-250kb	NM_177987
		# gene expression file (.cell file) : # 100730_s_4_export.txt.CEL.exp: # 200000_s_at     PRPF8   9.416984390952

		# class(subtype) merge
		python $bin_dir/merge_snp2.py `my_join ";" ${SNP_MGD_LIST[@]}` `my_join ";" ${GE_MGD_LIST[@]}` $corr_threshold $integ_result_dir/SNP_GE_MGD_${type_kind[$k]}.txt $integ_result_dir/SNP_GE_stat_${type_kind[$k]}.txt $integ_result_dir/SNP_GE_MAT_THRESH_${type_kind[$k]}.txt $integ_result_dir/SNP_GE_MAT_PNC_${type_kind[$k]}.txt &
		let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait

	# for all sample merge
	echo "[INFO] Merging all sample SNP data"
	python $bin_dir/merge_snp2.py `my_join ";" ${SNP_MGD_LIST_all[@]}` `my_join ";" ${GE_MGD_LIST_all[@]}` $corr_threshold $integ_result_dir/SNP_GE_MGD.txt $integ_result_dir/SNP_GE_stat.txt $integ_result_dir/SNP_GE_MAT_THRESH.txt $integ_result_dir/SNP_GE_MAT_PNC.txt
COMMENT

	# FOR PACKAGE, THIS IS NOT USEFUL, so commented
<<'COMMENT'
# 4.1.3 ENCODE/GEO data association
## 4.1.3.1 Intersect with DHS cluster
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
### 4.1.3.1.0 ME-GE bins 
		# 13763956
		#total_DMR_bin_num=`cut -f1,2,3 $ME_GE_COR_all_merged"_"${type_kind[$k]} | uniq | wc -l`
		bedtools intersect -wa -a $ME_GE_COR_all_merged"_"${type_kind[$k]} -b $DHS_CLUSTER_UCSC > $ME_GE_COR_all_merged"_"${type_kind[$k]}"_DHS"

		python $bin_dir/profile_ME_GE.py $ME_GE_COR_all_merged"_"${type_kind[$k]}"_DHS" $integ_result_dir"/ME_GE_stat_"${type_kind[$k]}"_DHS.txt" ${type_len[$k]} $corr_threshold $integ_result_dir/ME_GE_MAT_THRESH_${type_kind[$k]}_DHS.txt $integ_result_dir/ME_GE_MAT_PNC_${type_kind[$k]}_DHS.txt

		
		# draw plot
		echo "[INFO] Draw bar plots"
		Rscript $bin_dir/stat_bar.r $integ_result_dir/ME_GE_MAT_PNC_${type_kind[$k]}_DHS.txt $VIZ_DIR/BAR_ME_GE_MAT_PNC_${type_kind[$k]}"_DHS.jpg" "Among ME GE bins over threshold, PC/NC ratio with GE in ${type_kind[$k]}" 2
		Rscript $bin_dir/stat_bar.r $integ_result_dir/ME_GE_MAT_THRESH_${type_kind[$k]}_DHS.txt $VIZ_DIR/BAR_ME_GE_MAT_THRESH_${type_kind[$k]}"_DHS.jpg" "Ratio of ME GE intersected bins having CORR with GE over threshold in ${type_kind[$k]}" 1

		# copy to WEB_ACCESSIBLE_DIR
		cp $integ_result_dir/ME_GE_stat_${type_kind[$k]}"_DHS.txt" $VIZ_DIR/BAR_ME_GE_MAT_PNC_${type_kind[$k]}"_DHS.jpg" $VIZ_DIR/BAR_ME_GE_MAT_THRESH_${type_kind[$k]}"_DHS.jpg" $WEB_ACCESSIBLE_DIR

		# 3781450
		#DMR_bin_in_DHS_num=`cut -f1,2,3 $ME_GE_COR_all_merged"_"${type_kind[$k]}"_DHS" | uniq | wc -l`
### 4.1.3.1.1 ME-GE bins over threshold 
		# get ME-GE having corr over trheshold	
		#awk -v threshold=$corr_threshold '{if (($NF > threshold) || ($NF < -threshold)) print}' $ME_GE_COR_all_merged"_"${type_kind[$k]} | bedtools intersect -wa -a - -b $DHS_CLUSTER_UCSC > $ME_GE_COR_all_merged"_"${type_kind[$k]}"_OVT_DHS"
		# 2443065
		#OVT_DMR_bin_in_DHS_num=`cut -f1,2,3 $ME_GE_COR_all_merged"_"${type_kind[$k]}"_OVT_DHS" | uniq | wc -l`
	done
### 4.1.3.1.2 SNP-GE bins 
### 4.1.3.1.3 SNP-GE bins over threshold 

## 4.1.3.2 Intersect with Histone data (H3K4me3, H3K27me3)
### 4.1.3.2.0 ME-GE bins 
### 4.1.3.2.1 ME-GE bins over threshold 

### 4.1.3.2.2 SNP-GE bins 
### 4.1.3.2.3 SNP-GE bins over threshold 
COMMENT

#################################
# 4.2 Between Subtype
#################################

# 4.2.1 GDMR vs GDEG

	# 4.2.1.1 get GDMRs and GDEGs

	# example: 100730_s_7.fq_vs_100831_s_6.fq.mr.edgeR.s.gain.txt, use Q-value 0.01 for DMR
	#chr	start	stop	CF	X100730_s_7.fq.sorted.bam.counts	X100831_s_6.fq.sorted.bam.counts	X100730_s_7.fq.sorted.bam.rpkm	X100831_s_6.fq.sorted.bam.rpkm	X100730_s_7.fq.sorted.bam.rms	X100831_s_6.fq.sorted.bam.rms	X100730_s_7.fq.sorted.bam.prob	X100831_s_6.fq.sorted.bam.prob	X100730_s_7.fq.sorted.bam.counts.1	X100831_s_6.fq.sorted.bam.counts.1	X100730_s_7.fq.sorted.bam.rpkm.1	X100831_s_6.fq.sorted.bam.rpkm.1	MSets1.counts.mean	MSets1.rpkm.mean	MSets1.rms.mean	MSets1.prob.mean	MSets2.counts.mean	MSets2.rpkm.mean	MSets2.rms.mean	MSets2.prob.mean	edgeR.logFC	edgeR.logCPM	edgeR.p.value	edgeR.adj.p.value
	#chr1	63401	63500	1	39	5	12.2586761723585	1.29318962375682	4.64035822053151	0.449371294449568	1	0.927899578169853	39	5	12.2586761723585	1.29318962375682	39	12.2586761723585	4.64035822053151	1	5	1.29318962375682	0.449371294449568	0.927899578169853	3.20842545778672	-0.493264166192471	5.28317575583928e-09	0.0673129581546756 
	## rms per file


# all ME file
cut -f1-$((3+${#type_list[@]})) $ME_GE_all_merged".txt" | awk 'BEGIN{FS=OFS="\t"}{print $1"_"$2"_"$3,$0}' | cut --complement -f2-4 | grep -v "NA" | uniq > $integ_result_dir/ME_MAT.txt &

# all GE file
cut -f1-3,$((${#type_list[@]}+11))- $ME_GE_all_merged".txt" | awk 'BEGIN{FS=OFS="\t"}{print $1"_"$2"_"$3,$0}' | cut --complement -f2-4 | grep -v "n\/a" | uniq > $integ_result_dir/GE_MAT.txt &

wait

<<'COMMENT'
	paste_list=()
	paste_list2=()
	for (( i=0; i<${#me_list[@]}; i++ )); do
		work_file=$integ_result_dir/`basename \${me_list[$i]}`.methyl
		work_file2=$integ_result_dir/`basename \${ge_list[$i]}`.AVG_GE
		paste_list+=($work_file)
		paste_list2+=($work_file2)
	done 
	
	awk '{print $1"_"$2"_"$3}' $integ_result_dir/temp_prefix.txt | paste - ${paste_list[@]} | grep -v "NA" | uniq > $integ_result_dir/ME_MAT.txt &
	awk '{print $1"_"$2"_"$3}' $integ_result_dir/temp_prefix.txt | paste - ${paste_list2[@]} | grep -v "n\/a" | uniq > $integ_result_dir/GE_MAT.txt &
	wait

COMMENT
	###########################################
	# BY kruscal & and bonferroni
	###########################################

		python $bin_dir/kruskal_ver4_wo_rpy.py $integ_result_dir/ME_MAT.txt `my_join "," ${type_len[@]}` $kw_threshold $NUM_CPUS $integ_result_dir/ME_MAT.kw.txt $multiple_correction

#-a 0.05 -n $NUM_CPUS -c `my_join "," ${type_list[@]}` -m 0 --o1 $integ_result_dir/ME_MAT.kw.sigle_diff --o2 $integ_result_dir/ME_MAT.kw.all_diff

		python $bin_dir/kruskal_ver4_wo_rpy.py $integ_result_dir/GE_MAT.txt `my_join "," ${type_len[@]}` $kw_threshold $NUM_CPUS $integ_result_dir/GE_MAT.kw.txt $multiple_correction

#-a 0.05 -n $NUM_CPUS -c `my_join "," ${type_list[@]}` -m 0 --o1 $integ_result_dir/GE_MAT.kw.sigle_diff --o2 $integ_result_dir/GE_MAT.kw.all_diff
		
	# adjusted pvalue selection
		awk -v kw_threshold=$kw_threshold '{if(($(NF) < kw_threshold) && ($(NF) != "nan")){split($1, tokens, "_"); print tokens[1]"\t"tokens[2]"\t"tokens[3]}}' $integ_result_dir/ME_MAT.kw.txt  > $integ_result_dir/GDMR_kw &
		awk -v kw_threshold=$kw_threshold '{if(($(NF) < kw_threshold) && ($(NF) != "nan")){split($1, tokens, "_"); print tokens[1]"\t"tokens[2]"\t"tokens[3]}}' $integ_result_dir/GE_MAT.kw.txt  > $integ_result_dir/GDEG_kw &
		wait

		# too small
		#awk '{split($1, tokens, "_"); print tokens[1]"\t"tokens[2]"\t"tokens[3]}' $integ_result_dir/ME_matrix.txt.kw_test.signif.cut0.05.txt > $integ_result_dir/GDMR_kw
		#awk '{split($1, tokens, "_"); print tokens[1]"\t"tokens[2]"\t"tokens[3]}' $integ_result_dir/GE_matrix.txt.kw_test.signif.cut0.05.txt > $integ_result_dir/GDEG_kw

	###########################################
	# BY Entropy
	###########################################

	# Select GDEG, GDMR from various, THIS ONLY CONTAINS BED REGIONS	
	GDEG=$integ_result_dir/GDEG_kw
	GDMR=$integ_result_dir/GDMR_kw

	# compute % of overlap between GDMR and GDEG and their correlation by bedtools intersection with $ME_GE_COR_all_merged
	num_GDMR_bins=`wc $GDMR | awk '{print $1}'`
	# NOTE : name depandancy on "gene_symbol" in the header section , NOTE2 : -wa -wb is not working, so just use -wb -wa and cut -f4-
	temp_genesymbol_pos_file=$integ_result_dir/temp_chrom_pos_genesymbol.bed

	awk 'NR==1{split($0, tokens, "\t");for(i in tokens){if(tokens[i]=="gene_symbol") gene_symbol_index=i;};next;}{print $1"\t"$2"\t"$3"\t"$gene_symbol_index}' $integ_result_dir/ME_GE_MGD_${type_kind[0]} > $temp_genesymbol_pos_file
 
	num_GDEG=`bedtools intersect -wb -wa -a $GDEG -b $temp_genesymbol_pos_file | cut -f4-| tee $GDEG"_GSym" |cut -f4 | sort | uniq | tee $GDEG"_genelist" | wc | awk '{print $1}'`
	
	bedtools intersect -wa -wb -a $GDMR -b $GDEG"_GSym" | cut -f4- | sort -k1,1 -k2,2n -k4,4 | uniq > $integ_result_dir/GDMR_GDEG_w_GSym.txt

	GDEG_GDMR_MGD=$integ_result_dir/GDMR_GDEG_w_GSym.txt.all_merged_data.txt
	GDEG_GDMR_MGD_CORR=$GDEG_GDMR_MGD".corr.txt"
	bedtools intersect -wa -a $ME_GE_all_merged".txt" -b $integ_result_dir/GDMR_GDEG_w_GSym.txt > $GDEG_GDMR_MGD

	rm -rf $GDEG_GDMR_MGD_CORR

	# Compute correaltion
	cat $GDEG_GDMR_MGD |  parallel --no-notice --pipe -j"$NUM_CPUS" -L1000000 -k python $bin_dir/pearson.py - 3 $((9+${#type_list[@]})) ${#type_list[@]} $pear_sp_cor >> $GDEG_GDMR_MGD_CORR

	
	# compute % of over correlation threshold
	echo "[INFO] Compute % of correlation over threshold - GDMR & GDEG"
	python $bin_dir/profile_ME_GE.py $GDEG_GDMR_MGD_CORR $GDEG_GDMR_MGD_CORR"_stat.txt" ${#type_list[@]} $corr_threshold $integ_result_dir/GDMR_GDEG_MAT_THRESH.txt $integ_result_dir/GDMR_GDEG_MAT_PNC.txt

	
	# GDEG exp value table
	GDEG_exp_value_table=$integ_result_dir/ALL.DEG.table.txt

	ge_basename_list=()
	for ((i=0; i<${#ge_list[@]}; i++)); do
		ge_basename_list+=(`basename \${ge_list[\$i]}`)
	done

#	ge_basename_list_string=$(IFS=$'\t'; echo "${ge_basename_list[*]}")
	ge_basename_list_string=`my_join "/" ${ge_basename_list[@]}`

 	
	echo -e "GeneSymbol\t"$ge_basename_list_string"\tP.val\tadj.P.val" | tr '/' '\t' > $integ_result_dir/ALL.Total.table.txt

	awk 'NR>1{if($(NF) != "nan"){split($1, tokens, "_"); print tokens[1]"\t"tokens[2]"\t"tokens[3]"\t",$0}}' $integ_result_dir/GE_MAT.kw.txt | bedtools intersect -wb -wa -a $temp_genesymbol_pos_file -b - | cut --complement -f1-3,5-8 | sort -k1,1 | uniq | sort -k$((3+${#type_list[@]})),$((3+${#type_list[@]}))g | awk 'BEGIN{FS=OFS="\t"}{printf "%s", $1; for(i=2;i<=NF;i++){printf "\t%.4g", $(i)};printf "\n"}' >> $integ_result_dir/ALL.Total.table.txt
	
	awk -v kw_threshold=$kw_threshold 'NR==1{print $0;next}{if($(NF) < kw_threshold){print $0}}' $integ_result_dir/ALL.Total.table.txt > $GDEG_exp_value_table

 # top 30 gdeg for oncoprint
	awk 'BEGIN{FS=OFS="\t";i=0}NR>1{if(i<30){if(!($1 in arr) && ($1 != "NA")){i=i+1;arr[$1]=""}}else{exit}}END{for(x in arr){print x}}' $GDEG_exp_value_table > $integ_result_dir/ALL.DEG.Top30.genelist

	# top 100 gdegs
	head -n 101 $GDEG_exp_value_table > $integ_result_dir/ALL.DEG.Top100.table.txt 

	# cp to result dir
	cp $integ_result_dir/ALL.Total.table.txt $GDEG_exp_value_table $integ_result_dir/ALL.DEG.Top100.table.txt $final_result_dir_ge

<<'COMMENT'
	num_GDMR_bins_related_w_GDEG=`cut -f1,2,3 $integ_result_dir/GDMR_GDEG_w_GSym.txt | uniq  | wc | awk '{print $1}'`
	num_GDEGs_related_w_GDMR=`cut -f4 $integ_result_dir/GDMR_GDEG_w_GSym.txt | sort | uniq | tee $integ_result_dir/GDEGs_that_have_GDMR_in_genebody.txt | wc | awk '{print $1}'`

	{
	echo "Number of GDMR bins are : "$num_GDMR_bins;
	echo "Number of GDEGs are : "$num_GDEG;
	echo "Number of GDMR bins related with GDEG(genebody)s are : "$num_GDMR_bins_related_w_GDEG;
	echo "Number of GDEG(genebody)s related with GDMRs are : "$num_GDEGs_related_w_GDMR;

	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		num_GDEGs_that_have_GDMR_in_promoter=`grep promoter $ME_GE_COR_all_merged"_"${type_kind[$k]} | bedtools intersect -wa -wb -a $integ_result_dir/GDMR_GDEG_w_GSym.txt -b - | tee $integ_result_dir/GDEGs_that_have_GDMR_in_promoter"_"${type_kind[$k]}.txt | cut -f4 | sort | uniq | wc | awk '{print $1}'`
		echo "Number of GDEGs that has GDMR in promoter in ${type_kind[$k]}: "$num_GDEGs_that_have_GDMR_in_promoter;
	done
	echo "";
	} > $integ_result_dir/GDMR_GDEG_stat.txt
COMMENT

<<'COMMENT'
# 4.2.2 GDSNP vs GDEG 
	# 4.2.2.1 get GDSNPs

	# TODO : Performance improve. Currently GDSNP.txt contains everything -> reduce size to speed up.
	echo "[INFO] Get GDSNP"
	python $bin_dir/get_GDSNP.py $integ_result_dir/SNP_GE_MGD.txt ${#type_kind[@]} `my_join ";" ${type_list[@]}` $integ_result_dir/GDSNP.bed > $integ_result_dir/GDSNP.txt


	bedtools intersect -wa -wb -b $integ_result_dir/GDSNP.txt -a $GDEG"_GSym" > $integ_result_dir/GDEG_GDSNP.txt
	num_GDSNP=`cut -f1,2,3 $integ_result_dir/GDSNP.txt | sort -k1,1 -k2,2n | uniq | wc | awk '{print $1}'`
	num_GDSNP_bins_related_w_GDEG=`cut -f1,2,3 $integ_result_dir/GDEG_GDSNP.txt | sort -k1,1 -k2,2n | uniq | wc | awk '{print $1}'`
	num_GDEGs_related_w_GDSNP=`cut -f4 $integ_result_dir/GDEG_GDSNP.txt | sort | uniq | tee $integ_result_dir/GDEGs_related_GDSNP.txt | wc | awk '{print $1}'`
	num_GDEGs_that_have_GDSNP_in_promoter=`grep promoter $integ_result_dir/GDEG_GDSNP.txt | cut -f4 - | sort | uniq | tee $integ_result_dir/GDEGs_that_have_GDSNP_in_promoter.txt | wc | awk '{print $1}'`

	{
		echo "Number of GDSNP is : "$num_GDSNP;
		echo "Number of GDEGs are : "$num_GDEG;
		echo "Number of GDSNP bins related with GDEGs are : "$num_GDSNP_bins_related_w_GDEG;
		echo "Number of GDEGs related with GDSNPs are : "$num_GDEGs_related_w_GDSNP;
		echo "Number of GDEGs that has GDSNP in promoter : "$num_GDEGs_that_have_GDSNP_in_promoter;
	} > $integ_result_dir/GDEG_GDSNP_stat.txt

# 4.2.2 GDSNP vs GDEG vs GDMR
	bedtools intersect -wa -wb -b $integ_result_dir/GDEG_GDSNP.txt -a $integ_result_dir/GDMR_GDEG_w_GSym.txt > $integ_result_dir/GDEG_GDSNP_GDMR.txt

	#num_GDMR_bins=`wc $integ_result_dir/GDSNP.txt | awk '{print $1}'`
	num_GDMR_and_GDSNP_bins_related_w_GDEG=`cut -f1,2,3 $integ_result_dir/GDEG_GDSNP_GDMR.txt | uniq | wc | awk '{print $1}'`
	num_GDEGs_related_w_GDSNP_GDMR=`cut -f4 $integ_result_dir/GDEG_GDSNP_GDMR.txt | sort | uniq | tee $integ_result_dir/GDEGs_related_GDSNP_and_GDMR.txt | wc | awk '{print $1}'`
	num_GDEGs_that_have_GDSNP_GDMR_in_promoter=`grep promoter $integ_result_dir/GDEG_GDSNP_GDMR.txt | cut -f4 - | sort | uniq | tee $integ_result_dir/GDEGs_that_have_GDSNP_GDMR_in_promoter.txt | wc | awk '{print $1}'`
	{
		echo "Number of GDMR and GDSNP bins related with GDEGs are : "$num_GDMR_and_GDSNP_bins_related_w_GDEG;
		echo "Number of GDEGs related with GDSNP and GDMR is : "$num_GDEGs_related_w_GDSNP_GDMR;
		echo "Number of GDEGs that has GDSNP and GDMR in promoter : "$num_GDEGs_that_have_GDSNP_GDMR_in_promoter;
	} > $integ_result_dir/GDEG_GDSNP_GDMR_stat.txt
COMMENT


# 4.2.3 DEG vs DMR, pairwise subtype comparison


<<'COMMENT'
	for ((k=0;k<${#type_kind[@]};k++)); do
		for ((l=$k+1;l<${#type_kind[@]};l++)); do
			cut -f1-3,${sample_pos_list[$k]}-$((${sample_pos_list[$k]}+${type_len[$k]}-1)),${sample_pos_list[$l]}-$((${sample_pos_list[$l]}+${type_len[$l]}-1)),${sample_pos_list[${#type_kind[@]}]}- > 
COMMENT


	# Make SUBTYPE FILES
	for (( k=0; k<${#type_kind[@]}; k++)); do
		for (( l=$k+1; l<${#type_kind[@]}; l++)); do
			work_file=$ME_GE_all_merged"_"${type_kind[$k]}"_vs_"${type_kind[$l]}
			work_corr_file=$ME_GE_COR_all_merged"_"${type_kind[$k]}"_vs_"${type_kind[$l]}
		
			cut -f1-3,${sample_pos_list[$k]}-$((${sample_pos_list[$k]}+${type_len[$k]}-1)),${sample_pos_list[$l]}-$((${sample_pos_list[$l]}+${type_len[$l]}-1)),${sample_pos_list[${#type_kind[@]}]}-$((${sample_pos_list[${#type_kind[@]}]}+6)),$((${sample_pos_list[${#type_kind[@]}]}+3+${sample_pos_list[$k]}))-$((${sample_pos_list[${#type_kind[@]}]}+3+${sample_pos_list[$k]}+${type_len[$k]}-1)),$((${sample_pos_list[${#type_kind[@]}]}+3+${sample_pos_list[$l]}))-$((${sample_pos_list[${#type_kind[@]}]}+3+${sample_pos_list[$l]}+${type_len[$l]}-1)) $ME_GE_all_merged".txt" > $work_file
			
<<'COMMENT'
			first_sample_file=$ME_GE_all_merged"_"${type_kind[$k]}
			second_sample_file=$ME_GE_all_merged"_"${type_kind[$l]}

			# make SUBTYPE_DIFF base file
			# chr stat end (TYPE_1 ME) (TYPE_2 ME) TSS+-range/TSS TSS+-range/TSS strand GeneSymbol range_kind refseq (TYPE_1 GE) (TYPE_2 GE)
			paste <(cut -f1-$((${type_len[$k]}+3)) $first_sample_file) <(cut -f4-$((${type_len[$l]}+3)) $second_sample_file) $integ_result_dir/temp_pos.txt <(cut -f$((${type_len[$k]}+10))- $first_sample_file) <(cut -f$((${type_len[$l]}+10))- $second_sample_file) > $work_file

			rm -rf $work_corr_file
COMMENT
		
			echo -n "" > $work_corr_file
			# Compute correaltion
			awk 'NR==1{next;}{print $0}' $work_file | grep -v "n\/a" | parallel --no-notice --pipe -j"$NUM_CPUS" -L1000000 -k python $bin_dir/pearson.py - 3 $((10 + ${type_len[$k]} + ${type_len[$l]} )) $((${type_len[$k]}+${type_len[$l]})) $pear_sp_cor >> $work_corr_file

			# filter only DMR & DEG
			SUBTYPE_DMR=$me_result_dir"/"${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.bed"
			# chr start end
			SUBTYPE_DEG=$ge_result_dir"/"${type_kind[$k]}"_vs_"${type_kind[$l]}".DEG.list"
			# gene_symbol

			work_corr_file_Diff=$work_corr_file".DMR_DEG"
			bedtools intersect -wa -a $work_corr_file -b $SUBTYPE_DMR | awk -v gene_symbol_index=$((8+${type_len[$k]} + ${type_len[$l]})) 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""}FILENAME==ARGV[2]{if($(gene_symbol_index) in arr){print $0}}' $SUBTYPE_DEG - > $work_corr_file_Diff
#| grep -w -f $SUBTYPE_DEG > $work_corr_file_Diff		
	
			# compute % of over correlation threshold
			echo "[INFO] Compute % of correlation over threshold"
			python $bin_dir/profile_ME_GE.py $work_corr_file_Diff $integ_result_dir"/ME_GE_stat_"${type_kind[$k]}"_vs_"${type_kind[$l]}".txt" $((${type_len[$k]}+${type_len[$l]})) $corr_threshold $integ_result_dir/ME_GE_MAT_THRESH_${type_kind[$k]}"_vs_"${type_kind[$l]}.txt $integ_result_dir/ME_GE_MAT_PNC_${type_kind[$k]}"_vs_"${type_kind[$l]}.txt
		done
	done

<<'COMMENT'
# TODO : not finished
	# get DMR DEG intersection for each sample pairs
	for (( i=0; i<${#type_kind[@]}; i++ )); do 
		for (( j=$i+1; j<${#type_kind[@]}; j++ )); do 

			# dmr file 
			dmr_prefix=${type_kind[$i]}"_vs_"${type_kind[$j]}
			dmr_file=$me_result_dir/$dmr_prefix".mr.edgeR.s.txt" 

			# deg file 
			cel_prefix=${type_kind[$i]}"_vs_"${type_kind[$j]}
			deg_file=$cel_result_dir/$cel_prefix".deg"

			#echo $deg_file
			# convert deg file to bed 
			deg_header=`head -n 1 $deg_file`
			echo -e "#chr\tstart\tend\t$deg_header" > $integ_result_dir"/"$cel_prefix".deg.bed"
			echo -e "#chr\tpromoter_start\tpromoter_end\t$deg_header" > $integ_result_dir"/"$cel_prefix".deg_promoter.bed"
			grep promoter $REF_HUMAN_GENE_RANGE_INFO > $REF_HUMAN_GENE_RANGE_INFO"_promoter"
			awk -v input_file=$deg_file 'BEGIN { while ((getline < input_file ) > 0) data[$2] = $0 }{if (data[$5]) print $1"\t"$3"\t"$4"\t"data[$5] }' $WORK_DIR/lib/gene_symbol_chr_start_end.txt | sort -k1,1 -k2,2n >> $integ_result_dir"/"$cel_prefix".deg.bed"
			awk -v input_file=$deg_file 'BEGIN { while ((getline < input_file ) > 0) data[$2] = $0 }{if (data[$6]) print $1"\t"$2"\t"$3"\t"data[$5] }' $REF_HUMAN_GENE_RANGE_INFO"_promoter" | sort -k1,1 -k2,2n >> $integ_result_dir"/"$cel_prefix".deg_promoter.bed"

			# get intersection between DEG and DMR 
			# create intersection header
			header1=`head -n 1 $dmr_file | awk '{print $1"\t"$2"\t"$3"\t"$(NF-3)"\t"$(NF-1)"\t"$NF}'`
			header2=`head -n 1 $integ_result_dir"/"$cel_prefix".deg.bed"`
			echo -e "#$header1\t$header2" > $integ_result_dir/$dmr_prefix".dmr_deg"
			echo -e "#$header1\t$header2" > $integ_result_dir/$dmr_prefix".dmr_deg_promoter"

			awk 'NR==1{next;}{print $0}' $dmr_file| awk '{print $1"\t"$2"\t"$3"\t"$(NF-3)"\t"$(NF-1)"\t"$NF}'  | bedtools intersect -wa -wb -a - -b $integ_result_dir"/"$cel_prefix".deg.bed" >> $integ_result_dir/$dmr_prefix".dmr_deg"
			awk 'NR==1{next;}{print $0}' $dmr_file| awk '{print $1"\t"$2"\t"$3"\t"$(NF-3)"\t"$(NF-1)"\t"$NF}'  | bedtools intersect -wa -wb -a - -b $integ_result_dir"/"$cel_prefix".deg_promoter.bed" >> $integ_result_dir/$dmr_prefix".dmr_deg_promoter"

			# get % DMR related with DEGs(promoter(total, NC, PC), genebody(total, NC, PC)

			num_dmr_bins=`wc $dmr_file | awk '{print $1}'`
			num_deg=`more +2 $deg_file | cut -f2 |sort | uniq |wc | awk '{print $1}'`

			num_dmr_bins_that_intersect_w_deg=`cut -f1,2,3 $integ_result_dir/$dmr_prefix".dmr_deg" | uniq |wc -l`
			# genebody
			awk '{if ((($4>0) && $(NF-3)<0) || (($4<0) && $(NF-3)>0)) print;}' $integ_result_dir/$dmr_prefix".dmr_deg" > $integ_result_dir/$dmr_prefix".dmr_deg_NC"
			awk '{if ((($4<0) && $(NF-3)<0) || (($4>0) && $(NF-3)>0)) print;}' $integ_result_dir/$dmr_prefix".dmr_deg" > $integ_result_dir/$dmr_prefix".dmr_deg_PC"
			# promoter
			awk '{if ((($4>0) && $(NF-3)<0) || (($4<0) && $(NF-3)>0)) print;}' $integ_result_dir/$dmr_prefix".dmr_deg_promoter" > $integ_result_dir/$dmr_prefix".dmr_deg_promoter_NC"
			awk '{if ((($4<0) && $(NF-3)<0) || (($4>0) && $(NF-3)>0)) print;}' $integ_result_dir/$dmr_prefix".dmr_deg_promoter" > $integ_result_dir/$dmr_prefix".dmr_deg_promoter_PC"

			num_dmr_bins_that_have_NC_w_deg=`cut -f1,2,3 $integ_result_dir/$dmr_prefix".dmr_deg_NC" | uniq | wc -l`
			num_dmr_bins_that_have_PC_w_deg=`cut -f1,2,3 $integ_result_dir/$dmr_prefix".dmr_deg_PC" | uniq | wc -l`
			
			num_dmr_bins_that_intersect_w_deg_promoter=`cut -f1,2,3 $integ_result_dir/$dmr_prefix".dmr_deg_promoter" | uniq |wc | awk '{print $1}'`
			num_dmr_bins_that_have_NC_w_deg_promoter=`cut -f1,2,3 $integ_result_dir/$dmr_prefix".dmr_deg_promoter_NC" | uniq | wc -l`
			num_dmr_bins_that_have_PC_w_deg_promoter=`cut -f1,2,3 $integ_result_dir/$dmr_prefix".dmr_deg_promoter_PC" | uniq | wc -l`
			{
			echo "<<Pair : "$dmr_prefix">>";
			echo "Number of DMR bins : " $num_dmr_bins;
			echo "Number of DEG : " $num_deg;
			echo "Number of DMR bins in DEG genebody : "$num_dmr_bins_that_intersect_w_deg;
			echo "Number of DMR bins in DEG genebody has (PC:"$num_dmr_bins_that_have_PC_w_deg", NC:"$num_dmr_bins_that_have_NC_w_deg") rellation with DEG";
			echo "Number of DMR bins in DEG promoter : "$num_dmr_bins_that_intersect_w_deg_promoter;
			echo "Number of DMR bins in DEG promoter has (PC:"$num_dmr_bins_that_have_PC_w_deg_promoter", NC:"$num_dmr_bins_that_have_NC_w_deg_promoter") rellation with DEG";
			} > $integ_result_dir/$dmr_prefix"_DMR_DEG_stat.txt" 
		done
	done
COMMENT

# 4.2.4 DEG vs DSNP, pairwise comparison
	# get subtype snp, GE file list
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do 

			TYPE_PAIR=${type_kind[$k]}"_vs_"${type_kind[$l]}
			echo "[INFO] Currently processing $TYPE_PAIR"

			# deg file 
			deg_file=$ge_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DEG.list"

			cut -f1-3,${sample_pos_list[$k]}-$((${sample_pos_list[$k]}+${type_len[$k]}-1)),${sample_pos_list[$l]}-$((${sample_pos_list[$l]}+${type_len[$l]}-1)),$((${sample_pos_list[$k]}+${#type_list[@]}))-$((${sample_pos_list[$k]}+${type_len[$k]}-1+${#type_list[@]})),$((${sample_pos_list[$l]}+${#type_list[@]}))-$((${sample_pos_list[$l]}+${type_len[$l]}-1+${#type_list[@]})),$((${sample_pos_list[${#type_kind[@]}]}+${#type_list[@]}))-$((${sample_pos_list[${#type_kind[@]}]}+6+${#type_list[@]})),$((${sample_pos_list[${#type_kind[@]}]}+3+${sample_pos_list[$k]}+${#type_list[@]}))-$((${sample_pos_list[${#type_kind[@]}]}+3+${sample_pos_list[$k]}+${type_len[$k]}-1+${#type_list[@]})),$((${sample_pos_list[${#type_kind[@]}]}+3+${sample_pos_list[$l]}+${#type_list[@]}))-$((${sample_pos_list[${#type_kind[@]}]}+3+${sample_pos_list[$l]}+${type_len[$l]}-1+${#type_list[@]})) $integ_result_dir/SNP_GE_MGD.txt >	$integ_result_dir/SNP_GE_MGD_$TYPE_PAIR".txt"

<<'COMMENT'
			# init lists
			SNP_MGD_LIST=(); GE_MGD_LIST=(); temp_subtype_list=()
			for (( k=0; k<${#mu_list[@]}; k++ )); do 

				snp_file=$integ_result_dir"/"`basename \${mu_list[$k]}`".SNP_RNGs.MGDrefS"
				deg_only_snp_file=$snp_file"."$TYPE_PAIR"_deg_only"
				exp_file=""

				if [ $ge_data_type -eq 1 ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
					exp_file=$ge_result_dir"/"`basename \${ge_list[\$i]}`".exp.geneSymbol_avg_exp.txt"
				else #Li1_head_100000_1.htseq.geneSymbol
					file_extension=$(echo `basename \${ge_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
					exp_file=$ge_result_dir"/"`basename \${ge_list[0]} "."$file_extension`".htseq.geneSymbol"
				fi
				
			
				if [ ${type_list[$k]} == ${type_kind[$i]} ] || [ ${type_list[$k]} == ${type_kind[$j]} ] ; then
					# Filter SNP related with DEG
					awk -v input_file=$deg_file 'BEGIN{while((getline < input_file)>0) data[$1]=$0}{if (data[$10]) print}' $snp_file > $deg_only_snp_file

					# add to list
					SNP_MGD_LIST+=($deg_only_snp_file); GE_MGD_LIST+=($exp_file); temp_subtype_list+=(${type_list[$k]});
				fi
			done

			# merge DEG-only SNP with DEG's each gene expression and compute correlation
			# merge & corr	
			python $bin_dir/merge_snp3.py `my_join ";" ${SNP_MGD_LIST[@]}` `my_join ";" ${GE_MGD_LIST[@]}` $corr_threshold $integ_result_dir/SNP_DEG_AF_GE_CORR_$TYPE_PAIR".txt" $integ_result_dir/SNP_DEG_stat_$TYPE_PAIR".txt"
COMMENT
			temp_subtype_list=()

			for((i=0;i<${#mu_list[@]};i++)); do
				if [ ${type_list[$i]} == ${type_kind[$k]} ] || [ ${type_list[$i]} == ${type_kind[$l]} ] ; then
					temp_subtype_list+=(${type_list[$i]})
				fi
			done
	
		#	grep -f $deg_file $integ_result_dir/SNP_GE_MGD_$TYPE_PAIR".txt"	> $integ_result_dir/SNP_GE_MGD_$TYPE_PAIR"_CORR_DEG"
			awk -v gene_symbol_index=$((2*(${type_len[$K]}+${type_len[$l]}) + 8)) 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""}FILENAME==ARGV[2]{if(FNR==1){print $0}else{if($(gene_symbol_index) in arr){print $0}}}' $deg_file $integ_result_dir/SNP_GE_MGD_$TYPE_PAIR".txt" > $integ_result_dir/SNP_GE_MGD_$TYPE_PAIR"_CORR_DEG"

			# get DSNP 
			work_file=$integ_result_dir/DSNP_DEG_CORR_$TYPE_PAIR".txt"  
			python $bin_dir/get_GDSNP2.py  $integ_result_dir/SNP_GE_MGD_$TYPE_PAIR"_CORR_DEG" 2 `my_join ";" ${temp_subtype_list[@]}` > $work_file

			num_DSNP_related_w_DEG=`cut -f1,2,3 $work_file | sort | uniq | wc -l`
			num_DEGs_related_w_DSNP=`awk -v first_TYPE_LEN=${type_len[$k]} -v second_TYPE_LEN=${type_len[$j]} '{print $(3+(first_TYPE_LEN+second_TYPE_LEN)*2 + 4)}' $work_file | sort | uniq | tee "$integ_result_dir"/DEGs_related_DSNP_"$TYPE_PAIR"".txt" | wc -l`
			num_DEGs_that_have_DSNP_in_promoter=`grep promoter $work_file | awk -v first_TYPE_LEN=${type_len[$k]} -v second_TYPE_LEN=${type_len[$l]} '{print $(3+(first_TYPE_LEN+second_TYPE_LEN)*2 + 4)}' | sort | uniq | tee "$integ_result_dir"/DEGs_having_DSNP_in_promoter_$TYPE_PAIR".txt" | wc -l`
			{
				echo "Number of DSNP related with DEGs are : "$num_DSNP_related_w_DEG;
				echo "Number of DEGs related with DSNPs are : "$num_DEGs_related_w_DSNP;
				echo "Number of DEGs that has DSNP in promoter : "$num_DEGs_that_have_DSNP_in_promoter;
			} > $integ_result_dir/DSNP_DEG_stat_$TYPE_PAIR".txt"
		done
	done	


################################################
#### Cytoscape Network Data Generation #########
###############################################

candidate_degs_temp=$integ_result_dir/"cyto_input_deg.txt.temp"

echo -n "" > $candidate_degs_temp

for ((k=0; k<${#type_kind[@]}; k++)); do
	for ((l=$k+1; l<${#type_kind[@]}; l++)); do
		SUBTYPE_DEG=$ge_result_dir"/"${type_kind[$k]}"_vs_"${type_kind[$l]}".DEG.list"
		cat $SUBTYPE_DEG >> $candidate_degs_temp
	done
done

pair_wise_deg=$integ_result_dir/"pair_wise_degs.list"

sort $candidate_degs_temp | uniq > $pair_wise_deg

cat $GDEG"_genelist" >> $candidate_degs_temp

candidate_degs=$integ_result_dir/"cyto_input_deg.txt"

sort $candidate_degs_temp | uniq > $candidate_degs
# awk 'BEGIN{FS=OFS="\t"}{print $1, "NM_"}' > $candidate_degs

rm -rf $candidate_degs_temp

candidate_degs_tf_info=$candidate_degs".TF_info.txt"

#grep -f $candidate_degs $ENCODE_TF_GENE_INFO > $candidate_degs_tf_info
awk 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""}FILENAME==ARGV[2]{if($4 in arr){print $0}}' $candidate_degs $ENCODE_TF_GENE_INFO > $candidate_degs_tf_info


# filter deg-tf
candidate_degs_non_deg_tf_info=$candidate_degs".non_deg_TF_info.txt"
awk 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""}FILENAME==ARGV[2]{if(!($10 in arr)){print $0}}' $pair_wise_deg $candidate_degs_tf_info | bedtools groupby -i - -g 1,2,3,4,5,6 -c 10 -o collapse > $candidate_degs_non_deg_tf_info

tfbs_pos_file=$candidate_degs_non_deg_tf_info".pos.txt"
cut -f1-3 $candidate_degs_non_deg_tf_info | uniq > $tfbs_pos_file


# DO NOT USE met_level file...

# methylation file for tfbs-specific methyl level
	count=0
	for (( i=0; i<${#me_list[@]}; i++ )); do
		temp_filename_only=`basename \${me_list[$i]}`
		met_level_file=$me_result_dir/$temp_filename_only".met_level"
		met_level_file_for_tfbs=$integ_result_dir/$temp_filename_only".tfbs_overlap.met_level"

		bedtools intersect -wa -a <(tail -n+2 $ME_all_merged".txt" | cut -f1-3,$((4+$i))) -b $tfbs_pos_file | sort -k1,1 -k2,2n | uniq > $met_level_file_for_tfbs &

		let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait
	

 # run bedtools map
	tfbs_avg_methyl_file_list=()
	count=0
	for (( i=0; i<${#me_list[@]}; i++ )); do
		temp_filename_only=`basename \${me_list[$i]}`
		met_level_file_for_tfbs=$integ_result_dir/$temp_filename_only".tfbs_overlap.met_level"
		met_level_file_for_tfbs_avg=$integ_result_dir/$temp_filename_only".tfbs_avg.met_level"
		
		tfbs_avg_methyl_file_list+=($met_level_file_for_tfbs_avg)
		bedtools map -a $tfbs_pos_file -b $met_level_file_for_tfbs -c 4 -o mean | awk 'BEGIN{FS=OFS="\t"}{if($4=="."){print 0}else{print $4}}' > $met_level_file_for_tfbs_avg & 

		let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait

	tfbs_avg_methyl_file=$integ_result_dir/"TFBS_methyl_AVG_MGD"
	paste $tfbs_pos_file ${tfbs_avg_methyl_file_list[@]} | awk 'BEGIN{FS=OFS="\t"}{sum=0; for(i=4;i<=NF;i++){sum=sum+$(i)}; if(sum != 0){print $0}}' > $tfbs_avg_methyl_file
	# befor intersect - reduce size 
	

	tfbs_avg_methyl_TF_INFO=$integ_result_dir/"TFBS_methyl_AVG_TF_INFO_MGD"
	bedtools intersect -wa -wb -f 1.0 -a $candidate_degs_non_deg_tf_info -b $tfbs_avg_methyl_file > $tfbs_avg_methyl_TF_INFO

 # average promoter methylation
	tfbs_avg_methyl_TF_INFO_filtered=$tfbs_avg_methyl_TF_INFO".filtered"
	cut -f4,7,11- $tfbs_avg_methyl_TF_INFO | sort -k1,1 | uniq | awk 'BEGIN{FS=OFS="\t"; name="";num=0}{if(name != $1){if(name != ""){printf "%s\t", name; i=0;for(x in tf_arr){if(i==0){printf "%s", x;i=i+1}else{printf ",%s", x}}; for(i=3;i<=NF;i++){printf "\t%f", value_arr[i]/num}; printf "\n"}; name=$1; delete tf_arr; num=split($2,temp_arr,","); for(i=1;i<=num;i++){if(!(temp_arr[i] in tf_arr)){tf_arr[temp_arr[i]]=""}}; for(i=3;i<=NF;i++){value_arr[i] = $(i)}; num=1}else{num=split($2,temp_arr,","); for(i=1;i<=num;i++){if(!(temp_arr[i] in tf_arr)){tf_arr[temp_arr[i]]=""}};num=num+1; for(i=3;i<=NF;i++){value_arr[i] = value_arr[i] + $(i)}}}END{printf "%s\t", name; i=0;for(x in tf_arr){if(i==0){printf "%s", x;i=i+1}else{printf ",%s", x}}; for(i=3;i<=NF;i++){printf "\t%f", value_arr[i]/num}; printf "\n"}' > $tfbs_avg_methyl_TF_INFO_filtered

	ge_avg_list=()
	
	ge_header="Gene_Symbol"

		for (( i=0; i<${#ge_list[@]}; i=i+1 )); do

				temp_exp_input=""

				if [ $ge_data_type -eq 1 ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
					temp_exp_input=$ge_result_dir"/"`basename \${ge_list[\$i]}`".exp.geneSymbol_avg_exp.txt"
				else #Li1_head_100000_1.htseq.geneSymbol
					file_extension=$(echo `basename \${ge_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
					temp_exp_input=$ge_result_dir"/"`basename \${ge_list[0]} "."$file_extension`".htseq.geneSymbol"
				fi
				
				ge_avg_list+=($temp_exp_input)
				ge_header=$ge_header"\t"`basename \${ge_list[\$i]}`
		done


		ge_genesymbol_avg_exp_temp=$integ_result_dir/GE_genesymbol_avg_exp_temp

		paste ${ge_avg_list[@]} > $ge_genesymbol_avg_exp_temp

	# GE genesymbol
		ge_genesymbol_MGD=$integ_result_dir/GE_genesymbol_exp_MGD
		
		{
			echo -e $ge_header ;
			awk 'BEGIN{FS=OFS="\t"}{printf "%s", $1; for(i=2;i<=NF;i=i+2){printf "\t%f", $(i)}; printf "\n"}' $ge_genesymbol_avg_exp_temp ;
		} > $ge_genesymbol_MGD

#	cp $ge_genesymbol_MGD $final_result_dir_ge
 
# revise for average promoter methyl
	tfbs_avg_methyl_GE_COR=$integ_result_dir/"TFBS_methyl_AVG_TF_GE_MGD_COR"
	awk 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{if(NR!=1){exp_arr[$1]=$0}} FILENAME==ARGV[2]{if($1 in exp_arr){print $0,exp_arr[$1]}}' $ge_genesymbol_MGD $tfbs_avg_methyl_TF_INFO_filtered | python $bin_dir/pearson.py - 2 $((3+${#type_list[@]})) ${#type_list[@]} 0 > $tfbs_avg_methyl_GE_COR

	
#	bedtools intersect -wa -wb -a $tfbs_avg_methyl_file -b <(tail -n+2 $integ_result_dir/GE_MAT.txt | tr '_' '\t') | cut --complement -f$((11+${#type_list[@]}))-$((14+${#type_list[@]})) | sort -k1,1 -k2,2n | uniq | python $bin_dir pearson.py - 11 $((11+${#type_list[@]})) ${#type_list[@]} 0 > $tfbs_avg_methyl_GE_COR


#  if corr !='-' and corr !="NA" and corr !="nan":
  tfbs_avg_methyl_GE_COR_cut=$tfbs_avg_methyl_GE_COR".cut.txt"

#	awk -v cor_thres=$corr_threshold 'BEGIN{FS=OFS="\t"}{if( ($(NF) < (-1*cor_thres)) && ( ($(NF) != "-") && ($(NF) != "NA") && ($(NF) != "nan"))){print $4"_"$7, $0}}' $tfbs_avg_methyl_GE_COR | cut --complement -f2-11 > $tfbs_avg_methyl_GE_COR_cut
	awk -v cor_thres=$corr_threshold 'BEGIN{FS=OFS="\t"}{if( ($(NF) < (-1*cor_thres)) && ( ($(NF) != "-") && ($(NF) != "NA") && ($(NF) != "nan"))){printf "%s", $1"_"$2; for(i=3;i<=NF;i++){printf "\t%s", $(i)};printf "\n"}}' $tfbs_avg_methyl_GE_COR > $tfbs_avg_methyl_GE_COR_cut


# calculate DM-Promoter
  python $bin_dir/kruskal_ver4_wo_rpy_noheader.py $tfbs_avg_methyl_GE_COR_cut `my_join "," ${type_len[@]}` $kw_threshold $NUM_CPUS $tfbs_avg_methyl_GE_COR_cut".kw.txt" $multiple_correction

# cut DM-Promoter
 awk -v kw_threshold=$kw_threshold 'BEGIN{FS=OFS="\t"}{if(($(NF) < kw_threshold) && ($(NF) != "nan")){split($1, tokens, "_"); print tokens[1]"\t"tokens[2],$0}}' $tfbs_avg_methyl_GE_COR_cut".kw.txt" > $integ_result_dir/TFBS_DEG_DMR 


# methyl average for cytoscape
	awk -v sample_num_str=`my_join "," ${type_len[@]}` 'BEGIN{FS=OFS="\t"; num=split(sample_num_str, arr, ",")}{printf "%s\t%s\t", $1,$2; sum=0;start_index=4; for(i=1;i<=num;i++){sum=0;for(j=1;j<=arr[i];j++){sum=sum+$(start_index+j)};if(i!=num){printf "%f\t", sum/arr[i]}else{printf "%f", sum/arr[i]} ;start_index=start_index+arr[i]};printf "\n"}' $integ_result_dir/TFBS_DEG_DMR > $integ_result_dir/TFBS_DEG_DMR_me_avg

	# class average for GE
		ge_genesymbol_avg_exp=$integ_result_dir/GE_genesymbol_avg_exp

	  awk -v sample_num_str=`my_join "," ${type_len[@]}` 'BEGIN{FS=OFS="\t"; num=split(sample_num_str, arr, ",")}{printf "%s\t", $1; sum=0; start_index=0; for(i=1;i<=num;i++){sum=0;for(j=1;j<=arr[i];j++){sum=sum+$(start_index+(2*j))};printf "%f\t", sum/arr[i] ;start_index=start_index+(2*arr[i])};printf "\n"}' $ge_genesymbol_avg_exp_temp > $ge_genesymbol_avg_exp



	ALL_avg_for_cyto=$integ_result_dir/TF_GE_ME_avg
  awk 'BEGIN{FS=OFS="\t"} FILENAME==ARGV[1]{arr[$1]=$0;next}
			FILENAME==ARGV[2]{met_values=$3;for(i=4;i<=NF;i++){met_values=met_values"\t"$(i)};tf_num=split($2,tf_arr,",");for(i=1;i<=tf_num;i++){if( (tf_arr[i] in arr) && ($1 in arr)){printf "%s", arr[tf_arr[i]]; printf "%s", arr[$1]; printf "%s\t%s\t%s\n", tf_arr[i],$1,met_values;}}}' $ge_genesymbol_avg_exp $integ_result_dir/TFBS_DEG_DMR_me_avg | sort | uniq > $ALL_avg_for_cyto


	ALL_avg_for_cyto_100=$ALL_avg_for_cyto"_100_percent.txt"

	awk -v class_num=${#type_kind[@]} 'BEGIN{FS=OFS="\t"}{sum=0; start_index=1;index_arr[1]=1;index_arr[2]=(2+class_num);index_arr[3]=(4+(2*class_num)); for(j=1;j<=3;j++){start_index=index_arr[j]; sum=0;for(i=1;i<=class_num;i++){sum=sum+$(start_index+i)}; if(j<=2){printf "%s\t", $(start_index)}else{printf "%s\t", $(start_index)"_Promoter"}; per_sum=0;for(i=1;i<class_num;i++){temp=int($(start_index+i)/sum*100);per_sum=per_sum+temp; printf "%d\t", temp};printf "%d", (100-per_sum); if(j<=2){printf "\t"}}printf "\n"}' $ALL_avg_for_cyto > $ALL_avg_for_cyto_100

	ALL_avg_for_cyto_data=$ALL_avg_for_cyto".input_data.txt"
	
	awk -v class_num=${#type_kind[@]} -f $bin_dir/cytoscape_data_gen_ver2.awk $ALL_avg_for_cyto_100 > $ALL_avg_for_cyto_data

	cp $ALL_avg_for_cyto_data $final_result_dir_me/cyto.txt




############################
	# make cyto panel data#
###########################
	#template info
	template_TF_TG=$integ_result_dir"/Cyto_Template_TF_TG.txt"

	# TF TF_exp(N:type_kind) TG TF_exp(N) TF TG Pro_Methyl(N)
	
	cut --complement -f2-$(( ${#type_kind[@]} + 1 )),$(( $(( ${#type_kind[@]}*2 )) +3 )),$(( $(( ${#type_kind[@]}*2 )) + 4 )) $ALL_avg_for_cyto > $template_TF_TG

	panel_TF_TG=$integ_result_dir"/Cyto_panel_TF_TG.txt"
	echo -e "TF\tTG" > $panel_TF_TG
	cut -f1,2 $template_TF_TG >> $panel_TF_TG

	# cyto tg list
	template_TF_TG_temp_tg_list=$template_TF_TG".temp_TG.list"
	cut -f2 $template_TF_TG | sort | uniq > $template_TF_TG_temp_tg_list

	kind_header_list=()
	for((k=0;k<${#type_kind[@]};k++)); do
		kind_header_list+=(${type_kind[$k]})
	done

	# cyto tg exp
	panel_TG_exp=$integ_result_dir"/Cyto_panel_TG_exp.txt"
	echo -e "Gene\t"`my_join ";" ${kind_header_list[@]}` | tr ";" "\t" > $panel_TG_exp
	cut -f2-$((${#type_kind[@]}+2)) $template_TF_TG | sort -k1,1 | uniq >> $panel_TG_exp

	# DEG-Promoter Methyl Correlation
	DEG_Promoter_methyl_CORR=$integ_result_dir"/Cyto_Template_TG_Corr.txt"
	cut -f$((${#type_list[@]}+2)),$(( $(( ${#type_list[@]} * 2 )) + 3 )) $tfbs_avg_methyl_GE_COR_cut | awk 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]="";next}FILENAME==ARGV[2]{if($1 in arr){print $0}}' $template_TF_TG_temp_tg_list - > $DEG_Promoter_methyl_CORR


	# for promoter methyl
	panel_promoter=$integ_result_dir"/Cyto_panel_TG_promoter_methyl.txt"
	echo -e "Promoter\t"`my_join ";" ${kind_header_list[@]}` | tr ";" "\t" > $panel_promoter

	awk -v met_col=$((${#type_kind[@]}+3)) 'BEGIN{FS=OFS="\t"}{printf "%s", $2"_Promoter"; for(i=met_col;i<=NF;i++){printf "\t%f", $(i)}; printf "\n"}' $template_TF_TG | sort -k1,1 | uniq >> $panel_promoter

	# TG info file
	#<Gene_symbol, Promoter, corr>
	panel_TG_Info=$integ_result_dir"/Cyto_panel_TG_Info.txt"
	echo -e "Gene\tPromoter\tCorr" > $panel_TG_Info
	awk 'BEGIN{FS=OFS="\t"}{print $1,$1"_Promoter",$2}' $DEG_Promoter_methyl_CORR >> $panel_TG_Info



######################################################################################
# DEG & Mutation
########################################################################################
#for pair-wise
for((k=0; k<${#type_kind[@]}; k++)); do
  for((l=$k+1; l<${#type_kind[@]}; l++)); do
    work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}

		work_gene_list=$ge_result_dir/$work_file".DEG.list"

		snp_ge_MGD=$integ_result_dir/SNP_GE_MGD_$work_file".txt"

		DEG_with_mutation=$integ_result_dir"/"$work_file".DEG_with_mutation.genelist"
		tail -n+2 $snp_ge_MGD | grep -w "genebody" | awk -v genesymbol_index=$(( $(( ${type_len[$k]} + ${type_len[$l]} )) * 2 + 8 )) -v sample_nums=$((${type_len[$k]} + ${type_len[$l]})) 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""} FILENAME==ARGV[2]{if($(genesymbol_index) in arr){sum=0; for(i=4;i<(4+sample_nums);i++){if($(i) != "-/-"){sum=sum+1}};if(sum > 0){print $(genesymbol_index)}}}' $work_gene_list - | sort -k1,1 | uniq > $DEG_with_mutation

	#filter DEG table
	DEG_table=$ge_result_dir"/"$work_file".DEG.table.txt"
	DEG_with_mut_table=$integ_result_dir"/"$work_file".DEG_with_mutation.table.txt"
	DEG_with_mut_table_top100=$DEG_with_mut_table".top100.txt"

	temp_gs_index=1
	if [ "$ge_data_type" -eq "1" ]; then # array
		# $2 is gene symbol
		temp_gs_index=2
	else # rna-seq
		temp_gs_index=1
	fi

	awk -v genesymbol_index=$temp_gs_index 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""}FILENAME==ARGV[2]{if(FNR==1){print $0}else{if($(genesymbol_index) in arr){print $0}}}' $DEG_with_mutation $DEG_table > $DEG_with_mut_table

  head -n 101 $DEG_with_mut_table > $DEG_with_mut_table_top100

	# for onco-print
	DEG_mutation_onco_print_top30=$integ_result_dir"/"$work_file".DEG_w_mut_onco_top30.txt"

	awk -v genesymbol_index=$temp_gs_index 'BEGIN{FS=OFS="\t";i=0}NR>1{if(i<30){if(!($(genesymbol_index) in arr) && ($(genesymbol_index) != "NA")){i=i+1;arr[$(genesymbol_index)]=""}}else{exit}}END{for(x in arr){print x}}' $DEG_with_mut_table > $DEG_mutation_onco_print_top30	

	# run DAVID
	bash $bin_dir/gsea3.sh $DEG_with_mutation $integ_result_dir "GENE_SYMBOL" $work_file".DEG_with_mutation.genelist"


	done
done

# for ALL

	work_gene_list=$GDEG"_genelist"

	snp_ge_MGD=$integ_result_dir/"SNP_GE_MGD.txt"

	work_file="ALL"

	DEG_with_mutation=$integ_result_dir"/"$work_file".DEG_with_mutation.genelist"

  tail -n+2 $snp_ge_MGD | grep -w "genebody" | awk -v genesymbol_index=$((${#mu_list[@]}*2+8)) -v sample_nums=${#mu_list[@]} 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""} FILENAME==ARGV[2]{if($(genesymbol_index) in arr){sum=0; for(i=4;i<(4+sample_nums);i++){if($(i) != "-/-"){sum=sum+1}};if(sum > 0){print $(genesymbol_index)}}}' $work_gene_list - | sort -k1,1 | uniq > $DEG_with_mutation

	#filter DEG table
	DEG_table=$GDEG_exp_value_table
	DEG_with_mut_table=$integ_result_dir"/"$work_file".DEG_with_mutation.table.txt"
	DEG_with_mut_table_top100=$DEG_with_mut_table".top100.txt"

	# first column is gene symbol
	temp_gs_index=1

	awk -v genesymbol_index=$temp_gs_index 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""}FILENAME==ARGV[2]{if(FNR==1){print $0}else{if($(genesymbol_index) in arr){print $0}}}' $DEG_with_mutation $DEG_table > $DEG_with_mut_table

  head -n 101 $DEG_with_mut_table > $DEG_with_mut_table_top100


	# for onco-print
	DEG_mutation_onco_print_top30=$integ_result_dir"/"$work_file".DEG_w_mut_onco_top30.txt"

	awk -v genesymbol_index=$temp_gs_index 'BEGIN{FS=OFS="\t";i=0}NR>1{if(i<30){if(!($(genesymbol_index) in arr) && ($(genesymbol_index) != "NA")){i=i+1;arr[$(genesymbol_index)]=""}}else{exit}}END{for(x in arr){print x}}' $DEG_with_mut_table > $DEG_mutation_onco_print_top30	


	# run DAVID
	bash $bin_dir/gsea3.sh $DEG_with_mutation $integ_result_dir "GENE_SYMBOL" $work_file".DEG_with_mutation.genelist"

	# move results
	cp $integ_result_dir"/"*.DEG_with_mutation.table.txt $integ_result_dir"/"*.DEG_with_mutation.table.txt.top100.txt $integ_result_dir"/"*.DEG_with_mutation.genelist.DAVID_*.txt $final_result_dir_mu
 
#########################################################################################

############################################3
# DEG & DMR (Promoter)
###############################################
#for pair-wise
for (( k=0; k<${#type_kind[@]}; k++)); do
	for (( l=$k+1; l<${#type_kind[@]}; l++)); do
  	work_file=$ME_GE_all_merged"_"${type_kind[$k]}"_vs_"${type_kind[$l]}
    work_corr_file=$ME_GE_COR_all_merged"_"${type_kind[$k]}"_vs_"${type_kind[$l]}
		work_corr_file_Diff=$work_corr_file".DMR_DEG"

		work_file_name=${type_kind[$k]}"_vs_"${type_kind[$l]}
		DEG_DMR_promoter=$integ_result_dir"/"$work_file_name".DEG_with_DMR.genelist"

		grep -w "promoter" $work_corr_file_Diff | awk -v cor_thres=$corr_threshold -v genesymbol_index=$((${type_len[$k]}+${type_len[$l]}+8)) 'BEGIN{FS=OFS="\t"}{if($(NF) <= (-1*cor_thres)){print $(genesymbol_index)}}' | sort | uniq > $DEG_DMR_promoter

	#filter DEG table
	DEG_table=$ge_result_dir"/"$work_file_name".DEG.table.txt"
	DEG_with_DMR_table=$integ_result_dir"/"$work_file_name".DEG_with_DMR.table.txt"
	DEG_with_DMR_table_top100=$DEG_with_DMR_table".top100.txt"

	temp_gs_index=1
	if [ "$ge_data_type" -eq "1" ]; then # array
		# $2 is gene symbol
		temp_gs_index=2
	else # rna-seq
		temp_gs_index=1
	fi

	awk -v genesymbol_index=$temp_gs_index 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""}FILENAME==ARGV[2]{if(FNR==1){print $0}else{if($(genesymbol_index) in arr){print $0}}}' $DEG_DMR_promoter $DEG_table > $DEG_with_DMR_table

  head -n 101 $DEG_with_DMR_table > $DEG_with_DMR_table_top100


	# run DAVID
	bash $bin_dir/gsea3.sh $DEG_DMR_promoter $integ_result_dir "GENE_SYMBOL" $work_file_name".DEG_with_DMR.genelist"


	done
done

#for ALL
  GDEG_GDMR_MGD=$integ_result_dir/GDMR_GDEG_w_GSym.txt.all_merged_data.txt
  GDEG_GDMR_MGD_CORR=$GDEG_GDMR_MGD".corr.txt"

	work_file_name="ALL"

	DEG_DMR_promoter=$integ_result_dir"/"$work_file_name".DEG_with_DMR.genelist"

	grep -w "promoter" $GDEG_GDMR_MGD_CORR | awk -v cor_thres=$corr_threshold -v genesymbol_index=$((${#type_list[@]}+8)) 'BEGIN{FS=OFS="\t"}{if($(NF) <= (-1*cor_thres)){print $(genesymbol_index)}}' | sort | uniq > $DEG_DMR_promoter


	#filter DEG table
	DEG_table=$GDEG_exp_value_table
	DEG_with_DMR_table=$integ_result_dir"/"$work_file_name".DEG_with_DMR.table.txt"
	DEG_with_DMR_table_top100=$DEG_with_DMR_table".top100.txt"

	temp_gs_index=1

	awk -v genesymbol_index=$temp_gs_index 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""}FILENAME==ARGV[2]{if(FNR==1){print $0}else{if($(genesymbol_index) in arr){print $0}}}' $DEG_DMR_promoter $DEG_table > $DEG_with_DMR_table

  head -n 101 $DEG_with_DMR_table > $DEG_with_DMR_table_top100

	# run DAVID
	bash $bin_dir/gsea3.sh $DEG_DMR_promoter $integ_result_dir "GENE_SYMBOL" $work_file_name".DEG_with_DMR.genelist"


	# move results to final directory
	cp $integ_result_dir"/"*.DEG_with_DMR.table.txt $integ_result_dir"/"*.DEG_with_DMR.table.txt.top100.txt $integ_result_dir"/"*.DEG_with_DMR.genelist.DAVID_*.txt $final_result_dir_me
	
#######################################################################################

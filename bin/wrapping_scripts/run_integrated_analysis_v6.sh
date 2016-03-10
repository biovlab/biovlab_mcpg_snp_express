#!/bin/bash
source `dirname $0`/../../env.sh

# parse experiment information
source `dirname $0`/parse_exp_info.sh

# directories
bin_dir="$WORK_DIR/bin"
result_dir="$WORK_DIR/result"
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
kw_threshold=0.15
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
if [ $ge_type -eq 1 ]; then
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
TYPE_LEN=($(echo ${ge_type_list[@]} | tr ' ' '\n' | sort | uniq -c | awk '{print $1}'| tr '\n' ' '))

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
multiple_correction="fdr" #'bonferroni', 'fdr', 'holm' and so on. see R stats.p.adjust

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

				if [ $ge_type -eq 1 ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
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
		echo "[INFO] subtype length : "${TYPE_LEN[$k]}

		head -n1 $work_file | xargs -I {} echo -e "{}\tcorrelation" > $work_file2

		# Pearson
		# awk 'NR==1{next;}{print $0}' $work_file | parallel --no-notice --pipe -j"$NUM_CPUS" -L1000000 -k python $bin_dir/pearson.py - 3 $((9 + ${TYPE_LEN[$k]} )) ${TYPE_LEN[$k]} >> $work_file2

		# Spearman : Score is bad, so not use currently
		#Rscript p_spearman.r -i test.txt -n 3 --s1 2 --s2 8 -p 1 -o ./test_result.txt 

		# compute % of over correlation threshold
		#echo "[INFO] Compute % of correlation over threshold"
		#python $bin_dir/profile_ME_GE.py $work_file2 $integ_result_dir"/ME_GE_stat_"${type_kind[$k]}".txt" ${TYPE_LEN[$k]} $corr_threshold $integ_result_dir/ME_GE_MAT_THRESH_${type_kind[$k]}.txt $integ_result_dir/ME_GE_MAT_PNC_${type_kind[$k]}.txt
	done
COMMENT

# 4.1.2 SNP vs GE : TSS, TSE +- range, Correlation(Spearman's rho), % over threshold, % of NC and PC
# TODO: Change as ME vs GE to speed up
	# extract only snps in range
	count=0
	
	snp_file_extension=$(echo `basename \${mu_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
	for (( i=0; i<${#mu_list[@]}; i++ )); do
		temp_filename_only=`basename \${mu_list[\$i]}`
		temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`
		snp_level_file=""
		if [ "$mu_data_type" -eq "0" ]; then
			snp_level_file=$mu_result_dir/$temp_filename_only".vcf.bed"
		elif [ "$mu_data_type" -eq "1" ]; then
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

				if [ $ge_type -eq 1 ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
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

	# FOR PACKAGE, THIS IS NOT USEFUL, so commented
<<'COMMENT'
# 4.1.3 ENCODE/GEO data association
## 4.1.3.1 Intersect with DHS cluster
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
### 4.1.3.1.0 ME-GE bins 
		# 13763956
		#total_DMR_bin_num=`cut -f1,2,3 $ME_GE_COR_all_merged"_"${type_kind[$k]} | uniq | wc -l`
		bedtools intersect -wa -a $ME_GE_COR_all_merged"_"${type_kind[$k]} -b $DHS_CLUSTER_UCSC > $ME_GE_COR_all_merged"_"${type_kind[$k]}"_DHS"

		python $bin_dir/profile_ME_GE.py $ME_GE_COR_all_merged"_"${type_kind[$k]}"_DHS" $integ_result_dir"/ME_GE_stat_"${type_kind[$k]}"_DHS.txt" ${TYPE_LEN[$k]} $corr_threshold $integ_result_dir/ME_GE_MAT_THRESH_${type_kind[$k]}_DHS.txt $integ_result_dir/ME_GE_MAT_PNC_${type_kind[$k]}_DHS.txt

		
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
	###########################################
	# BY kruscal & and bonferroni
	###########################################

		python $bin_dir/kruskal_ver4.py $integ_result_dir/ME_MAT.txt `my_join "," ${TYPE_LEN[@]}` $kw_threshold $NUM_CPUS $integ_result_dir/ME_MAT.kw.txt $multiple_correction

#-a 0.05 -n $NUM_CPUS -c `my_join "," ${type_list[@]}` -m 0 --o1 $integ_result_dir/ME_MAT.kw.sigle_diff --o2 $integ_result_dir/ME_MAT.kw.all_diff

		python $bin_dir/kruskal_ver4.py $integ_result_dir/GE_MAT.txt `my_join "," ${TYPE_LEN[@]}` $kw_threshold $NUM_CPUS $integ_result_dir/GE_MAT.kw.txt $multiple_correction

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
	num_GDEG=`awk 'NR==1{split($0, tokens, "\t");for(i in tokens){if(tokens[i]=="gene_symbol") gene_symbol_index=i;};next;}{print $1"\t"$2"\t"$3"\t"$gene_symbol_index}' $integ_result_dir/ME_GE_MGD_${type_kind[0]} | bedtools intersect -wb -wa -a $GDEG -b - |cut -f4-| tee $GDEG"_GSym" |cut -f4 | sort | uniq | tee $GDEG"_genelist" | wc | awk '{print $1}'`
	
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

# 4.2.3 DEG vs DMR, pairwise subtype comparison


	# Make SUBTYPE FILES
	for (( k=0; k<${#type_kind[@]}; k++)); do
		for (( l=$k+1; l<${#type_kind[@]}; l++)); do
			work_file=$ME_GE_all_merged"_"${type_kind[$k]}"_vs_"${type_kind[$l]}
			work_corr_file=$ME_GE_COR_all_merged"_"${type_kind[$k]}"_vs_"${type_kind[$l]}
			

			first_sample_file=$ME_GE_all_merged"_"${type_kind[$k]}
			second_sample_file=$ME_GE_all_merged"_"${type_kind[$l]}

			# make SUBTYPE_DIFF base file
			# chr stat end (TYPE_1 ME) (TYPE_2 ME) TSS+-range/TSS TSS+-range/TSS strand GeneSymbol range_kind refseq (TYPE_1 GE) (TYPE_2 GE)
			paste <(cut -f1-$((${TYPE_LEN[$k]}+3)) $first_sample_file) <(cut -f4-$((${TYPE_LEN[$l]}+3)) $second_sample_file) $integ_result_dir/temp_pos.txt <(cut -f$((${TYPE_LEN[$k]}+10))- $first_sample_file) <(cut -f$((${TYPE_LEN[$l]}+10))- $second_sample_file) > $work_file

			rm -rf $work_corr_file

			# Compute correaltion
			awk 'NR==1{next;}{print $0}' $work_file | parallel --no-notice --pipe -j"$NUM_CPUS" -L1000000 -k python $bin_dir/pearson.py - 3 $((9 + ${TYPE_LEN[$k]} + ${TYPE_LEN[$l]} )) $((${TYPE_LEN[$k]}+${TYPE_LEN[$l]})) $pear_sp_cor >> $work_corr_file

			# filter only DMR & DEG
			SUBTYPE_DMR=$me_result_dir"/"${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.bed"
			# chr start end
			SUBTYPE_DEG=$ge_result_dir"/"${type_kind[$k]}"_vs_"${type_kind[$l]}".DEG.list"
			# gene_symbol

			work_corr_file_Diff=$work_corr_file".DMR_DEG"
			bedtools intersect -wa -a $work_corr_file -b $SUBTYPE_DMR | grep -w -f $SUBTYPE_DEG > $work_corr_file_Diff		
	
			# compute % of over correlation threshold
			echo "[INFO] Compute % of correlation over threshold"
			python $bin_dir/profile_ME_GE.py $work_corr_file_Diff $integ_result_dir"/ME_GE_stat_"${type_kind[$k]}"_vs_"${type_kind[$l]}".txt" $((${TYPE_LEN[$k]}+${TYPE_LEN[$l]})) $corr_threshold $integ_result_dir/ME_GE_MAT_THRESH_${type_kind[$k]}"_vs_"${type_kind[$l]}.txt $integ_result_dir/ME_GE_MAT_PNC_${type_kind[$k]}"_vs_"${type_kind[$l]}.txt
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
	for (( i=0; i<${#type_kind[@]}; i++ )); do 
		for (( j=$i+1; j<${#type_kind[@]}; j++ )); do 

			TYPE_PAIR=${type_kind[$i]}_${type_kind[$j]}
			echo "[INFO] Currently processing $TYPE_PAIR"

			# deg file 
			deg_file=$ge_result_dir/${type_kind[$i]}"_vs_"${type_kind[$j]}".deg"

			# init lists
			SNP_MGD_LIST=(); GE_MGD_LIST=(); temp_subtype_list=()
			for (( k=0; k<${#mu_list[@]}; k++ )); do 

				snp_file=$integ_result_dir"/"`basename \${mu_list[$k]}`".SNP_RNGs.MGDrefS"
				deg_only_snp_file=$snp_file"."$TYPE_PAIR"_deg_only"
				exp_file=""

				if [ $ge_type -eq 1 ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
					exp_file=$ge_result_dir"/"`basename \${ge_list[\$i]}`".exp.geneSymbol_avg_exp.txt"
				else #Li1_head_100000_1.htseq.geneSymbol
					file_extension=$(echo `basename \${ge_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
					exp_file=$ge_result_dir"/"`basename \${ge_list[0]} "."$file_extension`".htseq.geneSymbol"
				fi
				
			
				if [ ${type_list[$k]} == ${type_kind[$i]} ] || [ ${type_list[$k]} == ${type_kind[$j]} ] ; then
					# Filter SNP related with DEG
					awk -v input_file=$deg_file 'BEGIN{while((getline < input_file)>0) data[$2]=$0}{if (data[$10]) print}' $snp_file > $deg_only_snp_file

					# add to list
					SNP_MGD_LIST+=($deg_only_snp_file); GE_MGD_LIST+=($exp_file); temp_subtype_list+=(${type_list[$k]});
				fi
			done

			# merge DEG-only SNP with DEG's each gene expression and compute correlation
			# merge & corr	
			python $bin_dir/merge_snp3.py `my_join ";" ${SNP_MGD_LIST[@]}` `my_join ";" ${GE_MGD_LIST[@]}` $corr_threshold $integ_result_dir/SNP_DEG_AF_GE_CORR_$TYPE_PAIR".txt" $integ_result_dir/SNP_DEG_stat_$TYPE_PAIR".txt"

			# get DSNP 
			work_file=$integ_result_dir/DSNP_DEG_CORR_$TYPE_PAIR".txt"  
			python $bin_dir/get_GDSNP2.py $integ_result_dir/SNP_DEG_AF_GE_CORR_$TYPE_PAIR".txt" 2 `my_join ";" ${temp_subtype_list[@]}` > $work_file

			num_DSNP_related_w_DEG=`cut -f1,2,3 $work_file | sort | uniq | wc -l`
			num_DEGs_related_w_DSNP=`awk -v first_TYPE_LEN=${TYPE_LEN[$i]} -v second_TYPE_LEN=${TYPE_LEN[$j]} '{print $(3+(first_TYPE_LEN+second_TYPE_LEN)*2 + 4)}' $work_file | sort | uniq | tee "$integ_result_dir"/DEGs_related_DSNP_"$TYPE_PAIR"".txt" | wc -l`
			num_DEGs_that_have_DSNP_in_promoter=`grep promoter $work_file | awk -v first_TYPE_LEN=${TYPE_LEN[$i]} -v second_TYPE_LEN=${TYPE_LEN[$j]} '{print $(3+(first_TYPE_LEN+second_TYPE_LEN)*2 + 4)}' | sort | uniq | tee "$integ_result_dir"/DEGs_having_DSNP_in_promoter_$TYPE_PAIR".txt" | wc -l`
			{
				echo "Number of DSNP related with DEGs are : "$num_DSNP_related_w_DEG;
				echo "Number of DEGs related with DSNPs are : "$num_DEGs_related_w_DSNP;
				echo "Number of DEGs that has DSNP in promoter : "$num_DEGs_that_have_DSNP_in_promoter;
			} > $integ_result_dir/DSNP_DEG_stat_$TYPE_PAIR".txt"
		done
	done	

<<'COMMENT'
# 4.3 MBDseq, BS seq comparioson ICBP project specific
if [ $mbdseq_bs_corr -eq 1 ]; then
	# get methylation correlation between MBDseq and BSseq

	gene_5p_plus_10000=$WORK_DIR/lib/region_info/bs_seq_target_region_sorted.txt # This is 10000bp extended to 5' upstream to see influence from TF
	gene_only=$integ_result_dir/bs_seq_target_region_sorted_gene_only.txt 
	promoter_only=$integ_result_dir/bs_seq_target_region_sorted_promoter_only.txt
	awk '{if($6=="-") {print $1"\t"$2"\t"$3-10000"\t"$4-10000"\t"$5"\t"$6}else{print $1"\t"$2+10000"\t"$3"\t"$4-10000"\t"$5"\t"$6}}' $gene_5p_plus_10000 > $gene_only
	awk '{if($6=="-") {print $1"\t"$3-10000"\t"$3-10000+2000"\t2000\t"$5"\t"$6}else{print $1"\t"$2+10000-2000"\t"$2+10000"\t2000\t"$5"\t"$6}}' $gene_5p_plus_10000 > $promoter_only

	targeted_regions=$gene_5p_plus_10000
	count=0
	for (( i=0; i<${#bs_r1_list[@]}; i++ )); do 

		# get extention and filenale only w/o extention
		input_fq1_filename_only=`basename \${bs_r1_list[$i]}`

		# BSseq methylation file having coverage
		BS_methyl_file=$methylkit_from_bs_result_dir/$input_fq1_filename_only"_CpG.txt"
		BS_cov_file=$integ_result_dir/$input_fq1_filename_only".cov"
		BS_HML_file=$integ_result_dir/$input_fq1_filename_only".HML"
		MBD_methyl_file=$me_result_dir/`basename \${bs_paired_me_list[$i]}`".met_level";
		MBD_bam_file=$mbd_result_dir/`basename \${bs_paired_me_list[$i]}`".sorted.bam"


		# Compute genome coverage of mbdseq
		echo "[INFO] Compute correlation between $BS_methyl_file and $MBD_methyl_file"
	
		{ 

		# BS : convert BS converage file to bed format
		awk 'BEGIN{OFS="\t"}NR==1{next;}{if($4=="R"){print $2, $3-1, $3, $5*$6/100, "-"}else{print $2, $3-1, $3, $5*$6/100, "+"}}' $BS_methyl_file | sort -k1,1 -k2,2n > $BS_cov_file".bed";

		# THIS TAKE LONG, so TEMP COMMENTED
		#bedtools genomecov -d -ibam $MBD_bam_file -g $REF_HUMAN_CHR_SIZE | awk '{OFS="\t"; print $1, $2-1, $2, $3}'  | bedtools intersect -wa -wb -a - -b $BS_cov_file".bed" > $BS_cov_file".MBD_BS_intersect";

		#sort
		sort -k1,1 -k2,2n $BS_cov_file".MBD_BS_intersect" > $BS_cov_file".MBD_BS_intersect.sorted"

		# MBD : extract only targeted region from MBDseq result
		work_file=$integ_result_dir/bs_related_MBD_`basename \${bs_paired_me_list[$i]}`;
		awk 'NR==1{next;}{print $0}' $MBD_methyl_file | bedtools intersect -wa -wb -b $targeted_regions -a - |sort -k1,1 -k2,2n > $work_file;

		# MBD : merge near bins 
		awk -v merge_count=20 'BEGIN{OFS="\t"}{if(NR!=1 && count%merge_count!=0){split(prev, tokens,"\t"); if(tokens[1]!=$1){print tokens[3], tokens[4], tokens[5], tokens[6], tokens[7], tokens[8], tokens[9], tokens[10], tokens[11]; count=0;}};count=count+1;if (count%merge_count==1){printf "%s", $1"\t"$2"\t";temp_sum=$4;temp_sum2=$5;}else if(count%merge_count==0 ){print $3, temp_sum+$4, temp_sum2+$5, $6, $7, $8, $9, $10, $11}else{temp_sum=temp_sum+$4; temp_sum2=temp_sum2+$5;}; prev=$0;}END{if(count%merge_count!=0){print $3, temp_sum+$4, temp_sum2+$5, $6, $7, $8, $9, $10, $11}}' $work_file > $work_file"_merged";


		# MBD BS CORR, genome coverage based read count based correlation (read count on CpG from MBD, methylated read count on CpG site from BS)
		cut -f1,2,3 $work_file"_merged" | /packages/bedtools2-v2.20.1/bin/bedtools map -a - -b $BS_cov_file".MBD_BS_intersect.sorted" -c 4,8 -o sum,sum | cut -f4,5 | awk 'BEGIN{OFS="\t"}{if ($1=="."){if($2=="."){print "0", "0"}else{print "0", $2}}else{print}}' > $BS_cov_file".MBD_BS_intersect.sorted.val"

	       	awk -f $bin_dir/correlation.awk $BS_cov_file".MBD_BS_intersect.sorted.val" > $BS_cov_file".MBD_BS_intersect.sorted.val.BS_MBD_CORR"

		awk '{if($1>50){print $0"\tH"}else{print $0"\tL"}}' $BS_cov_file".MBD_BS_intersect.sorted.val" > $BS_HML_file










		#cp $work_file $work_file"_merged"
		# BS : binning BS seq based on MBDseq , DEDIPS result based (300bp extension on each read, NOT USED)
		/packages/bedtools2-v2.20.1/bin/bedtools map -a $work_file"_merged" -b $BS_cov_file".bed" -c 4 -o sum |tee $BS_cov_file".merged" | awk 'BEGIN{OFS="\t"}{if ($12=="."){print $5, "0"}else{print $5, $12}}' > $BS_cov_file".binned";

		# compute correlation: $5(mbd) $12(BS binned)
		# Pearson correlation
	       	awk -f $bin_dir/correlation.awk $BS_cov_file".binned" > $integ_result_dir/$input_fq1_filename_only".BS_MBD_CORR";

		# Spearman's rho (Rank correlation) : BAD RESULT
		#Rscript $bin_dir/p_spearman_2column.r -o $integ_result_dir/$input_fq1_filename_only".BS_MBD_SPEARMAN_CORR" -i $BS_cov_file".binned";

		# get H/M/L label on BS 

		# BS : convert BS converage file to bed format
		awk 'BEGIN{OFS="\t";}NR==1{next;}{if($4=="R"){print $2, $3-1, $3, $6, "-"}else{print $2, $3-1, $3, $6, "+"}}' $BS_methyl_file | sort -k1,1 -k2,2n > $BS_cov_file".level";

		# merge based on MBD bin
		cut -f1-3 $work_file"_merged" | bedtools map -a - -b $BS_cov_file".level" -c 4 -o mean | awk '{if ($4=="."){print $1"\t"$2"\t"$3"\t0"}else{print}}' |awk '{if($4>70){print $0"\tH"}else if($4 <30){print $0"\tL"}else{print $0"\tM"}}' > $BS_HML_file
		#cut -f1-3 $work_file"_merged" | /packages/bedtools2-v2.20.1/bin/bedtools map -a - -b $BS_cov_file".level" -c 4 -o mean | awk 'BEGIN{OFS="\t"}{if ($4=="."){print $1, $2, $3, "0"}else{print}}' |awk 'BEGIN{OFS="\t"}{if($4>50){print $0, "H"}else{print $0, "L"}}' > $BS_HML_file
		#awk 'NR==1{next;}{if($6>70){print $2"\t"$3-1"\t"$3"\t"$6"\tH"}else if($6 <30){print $2"\t"$3-1"\t"$3"\t"$6"\tL"}else{print $2"\t"$3-1"\t"$3"\t"$6"\tM"}}' $BS_methyl_file > $BS_HML_file;
		} &	
		let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait


		
	done
	wait
fi
exit
COMMENT
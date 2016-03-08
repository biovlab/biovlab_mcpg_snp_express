#!/bin/bash
source `dirname $0`/../../env.sh

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
mutation_result_dir=$result_dir"/snp/from_mbd"
viz_dir=$result_dir"/visualization"

# directory based variables
ME_all_merged=$integ_result_dir/ME_MGD		# prepare initiall prefix for all merged file
ME_GE_all_merged=$integ_result_dir/ME_GE_MGD
ME_GE_COR_all_merged=$integ_result_dir/ME_GE_COR_MGD
ME_GE_COR_all_merged_THRESH=$integ_result_dir/ME_GE_COR_MGD_THRESH.txt
#REF_HUMAN_GENE_RNAGE_INFO=$integ_result_dir/REF_HUMAN_GENE_RNAGE_INFO.txt

# for testing data
fastq_sample1="/data/project/mcpg/test_data/small_test_data/icbp/BrCa-02_head10000000.fastq"
fastq_sample2="/data/project/mcpg/test_data/small_test_data/icbp/BrCa-03_head10000000.fastq"
input_fq1_filename_only=`basename \$fastq_sample1`
input_fq2_filename_only=`basename \$fastq_sample2`
sample1_list=$result_dir/$input_fq1_filename_only".sorted.bam"
sample2_list=$result_dir/$input_fq2_filename_only".sorted.bam"

# variables
window_size=100 
TSS_range_for_methyl=2000
TSS_range_for_snp=250000
gene_exp_fold_chage=2
corr_threshold=0.3
kw_threshold=0.05
no_replicate=0

mbdseq_bs_corr=1


# functions
function my_join { local IFS="$1"; shift; echo "$*"; }



####################################################################################
# 0. set each of data type 
####################################################################################

# methyl
methyl_type=1			# 1:MBD, 2:BS, 3:infinium

# gene express
ge_type=1					# 1:array, 2:RNAseq

# mutation
mu_type=2					# 1:RNAseq, 2:DNAseq

# input class
type_kind=("lu" "baa" "bab")
type_list=("lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "baa" "baa" "baa" "baa" "baa" "baa" "baa" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab")
TYPE_LEN=(13 7 10)


#######################################################################################
# TEMP INPUTS
#######################################################################################

# MBD
		mbd_lu1=/data/project/mcpg/test_data/icbp/fastq/100730_s_1.fq.head1000000 
		mbd_lu2=/data/project/mcpg/test_data/icbp/fastq/100730_s_7.fq.head1000000
		mbd_lu3=/data/project/mcpg/test_data/icbp/fastq/100730_s_8.fq.head1000000 
		mbd_lu4=/data/project/mcpg/test_data/icbp/fastq/100803_s_5.fq.head1000000
		mbd_lu5=/data/project/mcpg/test_data/icbp/fastq/100803_s_6.fq.head1000000
		mbd_lu6=/data/project/mcpg/test_data/icbp/fastq/100812_s_3.fq.head1000000
		mbd_lu7=/data/project/mcpg/test_data/icbp/fastq/100812_s_4.fq.head1000000
		mbd_lu8=/data/project/mcpg/test_data/icbp/fastq/100824_s_2.fq.head1000000
		mbd_lu9=/data/project/mcpg/test_data/icbp/fastq/100831_s_5.fq.head1000000
		mbd_lu10=/data/project/mcpg/test_data/icbp/fastq/100902_s_7.fq.head1000000
		mbd_lu11=/data/project/mcpg/test_data/icbp/fastq/100908_s_3.fq.head1000000
		mbd_lu12=/data/project/mcpg/test_data/icbp/fastq/100908_s_4.fq.head1000000
		mbd_lu13=/data/project/mcpg/test_data/icbp/fastq/100908_s_8.fq.head1000000

		mbd_baa1=/data/project/mcpg/test_data/icbp/fastq/100730_s_5.fq.head1000000
		mbd_baa2=/data/project/mcpg/test_data/icbp/fastq/100831_s_6.fq.head1000000 
		mbd_baa3=/data/project/mcpg/test_data/icbp/fastq/100831_s_8.fq.head1000000
		mbd_baa4=/data/project/mcpg/test_data/icbp/fastq/100730_s_4.fq.head1000000
		mbd_baa5=/data/project/mcpg/test_data/icbp/fastq/100812_s_2.fq.head1000000
		mbd_baa6=/data/project/mcpg/test_data/icbp/fastq/100730_s_6.fq.head1000000
		mbd_baa7=/data/project/mcpg/test_data/icbp/fastq/100908_s_5.fq.head1000000

		mbd_bab1=/data/project/mcpg/test_data/icbp/fastq/100730_s_2.fq.head1000000
		mbd_bab2=/data/project/mcpg/test_data/icbp/fastq/100812_s_1.fq.head1000000
		mbd_bab3=/data/project/mcpg/test_data/icbp/fastq/100803_s_3.fq.head1000000
		mbd_bab4=/data/project/mcpg/test_data/icbp/fastq/100803_s_7.fq.head1000000
		mbd_bab5=/data/project/mcpg/test_data/icbp/fastq/100812_s_5.fq.head1000000
		mbd_bab6=/data/project/mcpg/test_data/icbp/fastq/100824_s_4.fq.head1000000
		mbd_bab7=/data/project/mcpg/test_data/icbp/fastq/100831_s_4.fq.head1000000 # Not in GEO, so created from ICBP data ELAND export format
		mbd_bab8=/data/project/mcpg/test_data/icbp/fastq/100902_s_6.fq.head1000000
		mbd_bab9=/data/project/mcpg/test_data/icbp/fastq/100908_s_6.fq.head1000000 # Not in GEO, so created from ICBP data ELAND export format
		mbd_bab10=/data/project/mcpg/test_data/icbp/fastq/100910_s_4.fq.head1000000


# set sample lists

mbd_list="..."
mbd_list=($mbd_lu1 $mbd_lu2 $mbd_lu3 $mbd_lu4 $mbd_lu5 $mbd_lu6 $mbd_lu7 $mbd_lu8 $mbd_lu9 $mbd_lu10 $mbd_lu11 $mbd_lu12 $mbd_lu13 $mbd_baa1 $mbd_baa2 $mbd_baa3 $mbd_baa4 $mbd_baa5 $mbd_baa6 $mbd_baa7 $mbd_bab1 $mbd_bab2 $mbd_bab3 $mbd_bab4 $mbd_bab5 $mbd_bab6 $mbd_bab7 $mbd_bab8 $mbd_bab9 $mbd_bab10)

#mbd_type_list=("lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "baa" "baa" "baa" "baa" "baa" "baa" "baa" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab")
#type_kind=("lu" "baa" "bab")
#sample_num_in_class=("13" "7" "8")
#fold_change=2


# CEL
		cel_lu1=/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_1_export.txt.CEL
		cel_lu2=/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_7_export.txt.CEL
		cel_lu3=/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_8_export.txt.CEL
		cel_lu4=/data/project/mcpg/test_data/icbp/gene_expression_array/100803_s_5_export.txt.CEL
		cel_lu5=/data/project/mcpg/test_data/icbp/gene_expression_array/100803_s_6_export.txt.CEL
		cel_lu6=/data/project/mcpg/test_data/icbp/gene_expression_array/100812_s_3_export.txt.CEL
		cel_lu7=/data/project/mcpg/test_data/icbp/gene_expression_array/100812_s_4_export.txt.CEL
		cel_lu8=/data/project/mcpg/test_data/icbp/gene_expression_array/100824_s_2_export.txt.CEL
		cel_lu9=/data/project/mcpg/test_data/icbp/gene_expression_array/100831_s_5_export.txt.CEL
		cel_lu10=/data/project/mcpg/test_data/icbp/gene_expression_array/100902_s_7_export.txt.CEL
		cel_lu11=/data/project/mcpg/test_data/icbp/gene_expression_array/100908_s_3_export.txt.CEL
		cel_lu12=/data/project/mcpg/test_data/icbp/gene_expression_array/100908_s_4_export.txt.CEL
		cel_lu13=/data/project/mcpg/test_data/icbp/gene_expression_array/100908_s_8_export.txt.CEL

		cel_baa1=/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_5_export.txt.CEL
		cel_baa2=/data/project/mcpg/test_data/icbp/gene_expression_array/100831_s_6_export.txt.CEL
		cel_baa3=/data/project/mcpg/test_data/icbp/gene_expression_array/100831_s_8_export.txt.CEL
		cel_baa4=/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_4_export.txt.CEL
		cel_baa5=/data/project/mcpg/test_data/icbp/gene_expression_array/100812_s_2_export.txt.CEL
		cel_baa6=/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_6_export.txt.CEL
		cel_baa7=/data/project/mcpg/test_data/icbp/gene_expression_array/100908_s_5_export.txt.CEL

		cel_bab1=/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_2_export.txt.CEL
		cel_bab2=/data/project/mcpg/test_data/icbp/gene_expression_array/100812_s_1_export.txt.CEL
		cel_bab3=/data/project/mcpg/test_data/icbp/gene_expression_array/100803_s_3_export.txt.CEL
		cel_bab4=/data/project/mcpg/test_data/icbp/gene_expression_array/100803_s_7_export.txt.CEL
		cel_bab5=/data/project/mcpg/test_data/icbp/gene_expression_array/100812_s_5_export.txt.CEL
		cel_bab6=/data/project/mcpg/test_data/icbp/gene_expression_array/100824_s_4_export.txt.CEL
		cel_bab7=/data/project/mcpg/test_data/icbp/gene_expression_array/100831_s_4_export.txt.CEL
		cel_bab8=/data/project/mcpg/test_data/icbp/gene_expression_array/100902_s_6_export.txt.CEL
		cel_bab9=/data/project/mcpg/test_data/icbp/gene_expression_array/100908_s_6_export.txt.CEL
		cel_bab10=/data/project/mcpg/test_data/icbp/gene_expression_array/100910_s_4_export.txt.CEL
		
		cel_list=($cel_lu1 $cel_lu2 $cel_lu3 $cel_lu4 $cel_lu5 $cel_lu6 $cel_lu7 $cel_lu8 $cel_lu9 $cel_lu10 $cel_lu11 $cel_lu12 $cel_lu13 $cel_baa1 $cel_baa2 $cel_baa3 $cel_baa4 $cel_baa5 $cel_baa6 $cel_baa7 $cel_bab1 $cel_bab2 $cel_bab3 $cel_bab4 $cel_bab5 $cel_bab6 $cel_bab7 $cel_bab8 $cel_bab9 $cel_bab10)
	
# SNP
		snp_lu1=/data/project/mcpg/test_data/icbp/fastq/100730_s_1.fq.vcf
		snp_lu2=/data/project/mcpg/test_data/icbp/fastq/100730_s_7.fq.vcf
		snp_lu3=/data/project/mcpg/test_data/icbp/fastq/100730_s_8.fq.vcf
		snp_lu4=/data/project/mcpg/test_data/icbp/fastq/100803_s_5.fq.vcf
		snp_lu5=/data/project/mcpg/test_data/icbp/fastq/100803_s_6.fq.vcf
		snp_lu6=/data/project/mcpg/test_data/icbp/fastq/100812_s_3.fq.vcf
		snp_lu7=/data/project/mcpg/test_data/icbp/fastq/100812_s_4.fq.vcf
		snp_lu8=/data/project/mcpg/test_data/icbp/fastq/100824_s_2.fq.vcf
		snp_lu9=/data/project/mcpg/test_data/icbp/fastq/100831_s_5.fq.vcf
		snp_lu10=/data/project/mcpg/test_data/icbp/fastq/100902_s_7.fq.vcf
		snp_lu11=/data/project/mcpg/test_data/icbp/fastq/100908_s_3.fq.vcf
		snp_lu12=/data/project/mcpg/test_data/icbp/fastq/100908_s_4.fq.vcf
		snp_lu13=/data/project/mcpg/test_data/icbp/fastq/100908_s_8.fq.vcf

		snp_baa1=/data/project/mcpg/test_data/icbp/fastq/100730_s_5.fq.vcf
		snp_baa2=/data/project/mcpg/test_data/icbp/fastq/100831_s_6.fq.vcf
		snp_baa3=/data/project/mcpg/test_data/icbp/fastq/100831_s_8.fq.vcf
		snp_baa4=/data/project/mcpg/test_data/icbp/fastq/100730_s_4.fq.vcf
		snp_baa5=/data/project/mcpg/test_data/icbp/fastq/100812_s_2.fq.vcf
		snp_baa6=/data/project/mcpg/test_data/icbp/fastq/100730_s_6.fq.vcf
		snp_baa7=/data/project/mcpg/test_data/icbp/fastq/100908_s_5.fq.vcf

		snp_bab1=/data/project/mcpg/test_data/icbp/fastq/100730_s_2.fq.vcf
		snp_bab2=/data/project/mcpg/test_data/icbp/fastq/100812_s_1.fq.vcf
		snp_bab3=/data/project/mcpg/test_data/icbp/fastq/100803_s_3.fq.vcf
		snp_bab4=/data/project/mcpg/test_data/icbp/fastq/100803_s_7.fq.vcf
		snp_bab5=/data/project/mcpg/test_data/icbp/fastq/100812_s_5.fq.vcf
		snp_bab6=/data/project/mcpg/test_data/icbp/fastq/100824_s_4.fq.vcf
		snp_bab7=/data/project/mcpg/test_data/icbp/fastq/100831_s_4.fq.vcf # Not in GEO, so created from ICBP data ELAND export format
		snp_bab8=/data/project/mcpg/test_data/icbp/fastq/100902_s_6.fq.vcf
		snp_bab9=/data/project/mcpg/test_data/icbp/fastq/100908_s_6.fq.vcf # Not in GEO, so created from ICBP data ELAND export format
		snp_bab10=/data/project/mcpg/test_data/icbp/fastq/100910_s_4.fq.vcf


# set sample lists

snp_list="..."
snp_list=($snp_lu1 $snp_lu2 $snp_lu3 $snp_lu4 $snp_lu5 $snp_lu6 $snp_lu7 $snp_lu8 $snp_lu9 $snp_lu10 $snp_lu11 $snp_lu12 $snp_lu13 $snp_baa1 $snp_baa2 $snp_baa3 $snp_baa4 $snp_baa5 $snp_baa6 $snp_baa7 $snp_bab1 $snp_bab2 $snp_bab3 $snp_bab4 $snp_bab5 $snp_bab6 $snp_bab7 $snp_bab8 $snp_bab9 $snp_bab10)


#######################################################################################


#########################################################
# Parameter settings
###################################################

# methylation data list
methyl_list=()
methyl_list=("${mbd_list[@]}")

# methylation data result directory
methyl_result_dir=""
if [ $methyl_type -eq 1 ]; then
	methyl_result_dir=$medips_from_mbd_result_dir
elif [ $methyl_type -eq 2 ]; then
	methyl_result_dir=$methylkit_from_bs_result_dir
else
	methyl_result_dir=""
fi

# gene expresison data list
rna_list=()
rna_list=("${cel_list[@]}")

# gene expression data result directory
rna_result_dir=""
if [ $ge_type -eq 1 ]; then
	rna_result_dir=$cel_result_dir
else
	rna_result_dir=$rnaseq_result_dir
fi

# mutation vcf list
mut_list=()
mut_list=("${snp_list[@]}")

snp_file_extension=$(echo `basename \${mut_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')

# mutation result directory
mutation_result_dir=""
mutation_result_dir="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/result/snp/from_mbd/"

if [ $mu_type -eq 1 ]; then # Human Map
	mutation_result_dir=""
elif [ $mu_type -eq 2 ]; then # RNA-seq
	mutation_result_dir=$result_dir"/mu/rna_seq/2pass/"
else
	mutation_result_dir=""
fi

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

	#################################################################
	# extract only bins in range & merge refseq and remove reseq with no gene symbol
	#################################################################
	count=0
	for (( i=0; i<${#methyl_list[@]}; i++ )); do
		temp_filename_only=`basename \${methyl_list[$i]}`
		met_level_file=$methyl_result_dir/$temp_filename_only".met_level"

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
	cut -f1-3 $integ_result_dir"/"`basename \${methyl_list[0]}`".RNGs.MGDrefS" > $integ_result_dir/temp_prefix.txt &
	cut -f6- $integ_result_dir"/"`basename \${methyl_list[0]}`".RNGs.MGDrefS" > $integ_result_dir/temp_pos.txt &
	wait


	paste_all_list=()
	for (( k=0; k<${#type_kind[@]}; k++ )); do 

		## rms per file
		paste_list=(); count=0
		for (( i=0; i<${#methyl_list[@]}; i++ )); do
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				work_file=$integ_result_dir/`basename \${methyl_list[$i]}`.methyl

				awk -v header=`basename \${methyl_list[$i]}`"-methyl" 'NR==1{print header;next;}{print $5}' $integ_result_dir"/"`basename \${methyl_list[$i]}`".RNGs.MGDrefS" > $work_file &
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
		for (( i=0; i<${#rna_list[@]}; i=i+1 )); do
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				work_file=$integ_result_dir/`basename \${rna_list[$i]}`".AVG_GE"

				echo `basename \${rna_list[$i]}`"-gene_exp" > $work_file
			
				temp_exp_input=""

				if [ $ge_type -eq 1 ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
					temp_exp_input=$rna_result_dir"/"`basename \${rna_list[\$i]}`".exp.geneSymbol_avg_exp.txt"
				else #Li1_head_100000_1.htseq.geneSymbol
					file_extension=$(echo `basename \${rna_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
					temp_exp_input=$rna_result_dir"/"`basename \${rna_list[0]} "."$file_extension`".htseq.geneSymbol"
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
	
	for (( i=0; i<${#mut_list[@]}; i++ )); do
		temp_filename_only=`basename \${mut_list[\$i]}`
		temp_filename_only_wo_extension=`basename \${mut_list[$i]} "."$snp_file_extension`
		snp_level_file=$mutation_result_dir/$temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf.bed"

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
		for (( i=0; i<${#mut_list[@]}; i++ )); do
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				SNP_MGD_LIST+=($integ_result_dir"/"`basename \${mut_list[$i]}`".SNP_RNGs.MGDrefS")
				SNP_MGD_LIST_all+=($integ_result_dir"/"`basename \${mut_list[$i]}`".SNP_RNGs.MGDrefS")
		
				temp_exp_input=""

				if [ $ge_type -eq 1 ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
					temp_exp_input=$rna_result_dir"/"`basename \${rna_list[\$i]}`".exp.geneSymbol_avg_exp.txt"
				else #Li1_head_100000_1.htseq.geneSymbol
					file_extension=$(echo `basename \${rna_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
					temp_exp_input=$rna_result_dir"/"`basename \${rna_list[0]} "."$file_extension`".htseq.geneSymbol"
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
	for (( i=0; i<${#methyl_list[@]}; i++ )); do
		work_file=$integ_result_dir/`basename \${methyl_list[$i]}`.methyl
		work_file2=$integ_result_dir/`basename \${rna_list[$i]}`.AVG_GE
		paste_list+=($work_file)
		paste_list2+=($work_file2)
	done 

	awk '{print $1"_"$2"_"$3}' $integ_result_dir/temp_prefix.txt | paste - ${paste_list[@]} | grep -v "NA" | uniq > $integ_result_dir/ME_MAT.txt &
	awk '{print $1"_"$2"_"$3}' $integ_result_dir/temp_prefix.txt | paste - ${paste_list2[@]} | grep -v "n\/a" | uniq > $integ_result_dir/GE_MAT.txt &
	wait
	###########################################
	# BY kruscal & and bonferroni
	###########################################

		python $bin_dir/kruskal_ver3.py $integ_result_dir/ME_MAT.txt `my_join "," ${TYPE_LEN[@]}` $kw_threshold $NUM_CPUS $integ_result_dir/ME_MAT.kw.txt $multiple_correction

#-a 0.05 -n $NUM_CPUS -c `my_join "," ${type_list[@]}` -m 0 --o1 $integ_result_dir/ME_MAT.kw.sigle_diff --o2 $integ_result_dir/ME_MAT.kw.all_diff

		python $bin_dir/kruskal_ver3.py $integ_result_dir/GE_MAT.txt `my_join "," ${TYPE_LEN[@]}` $kw_threshold $NUM_CPUS $integ_result_dir/GE_MAT.kw.txt $multiple_correction

#-a 0.05 -n $NUM_CPUS -c `my_join "," ${type_list[@]}` -m 0 --o1 $integ_result_dir/GE_MAT.kw.sigle_diff --o2 $integ_result_dir/GE_MAT.kw.all_diff
	
		# adjusted pvalue selection
		awk -v kw_threshold=$kw_threshold '{if(($(NF) < kw_threshold) && ($(NF) > 0)){split($1, tokens, "_"); print tokens[1]"\t"tokens[2]"\t"tokens[3]}}' $integ_result_dir/ME_MAT.kw.txt  > $integ_result_dir/GDMR_kw &
		awk -v kw_threshold=$kw_threshold '{if(($(NF) < kw_threshold) && ($(NF) > 0)){split($1, tokens, "_"); print tokens[1]"\t"tokens[2]"\t"tokens[3]}}' $integ_result_dir/GE_MAT.kw.txt  > $integ_result_dir/GDEG_kw &
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
	
	# compute % of over correlation threshold
	echo "[INFO] Compute % of correlation over threshold - GDMR & GDEG"
	python $bin_dir/profile_ME_GE.py $GDEG_GDMR_MGD $GDEG_GDMR_MGD_CORR ${#type_list[@]} $corr_threshold $integ_result_dir/GDMR_GDEG_MAT_THRESH.txt $integ_result_dir/GDMR_GDEG_MAT_PNC.txt
	

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
	num_GDSNP=`wc $integ_result_dir/GDSNP.txt | awk '{print $1}'`
	num_GDSNP_bins_related_w_GDEG=`cut -f1,2,3 $integ_result_dir/GDEG_GDSNP.txt | uniq | wc | awk '{print $1}'`
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


exit

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
			SUBTYPE_DMR=$methyl_result_dir"/"${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.bed"
			# chr start end
			SUBTYPE_DEG=$rna_result_dir"/"${type_kind[$k]}"_vs_"${type_kind[$l]}".DEG.list"
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
			dmr_file=$methyl_result_dir/$dmr_prefix".mr.edgeR.s.txt" 

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
			deg_file=$rna_result_dir/${type_kind[$i]}"_vs_"${type_kind[$j]}".deg"

			# init lists
			SNP_MGD_LIST=(); GE_MGD_LIST=(); temp_subtype_list=()
			for (( k=0; k<${#mut_list[@]}; k++ )); do 

				snp_file=$integ_result_dir"/"`basename \${mut_list[$k]}`".SNP_RNGs.MGDrefS"
				deg_only_snp_file=$snp_file"."$TYPE_PAIR"_deg_only"
				exp_file=""

				if [ $ge_type -eq 1 ]; then #100730_s_1_export.txt.CEL.exp.geneSymbol_avg_exp.txt
					temp_exp_input=$rna_result_dir"/"`basename \${rna_list[\$i]}`".exp.geneSymbol_avg_exp.txt"
				else #Li1_head_100000_1.htseq.geneSymbol
					file_extension=$(echo `basename \${rna_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
					temp_exp_input=$rna_result_dir"/"`basename \${rna_list[0]} "."$file_extension`".htseq.geneSymbol"
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
			} > $integ_result_dir/DSNP_DEG_stat_$TYPE_PAIR_PREF.txt
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
		MBD_methyl_file=$methyl_result_dir/`basename \${bs_paired_methyl_list[$i]}`".met_level";
		MBD_bam_file=$mbd_result_dir/`basename \${bs_paired_methyl_list[$i]}`".sorted.bam"


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
		work_file=$integ_result_dir/bs_related_MBD_`basename \${bs_paired_methyl_list[$i]}`;
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

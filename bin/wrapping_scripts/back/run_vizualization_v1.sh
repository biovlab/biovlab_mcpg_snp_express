#!/bin/bash
source `dirname $0`/../../env.sh

# directories
bin_dir="$WORK_DIR/bin"
result_dir="$WORK_DIR/result"
profile_result_dir="$WORK_DIR/profile"
viz_result_dir=$result_dir"/visualization"
methylkit_from_bs_result_dir=$result_dir"/methyl/bs/methylkit"
medips_from_mbd_result_dir=$result_dir"/methyl/mbd/medips"
snp_from_mbd_result_dir=$result_dir"/snp/from_mbd"
mbd_result_dir=$result_dir"/methyl/mbd"
circos_dir=$viz_result_dir"/circos"
rnaseq_result_dir=$result_dir"/gene_exp/rna_seq/deg"
deg_analysis_dir=$rnaseq_result_dir"/deg/"

# variables
ME_all_merged=$profile_result_dir/ME_MGD		# prepare initiall prefix for all merged file

# functions
function my_join { local IFS="$1"; shift; echo "$*"; }
####################################################################################
# 0. set parameters
####################################################################################

# TODO:set user-specific result dir to copy the result to user-specific web accessible dir

# MBD
		mbd_lu1=/data/project/mcpg/test_data/icbp/fastq/100730_s_1.fq 
		mbd_lu2=/data/project/mcpg/test_data/icbp/fastq/100730_s_7.fq 
		mbd_lu3=/data/project/mcpg/test_data/icbp/fastq/100730_s_8.fq 
		mbd_lu4=/data/project/mcpg/test_data/icbp/fastq/100803_s_5.fq
		mbd_lu5=/data/project/mcpg/test_data/icbp/fastq/100803_s_6.fq
		mbd_lu6=/data/project/mcpg/test_data/icbp/fastq/100812_s_3.fq
		mbd_lu7=/data/project/mcpg/test_data/icbp/fastq/100812_s_4.fq
		mbd_lu8=/data/project/mcpg/test_data/icbp/fastq/100824_s_2.fq
		mbd_lu9=/data/project/mcpg/test_data/icbp/fastq/100831_s_5.fq
		mbd_lu10=/data/project/mcpg/test_data/icbp/fastq/100902_s_7.fq
		mbd_lu11=/data/project/mcpg/test_data/icbp/fastq/100908_s_3.fq
		mbd_lu12=/data/project/mcpg/test_data/icbp/fastq/100908_s_4.fq
		mbd_lu13=/data/project/mcpg/test_data/icbp/fastq/100908_s_8.fq

		mbd_baa1=/data/project/mcpg/test_data/icbp/fastq/100730_s_5.fq 
		mbd_baa2=/data/project/mcpg/test_data/icbp/fastq/100831_s_6.fq 
		mbd_baa3=/data/project/mcpg/test_data/icbp/fastq/100831_s_8.fq 
		mbd_baa4=/data/project/mcpg/test_data/icbp/fastq/100730_s_4.fq
		mbd_baa5=/data/project/mcpg/test_data/icbp/fastq/100812_s_2.fq
		mbd_baa6=/data/project/mcpg/test_data/icbp/fastq/100730_s_6.fq
		mbd_baa7=/data/project/mcpg/test_data/icbp/fastq/100908_s_5.fq

		mbd_bab1=/data/project/mcpg/test_data/icbp/fastq/100730_s_2.fq 
		mbd_bab2=/data/project/mcpg/test_data/icbp/fastq/100812_s_1.fq 
		mbd_bab3=/data/project/mcpg/test_data/icbp/fastq/100803_s_3.fq
		mbd_bab4=/data/project/mcpg/test_data/icbp/fastq/100803_s_7.fq
		mbd_bab5=/data/project/mcpg/test_data/icbp/fastq/100812_s_5.fq
		mbd_bab6=/data/project/mcpg/test_data/icbp/fastq/100824_s_4.fq
		mbd_bab7=/data/project/mcpg/test_data/icbp/fastq/100831_s_4.fq # Not in GEO, so created from ICBP data ELAND export format
		mbd_bab8=/data/project/mcpg/test_data/icbp/fastq/100902_s_6.fq
		mbd_bab9=/data/project/mcpg/test_data/icbp/fastq/100908_s_6.fq # Not in GEO, so created from ICBP data ELAND export format
		mbd_bab10=/data/project/mcpg/test_data/icbp/fastq/100910_s_4.fq


# set sample lists

mbd_list="..."
mbd_list=($mbd_lu1 $mbd_lu2 $mbd_lu3 $mbd_lu4 $mbd_lu5 $mbd_lu6 $mbd_lu7 $mbd_lu8 $mbd_lu9 $mbd_lu10 $mbd_lu11 $mbd_lu12 $mbd_lu13 $mbd_baa1 $mbd_baa2 $mbd_baa3 $mbd_baa4 $mbd_baa5 $mbd_baa6 $mbd_baa7 $mbd_bab1 $mbd_bab2 $mbd_bab3 $mbd_bab4 $mbd_bab5 $mbd_bab6 $mbd_bab7 $mbd_bab8 $mbd_bab9 $mbd_bab10)

ME_input_type=1 #(1:MBD, 2:BS, 3:infinum...) 
GE_input_type=1 #(1:RNAseq, 2:microarray...) 
MU_input_type=1 #(1:RNAseq, 2:DNAseq, 3:humanhap...) 
type_list=("lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "lu" "baa" "baa" "baa" "baa" "baa" "baa" "baa" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab" "bab")
type_kind=("lu" "baa" "bab")
sample_num_in_class=("13" "7" "10")
fold_change=2


#########################################################################
mbd_result_dir="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/small_test_data/targeted_bs_breast_cancer_80genes/small_test/methyl/mbd/"
mbd_list=("BrCa-02_head10000000.fastq" "BrCa-03_head10000000.fastq")
type_kind=("A" "B")
type_list=("A" "B")


#########################################################
# Parameter settings
###################################################

# methylation data list
methyl_list=()
methyl_list=("${mbd_list[@]}")

# methylation data result directory
methyl_result_dir=""
if [ $ME_input_type -eq 1 ]; then
	methyl_result_dir=$medips_from_mbd_result_dir
elif [ $ME_input_type -eq 2 ]; then
	methyl_result_dir=$methylkit_from_bs_result_dir
else
	methyl_result_dir=""
fi

# gene expresison data list
rna_list=()
#rna_list=("${cel_list[@]}")

# gene expression data result directory
rna_result_dir=""
if [ $GE_input_type -eq 1 ]; then
	rna_result_dir=$cel_result_dir
else
	rna_result_dir=$rnaseq_result_dir
fi

# mutation vcf list
vcf_list=()
#vcf_list=("${snp_list[@]}")

# mutation result directory
mutation_result_dir=""
mutation_result_dir="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/result/snp/from_mbd/"



####################################################################################
# 6. Visualization 
####################################################################################

mkdir -p $viz_result_dir

echo "[INFO] Start visualizing results"
#################################
# 6.1 For Raw data. (NGSplot, UCSC links, ME-density)
#################################

############
# 6.1.1 GE 
############

<<'COMMENT_TEMP'
# TODO : set data kinds (RNAseq or array ..)
# UCSC links for each samples : BigWig
for (( i=0; i<${#rna_list[@]}; i++ )); do 
	# copy to web accessible directory
	# TODO :  this may be not necessary cuz already we are in the web accessible dir
	cp $rna_result_dir/`basename \${rna_list[$i]}`.bw $WEB_ACCESSIBLE_DIR
	# print UCSC url
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=test_SNP%20description=test2SNP%20visibility=dense%20bigDataUrl=$WEB_ACCESSIBLE_LOC/${rna_list[$i]}.bw"
done

COMMENT_TEMP


############
# 6.1.2 ME
############
<<'COMMENT_TEMP'
# UCSC links for each samples : BigWig
if [ $ME_input_type -eq 1 ]; then # MBD
	count=0
	for (( i=0; i<${#mbd_list[@]}; i++ )); do
		{
			met_level_file=$medips_from_mbd_result_dir/`basename \${mbd_list[$i]}`".met_level";

			# create bigWig file for UCSC
			awk 'NR==1{next;}{print}' $met_level_file | cut -f1,2,3,5 | bedtools intersect -a - -b $REF_HUMAN_GENOME_100_BIN_BED -sorted > $viz_result_dir/`basename \${mbd_list[$i]}`".met_level.bed" ;
			$bin_dir/bedGraphToBigWig $viz_result_dir/`basename \${mbd_list[$i]}`".met_level.bed" $REF_HUMAN_CHR_SIZE $viz_result_dir/`basename \${mbd_list[$i]}`".met_level.bw";

			cp $viz_result_dir/`basename \${mbd_list[$i]}`".met_level.bw" $WEB_ACCESSIBLE_DIR;
			# print UCSC url
			echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl=$WEB_ACCESSIBLE_LOC/`basename \${mbd_list[$i]}`.met_level.bw";
		} &
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait

elif [ $ME_input_type -eq 2 ]; then # BS
	ls $methylkit_from_bs_result_dir/*.sorted.bw | xargs -I {} cp {} $WEB_ACCESSIBLE_DIR
	ls $methylkit_from_bs_result_dir/*.sorted.bw | xargs -I {} echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl="$WEB_ACCESSIBLE_LOC/{};
fi

COMMENT_TEMP


# NGS plot
if [ $ME_input_type -eq 1 ] ; then # MBD
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		echo -n "" > $viz_result_dir/ngsplot_config_ME.txt
		echo -n "" > $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt


		for (( i=0; i<${#mbd_list[@]}; i++ )); do 
			temp_filename_only=`basename \${mbd_list[$i]}`

			# create configure file for all multiple file plotting
			echo -e "$mbd_result_dir/$temp_filename_only.sorted.bam\t-1\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_all.txt
			#echo -e "$mbd_result_dir/$temp_filename_only.sorted.bam\tintron.bed\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_all.txt

			# create configure file for subtype multiple file plotting
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				echo -e "$mbd_result_dir/$temp_filename_only.sorted.bam\t-1\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt
				#echo -e "$mbd_result_dir/$temp_filename_only.sorted.bam\tintron.bed\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt
			fi
		done

		# NGS plot itself is parallelized, so do not need to be parallelzed
		for region in "cgi" "genebody" "enhancer" "dhs" "cgi" "exon"; do 
<<'COMMENT'
			# THIS IS HARDCODED FOR ICBP PROJECT!!!
			if [ "$region" == "genebody" ]; then
				YAS="0.02,0.13"
				SC="0,0.5" # NOTE : NOT YET TESTED
			elif [ "$region" == "exon" ]; then
				YAS="0.04,0.21"
				SC="0,0.5" # NOTE : NOT YET TESTED
			elif [ "$region" == "cgi" ]; then
				YAS="0.01,0.15"
				SC="0,0.5" # NOTE : NGSPLOT HEATMAP NEED TO BE RESCALE TO FIT ALL THRE FIGURE HAVING SAME SCALE!!!
			elif [ "$region" == "dhs" ]; then
				YAS="0.02,0.15"
				SC="0,0.5" # NOTE : NOT YET TESTED
			else
				YAS="0.00,0.20"
				SC="0,0.5" # NOTE : NOT YET TESTED
			fi
COMMENT			
			ngs.plot.r -G hg19 -R $region -C $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt -O $viz_result_dir/ME_${type_kind[$k]}_$region 
		done
	done

	# for all sample
	for region in "bed" ;do #"genebody" "enhancer" "dhs" "cgi"; do
		/packages/test2/ngsplot-develop/bin/ngs.plot.r -SC $SC -G hg19 -R $region -C $viz_result_dir/ngsplot_config_ME_all.txt -O $viz_result_dir/ME_all_$region -E /data/project/mcpg/lib/region_info/intron.bed
	done

elif [ $ME_input_type -eq 2 ] ; then # BS
	# TODO : no figure for BS?
fi

exit

# TODO : ME density plot and comparison



############
# 6.1.3 MU
############

# UCSC links
for (( i=0; i<${#mbd_list[@]}; i++ )); do
	cp $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.gz" $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.gz.tbi" $WEB_ACCESSIBLE_DIR
	# print UCSC url
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=vcfTabix%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl=$WEB_ACCESSIBLE_LOC/`basename \${mbd_list[$i]}`.vcf.gz"
done




#################################
# 6.2 For DEG, DMR, DSNP (Circos, Heatmap, region plot, Mu-Ge plot)
#################################

<<'COMMENT'

############
# 6.2.1 GE 
############

# depic in the circos as DEGs

#########################################
# DEG : gene list -> info link
#########################################
	for (( k=0; k<${#CEL_TYPE_KIND[@]}; k++ )); do 
		for (( l=$k+1; l<${#CEL_TYPE_KIND[@]}; l++ )); do 
			work_file=${CEL_TYPE_KIND[$k]}"_vs_"${CEL_TYPE_KIND[$l]}".deg"
			
		
			# get deg gene list only
			cut -f2 $cel_result_dir/$work_file > $viz_result_dir/$work_file".gene_only"
			sh $bin_dir/create_annotation_link_html_based_on_gene_list.sh $viz_result_dir/$work_file".gene_only" > $viz_result_dir/$work_file".gene_only.annot.html"

			# create annotation page for deg
			cp $cel_result_dir/$work_file $viz_result_dir/$work_file".gene_only" $viz_result_dir/$work_file".gene_only.annot.html" $WEB_ACCESSIBLE_DIR
		done
	done



############
# 6.2.2 ME 
############

#########################################
# circos
#########################################
mkdir -p $circos_dir
echo "generate circos for subtype-specific avg methylation with difference and DEGs"

# create subtype avg methyl file from merged file
for (( i=0; i<${#type_kind[@]}; i++ )); do 
	#subtype_ME_merged_file=$ME_all_merged"_"${type_kind[$i]}
	subtype_ME_merged_file="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/temp_result_files/ME_MGD_"${type_kind[$i]}".head100" # temp test with samll datda
	subtype_ME_merged_avg_file=$viz_result_dir/`basename \$subtype_ME_merged_file`"_avg"
	subtype_ME_merged_avg_10Mb_file=$subtype_ME_merged_avg_file"_10Mb"


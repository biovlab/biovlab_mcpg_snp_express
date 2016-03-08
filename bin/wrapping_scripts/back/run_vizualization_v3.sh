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


	# file example
	#chr	bin_start	bin_end	100730_s_1.fq-rms	100730_s_7.fq-rms...	range_start	range_end	strand	gene_symbol	range_kind	refseq
	#chr1	9801	9900	0	0	0	0	0	0	0	0	0	0	0	0	0	9873	11873	+	DDX11L1	promoter	NR_046018
	#chr1	9901	10000	0	0	0	0	0	0	0	0	0	0	0	0	0	9873	11873	+	DDX11L1	promoter	NR_046018

	# TODO : check ME input type. file format should be same for all input type as above
	# compute avg methyl value
	temp_sample_num_in_class=${#sample_num_in_class[$i]}
	awk -v sample_num_in_class=$temp_sample_num_in_class 'NR==1{print "#chr","start","end","avg_value";next;}{avg_value=0; for (i=1+3; i<=sample_num_in_class+3; i++){avg_value+=$i;}; avg_value=avg_value/sample_num_in_class; print $1,$2,$3,avg_value;}' OFS='\t' $subtype_ME_merged_file > $subtype_ME_merged_avg_file

	# merge value the range in 10Mb for circos input
	bedtools map -header -c 4 -o mean -a $REF_HUMAN_GENOME_10M_BIN_BED -b $subtype_ME_merged_avg_file | sed -e 's/chr/hs/g' | sed -e 's/\t\./\t0/g' > $subtype_ME_merged_avg_10Mb_file
done



# crate circos for pairwise & all class

all_type_file_list=()
for (( i=0; i<${#type_kind[@]}; i++ )); do 
	for (( j=$i+1; j<${#type_kind[@]}; j++ )); do 

		# met level file
		# generate circos input list deliminated by ';'	
		subtype_ME_merged_file="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/temp_result_files/ME_MGD_"${type_kind[$i]}".head100" # temp test with samll datda
		subtype_ME_merged_file_2="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/temp_result_files/ME_MGD_"${type_kind[$j]}".head100" # temp test with samll datda
		subtype_ME_merged_avg_file=$viz_result_dir/`basename \$subtype_ME_merged_file`"_avg"
		subtype_ME_merged_avg_file_2=$viz_result_dir/`basename \$subtype_ME_merged_file_2`"_avg"
		subtype_ME_merged_avg_10Mb_file=$subtype_ME_merged_avg_file"_10Mb"
		subtype_ME_merged_avg_10Mb_file_2=$subtype_ME_merged_avg_file_2"_10Mb"


		# for deg
	  work_file=${type_kind[$i]}"_vs_"${type_kind[$j]}
    #DEG_file_name=$deg_analysis_dir/$work_file".FC"$fold_change".DEGs"
    DEG_file_name=/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/rna_paired_test_data/gene_exp/rna_seq/deg/Li_vs_Th.FC2.DEGs
		DEG_for_circos=$circos_dir/$work_file".FC"$fold_change".DEGs.circos"
		diff_file=$circos_dir"/circos_${type_kind[$i]}_${type_kind[$j]}.diff"
		file_list=$subtype_ME_merged_avg_10Mb_file";"$subtype_ME_merged_avg_10Mb_file_2


		# TODO : check the cut off for degs file name

		# convert deg file to circos file (add position information and change chr to hs)

		# genebody info
		#	chr1	11873	14409	DDX11L1	NR_046018	+
		# chr1	14361	29370	WASH7P	NR_024540	-

		join -1 1 -2 4 <(sort $DEG_file_name) <(sort -k4 $REF_HUMAN_GENE_GENEBODY_INFO) | awk '{print $2,$3,$4,$1}' OFS='\t' | sed -e 's/chr/hs/g' > $DEG_for_circos

		# SLCO1B1 chr12 21284127 21392730 NM_006446 +
		# SLCO1B3 chr12 20963637 21069843 NM_019844 +kI

		# create diff file
		# chr1  0 10000000  5.532632995
		paste $subtype_ME_merged_avg_10Mb_file $subtype_ME_merged_avg_10Mb_file_2 | awk '{print $1,$2,$3,$4-$8}' OFS='\t' > $diff_file
		
		# create subtype pair-wise circos
		# TODO : set min max value for circos
		circos_min=0
		circos_max=15
		temp_circos_file=$circos_dir"/circos_${type_kind[$i]}_${type_kind[$j]}.conf"
		sh $bin_dir/generate_circos_conf.sh 2 $file_list $diff_file $DEG_for_circos $CIRCOS_CONFIG_LOC $circos_min $circos_max > $temp_circos_file
		circos -conf $temp_circos_file -outputdir $circos_dir -outputfile $work_file
	done
	all_type_file_list+=($subtype_ME_merged_avg_10Mb_file)
done
all_type_file_list+=($subtype_ME_merged_avg_10Mb_file_2)

# TODO : snagseon
# measure entropy for all diff file, set proper min mac for diff as well... maybe it should be set seperately in generate_circos_conf.
diff_all_file=
# create all subtype circos
circos_min=0
circos_max=15
temp_circos_file=$circos_dir"/circos_all.conf"
sh $bin_dir/generate_circos_conf.sh ${#type_kind[@]} $all_type_file_list $diff_all_file $CIRCOS_CONFIG_LOC $circos_min $circos_max > $temp_circos_file
circos -conf $temp_circos_file -outputdir $circos_dir -outputfile "all"




########################################
# DMR : Fold change -> bedgraph -> bw
########################################

for (( k=0; k<${#type_kind[@]}; k++ )); do 
	for (( l=$k+1; l<${#type_kind[@]}; l++ )); do 
		# create subtype fold change DMR region based on 'subtype1_vs_subtype2.mr.edgeR.s.txt'
		awk 'NR==1{next;}{print $1"\t"$2"\t"$3"\t"$(NF-3)}' $medips_from_mbd_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".mr.edgeR.p0.05.txt"  > $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.FC"
		bedGraphToBigWig $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.FC" $REF_HUMAN_CHR_SIZE $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.FC.bw"

		# copy to web accessible directory
		cp $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.FC.bw" $WEB_ACCESSIBLE_DIR

		# print UCSC url
		echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=${type_kind[$k]}_${type_kind[$l]}_DMR%20description=${type_kind[$k]}_${type_kind[$l]}_DMR%20visibility=dense%20bigDataUrl="$WEB_ACCESSIBLE_LOC/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.FC.bw"
	done
done

########################################
# GDMR : bed -> bb
########################################
	bedToBigBed $PROF_DIR/GDMR_kw $REF_HUMAN_CHR_SIZE $viz_result_dir/GDMR_by_kw.bb
	
	# copy to web accessible directory
	cp $viz_result_dir/GDMR_by_kw.bb $WEB_ACCESSIBLE_DIR

	# print UCSC url
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigBed%20name=GDMR_kw%20description=GDMR_kw%20visibility=dense%20bigDataUrl="$WEB_ACCESSIBLE_LOC/GDMR_by_kw.bb 


########################################
# DMR all : heatmap
########################################
	# get top 10% by pvalue (NF -1). NF:post hoc test(0:fail 1:pass)
	pvalue_index=`awk '{print NF-1;exit}' /data/project/mcpg/result/profile/ME_MAT.kw.sigle_diff`
	line_num=`grep -v NaN /data/project/mcpg/result/profile/ME_MAT.kw.sigle_diff | sort -k"$pvalue_index","$pvalue_index"n | tee $viz_result_dir/ME_MAT.kw.sigle_diff.sorted | wc -l`
	top_10_line_num=$(( $line_num / 10 ))

	echo "total line : " $line_num
	echo "10% line : " $top_10_line_num

	head -n $top_10_line_num $viz_result_dir/ME_MAT.kw.single_diff.sorted > $viz_result_dir/ME_MAT.kw.single_diff.sorted.top10
	#grep -v NaN /data/project/mcpg/result/profile/ME_matrix.kw.sigle_diff | awk '{print NF-1;exit}' | xargs -I {} sort -k{},{}n /data/project/mcpg/result/profile/ME_matrix.kw.sigle_diff
	# get

	# create heatmap
	Rscript $bin_dir/make_heatmap_flex.r "$viz_result_dir/ME_matrix.kw.sigle_diff.sorted.top10" `my_join "," ${mbd_type_list[@]}` $sample_num `my_join "," ${TYPE_LEN[@]}` 100 "$viz_result_dir/DMR_heatmap.png" "Top_10_DMRs_by_Kruscal_Wallis"

	# copy to WEB_ACCESSIBLE_DIR
	cp $viz_result_dir/DMR_heatmap.png $WEB_ACCESSIBLE_DIR



########################################
# DMR region stat : pie chart
########################################
	echo "[INFO] DMR region stat on bar chart"

	echo -e "Region\tcounts\tsubtype_pair" > $viz_result_dir/DMR_region_stat_matrix.txt
	echo -e "Region\tcounts\tsubtype_pair" > $viz_result_dir/DMR_cpgi_stat_matrix.txt



	echo "[INFO] count number of DMR in each genomic region"
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do 
			subtype_pair="${type_kind[$k]}_vs_${type_kind[$l]}"
			#for region in "3 prime UTR" "5 prime UTR" "Promoter" "Exon" "Intron"; do
			for region in "3utr" "5utr" "promoter" "exon" "intron"; do
				temp_region_dmr_count=`bedtools intersect -wa -wb -a $viz_result_dir/$subtype_pair.DMR.FC -b $human_refseq_various_ranges | grep -w $region | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`

				echo -e $region"\t"$temp_region_dmr_count"\t"$subtype_pair >> $viz_result_dir/DMR_region_stat_matrix.txt
			done
				# create cpg_island, cpg_shore, cpg_shelf
			for region in "CpG Island" "CpGI Shore" "CpGI Shelf"; do
				temp_region_dmr_count=`bedtools intersect -wa -wb -a $viz_result_dir/$subtype_pair.DMR.FC -b $human_refseq_various_ranges | grep -w $region | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`
				echo -e $region"\t"$temp_region_dmr_count"\t"$subtype_pair >> $viz_result_dir/DMR_cpgi_stat_matrix.txt
			done
		done
	done

	# draw pie chart
	echo "[INFO] Draw bar chart"
	Rscript $bin_dir/piechart.r $viz_result_dir/DMR_region_stat_matrix.txt $viz_result_dir/DMR_region_stat_barplot.jpg "Number of differentially methylated bins"
	Rscript $bin_dir/piechart.r $viz_result_dir/DMR_cpgi_stat_matrix.txt $viz_result_dir/DMR_cpgi_stat_barplot.jpg "Number of differentially methylated bins"

	cp $viz_result_dir/DMR_region_stat_barplot.jpg $viz_result_dir/DMR_cpgi_stat_barplot.jpg $WEB_ACCESSIBLE_DIR




############
# 6.2.3 MU 
############

													# SNP Region stat NO need?
													#NUM_CPUS=5
													#count=0 
													#for (( i=0; i<${#mbd_list[@]}; i++ )); do 
													#	bedtools intersect -wa -wb -a $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.bed.all.cov" -b $human_refseq_various_ranges > $viz_result_dir/`basename \${mbd_list[$i]}`".SNP_bin_ranges" &
													#	let count+=1
													# 	[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
													#done

# create prefix ID for stat matrix
{
echo -e "relative_range\n5p"; 
for ((i=1; i<99; i++)); do 
	echo $i"%";
done; 
echo "3p";
} > $viz_result_dir/prefix_percentage_ids

for range in "genebody" "cpgIsland" "exon"; do 
	grep $range $REF_HUMAN_GENE_RANGE_INFO_DIV_100 > $viz_result_dir/temp_range_div_100 # get range bins(1-100)

	# For GDSNP
	work_file_for_GDSNP=$viz_result_dir/GSNP_$range"_MAT"
	echo "GDSNP_in_$range" > $work_file_for_GDSNP

	grep -w $range $PROF_DIR/GDSNP.txt | cut -f1,2,3 - | bedtools coverage -b $viz_result_dir/temp_range_div_100 -a - -counts | cut -f4,6  | sort -k1,1n | bedtools groupby -g 1 -c 2 -o sum | cut -f2 >> $work_file_for_GDSNP 

	paste $viz_result_dir/prefix_percentage_ids $work_file > $viz_result_dir/GSNP_$range"_matrix4plot_merged"
	Rscript $bin_dir/plot.r $viz_result_dir/GSNP_$range"_matrix4plot_merged" $viz_result_dir/GSNP_$range"_matrix4plot_merged.jpg"


	# for DSNP
	paste_list_all=()
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		paste_list_sub_type=()
		count=0
		for (( i=0; i<${#mbd_list[@]}; i++ )); do 

			# create range stat file per sample
			work_file=$viz_result_dir/SNP_$range"_MAT_"`basename \${mbd_list[$i]}`

			if [ ${mbd_type_list[$i]} == ${type_kind[$k]} ]; then
				
				echo "[INFO] Estimate relative coverage in range of all $range region in `basename \${mbd_list[$i]}`"

				paste_list_sub_type+=($work_file)
				paste_list_all+=($work_file)
				echo "`basename \${mbd_list[$i]}`" > $work_file

				# parallization
				cut -f1,2,3 $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.bed" | bedtools coverage -b $viz_result_dir/temp_range_div_100 -a - -counts | cut -f4,6  | sort -k1,1n | bedtools groupby -g 1 -c 2 -o sum | cut -f2 >> $work_file &

				# NOTE : MEMORY INTENSE, so ONLY USE 5 CPUS			
				let count+=1; [[ $((count%5)) -eq 0 ]] && wait
			fi
		done; wait
		paste $viz_result_dir/prefix_percentage_ids ${paste_list_sub_type[@]} > $viz_result_dir/SNP_$range"_MAT_"${type_kind[$k]}

		Rscript $bin_dir/plot.r $viz_result_dir/SNP_$range"_MAT_"${type_kind[$k]} $viz_result_dir/SNP_$range"_"${type_kind[$k]}".jpg"
	done

	# merge all stat files
	paste $viz_result_dir/prefix_percentage_ids ${paste_list_all[@]} > $viz_result_dir/SNP_$range"_MAT_all"

	Rscript $bin_dir/plot.r $viz_result_dir/SNP_$range"_matrix4plot_all" $viz_result_dir/SNP_$range"_all.jpg"

	cp $viz_result_dir/SNP_*.jpg $WEB_ACCESSIBLE_DIR
done


# GDSNP all :  bed w/ label -> bb
	sort -k1,1 -k2,2n $PROF_DIR/GDSNP.bed | uniq > $viz_result_dir/GDSNP_sorted.bed
	bedToBigBed $viz_result_dir/GDSNP_sorted.bed $REF_HUMAN_CHR_SIZE $viz_result_dir/GDSNP_sorted.bb

	# copy BigBed to web accessible dir
	cp $viz_result_dir/GDSNP_sorted.bb $WEB_ACCESSIBLE_DIR

	# create UCSC link
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigBed%20name=GDSNP%20description=GDSNP%20visibility=pack%20bigDataUrl=$WEB_ACCESSIBLE_LOC/GDSNP_sorted.bb"

# GDSNP region stat : barplot
	echo "[INFO] GDSNP region stat on bar chart"
	echo -e "Region\tcounts\tSubtype" > $viz_result_dir/GDSNP_RNG_STAT_MAT.txt
	echo -e "Region\tcounts\tSubtype" > $viz_result_dir/GDSNP_CPGs_STAT_MAT.txt
	echo -e "Region\tcounts\tSubtype" > $viz_result_dir/GDSNP_ALL_STAT_MAT.txt

	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for region in "3utr" "5utr" "promoter" "exon" "intron" "genebody"; do
			temp_region_gdsnp_count=`bedtools intersect -wa -wb -a $viz_result_dir/GDSNP_sorted.bed -b $human_refseq_various_ranges | grep -w $region | grep -w ${type_kind[$k]} | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`

			echo -e $region"\t"$temp_region_gdsnp_count"\t"${type_kind[$k]} >> $viz_result_dir/GDSNP_RNG_STAT_MAT.txt
			echo -e $region"\t"$temp_region_gdsnp_count"\t"${type_kind[$k]} >> $viz_result_dir/GDSNP_ALL_STAT_MAT.txt
		done

		for region in "CpG Island" "CpGI Shelf" "CpGI Shore"; do 
			temp_region_gdsnp_count=`bedtools intersect -wa -wb -a $viz_result_dir/GDSNP_sorted.bed -b $human_refseq_various_ranges | grep -w $region | grep -w ${type_kind[$k]} | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`
			echo -e $region"\t"$temp_region_gdsnp_count"\t"${type_kind[$k]} >> $viz_result_dir/GDSNP_CPGs_STAT_MAT.txt
			echo -e $region"\t"$temp_region_gdsnp_count"\t"${type_kind[$k]} >> $viz_result_dir/GDSNP_ALL_STAT_MAT.txt
		done
	done

	# draw pie chart
	echo "[INFO] Draw bar chart"
	Rscript $bin_dir/piechart.r $viz_result_dir/GDSNP_RNG_STAT_MAT.txt $viz_result_dir/BAR_GDSNP_RNG.pdf "Number of tumor subtype specific mutation"
	Rscript $bin_dir/piechart.r $viz_result_dir/GDSNP_CPGs_STAT_MAT.txt $viz_result_dir/BAR_GDSNP_CPGs.pdf "Number of tumor subtype specific mutation"
	Rscript $bin_dir/barplot.r $viz_result_dir/GDSNP_ALL_STAT_MAT.txt $viz_result_dir/BAR_GDSNP_ALL.pdf "Number of tumor subtype specific mutation"

	cp $viz_result_dir/BAR_GDSNP_RNG.jpg $viz_result_dir/GDSNP_RNG_STAT_MAT.txt $viz_result_dir/GDSNP_CPGs_STAT_MAT.txt $viz_result_dir/GDSNP_ALL_STAT_MAT.txt $viz_result_dir/BAR_GDSNP_ALL.jpg $viz_result_dir/BAR_GDSNP_CPGs.jpg $WEB_ACCESSIBLE_DIR


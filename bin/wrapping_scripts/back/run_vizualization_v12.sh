#!/bin/bash
source `dirname $0`/../../env.sh

# parse experiment information
source `dirname $0`/parse_exp_info.sh

# directories
bin_dir="$WORK_DIR/bin"
result_dir="$WORK_DIR/result"
#profile_result_dir="$WORK_DIR/profile"
integ_result_dir=$result_dir"/integrated"
viz_result_dir=$result_dir"/visualization"
mbd_result_dir=$result_dir"/methyl/mbd"
bs_result_dir=$result_dir"/methyl/bs"
methylkit_from_bs_result_dir=$result_dir"/methyl/bs/methylkit"
medips_from_mbd_result_dir=$result_dir"/methyl/mbd/medips"
cel_result_dir=$result_dir"/gene_exp/cel"
rnaseq_result_dir=$result_dir"/gene_exp/rna_seq/deg/"
mu_dnaseq_result_dir=$result_dir"/mu/dna_seq/"
mu_rnaseq_result_dir=$result_dir"/mu/rna_seq/"
runDir_2pass=$rnaseq_result_dir"/2pass"


# directory based variables
ME_all_merged=$integ_result_dir/ME_MGD    # prepare initiall prefix for all merged file
ME_GE_all_merged=$integ_result_dir/ME_GE_MGD
ME_GE_COR_all_merged=$integ_result_dir/ME_GE_COR_MGD
ME_GE_COR_all_merged_THRESH=$integ_result_dir/ME_GE_COR_MGD_THRESH.txt


# functions
function my_join { local IFS="$1"; shift; echo "$*"; }



# functions
function binning_value {
		local local_interval=$1
		local local_max=$2
		local local_value_column=$3
		local level_file=$4

		awk -v interval=$local_interval -v max_range=$local_max -v value_column=$local_value_column '{ ind=int($value_column/interval)*interval;
								if($value_column>max_range){ind=max_range;}; value[ind] = value[ind] + 1;}
								END{for (v in value) {print v,value[v]}
							}' $level_file  | sort -k1,1g 
}



####################################################################################
# 0. set parameters
####################################################################################


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
#type_len=(13 7 10)

# multiple correciton
multiple_correction="fdr" #'bonferroni', 'fdr', 'holm' and so on. see R stats.p.adjust

# pearson or spearman
pear_sp_cor=0 # 0 : pearson / 1: spearman



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

# TODO : set data kinds (RNAseq or array ..)
# UCSC links for each samples : BigWig
<<'COMMENT_TEMP'
for (( i=0; i<${#ge_list[@]}; i++ )); do 
	# copy to web accessible directory
	# TODO :  this may be not necessary cuz already we are in the web accessible dir
	cp $ge_result_dir/`basename \${ge_list[$i]}`.bw $WEB_ACCESSIBLE_DIR
	# print UCSC url
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=test_SNP%20description=test2SNP%20visibility=dense%20bigDataUrl=$WEB_ACCESSIBLE_LOC/${ge_list[$i]}.bw"
done
COMMENT_TEMP


############
# 6.1.2 ME
############
# UCSC links for each samples : BigWig
<<'COMMENT_TEMP'
echo "[INFO] Methyl UCSC Bigwig"

if [ "$me_data_type" -eq "1" ]; then # MBD
	count=0
	for (( i=0; i<${#me_list[@]}; i++ )); do
		{
			temp_filename_only=`basename \${me_list[$i]}`
			met_level_file=$me_result_dir/$temp_filename_only".met_level";

			# create bigWig file for UCSC
			awk 'NR==1{next;}{print}' $met_level_file | cut -f1,2,3,5 | bedtools intersect -a - -b $REF_HUMAN_GENOME_100_BIN_BED -sorted > $viz_result_dir/$temp_filename_only".met_level.bed" ;
			$bin_dir/bedGraphToBigWig $viz_result_dir/$temp_filename_only".met_level.bed" $REF_HUMAN_CHR_SIZE $viz_result_dir/$temp_filename_only".met_level.bw";

			cp $viz_result_dir/$temp_filename_only".met_level.bw" $WEB_ACCESSIBLE_DIR;
			# print UCSC url
			echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl=$WEB_ACCESSIBLE_LOC/"$temp_filename_only".met_level.bw";
		} &
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait

elif [ "$me_data_type" -eq "2" ]; then # BS
	ls $methylkit_from_bs_result_dir/*.sorted.bw | xargs -I {} cp {} $WEB_ACCESSIBLE_DIR
	ls $methylkit_from_bs_result_dir/*.sorted.bw | xargs -I {} echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl="$WEB_ACCESSIBLE_LOC/{};
fi
COMMENT_TEMP

# NGS plot
<<'COMMENT_TEMP'
echo "[INFO] Methyl NGS plots"

if [ "$me_data_type" -eq "1" ] ; then # MBD
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		echo -n "" > $viz_result_dir/ngsplot_config_ME.txt
		echo -n "" > $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt


		for (( i=0; i<${#me_list[@]}; i++ )); do 
			temp_filename_only=`basename \${me_list[$i]}`

			# create configure file for all multiple file plotting
			echo -e "$me_result_dir/$temp_filename_only.sorted.bam\t-1\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_all.txt
			#echo -e "$mbd_result_dir/$temp_filename_only.sorted.bam\tintron.bed\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_all.txt

			# create configure file for subtype multiple file plotting
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				echo -e "$me_result_dir/$temp_filename_only.sorted.bam\t-1\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt
				#echo -e "$mbd_result_dir/$temp_filename_only.sorted.bam\tintron.bed\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt
			fi
		done

		# NGS plot itself is parallelized, so do not need to be parallelzed
		for region in "cgi" "genebody" "exon"; do #"enhancer" "dhs" "cgi" "exon"; do 
			ngs.plot.r -G hg19 -R $region -C $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt -O $viz_result_dir/ME_${type_kind[$k]}_$region 
		done
	done

	# for all sample
	for region in "bed" ;do #"genebody" "enhancer" "dhs" "cgi"; do
		/packages/test2/ngsplot-develop/bin/ngs.plot.r -G hg19 -R $region -C $viz_result_dir/ngsplot_config_ME_all.txt -O $viz_result_dir/ME_all_$region -E /data/project/mcpg/lib/region_info/intron.bed
	done

elif [ "$me_data_type" -eq "3" ] ; then # BS
	# TODO : no figure for BS?
	echo "no figure"
fi
COMMENT_TEMP
# TODO : ME density plot and comparison

<<'COMMENT_TEMP'
echo "[INFO] Drawing ME density plots"
	# make subtype avg met level
	for (( k=0; k<${#type_kind[@]}; k++)); do
	{	
		work_file=$ME_all_merged"_"${type_kind[$k]}
		filename_only=`basename $work_file`
		avg_met_file=$viz_result_dir/$filename_only"_avg"

		awk -v sample_num=${type_len[$k]} 'BEGIN{FS=OFS="\t";print "#chr","start","end","avg_value";}NR > 1{sum=0.0;for(i=4;i<=sample_num+3;i++){sum=sum+$(i)};print $1,$2,$3,sum/sample_num}' $work_file > $avg_met_file
	} &
	done
	wait
COMMENT_TEMP

<<'COMMENT_TEMP'
	count=0
	# draw ME density plot
	for (( k=0; k<${#type_kind[@]}; k++)); do
		for (( l=$k+1 ; l<${#type_kind[@]}; l++)); do
		{	
			work_file1=$ME_all_merged"_"${type_kind[$k]}
			filename_only1=`basename $work_file1`
			avg_met_file1=$viz_result_dir/$filename_only1"_avg"

			work_file2=$ME_all_merged"_"${type_kind[$l]}
			filename_only2=`basename $work_file2`
			avg_met_file2=$viz_result_dir/$filename_only2"_avg"

			output_prefix=$viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}"_ME_density_plot"			

			pair_DMR=$me_result_dir"/"${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.bed"

#			bedtools intersect -wa -a $avg_met_file1 -b $pair_DMR | cut -f4 > $output_prefix"_DMR_input1.txt"
#			bedtools intersect -wa -a $avg_met_file2 -b $pair_DMR | cut -f4 > $output_prefix"_DMR_input2.txt"

			max_val=100
			step_val=20
			
			if [ "$me_data_type" -eq "0" ]; then
				max_val=1.5
				step_val=0.3
			elif [ "$me_data_type" -eq "1" ]; then
				max_val=1
				step_val=0.2
			elif [ "$me_data_type" -eq "2" ]; then
				max_val=1
				step_val=0.2
			elif [ "$me_data_type" -eq "3" ]; then
				max_val=100
				step_val=20
			fi
			
			echo "[TEST] : "$output_prefix	
			$NEW_R_DIR/Rscript $bin_dir/me_density_v2.R $output_prefix"_DMR_input1.txt" $output_prefix"_DMR_input2.txt" ${type_kind[$k]} ${type_kind[$l]} $max_val $step_val $output_prefix".png"	
		} &
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
		done; wait
	done
	wait
COMMENT_TEMP

############
# 6.1.3 MU
############

<<'COMMENT_TEMP'
# UCSC links
for (( i=0; i<${#mu_list[@]}; i++ )); do
	#already moved in run_mutation code
#	cp $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.gz" $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.gz.tbi" $WEB_ACCESSIBLE_DIR
	# print UCSC url
	temp_filename_only=`basename \${mu_list[\$i]}`
	temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`
		
	vcf_filename=""

	if [ "$mu_data_type" -eq "0" ]; then
		vcf_filename=$temp_filename_only".vcf"
	elif [ "$mu_data_type" -eq "1" ]; then
		vcf_filename=$temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf"
	else
		vcf_filename=""
	fi
	
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=vcfTabix%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl=$WEB_ACCESSIBLE_LOC/$vcf_filename.gz"
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
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do 
			work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}".DEG.list"
		
			bash $bin_dir/create_annotation_link_html_based_on_gene_list.sh $ge_result_dir/$work_file > $viz_result_dir/$work_file".gene_only.annot.html"

			# create annotation page for deg
			cp $ge_result_dir/$work_file $viz_result_dir/$work_file".gene_only.annot.html" $WEB_ACCESSIBLE_DIR
		done
	done
COMMENT_TEMP

<<COMMENT
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do 
			work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}".deg"
			
		
			# get deg gene list only
			cut -f2 $rna_result_dir/$work_file > $viz_result_dir/$work_file".gene_only"
			sh $bin_dir/create_annotation_link_html_based_on_gene_list.sh $viz_result_dir/$work_file".gene_only" > $viz_result_dir/$work_file".gene_only.annot.html"

			# create annotation page for deg
			cp $rna_result_dir/$work_file $viz_result_dir/$work_file".gene_only" $viz_result_dir/$work_file".gene_only.annot.html" $WEB_ACCESSIBLE_DIR
		done
	done
COMMENT

############
# 6.2.2 ME 
############
#########################################
# circos
#########################################
circos_dir=$viz_result_dir/circos
mkdir -p $circos_dir
<<'COMMENT_TEMP'
echo "generate circos for subtype-specific avg methylation with difference and DEGs"

# create subtype avg methyl file from merged file
for (( i=0; i<${#type_kind[@]}; i++ )); do 
	subtype_ME_merged_file=$ME_all_merged"_"${type_kind[$i]}
	subtype_ME_merged_avg_file=$viz_result_dir/`basename \$subtype_ME_merged_file`"_avg"
	subtype_ME_merged_avg_10Mb_file=$subtype_ME_merged_avg_file"_10Mb"


	# file example
	#chr	bin_start	bin_end	100730_s_1.fq-rms	100730_s_7.fq-rms...	range_start	range_end	strand	gene_symbol	range_kind	refseq
	#chr1	9801	9900	0	0	0	0	0	0	0	0	0	0	0	0	0	9873	11873	+	DDX11L1	promoter	NR_046018
	#chr1	9901	10000	0	0	0	0	0	0	0	0	0	0	0	0	0	9873	11873	+	DDX11L1	promoter	NR_046018

	# TODO : check ME input type. file format should be same for all input type as above
	# compute avg methyl value

#	temp_sample_num_in_class=${#sample_num_in_class[$i]}
#	awk -v sample_num_in_class=$temp_sample_num_in_class 'NR==1{print "#chr","start","end","avg_value";next;}{avg_value=0; for (i=1+3; i<=sample_num_in_class+3; i++){avg_value+=$i;}; avg_value=avg_value/sample_num_in_class; print $1,$2,$3,avg_value;}' OFS='\t' $subtype_ME_merged_file > $subtype_ME_merged_avg_file

	# merge value the range in 10Mb for circos input
	bedtools map -header -c 4 -o mean -a $REF_HUMAN_GENOME_10M_BIN_BED -b $subtype_ME_merged_avg_file | sed -e 's/chr/hs/g' | sed -e 's/\t\./\t0/g' > $subtype_ME_merged_avg_10Mb_file
done



# crate circos for pairwise & all class

all_type_file_list=()
for (( i=0; i<${#type_kind[@]}; i++ )); do 

		subtype_ME_merged_file=$ME_all_merged"_"${type_kind[$i]}
		subtype_ME_merged_avg_file=$viz_result_dir/`basename \$subtype_ME_merged_file`"_avg"
		subtype_ME_merged_avg_10Mb_file=$subtype_ME_merged_avg_file"_10Mb"
	for (( j=$i+1; j<${#type_kind[@]}; j++ )); do 

		# met level file
		# generate circos input list deliminated by ';'	
		subtype_ME_merged_file=$ME_all_merged"_"${type_kind[$i]}
		subtype_ME_merged_file_2=$ME_all_merged"_"${type_kind[$j]}
	#	subtype_ME_merged_file="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/temp_result_files/ME_MGD_"${type_kind[$i]}".head100" # temp test with samll datda
	#	subtype_ME_merged_file_2="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/temp_result_files/ME_MGD_"${type_kind[$j]}".head100" # temp test with samll datda
		subtype_ME_merged_avg_file=$viz_result_dir/`basename \$subtype_ME_merged_file`"_avg"
		subtype_ME_merged_avg_file_2=$viz_result_dir/`basename \$subtype_ME_merged_file_2`"_avg"
		subtype_ME_merged_avg_10Mb_file=$subtype_ME_merged_avg_file"_10Mb"
		subtype_ME_merged_avg_10Mb_file_2=$subtype_ME_merged_avg_file_2"_10Mb"


		# for deg
	  work_file=${type_kind[$i]}"_vs_"${type_kind[$j]}
    DEG_file_name=$ge_result_dir"/"$work_file".DEG.list"
 #   DEG_file_name=/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/rna_paired_test_data/gene_exp/rna_seq/deg/Li_vs_Th.FC2.DEGs
		DEG_for_circos=$circos_dir/$work_file".DEGs.circos"
	#	DEG_for_circos=$circos_dir/$work_file".FC"$fold_change".DEGs.circos"
		diff_file=$circos_dir"/circos_${type_kind[$i]}_${type_kind[$j]}.diff"
		file_list=$subtype_ME_merged_avg_10Mb_file";"$subtype_ME_merged_avg_10Mb_file_2


		# TODO : check the cut off for degs file name

		# convert deg file to circos file (add position information and change chr to hs)

		# genebody info
		#	chr1	11873	14409	DDX11L1	NR_046018	+
		# chr1	14361	29370	WASH7P	NR_024540	-

		join -1 1 -2 4 <(sort $DEG_file_name) $REF_HUMAN_GENE_GENEBODY_INFO_SORTED | awk '{print $2,$3,$4,$1}' OFS='\t' | sed -e 's/chr/hs/g' > $DEG_for_circos

		# SLCO1B1 chr12 21284127 21392730 NM_006446 +
		# SLCO1B3 chr12 20963637 21069843 NM_019844 +kI

		# create diff file
		# chr1  0 10000000  5.532632995
		paste $subtype_ME_merged_avg_10Mb_file $subtype_ME_merged_avg_10Mb_file_2 | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1,$2,$3,abs($4-$8)}' OFS='\t' > $diff_file
		
		# create subtype pair-wise circos
		# TODO : set min max value for circos
		circos_min=0
		circos_max=15

		circos_diff_min="0"
		circos_diff_max="0.05"
		
		if [ "$ME_input_type" -eq "0" ]; then
				circos_diff_min="0"				
				circos_diff_max="2"
		elif [ "$ME_input_type" -eq "3" ]; then
				circos_diff_min="0"
				circos_diff_max="50"
		else
				circos_diff_min="0"
				circos_diff_max="0.05"
		fi

		temp_circos_file=$circos_dir"/circos_${type_kind[$i]}_${type_kind[$j]}.conf"
		bash $bin_dir/generate_circos_conf.sh 2 $file_list $diff_file $DEG_for_circos $CIRCOS_CONFIG_LOC $circos_min $circos_max $circos_diff_min $circos_diff_max > $temp_circos_file
		circos -conf $temp_circos_file -outputdir $circos_dir -outputfile $work_file
	done
	all_type_file_list+=($subtype_ME_merged_avg_10Mb_file)
done
#all_type_file_list+=($subtype_ME_merged_avg_10Mb_file_2)


# TODO : snagseon
# measure entropy for all diff file, set proper min mac for diff as well... maybe it should be set seperately in generate_circos_conf.

	diff_all_file=$circos_dir"/circos_all.diff"
	
	circos_diff_max=5
	paste ${all_type_file_list[@]} | awk -v class_num=${#type_kind[@]} -v max_val=$((circos_diff_max + 1)) 'BEGIN{FS=OFS="\t"}{sum=0.0;for(i=1;i<=class_num;i++){met_arr[i]= 0.1 + $(4*i); sum=sum+met_arr[i]}; ent = 0.0; for(i=1;i<=class_num;i++){ent=ent-((met_arr[i]/sum)*(log(met_arr[i]/sum)/log(2)))}; print $1,$2,$3, max_val*(1-ent/(log(class_num)/log(2)))}' > $diff_all_file

# create all subtype circos

# TODO : comment for test
GDEG=$integ_result_dir/GDEG_kw
GDEG_for_circos=$circos_dir"/GDEG_for_circos"
join -1 1 -2 4 <(sort $GDEG"_GSym") $REF_HUMAN_GENE_GENEBODY_INFO_SORTED | awk '{print $2,$3,$4,$1}' OFS='\t' | sed -e 's/chr/hs/g' > $GDEG_for_circos

circos_min=0
circos_max=15
circos_diff_min=0
temp_circos_file=$circos_dir"/circos_all.conf"
bash $bin_dir/generate_circos_conf.sh ${#type_kind[@]} `my_join ";" ${all_type_file_list[@]}` $diff_all_file $GDEG_for_circos $CIRCOS_CONFIG_LOC $circos_min $circos_max $circos_diff_min $circos_diff_max > $temp_circos_file
circos -conf $temp_circos_file -outputdir $circos_dir -outputfile "all"

COMMENT_TEMP

<<'COMMENT_TEMP'
########################################
# DMR : Fold change -> bedgraph -> bw
########################################
if [ "$me_data_type" -eq "0" ]; then
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do 
			# create subtype fold change DMR region based on 'subtype1_vs_subtype2.mr.edgeR.s.txt'
			awk 'NR==1{next;}{print $1"\t"$2"\t"$3"\t"$(NF-3)}' $medips_from_mbd_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".mr.edgeR.s.txt"  > $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff"
			$bin_dir/bedGraphToBigWig $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff" $REF_HUMAN_CHR_SIZE $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw"

			# copy to web accessible directory
			cp $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw" $WEB_ACCESSIBLE_DIR

			# print UCSC url
			echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=${type_kind[$k]}_${type_kind[$l]}_DMR%20description=${type_kind[$k]}_${type_kind[$l]}_DMR%20visibility=dense%20bigDataUrl="$WEB_ACCESSIBLE_LOC/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw"
		done
	done
elif [ "$me_data_type" -eq "3" ]; then

	for (( k=0; k<${#type_kind[@]}; k++ )); do
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do
			awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$7}' $methylkit_from_bs_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR_list.txt" > $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff"

			$bin_dir/bedGraphToBigWig $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff" $REF_HUMAN_CHR_SIZE $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw"

			# copy to web accessible directory
			cp $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw" $WEB_ACCESSIBLE_DIR

			# print UCSC url
			echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=${type_kind[$k]}_${type_kind[$l]}_DMR%20description=${type_kind[$k]}_${type_kind[$l]}_DMR%20visibility=dense%20bigDataUrl="$WEB_ACCESSIBLE_LOC/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw"


		done
	done

else

	#TODO infinium
	echo "infinium"

fi
COMMENT_TEMP
########################################
# GDMR : bed -> bb
########################################
<<'COMMENT_TEMP'
	$bin_dir/bedToBigBed $integ_result_dir/GDMR_kw $REF_HUMAN_CHR_SIZE $viz_result_dir/GDMR_by_kw.bb
	
	# copy to web accessible directory
	cp $viz_result_dir/GDMR_by_kw.bb $WEB_ACCESSIBLE_DIR

	# print UCSC url
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigBed%20name=GDMR_kw%20description=GDMR_kw%20visibility=dense%20bigDataUrl="$WEB_ACCESSIBLE_LOC/GDMR_by_kw.bb 

########################################
# DMR all : heatmap
########################################
	# get top 10% by adj pavlue (NF)
	pvalue_index=`awk '{print NF;exit}' \$integ_result_dir/ME_MAT.kw.txt`
	line_num=`grep -v nan \$integ_result_dir/ME_MAT.kw.txt | sort -k"$pvalue_index","$pvalue_index"g | tee $viz_result_dir/ME_MAT.kw.txt.sorted | wc -l`
	top_10_line_num=$(( $line_num / 10 ))

	echo "total line : " $line_num
	echo "10% line : " $top_10_line_num

	head -n $top_10_line_num $viz_result_dir/ME_MAT.kw.txt.sorted | cut -f1-$((${#type_list[@]}+1))> $viz_result_dir/ME_MAT.kw.txt.sorted.top10
	#grep -v NaN /data/project/mcpg/result/profile/ME_matrix.kw.sigle_diff | awk '{print NF-1;exit}' | xargs -I {} sort -k{},{}n /data/project/mcpg/result/profile/ME_matrix.kw.sigle_diff
	# get

	#base name list
	basename_list=()
	for((i=0;i<${#me_list[@]};i++)); do
		filename_only=`basename ${me_list[$i]}`
		basename_list+=($filename_only)
	done

	echo `my_join "," ${basename_list[@]}`
	echo `my_join "," ${type_kind[@]}`	
	# create heatmap
	$R_DIR/Rscript $bin_dir/make_heatmap3_ver2.r "$viz_result_dir/ME_MAT.kw.txt.sorted.top10" `my_join "," ${basename_list[@]}` "Class" `my_join "," ${type_kind[@]}` "$viz_result_dir/DMR_TOP10_heatmap.png" "Top_10_DMRs_by_Kruscal_Wallis"

	# copy to WEB_ACCESSIBLE_DIR
	cp $viz_result_dir/DMR_TOP10_heatmap.png $WEB_ACCESSIBLE_DIR

COMMENT_TEMP

########################################
# DMR region stat : bar chart
########################################
<<'COMMENT_TEMP'
	echo "[INFO] DMR region stat on bar chart"

	echo -e "Region\tcounts\tsubtype_pair" > $viz_result_dir/DMR_region_stat_matrix.txt
	echo -e "Region\tcounts\tsubtype_pair" > $viz_result_dir/DMR_cpgi_stat_matrix.txt



	echo "[INFO] count number of DMR in each genomic region"
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do 
			subtype_pair="${type_kind[$k]}_vs_${type_kind[$l]}"
			#for region in "3 prime UTR" "5 prime UTR" "Promoter" "Exon" "Intron"; do
			for region in "3utr" "5utr" "promoter" "exon" "intron"; do
				temp_region_dmr_count=`bedtools intersect -wa -wb -a $viz_result_dir/$subtype_pair.DMR.diff -b $REF_HUMAN_GENE_RANGE_INFO | grep -w $region | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`

				echo -e $region"\t"$temp_region_dmr_count"\t"$subtype_pair >> $viz_result_dir/DMR_region_stat_matrix.txt
			done
				# create cpg_island, cpg_shore, cpg_shelf
			for region in "cpgIsland" "cpgShore" "cpgShelf"; do
				temp_region_dmr_count=`bedtools intersect -wa -wb -a $viz_result_dir/$subtype_pair.DMR.diff -b $REF_HUMAN_GENE_RANGE_INFO | grep -w $region | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`
				echo -e $region"\t"$temp_region_dmr_count"\t"$subtype_pair >> $viz_result_dir/DMR_cpgi_stat_matrix.txt
			done
		done
	done

	# draw bar chart
	echo "[INFO] Draw bar chart"
	$NEW_R_DIR/Rscript $bin_dir/barplot_ver2.r $viz_result_dir/DMR_region_stat_matrix.txt $viz_result_dir/DMR_region_stat_barplot.png "region"
	$NEW_R_DIR/Rscript $bin_dir/barplot_ver2.r $viz_result_dir/DMR_cpgi_stat_matrix.txt $viz_result_dir/DMR_cpgi_stat_barplot.png "cgi"

	cp $viz_result_dir/DMR_region_stat_barplot.png $viz_result_dir/DMR_cpgi_stat_barplot.png $WEB_ACCESSIBLE_DIR
COMMENT_TEMP

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
<<'COMMENT_TEMP'
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
	work_file_for_GDSNP=$viz_result_dir/GDSNP_$range"_MAT"
	echo "GDSNP_in_$range" > $work_file_for_GDSNP

	grep -w $range $integ_result_dir/GDSNP.txt | cut -f1,2,3 - | bedtools coverage -a $viz_result_dir/temp_range_div_100 -b - -counts | cut -f4,6  | sort -k1,1n | bedtools groupby -g 1 -c 2 -o sum | cut -f2 >> $work_file_for_GDSNP 

	paste $viz_result_dir/prefix_percentage_ids $work_file_for_GDSNP > $viz_result_dir/GDSNP_$range"_matrix4plot_merged"
	# line plot
	$NEW_R_DIR/Rscript $bin_dir/plot.r $viz_result_dir/GDSNP_$range"_matrix4plot_merged" $viz_result_dir/GDSNP_$range"_matrix4plot_merged.png" 
COMMENT_TEMP
	snp_file_extension=$(echo `basename \${mu_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
<<'COMMENT_TEMP'
	# for DSNP
	paste_list_all=()
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		paste_list_sub_type=()
		count=0
		for (( i=0; i<${#mu_list[@]}; i++ )); do 
			temp_filename_only=`basename \${mu_list[\$i]}`
			temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`
			vcf_bed_file=""
			if [ "$mu_data_type" -eq "0" ]; then # DNR-seq
				vcf_bed_file=$mu_result_dir/$temp_filename_only".vcf.bed"
			elif [ "$mu_data_type" -eq "1" ]; then # RNA-seq
				vcf_bed_file=$mu_result_dir/$temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf.bed"
			else
				vcf_bed_file=""
			fi
			
			# create range stat file per sample
			work_file=$viz_result_dir/SNP_$range"_MAT_"`basename \${mu_list[$i]}`

			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				
				echo "[INFO] Estimate relative coverage in range of all $range region in `basename \${mu_list[$i]}`"

				paste_list_sub_type+=($work_file)
				paste_list_all+=($work_file)
				echo "`basename \${mu_list[$i]}`" > $work_file

				# parallization
				cut -f1,2,3 $vcf_bed_file | bedtools coverage -b $viz_result_dir/temp_range_div_100 -a - -counts | cut -f4,6  | sort -k1,1n | bedtools groupby -g 1 -c 2 -o sum | cut -f2 >> $work_file &

				# NOTE : MEMORY INTENSE, so ONLY USE 5 CPUS			
				let count+=1; [[ $((count%5)) -eq 0 ]] && wait
			fi
		done; wait
		
		paste $viz_result_dir/prefix_percentage_ids ${paste_list_sub_type[@]} > $viz_result_dir/SNP_$range"_MAT_"${type_kind[$k]}

		$NEW_R_DIR/Rscript $bin_dir/plot.r $viz_result_dir/SNP_$range"_MAT_"${type_kind[$k]} $viz_result_dir/SNP_$range"_"${type_kind[$k]}".png"
	done

	# merge all stat files
	paste $viz_result_dir/prefix_percentage_ids ${paste_list_all[@]} > $viz_result_dir/SNP_$range"_MAT_all"

	$NEW_R_DIR/Rscript $bin_dir/plot.r $viz_result_dir/SNP_$range"_matrix4plot_all" $viz_result_dir/SNP_$range"_all.png"

	cp $viz_result_dir/SNP_*.png $WEB_ACCESSIBLE_DIR
done

# GDSNP all :  bed w/ label -> bb
	sort -k1,1 -k2,2n $integ_result_dir/GDSNP.bed | uniq > $viz_result_dir/GDSNP_sorted.bed
	$bin_dir/bedToBigBed $viz_result_dir/GDSNP_sorted.bed $REF_HUMAN_CHR_SIZE $viz_result_dir/GDSNP_sorted.bb

	# copy BigBed to web accessible dir
	cp $viz_result_dir/GDSNP_sorted.bb $WEB_ACCESSIBLE_DIR

	# create UCSC link
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigBed%20name=GDSNP%20description=GDSNP%20visibility=pack%20bigDataUrl=$WEB_ACCESSIBLE_LOC/GDSNP_sorted.bb"


# GDSNP per sample : boxplot
	GDSNP_cnt=$viz_result_dir/"GDSNP_region_count_per_sample"
	
	paste $integ_result_dir/GDSNP.txt <(awk '{print $NF}' $integ_result_dir/GDSNP.bed) | cut -f1-$((${#type_list[@]}+3)),$((${#type_list[@]}*2+8)),$((${#type_list[@]}*3+11)) | grep -v "SNP_TSS_TSE_flanking_range+-250kb" | sort -k1,1 -k2,2n | uniq > $GDSNP_cnt
	
	sample_pos_list=()
	temp_sum=3
	for ((i=0;i<${#type_len[@]};i++)); do
		let temp_sum+=${type_len[$i]}
		sample_pos_list+=($temp_sum)
	done
	
	sample_pos_string=`my_join "," ${sample_pos_list[@]}`


	for region in "3utr" "5utr" "promoter" "exon" "intron" "cpgIsland" "cpgShelf" "cpgShore"; do
		grep -w $region $GDSNP_cnt | sort -k1,1 -k2,2n > $GDSNP_cnt"_"$region

		awk -v class_type=`my_join "," ${type_kind[@]}` -v sample_num_class=$sample_pos_string -v sample_num=${#type_list[@]} 'BEGIN{type_num=split(class_type, class_arr, ","); num=split(sample_num_class, sample_num_arr, ",");sample_num_arr[0]=0;for(i=4; i<=(3+sample_num); i++){a[i]=0;}} {for(i=4; i<=(3+sample_num); i++){if($i!="-/-") { for(j=1;j<=type_num;j++){ if($NF==class_arr[j] && (i > sample_num_arr[j-1] && i <= sample_num_arr[j])){a[i]=a[i]+1;break}}}}}END{for(i=4; i<=(3+sample_num); i++){print a[i]};}' $GDSNP_cnt"_"$region > $GDSNP_cnt"_"$region"_per_sample_cnt_w_zero_for_others"
	
	done

	# draw boxplot

	$NEW_R_DIR/Rscript $bin_dir/boxplot.r `ls $viz_result_dir/*sample_cnt_w_zero_for_others | tr "\n" ";"` `my_join ";" ${type_len[@]}` `my_join ";" ${type_kind[@]}` $viz_result_dir/GDSNP_boxplot_genomic.png $viz_result_dir/GDSNP_boxplot_cpgi.png

	cp $vis_result_dir/GDSNP_boxplot_genomic.png $viz_result_dir/GDSNP_boxplot_cpgi.png $WEB_ACCESSIBLE_DIR

COMMENT_TEMP

<<'COMMENT'
# GDSNP region stat : barplot
	echo "[INFO] GDSNP region stat on bar chart"
	echo -e "Region\tcounts\tSubtype" > $viz_result_dir/GDSNP_RNG_STAT_MAT.txt
	echo -e "Region\tcounts\tSubtype" > $viz_result_dir/GDSNP_CPGs_STAT_MAT.txt
#	echo -e "Region\tcounts\tSubtype" > $viz_result_dir/GDSNP_ALL_STAT_MAT.txt

	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for region in "3utr" "5utr" "promoter" "exon" "intron" "genebody"; do
			temp_region_gdsnp_count=`bedtools intersect -wa -wb -a $viz_result_dir/GDSNP_sorted.bed -b $REF_HUMAN_GENE_RANGE_INFO | grep -w $region | grep -w ${type_kind[$k]} | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`

			echo -e $region"\t"$temp_region_gdsnp_count"\t"${type_kind[$k]} >> $viz_result_dir/GDSNP_RNG_STAT_MAT.txt
#			echo -e $region"\t"$temp_region_gdsnp_count"\t"${type_kind[$k]} >> $viz_result_dir/GDSNP_ALL_STAT_MAT.txt
		done

		for region in "CpG Island" "CpGI Shelf" "CpGI Shore"; do 
			temp_region_gdsnp_count=`bedtools intersect -wa -wb -a $viz_result_dir/GDSNP_sorted.bed -b $REF_HUMAN_GENE_RANGE_INFO | grep -w $region | grep -w ${type_kind[$k]} | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`
			echo -e $region"\t"$temp_region_gdsnp_count"\t"${type_kind[$k]} >> $viz_result_dir/GDSNP_CPGs_STAT_MAT.txt
#			echo -e $region"\t"$temp_region_gdsnp_count"\t"${type_kind[$k]} >> $viz_result_dir/GDSNP_ALL_STAT_MAT.txt
		done
	done

	# draw pie chart
	echo "[INFO] Draw bar chart"
	/home/airavata/packages/R-3.2.3/bin/Rscript $bin_dir/boxplot.r $viz_result_dir/GDSNP_RNG_STAT_MAT.txt $viz_result_dir/BAR_GDSNP_RNG.png "Number of tumor subtype specific mutation"
	/home/airavata/packages/R-3.2.3/bin/Rscript $bin_dir/boxplot.r $viz_result_dir/GDSNP_CPGs_STAT_MAT.txt $viz_result_dir/BAR_GDSNP_CPGs.png "Number of tumor subtype specific mutation"
#	Rscript $bin_dir/barplot.r $viz_result_dir/GDSNP_ALL_STAT_MAT.txt $viz_result_dir/BAR_GDSNP_ALL.png "Number of tumor subtype specific mutation"

	cp $viz_result_dir/BAR_GDSNP_RNG.png $viz_result_dir/GDSNP_RNG_STAT_MAT.txt $viz_result_dir/GDSNP_CPGs_STAT_MAT.txt $viz_result_dir/BAR_GDSNP_CPGs.png $WEB_ACCESSIBLE_DIR
COMMENT

## case2
###################################################################################
### Extract single cpg site specific normalized methylation level
###################################################################################
					#cat $hg19_reference_sequence | $bin_dir/extract_cpg_site_from_fasta > $region_info_dir/hg19_cpg_site.bed

NUM_CPUS=3 # NOTE: below takes 10GB per CPU.
# create cpgsite me level & snp exist cpgsite me level 
for deg_pair in "all_gene"; do # "p_cut_genes" "lu_vs_baa.deg" "lu_vs_bab.deg" "baa_vs_bab.deg" "GDEG"; do 
	for region in  "whole_genome" ; do #cpgShore" "cpgi"; do # "whole_genome" "Intron" "Exon" "Promoter" "Genebody" "TFBS" "Promoter_Genebody"; do #
		# set region information
		# NOTE : FILE NAME DEPENDENCY!!!
		# promotor : lu_vs_baa.deg.list_Promoter_Refseq.bed
		# tbfs in promotor : lu_vs_baa.deg.list_Promoter_Refseq.bed_sequence.fa_match_minSUM_good.txt
		# promoter and genebody : lu_vs_baa.deg.list_Promoter_Genebody_Refseq.bed
<<'COMMENT'
		if [ $deg_pair == "all_gene" ]; then # not consider GE, just SNP & ME
			if [ $region == "whole_genome" ]; then
				region_file=""
			elif [ $region == "TFBS" ]; then
				continue
			elif [ $region == "cpgShore" ] || [ $region == "cpgi" ]; then
				region_file=$region_info_dir"$region"".bed"
			else
				region_file=$region_info_dir"$region""_Refseq.bed"
			fi
		elif [ $deg_pair == "p_cut_genes" ]; then 
			if [ $region == "whole_genome" ]; then
				continue
			elif [ $region == "cpgShore" ] || [ $region == "cpgi" ] || [ $region == "TFBS" ]; then
				continue
			else
				region_file=$region_info_dir"/p_cut_gene_related/p_cut_genes_"$region"_Refseq.bed"
			fi
		else # consider SNP, ME, GE together, so consider degs
			if [ $region == "cpgShore" ] || [ $region == "cpgi" ]; then
				# These regions are not related with degs, only consider when deg_pair is all_gene and skip for others
				continue
			elif [ $region == "TFBS" ]; then
				echo "[INFO] This is TFBS, so need to merge the region information"
				info_file=$deg_with_regions_dir"$deg_pair"".list_Promoter_Refseq.bed_sequence.fa_match_minSUM_good.txt.bed.sorted"
				TFBS_info_file_merged=$info_file".merged"

				# get TFBS merged region since TFBS overlaps themselves
				# NOTE : !!! THIS IS ONLY NEEDE ONE TIME, SO COMMENTED !!!
#				cut -f1,2,3,4 $info_file | sort -k1,1 -k2,2n | uniq | bedtools merge -c 4 -o distinct -i - > $TFBS_info_file_merged
				region_file=$TFBS_info_file_merged
			else
				region_file=$deg_with_regions_dir"$deg_pair"".list_"$region"_Refseq.bed"
			fi 
		fi
COMMENT
		region_file=""
		echo "[INFO] Region file : $region_file"
		for metric in "rms"; do # "rms"; do #rms

			# variables

			if [ "$me_data_type" -eq "0" ]; then	# MeDIP/MBD
				if [ $metric == "cnt" ]; then 
					value_column=4; interval=1; max=30
				elif [ $metric == "rms" ]; then
					value_column=5; interval=1; max=20
				fi
			elif [ "$me_data_type" -eq "3" ]; then # BS-seq
				value_column=5; interval=5; max=100
			fi


			category_pair=$deg_pair"_"$region"_"$metric"_"$interval"_"$max
		
			table_for_boxplot=$viz_result_dir/"table_for_boxplot_"$category_pair
		
			# create table for box plot
			echo -e "Subtype\tDNA_methylation_level\tMutation_rate" > $table_for_boxplot

<<'COMMENT'
			class_mutation_rate_table_list=()
			for ((k=0;k<${#type_kind[@]};k++)); do
				temp_calss_mutation_rate_table=$viz_result_dir"/mutation_rate_table_"${type_kind[$k]}"_"$category_pair
				echo -n "" > $temp_calss_mutation_rate_table
				class_mutation_rate_table_list+=($temp_calss_mutation_rate_table)
			done
COMMENT

<<'COMMENT'
			lu_mutation_rate_table=$vizs_result_dir/"mutation_rate_table_lu_"$category_pair
			baa_mutation_rate_table=$viz_result_dir/"mutation_rate_table_baa_"$category_pair
			bab_mutation_rate_table=$viz_result_dir/"mutation_rate_table_bab_"$category_pair
COMMENT
		<<'COMMENT'
			# create empty subtype mutation rate table
			echo -n "" > $lu_mutation_rate_table
			echo -n "" > $baa_mutation_rate_table
			echo -n "" > $bab_mutation_rate_table
COMMENT

			# create cpgsite me level & snp exist cpgsite me level 
			count=0
			echo "[INFO] Create cpgsite me level file & snp exist cpgsite me level file"
			for((i=0;i<${#me_list[@]};i++)); do 

				echo "[INFO] Currently processing $i sample"
				temp_filename_only=`basename \${me_list[$i]}`
    		met_level_file=$me_result_dir/$temp_filename_only".met_level"
				met_level_file_only=$temp_filename_only".met_level"


    		snp_temp_filename_only=`basename \${mu_list[\$i]}`
    		snp_temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`
    
				if [ "$mu_data_type" -eq "0" ]; then
      		snp_file=$mu_result_dir/$snp_temp_filename_only".vcf.bed"
    		elif [ "$mu_data_type" -eq "1" ]; then
      		snp_file=$mu_result_dir/$snp_temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf.bed"
    		else
      		snp_file=""
    		fi

				snp_filename_only=`basename $snp_file`


				cpgsite_me_level=$viz_result_dir/$met_level_file_only".cpg_site_me_level_"$category_pair
				snp_exist_cpgsite_me_level=$viz_result_dir/$met_level_file_only".snp_exist_cpg_site_me_level_"$category_pair
				subtype=${type_list[$i]}

<<'COMMENT'
				# skip samples depends on deg pair
				if [ $(check_sample_in_list $deg_pair $subtype) == "continue" ]; then 
					continue 
				fi
COMMENT
		
				# extract cpg site ME level, extract snp exist cpg site me level 
				if [ $region == "whole_genome" ]; then
					echo "here"
					tail -n+2 $met_level_file | bedtools intersect -wb -wa -a - -b $HUMAN_HG19_CPGSITE | awk -v value_column=$value_column '{OFS="\t"; print $6,$7,$8,$value_column}' | tee $cpgsite_me_level | bedtools intersect -wa -wb -a - -b $snp_file > $snp_exist_cpgsite_me_level &
				else
					tail -n+2 $met_level_file | bedtools intersect -wb -wa -a - -b $HUMAN_HG19_CPGSITE | awk -v value_column=$value_column '{OFS="\t"; print $6,$7,$8,$value_column}'| bedtools intersect -wa -wb -a - -b $region_file |tee $cpgsite_me_level | bedtools intersect -wa -wb -a - -b $snp_file > $snp_exist_cpgsite_me_level &
				fi
				
				let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
			done; wait

			echo "[INFO] Compute mutation rate and plot"
			count=0
			for((i=0;i<${#me_list[@]};i++)); do 

				echo "[INFO] Currently processing $i sample"
				# filter outlier samples # NOTE: HERE JUST FOR TESTING!!!!!!!!!!!!!!!!!!!
				#if [ $i -eq 5 ] || [ $i -eq 6 ] || [ $i -eq 7 ] || [ $i -eq 17 ] || [ $i -eq 21 ] || [ $i -eq 24 ]; then
				#	continue
				#fi


				temp_filename_only=`basename \${me_list[$i]}`
    		met_level_file=$me_result_dir/$temp_filename_only".met_level"
				met_level_file_only=$temp_filename_only".met_level"


    		snp_temp_filename_only=`basename \${mu_list[\$i]}`
    		snp_temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`
    
				if [ "$mu_data_type" -eq "0" ]; then
      		snp_file=$mu_result_dir/$snp_temp_filename_only".vcf.bed"
    		elif [ "$mu_data_type" -eq "1" ]; then
      		snp_file=$mu_result_dir/$snp_temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf.bed"
    		else
      		snp_file=""
    		fi

				snp_filename_only=`basename $snp_file`


				subtype=${type_list[$i]}

<<'COMMENT'
				# skip samples depends on deg pair
				if [ $(check_sample_in_list $deg_pair $subtype) == "continue" ]; then 
					continue 
			        fi
COMMENT

				# for debug
				echo "[INFO] Binning"
				echo "[INFO] ME file : " $met_level_file
				echo "[INFO] SNP file : " $snp_file
				echo "[INFO] Subtype : " $subtype
				echo "[INFO] Region : " $region
				echo "[INFO] Metric : " $metric
				echo "[INFO] DEG pair : " $deg_pair

				cpgsite_me_level=$viz_result_dir/$met_level_file_only".cpg_site_me_level_"$category_pair
				snp_exist_cpgsite_me_level=$viz_result_dir/$met_level_file_only".snp_exist_cpg_site_me_level_"$category_pair


				number_of_cpgsite_per_methylation_range=$viz_result_dir/$met_level_file_only".number_of_cpgsite_per_methylation_range_"$category_pair
				number_of_snp_exist_cpgsite_per_methylation_range=$viz_result_dir/$met_level_file_only".number_of_snp_exist_cpgsite_per_methylation_range_"$category_pair

				mutation_rate_per_methylation_range=$viz_result_dir/$met_level_file_only".mutation_rate_per_methylation_range_"$subtype"_"$category_pair
							#mutation_rate_table=$viz_result_dir/"mutation_rate_table_"$subtype"_"$category_pair

				# compute mutation rate of each methylation range
				## extract number of cpgsite per multation range
				binning_value $interval $max 4 $cpgsite_me_level > $number_of_cpgsite_per_methylation_range &
				binning_value $interval $max 4 $snp_exist_cpgsite_me_level  > $number_of_snp_exist_cpgsite_per_methylation_range &
				let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
			done;wait

			for((i=0;i<${#me_list[@]};i++)); do 

				echo "[INFO] Currently processing $i sample"
				# filter outlier samples # NOTE: HERE JUST FOR TESTING!!!!!!!!!!!!!!!!!!!
				#if [ $i -eq 5 ] || [ $i -eq 6 ] || [ $i -eq 7 ] || [ $i -eq 17 ] || [ $i -eq 21 ] || [ $i -eq 24 ]; then
				#	continue
				#fi

				temp_filename_only=`basename \${me_list[$i]}`
    		met_level_file=$me_result_dir/$temp_filename_only".met_level"
				met_level_file_only=$temp_filename_only".met_level"


    		snp_temp_filename_only=`basename \${mu_list[\$i]}`
    		snp_temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`
    
				if [ "$mu_data_type" -eq "0" ]; then
      		snp_file=$mu_result_dir/$snp_temp_filename_only".vcf.bed"
    		elif [ "$mu_data_type" -eq "1" ]; then
      		snp_file=$mu_result_dir/$snp_temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf.bed"
    		else
      		snp_file=""
    		fi

				snp_filename_only=`basename $snp_file`


				subtype=${type_list[$i]}


<<'COMMENT'
				# skip samples depends on deg pair
				if [ $(check_sample_in_list $deg_pair $subtype) == "continue" ]; then 
					continue 
			        fi
COMMENT
				# for debug
				echo "[INFO] Mutation rate computing"
				echo "[INFO] ME file : " $met_level_file
				echo "[INFO] SNP file : " $snp_file
				echo "[INFO] Subtype : " $subtype
				echo "[INFO] Region : " $region
				echo "[INFO] Metric : " $metric
				echo "[INFO] DEG pair : " $deg_pair


				number_of_cpgsite_per_methylation_range=$viz_result_dir/$met_level_file_only".number_of_cpgsite_per_methylation_range_"$category_pair
				number_of_snp_exist_cpgsite_per_methylation_range=$viz_result_dir/$met_level_file_only".number_of_snp_exist_cpgsite_per_methylation_range_"$category_pair

				mutation_rate_per_methylation_range=$viz_result_dir/$met_level_file_only".mutation_rate_per_methylation_range_"$subtype"_"$category_pair
							#mutation_rate_table=$viz_result_dir/"mutation_rate_table_"$subtype"_"$category_pair



				## compute mutation rate
				join -a1 -a2 <(sort -k1,1 $number_of_cpgsite_per_methylation_range) <(sort -k1,1 $number_of_snp_exist_cpgsite_per_methylation_range) | sort -k1,1n | awk '{OFS="\t"; if($3=="") $4=0; print $1,$3/$2 }' > $mutation_rate_per_methylation_range

				## plot mutation rate
							#echo "[INFO] out file" $mutation_rate_per_methylation_range".jpg"
							#Rscript $bin_dir/line_plot_with_two_inputs.r $mutation_rate_per_methylation_range $mutation_rate_per_methylation_range".jpg"

				# create table for ttest	
							#cut -f2 $mutation_rate_per_methylation_range | awk -f $bin_dir/transpose.awk >> $mutation_rate_table
							#for(i in 1:ncol(lu_data)){print(i);print(t.test(lu_data[,i],baa_data[,i]))}
		
				# create boxplot data	
				awk -v subtype=$subtype '{OFS="\t";print subtype, $1, $2}' $mutation_rate_per_methylation_range >> $table_for_boxplot
			done
			# draw boxplot
			$NEW_R_DIR/Rscript $bin_dir/mu_rate_boxplot.r $table_for_boxplot $table_for_boxplot".png"

		done
	done
done


<<'COMMENT'
################################################################
# mutation rate from bs 
################################################################

# all cpg site info file
# HCC70_TGACCA_L001_R1_001_val_1.fq_bismark_pecpg.raw.sorted.filtered..coverage.bedgraph

bs_region=/data/project/mcpg/lib/region_info/bs_seq_target_region_sorted.txt

# bs_me file example
# mDAmB_231_GGCTAC_L002_R1_001_val_1.fq_bismark_pe.sorted.sam.bed


NUM_CPUS=8
for region in "whole_bs_region"; do #"cpgShore" "cpgi" "Intron" "Promoter" "Genebody" "Exon" "Promoter_Genebody"; do
	echo "[INFO] Create cpgsite me level file & snp exist cpgsite me level file"


	count=0
	for((i=0;i<${#bs_r1_list[@]};i++)); do 

		# variables

		# bs_me
		input_fq1_filename_only=`basename \${bs_r1_list[$i]}`
		file_extention=$(echo $input_fq1_filename_only |awk -F . '{if (NF>1) {print $NF}}')
		filename1_wo_ext=`basename $input_fq1_filename_only "."$file_extention`

		bs_me_level_file=$bs_me_result_dir/$filename1_wo_ext"_val_1.fq_bismark_pe.sorted.sam.bed"
		bs_me_filename_only=`basename \$bs_me_level_file`

		# snp
		snp_from_bs_result_dir="/data/project/breast_TF_methyl/30celline_SNP_TFBS_ME_GE/result/snp_from_bs/"
		bam_prefix=$filename1_wo_ext"_val_1.fq_bismark_pe"

		bs_snp_file=$snp_from_bs_result_dir/$bam_prefix"snp.raw.sorted.filtered.vcf"
		bs_snp_filename_only=`basename \$bs_snp_file`

		bs_cpgsite_me_level=$viz_result_dir/$bs_me_filename_only".cpg_site_me_level_"$region
		bs_snp_exist_cpgsite_me_level=$viz_result_dir/$bs_me_filename_only".snp_exist_cpg_site_me_level_"$region
		bs_subtype=${bs_subtype_list[$i]}


		if [ $region == "whole_bs_region" ]; then
			region_file=""
		elif [ $region == "cpgShore" ] || [ $region == "cpgi" ]; then
			region_file=$region_info_dir"$region"".bed"
		else
			region_file=$region_info_dir"$region""_Refseq.bed"
		fi



		# extract cpg site ME level, extract snp exist cpg site me level 
		# NOTE : !!! THIS IS ONLY NEEDE ONE TIME, SO COMMENTED !!!

		# 1. intersect bs_region + hg19 cpgsite = cpgsite in bs region [all sample same]
		# 2. intersect bs snp file + cpgsite in bs region = snp exist cpgsite in bs region [each sample]
		# 3. intersect 

		if [ $region == "whole_bs_region" ]; then
			more +2 $bs_me_level_file | bedtools intersect -wa -wb -a - -b $hg_cpgsite_bed | awk -v value_column=$value_column '{OFS="\t"; print $1,$2,$3,$4}' |tee $bs_cpgsite_me_level | bedtools intersect -wa -wb -a - -b $bs_snp_file > $bs_snp_exist_cpgsite_me_level &
		else
			more +2 $bs_me_level_file | bedtools intersect -wa -wb -a - -b $hg_cpgsite_bed | awk -v value_column=$value_column '{OFS="\t"; print $1,$2,$3,$4}'| bedtools intersect -wa -wb -a - -b $region_file |tee $bs_cpgsite_me_level | bedtools intersect -wa -wb -a - -b $bs_snp_file > $bs_snp_exist_cpgsite_me_level &
		fi 
		let count+=1; [[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait

	# create table for box plot
	table_for_boxplot=$viz_result_dir/"table_for_boxplot_"$region
	echo -e "Subtype\tDNA_methylation_level\tMutation_rate" > $table_for_boxplot

	echo "[INFO] Compute mutation rate and plot"
	for((i=0;i<${#bs_r1_list[@]};i++)); do 

		# variables


		# bs_me
		input_fq1_filename_only=`basename \${bs_r1_list[$i]}`
		file_extention=$(echo $input_fq1_filename_only |awk -F . '{if (NF>1) {print $NF}}')
		filename1_wo_ext=`basename $input_fq1_filename_only "."$file_extention`

		bs_me_level_file=$bs_me_result_dir/$filename1_wo_ext"_val_1.fq_bismark_pe.sorted.sam.bed"
		bs_me_filename_only=`basename \$bs_me_level_file`

		# snp
		snp_from_bs_result_dir="/data/project/breast_TF_methyl/30celline_SNP_TFBS_ME_GE/result/snp_from_bs/"
		bam_prefix=$filename1_wo_ext"_val_1.fq_bismark_pe"

		bs_snp_file=$snp_from_bs_result_dir/$bam_prefix"snp.raw.sorted.filtered.vcf"
		bs_snp_filename_only=`basename \$bs_snp_file`

		bs_cpgsite_me_level=$viz_result_dir/$bs_me_filename_only".cpg_site_me_level_"$region
		bs_snp_exist_cpgsite_me_level=$viz_result_dir/$bs_me_filename_only".snp_exist_cpg_site_me_level_"$region
		bs_subtype=${bs_subtype_list[$i]}

		interval=5; max=100

		# for debug
		echo "[INFO] ME file : " $bs_me_level_file
		echo "[INFO] SNP file : " $bs_snp_file
		echo "[INFO] Subtype : " $bs_subtype
		#echo "[INFO] Region : " $region

		number_of_cpgsite_per_methylation_range=$viz_result_dir/$bs_me_filename_only".number_of_cpgsite_per_methylation_range_"$region
		number_of_snp_exist_cpgsite_per_methylation_range=$viz_result_dir/$bs_me_filename_only".number_of_snp_exist_cpgsite_per_methylation_range_"$region

		mutation_rate_per_methylation_range=$viz_result_dir/$bs_me_filename_only".mutation_rate_per_methylation_range_"$bs_subtype"_"$region

		# compute mutation rate of each methylation range
		## extract number of cpgsite per multation range
		binning_value $interval $max 4 $bs_cpgsite_me_level > $number_of_cpgsite_per_methylation_range &
		binning_value $interval $max 4 $bs_snp_exist_cpgsite_me_level  > $number_of_snp_exist_cpgsite_per_methylation_range &
		wait

		## compute mutation rate
		join -a1 -a2 <(sort -k1,1 $number_of_cpgsite_per_methylation_range) <(sort -k1,1 $number_of_snp_exist_cpgsite_per_methylation_range) | sort -k1,1n | awk '{OFS="\t"; if($3=="") $4=0; print $1,$3/$2 }' > $mutation_rate_per_methylation_range

		## plot mutation rate
					#echo "[INFO] out file" $mutation_rate_per_methylation_range".jpg"
					#Rscript $bin_dir/line_plot_with_two_inputs.r $mutation_rate_per_methylation_range $mutation_rate_per_methylation_range".jpg"

		# create table for ttest	
					#cut -f2 $mutation_rate_per_methylation_range | awk -f $bin_dir/transpose.awk >> $mutation_rate_table
					#for(i in 1:ncol(lu_data)){print(i);print(t.test(lu_data[,i],baa_data[,i]))}

		# create boxplot data	
		awk -v subtype=$bs_subtype '{OFS="\t";print subtype, $1, $2}' $mutation_rate_per_methylation_range >> $table_for_boxplot
	done
	# draw boxplot
	Rscript $bin_dir/draw_boxplot.R $table_for_boxplot $table_for_boxplot".jpg"

done
COMMENT

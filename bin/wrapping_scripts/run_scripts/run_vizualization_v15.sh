#!/bin/bash
source `dirname $0`/../../env.sh

# parse experiment information
source `dirname $0`/parse_exp_info.sh


NUM_CPUS=$SYS_NUM_CPUS

# directories
bin_dir="$WORK_DIR/bin"
<<'COMMENT'
result_dir="$WORK_DIR/result"
COMMENT
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

final_result_dir_ge=$final_result_root_dir/gene_expression
final_result_dir_me=$final_result_root_dir/methylation
final_result_dir_mu=$final_result_root_dir/mutation
final_result_dir_inte=$final_result_root_dir/integrated
UCSC_hub_dir_name="UCSC_hub/"
#UCSC_hub_result_dir="$final_result_dir_inte/$UCSC_hub_dir_name"
UCSC_hub_result_dir=$final_result_dir_inte

final_result_url="$WEB_ACCESSIBLE_LOC/images/$uid"
base_url_for_UCSC_hub="$final_result_url/integrated/"
delimeter_for_UCSC_hub=";"


mkdir -p $final_result_dir_ge
mkdir -p $final_result_dir_me
mkdir -p $final_result_dir_mu
mkdir -p $final_result_dir_inte
#mkdir -p $UCSC_hub_result_dir

final_result_url_ge=$final_result_url/gene_expression
final_result_url_me=$final_result_url/methylation
final_result_url_mu=$final_result_url/mutation
final_result_url_inte=$final_result_url/integrated



# directory based variables
ME_all_merged=$integ_result_dir/ME_MGD    # prepare initiall prefix for all merged file
ME_GE_all_merged=$integ_result_dir/ME_GE_MGD
ME_GE_COR_all_merged=$integ_result_dir/ME_GE_COR_MGD
ME_GE_COR_all_merged_THRESH=$integ_result_dir/ME_GE_COR_MGD_THRESH.txt


# KEGG color code based on log2 foldchange
red4="#ff6666" # 5<x : red4
red3="#ff8888" # 5>x>2 : red3
red2="#ffaaaa" # 2>x>1 : red2
red1="#ffcccc" # 1>x>0 : red1
white="#ffffff" # x=0 : white
blue4="#6666ff" # -5>x : blue4
blue3="#8888ff" # -5<x<-2 : blue3
blue2="#aaaaff" # -2<x<-1 : blue2





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


snp_file_extension=$(echo `basename \${mu_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')

####################################################################################
# 6. Visualization 
####################################################################################

mkdir -p $viz_result_dir
mkdir -p $UCSC_hub_result_dir

echo "[INFO] Start visualizing results"
#################################
# 6.1 For Raw data. (NGSplot, UCSC links, ME-density)
#################################

############
# 6.1.1 GE 
############

# TODO : set data kinds (RNAseq or array ..)
# UCSC links for each samples : BigWig
for (( i=0; i<${#ge_list[@]}; i++ )); do 
	# copy to web accessible directory
	# TODO :  this may be not necessary cuz already we are in the web accessible dir
	cp $ge_result_dir/`basename \${ge_list[$i]}`.bw $final_result_dir_ge
	# print UCSC url
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=test_SNP%20description=test2SNP%20visibility=dense%20bigDataUrl=$final_result_url_ge/${ge_list[$i]}.bw"
done


# GDEG heatmap & pca
#TODO
	GDEG_exp_value_table=$integ_result_dir/ALL.DEG.table.txt
	GDEG_exp_value_table_top100=$integ_result_dir/ALL.DEG.Top100.table.txt

	GDEG_top100_heatmap=$viz_result_dir/ALL.DEG.MGD.heatmap.png
	GDEG_PCA=$viz_result_dir/ALL.DEG.PCA.png

	# draw heatmap
	GDEG_top100_heatmap_input=$viz_result_dir"/ALL.DEG.Top100.table.txt.heatmap.input"

#type_len

	{
		echo -e "ID\t"`my_join "/" ${type_kind[@]}` | tr '/' '\t' ;

		awk -v sample_nums=`my_join "," ${type_len[@]}` 'BEGIN{FS=OFS="\t";num=split(sample_nums,arr,",")}NR>1{printf "%s", $1; my_index=1;for(i=1;i<=num;i++){sum=0;for(j=1;j<=arr[i];j++){sum=sum+$(my_index+j)}; printf "\t%f", sum/arr[i]; my_index=my_index+arr[i]};printf "\n"}' $GDEG_exp_value_table_top100 ;
	} > $GDEG_top100_heatmap_input 
	
	$R_DIR/Rscript $bin_dir/make_heatmap3_ver2_ge_v2.r $GDEG_top100_heatmap_input `my_join "," ${type_kind[@]}` "Class" `my_join "," ${type_kind[@]}` $GDEG_top100_heatmap "Top100_ALL_class_DEG"
				
	# draw PCA
	GDEG_PCA_input=$viz_result_dir"/ALL.DEG.table.txt.pca.input"
	cut -f1-$((${#type_list[@]}+1)) $GDEG_exp_value_table > $GDEG_PCA_input

	  ge_basename_list=()
  for ((i=0; i<${#ge_list[@]}; i++)); do
    ge_basename_list+=(`basename \${ge_list[\$i]}`)
  done

	$NEW_R_DIR/Rscript $bin_dir/pca.r $GDEG_PCA_input `my_join "," ${ge_basename_list[@]}` `my_join "," ${type_list[@]}` $viz_result_dir"/ALL.DEG.PCA.stat.txt" $viz_result_dir"/ALL.DEG_PCA.var.png" $GDEG_PCA
 	
	cp $GDEG_top100_heatmap $GDEG_PCA $final_result_dir_ge  

# TODO : DAVID & KEGG
	GDEG_gene_list=$integ_result_dir/GDEG_kw_genelist
	
	#DAVID
	bash $bin_dir/gsea3.sh $GDEG_gene_list $viz_result_dir "GENE_SYMBOL" "ALL.DEG.list"

	#KEGG

	# create KEGG input gene_symbol, color for up & down regulation
	# This will be depricated since we will generate KEGG figure locally below
			awk -v pval_threadhold=$P_VALUE_CUT -v red4=$red4 -v red3=$red3 -v red2=$red2 -v red1=$red1 \
				'NR==1{next;}{OFS="\t"; if ($(NF) <= (pval_threadhold*1/4)) print $1, red4;
					else if (($(NF) <= (pval_threadhold*2/4)) && ($(NF) > (pval_threadhold*1/4))) print $1, red3;
					else if (($(NF) <= (pval_threadhold*3/4) ) && ((pval_threadhold*2/4) < $(NF))) print $1, red2;
					else if (($(NF) <= pval_threadhold) && ( (pval_threadhold*3/4)< $(NF))) print $1, red1;}' $GDEG_exp_value_table > $viz_result_dir/"ALL.DEG.KEGG.txt"


	# generate KEGG figures based on DEG and their fold change or pvalue
	# TODO : kegg figure for GDEG since no FC

	cp $viz_result_dir"/"ALL.DEG.list.DAVID_*.txt $viz_result_dir"/"ALL.DEG.KEGG.txt $final_result_dir_ge

	

############
# 6.1.2 ME
############
# UCSC links for each samples : BigWig
echo "[INFO] Methyl UCSC Bigwig"

if [ "$me_data_type" -eq "0" ]; then # MBD
	count=0
	for (( i=0; i<${#me_list[@]}; i++ )); do
		{
			temp_filename_only=`basename \${me_list[$i]}`
			met_level_file=$me_result_dir/$temp_filename_only".met_level";

			# create bigWig file for UCSC
			awk 'NR==1{next;}{print}' $met_level_file | cut -f1,2,3,5 $viz_result_dir/$temp_filename_only".met_level.bed" ;
#| bedtools intersect -a - -b $REF_HUMAN_GENOME_100_BIN_BED -sorted > $viz_result_dir/$temp_filename_only".met_level.bed" ;
			$bin_dir/bedGraphToBigWig $viz_result_dir/$temp_filename_only".met_level.bed" $REF_HUMAN_CHR_SIZE $viz_result_dir/$temp_filename_only".met_level.bw";

			cp $viz_result_dir/$temp_filename_only".met_level.bw" $final_result_dir_me;
			# print UCSC url
			echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl=$final_result_url_me/"$temp_filename_only".met_level.bw";
		} &
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait

elif [ "$me_data_type" -eq "3" ]; then # BS
	ls $methylkit_from_bs_result_dir/*.sorted.bw | xargs -I {} cp {} $final_result_dir_me
	ls $methylkit_from_bs_result_dir/*.sorted.bw | xargs -I {} echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl="$final_result_url_me/{};
fi

# NGS plot
echo "[INFO] Methyl NGS plots"

if [ "$me_data_type" -eq "0" ] ; then # MBD
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
# TODO : ME density plot and comparison

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

			output_prefix=$viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.density"			

			pair_DMR=$me_result_dir"/"${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.bed"

			bedtools intersect -wa -a $avg_met_file1 -b $pair_DMR | cut -f4 > $output_prefix"_DMR_input1.txt"
			bedtools intersect -wa -a $avg_met_file2 -b $pair_DMR | cut -f4 > $output_prefix"_DMR_input2.txt"

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

	cp $viz_result_dir/*.DMR.density.png $final_result_dir_me

############
# 6.1.3 MU
############

# UCSC links
for (( i=0; i<${#mu_list[@]}; i++ )); do
	#already moved in run_mutation code
#	cp $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.gz" $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.gz.tbi" $final_result_dir
	# print UCSC url
	temp_filename_only=`basename \${mu_list[\$i]}`
	temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`
		
	vcf_filename=""

	if [ "$mu_data_type" -eq "0" ] || [ "$mu_data_type" -eq "1" ]; then
		vcf_filename=$temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf"
	else
		vcf_filename=""
	fi
	
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=vcfTabix%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl=$final_result_url_mu/$vcf_filename.gz"
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
			cp $ge_result_dir/$work_file $viz_result_dir/$work_file".gene_only.annot.html" $final_result_dir_inte
		done
	done

<<'COMMENT'
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do 
			work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}".deg"
			
		
			# get deg gene list only
			cut -f2 $rna_result_dir/$work_file > $viz_result_dir/$work_file".gene_only"
			sh $bin_dir/create_annotation_link_html_based_on_gene_list.sh $viz_result_dir/$work_file".gene_only" > $viz_result_dir/$work_file".gene_only.annot.html"

			# create annotation page for deg
			cp $rna_result_dir/$work_file $viz_result_dir/$work_file".gene_only" $viz_result_dir/$work_file".gene_only.annot.html" $final_result_dir
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
echo "generate circos for subtype-specific avg methylation with difference and DEGs"

# create subtype avg methyl file from merged file
for (( i=0; i<${#type_kind[@]}; i++ )); do 
	subtype_ME_merged_file=$ME_all_merged"_"${type_kind[$i]}
	subtype_ME_merged_avg_file=$viz_result_dir/`basename \$subtype_ME_merged_file`"_avg"
	subtype_ME_merged_avg_10Mb_file=$subtype_ME_merged_avg_file"_10Mb"
	subtype_ME_merged_avg_1Mb_file=$subtype_ME_merged_avg_file"_1Mb"

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
	bedtools map -header -c 4 -o mean -a $REF_HUMAN_GENOME_1M_BIN_BED -b $subtype_ME_merged_avg_file | sed -e 's/chr/hs/g' | sed -e 's/\t\./\t0/g' > $subtype_ME_merged_avg_1Mb_file
done


# TODO : Divide 1MB merged file to each chromosome to generated zoom ined each circos


# crate circos for pairwise & all class

all_type_file_list=()
for (( i=0; i<${#type_kind[@]}; i++ )); do 
	for (( j=$i+1; j<${#type_kind[@]}; j++ )); do 
		# met level file
		# generate circos input list deliminated by ';'	
		subtype_ME_merged_file=$ME_all_merged"_"${type_kind[$i]}
		subtype_ME_merged_file_2=$ME_all_merged"_"${type_kind[$j]}
		subtype_ME_merged_avg_file=$viz_result_dir/`basename \$subtype_ME_merged_file`"_avg"
		subtype_ME_merged_avg_file_2=$viz_result_dir/`basename \$subtype_ME_merged_file_2`"_avg"

		subtype_ME_merged_avg_10Mb_file=$subtype_ME_merged_avg_file"_10Mb"
		subtype_ME_merged_avg_10Mb_file_2=$subtype_ME_merged_avg_file_2"_10Mb"

		subtype_ME_merged_avg_1Mb_file=$subtype_ME_merged_avg_file"_1Mb"
		subtype_ME_merged_avg_1Mb_file_2=$subtype_ME_merged_avg_file_2"_1Mb"

	#	subtype_ME_merged_file="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/temp_result_files/ME_MGD_"${type_kind[$i]}".head100" # temp test with samll datda
	#	subtype_ME_merged_file_2="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/temp_result_files/ME_MGD_"${type_kind[$j]}".head100" # temp test with samll datda

		# for deg
	  work_file=${type_kind[$i]}"_vs_"${type_kind[$j]}
    DEG_file_name=$ge_result_dir"/"$work_file".DEG.list"
 #  DEG_file_name=/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/rna_paired_test_data/gene_exp/rna_seq/deg/Li_vs_Th.FC2.DEGs
		DEG_for_circos=$circos_dir/$work_file".DEGs.circos"
	#	DEG_for_circos=$circos_dir/$work_file".FC"$fold_change".DEGs.circos"
		diff_file_10Mb=$circos_dir"/circos_${type_kind[$i]}_${type_kind[$j]}_10Mb.diff"
		diff_file_1Mb=$circos_dir"/circos_${type_kind[$i]}_${type_kind[$j]}_1Mb.diff"
		file_list_10Mb=$subtype_ME_merged_avg_10Mb_file";"$subtype_ME_merged_avg_10Mb_file_2


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
		paste $subtype_ME_merged_avg_10Mb_file $subtype_ME_merged_avg_10Mb_file_2 | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1,$2,$3,abs($4-$8)}' OFS='\t' > $diff_file_10Mb
		paste $subtype_ME_merged_avg_1Mb_file $subtype_ME_merged_avg_1Mb_file_2 | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1,$2,$3,abs($4-$8)}' OFS='\t' > $diff_file_1Mb
		
		# create subtype pair-wise circos
		# TODO : set min max value for circos
		circos_min=0
		circos_max=15

		circos_diff_min="0"
		circos_diff_max="0.05"
		
		if [ "$me_data_type" -eq "0" ]; then # MBD
				circos_diff_min="0"				
				circos_diff_max="2"
				circos_min=0
				circos_max=15
		elif [ "$me_data_type" -eq "3" ]; then # BS
				circos_diff_min="0"
				circos_diff_max="50"
				circos_min=0
				circos_max=100
		else																	# array, TODO:Need to adjust this value for array
				circos_diff_min="0"
				circos_diff_max="0.05"
		fi

		temp_circos_file=$circos_dir"/circos_${type_kind[$i]}_${type_kind[$j]}.conf"
		bash $bin_dir/generate_circos_conf.sh 2 $file_list_10Mb $diff_file_10Mb $DEG_for_circos $CIRCOS_CONFIG_LOC $circos_min $circos_max $circos_diff_min $circos_diff_max $final_result_dir_inte > $temp_circos_file
		circos -conf $temp_circos_file -outputdir $circos_dir -outputfile $work_file".circos"

		# generate 1Mb circos per chromosome
		for hs in `cut -f1 REF_HUMAN_GENOME_1M_BIN_BED | sort -k1,1 | uniq`; do
			file_list_1Mb=$subtype_ME_merged_avg_1Mb_file"_$hs;"$subtype_ME_merged_avg_1Mb_file_2"_$hs"

			# get hs based data
			grep -w $hs $subtype_ME_merged_avg_1Mb_file > $subtype_ME_merged_avg_1Mb_file"_$hs"
			grep -w $hs $subtype_ME_merged_avg_1Mb_file_2 > $subtype_ME_merged_avg_1Mb_file_2"_$hs"

			temp_circos_file=$circos_dir"/circos_${type_kind[$i]}_${type_kind[$j]}_$hs.conf"
			bash $bin_dir/generate_circos_conf.sh 2 $file_list_1Mb $diff_file_1Mb $DEG_for_circos $CIRCOS_CONFIG_LOC $circos_min $circos_max $circos_diff_min $circos_diff_max $final_result_dir_inte > $temp_circos_file
			circos -conf $temp_circos_file -outputdir $circos_dir -outputfile $work_file".circos"
		done
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
join -1 1 -2 4 <(sort $GDEG"_genelist") $REF_HUMAN_GENE_GENEBODY_INFO_SORTED | awk '{print $2,$3,$4,$1}' OFS='\t' | sed -e 's/chr/hs/g' > $GDEG_for_circos

circos_min=0
circos_max=15
circos_diff_min=0

if [ "$me_data_type" -eq "0" ]; then # MBD
		circos_min=0
		circos_max=15
elif [ "$me_data_type" -eq "3" ]; then # BS
		circos_min=0
		circos_max=100
else																	# array, TODO:Need to adjust this value for array
	pass;
fi

# circos for all class
temp_circos_file=$circos_dir"/circos_all.conf"
bash $bin_dir/generate_circos_conf.sh ${#type_kind[@]} `my_join ";" ${all_type_file_list[@]}` $diff_all_file $GDEG_for_circos $CIRCOS_CONFIG_LOC $circos_min $circos_max $circos_diff_min $circos_diff_max $final_result_dir_inte > $temp_circos_file
circos -conf $temp_circos_file -outputdir $circos_dir -outputfile "ALL.circos"

# generate 1Mb circos per chromosome for all class
for hs in `cut -f1 REF_HUMAN_GENOME_1M_BIN_BED | sort -k1,1 | uniq`; do

	# get chromosomal file list	
	all_type_file_list_1Mb_hs=()
	for (( i=0; i<${#type_kind[@]}; i++ )); do
		# variables
		subtype_ME_merged_file=$ME_all_merged"_"${type_kind[$i]}
		subtype_ME_merged_avg_file=$viz_result_dir/`basename \$subtype_ME_merged_file`"_avg"
		subtype_ME_merged_avg_1Mb_file=$subtype_ME_merged_avg_file"_1Mb"
		all_type_file_list_1Mb_hs+=($subtype_ME_merged_avg_1Mb_file"_$hs")
	done

	diff_all_file_hs=$circos_dir"/circos_all_"$hs".diff"
	paste ${all_type_file_list_1Mb_hs[@]} | awk -v class_num=${#type_kind[@]} -v max_val=$((circos_diff_max + 1)) 'BEGIN{FS=OFS="\t"}{sum=0.0;for(i=1;i<=class_num;i++){met_arr[i]= 0.1 + $(4*i); sum=sum+met_arr[i]}; ent = 0.0; for(i=1;i<=class_num;i++){ent=ent-((met_arr[i]/sum)*(log(met_arr[i]/sum)/log(2)))}; print $1,$2,$3, max_val*(1-ent/(log(class_num)/log(2)))}' > $diff_all_file_hs

	temp_circos_file=$circos_dir"/circos_all_"$hs".conf"
	bash $bin_dir/generate_circos_conf.sh ${#type_kind[@]} `my_join ";" ${all_type_file_list_1Mb_hs[@]}` $diff_all_file_hs $GDEG_for_circos $CIRCOS_CONFIG_LOC $circos_min $circos_max $circos_diff_min $circos_diff_max $final_result_dir_inte > $temp_circos_file
	circos -conf $temp_circos_file -outputdir $circos_dir -outputfile "ALL_$hs"".circos"
done

cp $circos_dir/*.png $circos_dir/*.html $final_result_dir_inte

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
			cp $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw" $final_result_dir_me

			# print UCSC url
			echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=${type_kind[$k]}_${type_kind[$l]}_DMR%20description=${type_kind[$k]}_${type_kind[$l]}_DMR%20visibility=dense%20bigDataUrl="$final_result_url_me/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw"
		done
	done
elif [ "$me_data_type" -eq "3" ]; then

	for (( k=0; k<${#type_kind[@]}; k++ )); do
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do
			awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$7}' $methylkit_from_bs_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR_list.txt" > $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff"

			$bin_dir/bedGraphToBigWig $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff" $REF_HUMAN_CHR_SIZE $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw"

			# copy to web accessible directory
			cp $viz_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw" $final_result_dir_me

			# print UCSC url
			echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=${type_kind[$k]}_${type_kind[$l]}_DMR%20description=${type_kind[$k]}_${type_kind[$l]}_DMR%20visibility=dense%20bigDataUrl="$final_result_url_me/${type_kind[$k]}"_vs_"${type_kind[$l]}".DMR.diff.bw"


		done
	done

else

	#TODO infinium
	echo "infinium"

fi

########################################
# GDMR : bed -> bb
########################################
<<'COMMENT'
	$bin_dir/bedToBigBed $integ_result_dir/GDMR_kw $REF_HUMAN_CHR_SIZE $viz_result_dir/GDMR_by_kw.bb
	
	# copy to web accessible directory
	cp $viz_result_dir/GDMR_by_kw.bb $final_result_dir_me/ALL.DMR.bb
COMMENT
	# print UCSC url
#	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigBed%20name=GDMR_kw%20description=GDMR_kw%20visibility=dense%20bigDataUrl="$final_result_url_me/GDMR_by_kw.bb 

########################################
# DMR all : heatmap & pca
########################################
<<'COMMENT' #TOP 10%
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

	# copy to final_result_dir
	cp $viz_result_dir/DMR_TOP10_heatmap.png $final_result_dir_me
COMMENT

	#format
	##chr_bin_start_bin_end  HCC70_test_1.fq-methyl  HCC1954_test_1.fq-methyl        BT549_test_1.fq-methyl  mDAmB_231_test_1.fq-methyl  AU565_test_1.fq-methyl  HCC202_test_1.fq-methyl
	#chr10_8105400_8105500   96.126761       100.000000      87.837838       100.000000      45.744681       31.395349
	#chr10_8105800_8105900   94.214876       100.000000      85.714286       96.712329       29.545455       32.941176
	
	GDMR_with_value=$viz_result_dir/GDMR_methyl.txt
	awk -v kw_threshold=$P_VALUE_CUT 'NR==1{printf "%s", "chr_bin_start_bin_end";for(i=2;i<(NF-1);i++){printf "\t%s", $(i)};printf "\n";next}{if(($(NF) < kw_threshold) && ($(NF) != "nan")){printf "%s", $1;for(i=2;i<(NF-1);i++){printf "\t%f", $(i)};printf "\n"}}' $integ_result_dir/ME_MAT.kw.txt > $GDMR_with_value



  #base name list
  basename_list=()
  for((i=0;i<${#me_list[@]};i++)); do
    filename_only=`basename ${me_list[$i]}`
    basename_list+=($filename_only)
  done


	#drwa pca
		
	$NEW_R_DIR/Rscript $bin_dir/pca.r $GDMR_with_value `my_join "," ${basename_list[@]}`  `my_join "," ${type_list[@]}` $viz_result_dir"/ALL.DMR.PCA.stat.txt" $viz_result_dir"/ALL.DMR.PCA.var.png" $viz_result_dir"/ALL.DMR.PCA.png"

	# make heatmap input
	GDMR_with_value_promoter_bed=$GDMR_with_value".pro.bed"
	
	tail -n+2 $GDMR_with_value | tr '_' '\t' | bedtools intersect -wa -a - -b $REF_HUMAN_PROMOTER | uniq > $GDMR_with_value_promoter_bed

	GDMR_with_value_promoter_bed_heatmap_input=$GDMR_with_value_promoter_bed".heatmap.input"

	awk 'BEGIN{FS=OFS="\t"}{printf "%s_%s_%s", $1,$2,$3; for(i=4;i<=NF;i++){printf "\t%f", $(i)};printf "\n"}' $GDMR_with_value_promoter_bed > $GDMR_with_value_promoter_bed_heatmap_input

	# draw heatmap
	$R_DIR/Rscript $bin_dir/make_heatmap3_ver2.r $GDMR_with_value_promoter_bed_heatmap_input `my_join "," ${basename_list[@]}` "Class" `my_join "," ${type_kind[@]}` $viz_result_dir"/ALL.DMR_pro.heatmap.png" "ALL_class_DMR_overlap_with_Promoter"

##################################################
# DMR pairwise : heatmap & pca
#################################################
 for(( k=0; k<${#type_kind[@]}; k++)); do
	for(( l=$k+1; l<${#type_kind[@]}; l++)); do
		first_class_file=()
		second_class_file=()
		
		first_class_file_header=()
		second_class_file_header=()

		first_class=()
		second_class=()

		for(( i=0; i<${#me_list[@]}; i++)); do
			temp_filename_only=`basename \${me_list[$i]}`
			met_level_file_value_only=$integ_result_dir//$temp_filename_only".only.met_level"

			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				first_class_file+=($met_level_file_value_only)
				first_class_file_header+=($temp_filename_only)
				first_class+=(${type_kind[$k]})

			elif [ ${type_list[$i]} == ${type_kind[$l]} ]; then
				second_class_file+=($met_level_file_value_only)
				second_class_file_header+=($temp_filename_only)
				second_class+=(${type_kind[$l]})
			fi
		done

		work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}

		subtype_dmr=$me_result_dir"/"$work_file".DMR.bed"

		subtype_dmr_with_value=$viz_result_dir"/"$work_file".DMR.methyl.txt"
	
		echo "[DEBUG] : ${first_class_file[@]}"
		echo "[DEBUG] : ${second_class_file[@]}"
	
# TODO
	# echo header
#		class1_header_string=$(IFS=$'\t'; echo "${first_class_file_header[*]}")
#		class2_header_string=$(IFS=$'\t'; echo "${second_class_file_header[*]}")
		class1_header_string=`my_join "/" ${first_class_file_header[@]}`
		class2_header_string=`my_join "/" ${second_class_file_header[@]}`		

		bedtools unionbedg -i `my_join " " ${first_class_file[@]}` `my_join " " ${second_class_file[@]}` | bedtools intersect -wa -a - -b $subtype_dmr > $subtype_dmr_with_value
	
	
		{
			echo -e "ID\t"$class1_header_string"\t"$class2_header_string | tr '/' '\t' ;
			awk 'BEGIN{FS=OFS="\t"}{printf "%s", $1"_"$2"_"$3; for(i=4;i<=NF;i++){printf "\t%f", $(i)};printf "\n"}' $subtype_dmr_with_value ;
		} > $subtype_dmr_with_value".PCA.input"

		
		{
			echo -e "#ID\t"$class1_header_string"\t"$class2_header_string | tr '/' '\t' ;
			bedtools intersect -wa -a $subtype_dmr_with_value -b $REF_HUMAN_PROMOTER | uniq | awk 'BEGIN{FS=OFS="\t"}{printf "%s", $1"_"$2"_"$3; for(i=4;i<=NF;i++){printf "\t%f", $(i)};printf "\n"}' ;
		} > $subtype_dmr_with_value".promoter.heatmap.input"

			
		# draw heatmap
			
		$R_DIR/Rscript $bin_dir/make_heatmap3_ver2.r $subtype_dmr_with_value".promoter.heatmap.input" `my_join "," ${first_class_file_header[@]}`","`my_join "," ${second_class_file_header[@]}` "Class" ${type_kind[$k]}","${type_kind[$l]} $viz_result_dir"/"$work_file".DMR_pro.heatmap.png" "DMR_overlap_with_Promoter"

		# draw pca
		$NEW_R_DIR/Rscript $bin_dir/pca.r $subtype_dmr_with_value".PCA.input" `my_join "," ${first_class_file_header[@]}`","`my_join "," ${second_class_file_header[@]}`  `my_join "," ${first_class[@]}`","`my_join "," ${second_class[@]}` $viz_result_dir"/"$work_file".DMR.PCA.stat.txt" $viz_result_dir"/"$work_file".DMR.PCA.var.png" $viz_result_dir"/"$work_file".DMR.PCA.png"

	done
 done

	#move heatmap & pca plot
	cp $viz_result_dir/*.DMR.PCA.png $viz_result_dir/*.DMR_pro.heatmap.png $final_result_dir_me

########################################
# DMR region stat : bar chart
########################################
	echo "[INFO] DMR region stat on bar chart"

	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do 
		
			subtype_pair="${type_kind[$k]}_vs_${type_kind[$l]}"
			
			echo -e "Region\tcounts\tsubtype_pair" > $viz_result_dir/$subtype_pair".DMR_region_stat_matrix.txt"
			echo -e "Region\tcounts\tsubtype_pair" > $viz_result_dir/$subtype_pair".DMR_cpgi_stat_matrix.txt"

			temp_region_dmr=$viz_result_dir/$subtype_pair".DMR_region.temp"
	
			bedtools intersect -wa -wb -a $viz_result_dir/$subtype_pair.DMR.diff -b $REF_HUMAN_GENE_RANGE_INFO > $temp_region_dmr

			#for region in "3 prime UTR" "5 prime UTR" "Promoter" "Exon" "Intron"; do
			for region in "3utr" "5utr" "promoter" "exon" "intron"; do
				temp_region_dmr_count=`grep -w $region $temp_region_dmr | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`

				echo -e $region"\t"$temp_region_dmr_count"\t"$subtype_pair >> $viz_result_dir/$subtype_pair".DMR_region_stat_matrix.txt"
			done
				# create cpg_island, cpg_shore, cpg_shelf
			for region in "cpgIsland" "cpgShore" "cpgShelf"; do
				temp_region_dmr_count=`grep -w $region $temp_region_dmr | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`
				echo -e $region"\t"$temp_region_dmr_count"\t"$subtype_pair >> $viz_result_dir/$subtype_pair".DMR_cpgi_stat_matrix.txt"
			done
	
			# draw bar chart
		echo "[INFO] Draw bar chart : "$subtype_pair
	
		$NEW_R_DIR/Rscript $bin_dir/barplot_ver3.r $viz_result_dir/$subtype_pair".DMR_region_stat_matrix.txt" $viz_result_dir/$subtype_pair".DMR_cpgi_stat_matrix.txt" $viz_result_dir/$subtype_pair".DMR.barplot.png"

		
		done
	done

# for GDMR
	echo -e "Region\tcounts\tsubtype_pair" > $viz_result_dir/"ALL.DMR_region_stat_matrix.txt"
	echo -e "Region\tcounts\tsubtype_pair" > $viz_result_dir/"ALL.DMR_cpgi_stat_matrix.txt"

	temp_region_dmr=$viz_result_dir/"ALL.DMR_region.temp"
	
	bedtools intersect -wa -wb -a $integ_result_dir/GDMR_kw -b $REF_HUMAN_GENE_RANGE_INFO  > $temp_region_dmr
	
	for region in "3utr" "5utr" "promoter" "exon" "intron"; do
		temp_region_dmr_count=`grep -w $region $temp_region_dmr | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`

		echo -e $region"\t"$temp_region_dmr_count"\tALL_class" >> $viz_result_dir/ALL.DMR_region_stat_matrix.txt
	done

	for region in "cpgIsland" "cpgShore" "cpgShelf"; do
		temp_region_dmr_count=`grep -w $region $temp_region_dmr | cut -f1,2,3 | sort -k1,1 -k2,2n | uniq | wc -l`
		
		echo -e $region"\t"$temp_region_dmr_count"\tALL_class" >> $viz_result_dir/ALL.DMR_cpgi_stat_matrix.txt
	done

	echo "[INFO] Draw bar chart : ALL"
	
	$NEW_R_DIR/Rscript $bin_dir/barplot_ver3.r $viz_result_dir/"ALL.DMR_region_stat_matrix.txt" $viz_result_dir/"ALL.DMR_cpgi_stat_matrix.txt" $viz_result_dir/"ALL.DMR.barplot.png"

	
	# move DMR barplots to final result directory
	cp $viz_result_dir/*.DMR.barplot.png $final_result_dir_me

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
<<'COMMENT'
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
	snp_file_extension=$(echo `basename \${mu_list[0]}` | awk -F . '{if(NF>1) {print $NF}}')
	# for DSNP
	paste_list_all=()
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		paste_list_sub_type=()
		count=0
		for (( i=0; i<${#mu_list[@]}; i++ )); do 
			temp_filename_only=`basename \${mu_list[\$i]}`
			temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`
			vcf_bed_file=""
			if [ "$mu_data_type" -eq "0" ] || [ "$mu_data_type" -eq "1" ]; then # DNR-seq
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

	cp $viz_result_dir/SNP_*.png $final_result_dir_mu
done
COMMENT

# GDSNP all :  bed w/ label -> bb
	sort -k1,1 -k2,2n $integ_result_dir/GDSNP.bed | uniq > $viz_result_dir/GDSNP_sorted.bed
	$bin_dir/bedToBigBed $viz_result_dir/GDSNP_sorted.bed $REF_HUMAN_CHR_SIZE $viz_result_dir/GDSNP_sorted.bb

	# copy BigBed to web accessible dir
	cp $viz_result_dir/GDSNP_sorted.bb $final_result_dir_mu/ALL_SNP.bb

	# create UCSC link
#	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigBed%20name=GDSNP%20description=GDSNP%20visibility=pack%20bigDataUrl=$final_result_url_mu/GDSNP_sorted.bb"


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

	$NEW_R_DIR/Rscript $bin_dir/boxplot.r `ls $viz_result_dir/*sample_cnt_w_zero_for_others | tr "\n" ";"` `my_join ";" ${type_len[@]}` `my_join ";" ${type_kind[@]}` $viz_result_dir/ALL_SNP_boxplot_genomic.png $viz_result_dir/ALL_SNP_boxplot_cpgi.png

	cp $vis_result_dir/ALL_SNP_boxplot_genomic.png $viz_result_dir/ALL_SNP_boxplot_cpgi.png $final_result_dir_mu


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

	cp $viz_result_dir/BAR_GDSNP_RNG.png $viz_result_dir/GDSNP_RNG_STAT_MAT.txt $viz_result_dir/GDSNP_CPGs_STAT_MAT.txt $viz_result_dir/BAR_GDSNP_CPGs.png $final_result_dir
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
    
				if [ "$mu_data_type" -eq "0" ] || [ "$mu_data_type" -eq "1" ] ; then
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
    
				if [ "$mu_data_type" -eq "0" ] || [ "$mu_data_type" -eq "1" ] ; then
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
    
				if [ "$mu_data_type" -eq "0" ] || [ "$mu_data_type" -eq "1" ] ; then
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
#			$NEW_R_DIR/Rscript $bin_dir/mu_rate_boxplot.r $table_for_boxplot $table_for_boxplot"_mu_rate.png"			
			$NEW_R_DIR/Rscript $bin_dir/mu_rate_boxplot.r $table_for_boxplot $viz_result_dir/"mu_rate.png"
			
		done
	done
done

cp $viz_result_dir/"mu_rate.png" $final_result_dir_inte


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
#################################
# 6.3 Build UCSC hub for all UCSC browser data
#################################
<<'COMMENT_TEMP'
UCSC_result_file_list_with_path=()
result_file_type_list=()
UCSC_track_name_list=()
genome_version="hg19"
base_url_for_UCSC_hub="http://bhi2.snu.ac.kr:3000/UCSC_hub/"
bash $bin_dir/create_UCSC_hub.sh $UCSC_result_file_list_with_path $result_file_type_list $UCSC_track_name_list $UCSC_hub_result_dir $genome_version $base_url_for_UCSC_hub > .UCSC_hub.txt
COMMENT_TEMP

me_file_extension=$(echo ${me_list[0]} |awk -F . '{if (NF>1) {print $NF}}')

###################
# pairwise
####################

for((k=0; k<${#type_kind[@]}; k++)); do
  for((l=$k+1; l<${#type_kind[@]}; l++)); do
    work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}

		#TODO
		deg_file=$ge_result_dir"/"$work_file".DEG.list"
		deg_bed=$viz_result_dir"/"$work_file".DEG.bed"
		deg_bb=$viz_result_dir"/"$work_file".DEG.bb"
		#file1 : only DEG genesymbol
		#file2 : genesymbol genebody info
		# chr1    69090   70008   OR4F5
		
		awk 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""}FILENAME==ARGV[2]{if($4 in arr){print $0}}' $deg_file $GENESYMBOL_GENEBODY_INFO > $deg_bed

		$bin_dir/bedToBigBed $deg_bed $REF_HUMAN_CHR_SIZE $deg_bb
		
		dmr_bed=$me_result_dir/$work_file".DMR.bed"
		dmr_bb=$viz_result_dir/$work_file".DMR.bb"
		dmr_bed_sorted=$viz_result_dir/$work_file".DMR.bed.sorted"
	#TODO : need to sort. other methods are required to mimic sort	
		sort -k1,1 -k2,2n $dmr_bed > $dmr_bed_sorted
		$bin_dir/bedToBigBed $dmr_bed_sorted $REF_HUMAN_CHR_SIZE $dmr_bb

	#for gene expression
		first_class_exp_list=()
		first_class_exp_filename=()
		second_class_exp_list=()
		second_class_exp_filename=()
		exp_type_list=()

		for (( i=0; i<${#ge_list[@]}; i++ )); do
			temp_filename_only=`basename \${ge_list[\$i]}`  
			
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				first_class_exp_list+=($final_result_dir_ge/$temp_filename_only".bw")
				first_class_exp_filename+=($temp_filename_only)
				exp_type_list+=("bigWig")
			elif [ ${type_list[$i]} == ${type_kind[$l]} ]; then
				second_class_exp_list+=($final_result_dir_ge/$temp_filename_only".bw")
				second_class_exp_filename+=($temp_filename_only)
				exp_type_list+=("bigWig")
			fi
		done

	#for methylation
		first_class_methyl_list=()
    first_class_methyl_filename=()
    second_class_methyl_list=()
    second_class_methyl_filename=()
    methyl_type_list=()

    for (( i=0; i<${#me_list[@]}; i++ )); do
      temp_filename_only=`basename \${me_list[\$i]}`

			met_sample_bw=""

			temp_filename_only_wo_ext=`basename \${me_list[\$i]} "."$me_file_extension`
			
			if [ "$me_data_type" -eq "0" ]; then # MBD
				met_sample_bw=$final_result_dir_me/$temp_filename_only".met_level.bw"
			elif [ "$me_data_type" -eq "3" ]; then # BS
				met_sample_bw=$final_result_dir_me/$temp_filename_only_wo_ext"_bismark_bt2_pe.sorted.sam.bed.sorted.bw"
			fi
      
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
        first_class_methyl_list+=($met_sample_bw)
        first_class_methyl_filename+=($temp_filename_only)
        methyl_type_list+=("bigWig")
      elif [ ${type_list[$i]} == ${type_kind[$l]} ]; then
        second_class_methyl_list+=($met_sample_bw)
        second_class_methyl_filename+=($temp_filename_only)
        methyl_type_list+=("bigWig")
      fi
    done


	#for mutaiton

		first_class_vcf_list=()
		first_class_vcf_filename=()
		second_class_vcf_list=()
		second_class_vcf_filename=()
		vcf_type_list=()

		for((i=0;i<${#mu_list[@]};i++)); do
			temp_filename_only=`basename \${mu_list[\$i]}`
	  	temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`

	  	vcf_filename=""
			
	  	if [ "$mu_data_type" -eq "0" ] || [ "$mu_data_type" -eq "1" ]; then
  		  vcf_filename=$temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf"
			else
  		  vcf_filename=""
			
 			fi


			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				first_class_vcf_list+=($final_result_dir_mu/$vcf_filename".gz")	
 				first_class_vcf_filename+=($temp_filename_only)
	 			vcf_type_list+=("vcfTabix")
			elif [ ${type_list[$i]} == ${type_kind[$l]} ]; then
				second_class_vcf_list+=($final_result_dir_mu/$vcf_filename".gz")
  			second_class_vcf_filename+=($temp_filename_only)
				vcf_type_list+=("vcfTabix")
			fi
		done

		#UCSC result file path config
		UCSC_result_file_list_with_ge_path=$deg_bb";"`my_join ";" ${first_class_exp_list[@]}`";"`my_join ";" ${second_class_exp_list[@]}`
		UCSC_result_file_list_with_me_path=$dmr_bb";"`my_join ";" ${first_class_methyl_list[@]}`";"`my_join ";" ${second_class_methyl_list[@]}`
		UCSC_result_file_list_with_mu_path=`my_join ";" ${first_class_vcf_list[@]}`";"`my_join ";" ${second_class_vcf_list[@]}`
		UCSC_result_file_list_with_path=$UCSC_result_file_list_with_ge_path";"$UCSC_result_file_list_with_me_path";"$UCSC_result_file_list_with_mu_path

		# result file type config
		result_file_type_list="bigBed;"`my_join ";" ${exp_type_list[@]}`";bigBed;"`my_join ";" ${methyl_type_list[@]}`";"`my_join ";" ${vcf_type_list[@]}`

		# track name config
		UCSC_track_name_ge_list=$work_file".DEG;"`my_join ";" ${first_class_exp_filename[@]}`";"`my_join ";" ${second_class_exp_filename[@]}`
		UCSC_track_name_me_list=$work_file".DMR;"`my_join ";" ${first_class_methyl_filename[@]}`";"`my_join ";" ${second_class_methyl_filename[@]}`
		UCSC_track_name_mu_list=`my_join ";" ${first_class_vcf_filename[@]}`";"`my_join ";" ${second_class_vcf_filename[@]}`
		UCSC_track_name_list=$UCSC_track_name_ge_list";"$UCSC_track_name_me_list";"$UCSC_track_name_mu_list

		bash $bin_dir/create_UCSC_hub.sh $UCSC_result_file_list_with_path $result_file_type_list $UCSC_track_name_list $UCSC_hub_result_dir $GENOME_VERSION $delimeter_for_UCSC_hub $base_url_for_UCSC_hub $work_file > $UCSC_hub_result_dir/$work_file".UCSC_hub_link.txt"
	
	done
done

###################
# ALL class
####################
		#TODO
		deg_file=$integ_result_dir"/GDEG_kw_genelist"
		deg_bed=$viz_result_dir"/ALL.DEG.bed"
		deg_bb=$viz_result_dir"/ALL.DEG.bb"
			
		awk 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""}FILENAME==ARGV[2]{if($4 in arr){print $0}}' $deg_file $GENESYMBOL_GENEBODY_INFO > $deg_bed

		$bin_dir/bedToBigBed $deg_bed $REF_HUMAN_CHR_SIZE $deg_bb
		
		dmr_bed=$integ_result_dir/GDMR_kw
		dmr_bb=$viz_result_dir"/ALL.DMR.bb"
		dmr_bed_sorted=$viz_result_dir/"ALL.DMR.bed.sorted"
	#TODO : need to sort. other methods are required to mimic sort	
		sort -k1,1 -k2,2n $dmr_bed > $dmr_bed_sorted
		$bin_dir/bedToBigBed $dmr_bed_sorted $REF_HUMAN_CHR_SIZE $dmr_bb

	#for gene expression
		first_class_exp_list=()
		first_class_exp_filename=()
		exp_type_list=()

		for (( i=0; i<${#ge_list[@]}; i++ )); do
			temp_filename_only=`basename \${ge_list[\$i]}`  
			
			first_class_exp_list+=($final_result_dir_ge/$temp_filename_only".bw")
			first_class_exp_filename+=($temp_filename_only)
			exp_type_list+=("bigWig")
		done

	#for methylation
		first_class_methyl_list=()
    first_class_methyl_filename=()
    methyl_type_list=()

    for (( i=0; i<${#me_list[@]}; i++ )); do
      temp_filename_only=`basename \${me_list[\$i]}`

			met_sample_bw=""

			temp_filename_only_wo_ext=`basename \${me_list[\$i]} "."$me_file_extension`
			
			if [ "$me_data_type" -eq "0" ]; then # MBD
				met_sample_bw=$final_result_dir_me/$temp_filename_only".met_level.bw"
			elif [ "$me_data_type" -eq "3" ]; then # BS
				met_sample_bw=$final_result_dir_me/$temp_filename_only_wo_ext"_bismark_bt2_pe.sorted.sam.bed.sorted.bw"
			fi
      
       first_class_methyl_list+=($met_sample_bw)
       first_class_methyl_filename+=($temp_filename_only)
       methyl_type_list+=("bigWig")
    done


	#for mutaiton

		first_class_vcf_list=()
		first_class_vcf_filename=()
		vcf_type_list=()

		for((i=0;i<${#mu_list[@]};i++)); do
			temp_filename_only=`basename \${mu_list[\$i]}`
	  	temp_filename_only_wo_extension=`basename \${mu_list[$i]} "."$snp_file_extension`

	  	vcf_filename=""
			
	  	if [ "$mu_data_type" -eq "0" ] || [ "$mu_data_type" -eq "1" ]; then
  		  vcf_filename=$temp_filename_only_wo_extension"_Aligned.out.split.filtered.vcf"
			else
  		  vcf_filename=""
			
 			fi


			first_class_vcf_list+=($final_result_dir_mu/$vcf_filename".gz")	
 			first_class_vcf_filename+=($temp_filename_only)
	 		vcf_type_list+=("vcfTabix")
		done

		echo "[INFO] Create UCSC hub file"

		#UCSC result file path config
		UCSC_result_file_list_with_ge_path=$deg_bb";"`my_join ";" ${first_class_exp_list[@]}`
		UCSC_result_file_list_with_me_path=$dmr_bb";"`my_join ";" ${first_class_methyl_list[@]}`
		UCSC_result_file_list_with_mu_path=`my_join ";" ${first_class_vcf_list[@]}`
		UCSC_result_file_list_with_path=$UCSC_result_file_list_with_ge_path";"$UCSC_result_file_list_with_me_path";"$UCSC_result_file_list_with_mu_path

		# result file type config
		result_file_type_list="bigBed;"`my_join ";" ${exp_type_list[@]}`";bigBed;"`my_join ";" ${methyl_type_list[@]}`";"`my_join ";" ${vcf_type_list[@]}`


		# track name config
		UCSC_track_name_ge_list="ALL_class.DEG;"`my_join ";" ${first_class_exp_filename[@]}`
		UCSC_track_name_me_list="ALL_class.DMR;"`my_join ";" ${first_class_methyl_filename[@]}`
		UCSC_track_name_mu_list=`my_join ";" ${first_class_vcf_filename[@]}`
		UCSC_track_name_list=$UCSC_track_name_ge_list";"$UCSC_track_name_me_list";"$UCSC_track_name_mu_list

		bash $bin_dir/create_UCSC_hub.sh $UCSC_result_file_list_with_path $result_file_type_list $UCSC_track_name_list $UCSC_hub_result_dir $GENOME_VERSION $delimeter_for_UCSC_hub $base_url_for_UCSC_hub "ALL" > $UCSC_hub_result_dir/"ALL.UCSC_hub_link.txt"

	




####################################
# 6.4 GE - MU / Oncoprint
###################################
# FOR pair-wise

for((k=0; k<${#type_kind[@]}; k++)); do
	for((l=$k+1; l<${#type_kind[@]}; l++)); do
		work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}

		first_class_filename=()
		second_class_filename=()
		
		for((i=0; i<${#mu_list[@]};i++)); do
			basename_only=`basename \${mu_list[\$i]}`
			
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				first_class_filename+=($basename_only)
			elif [ ${type_list[$i]} == ${type_kind[$l]} ]; then
				second_class_filename+=($basename_only)
			fi
		done

		top30_degs_genelist=$ge_result_dir/$work_file".DEG.Top30.genelist"
		
		snp_ge_MGD=$integ_result_dir/SNP_GE_MGD_$work_file".txt"

		oncoprint_input_count=$viz_result_dir/SNP_GE_MGD_$work_file".oncoprint.input_count.txt"
		oncoprint_input=$viz_result_dir/SNP_GE_MGD_$work_file".oncoprint.input.txt"

		tail -n+2 $snp_ge_MGD | grep -w "genebody" | awk -v genesymbol_index=$(((${type_len[$k]}+${type_len[$l]})*2+8)) -v sample_nums=$((${type_len[$k]} + ${type_len[$l]})) 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""} FILENAME==ARGV[2]{if($(genesymbol_index) in arr){printf "%s", $(genesymbol_index); for(i=4;i<(4+sample_nums);i++){if($(i) == "-/-"){printf "\t0"}else{printf "\t1"}};printf "\n"}}' $top30_degs_genelist - | sort -k1,1 | awk 'BEGIN{FS=OFS="\t";sym=""}{if(sym != $1){if(sym != ""){printf "%s", sym; for(i=2;i<=NF;i++){printf "\t%d", arr[i]}; printf "\n"};sym=$1; for(i=2;i<=NF;i++){arr[i]=$(i)};}else{ for(i=2;i<=NF;i++){arr[i]=arr[i]+$(i)}}}END{printf "%s", sym; for(i=2;i<=NF;i++){printf "\t%d", arr[i]}; printf "\n"}' > $oncoprint_input_count
		
		awk -v class_string=${type_kind[$k]}","${type_kind[$l]} -v sample_nums=${type_len[$k]}","${type_len[$l]} 'BEGIN{FS=OFS="\t"; class_num=split(class_string, class_arr, ","); sample_num=split(sample_nums, sample_arr, ",")}{sum=0; for(i=2;i<=NF;i++){sum=sum+$(i)}; if(sum==0){next;}; printf "%s", $1; my_index=1; for(i=1;i<=class_num;i++){for(j=1;j<=sample_arr[i];j++){if($(my_index+j) > 0){printf "\t%s", class_arr[i]}else{printf "\t"}}; my_index = my_index + sample_arr[i]};printf "\n"}' $oncoprint_input_count > $oncoprint_input
	
		# make oncoprint R code
		class_list=${type_kind[$k]}";"${type_kind[$l]}
		color_list=${TEN_COLOR_CODE[$k]}";"${TEN_COLOR_CODE[$l]}
		bash $bin_dir/gen_oncoprint_R.ver3.sh $class_list $color_list > $viz_result_dir/$work_file".oncoprint.R"
		
		# draw oncoprint plot
		$NEW_R_DIR/Rscript $viz_result_dir/$work_file".oncoprint.R" $oncoprint_input $class_list `my_join ";" ${first_class_filename[@]}`";"`my_join ";" ${second_class_filename[@]}` $work_file" SNP" $viz_result_dir/$work_file".onco.png" $bin_dir/oncoprint_1.6_bugfix.R

	done
done

# FOR ALL
	header_string=()
	for((i=0;i<${#mu_list[@]};i++)); do
		basename_only=`basename \${mu_list[\$i]}`
		
		header_string+=($basename_only)
	done
	
	top30_degs_genelist=$integ_result_dir/ALL.DEG.Top30.genelist
	
	snp_ge_MGD=$integ_result_dir/"SNP_GE_MGD.txt"

	work_file="ALL"

  oncoprint_input_count=$viz_result_dir/SNP_GE_MGD_$work_file".oncoprint.input_count.txt"
  oncoprint_input=$viz_result_dir/SNP_GE_MGD_$work_file".oncoprint.input.txt"


	tail -n+2 $snp_ge_MGD | grep -w "genebody" | awk -v genesymbol_index=$((${#mu_list[@]}*2+8)) -v sample_nums=${#mu_list[@]} 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=""} FILENAME==ARGV[2]{if($(genesymbol_index) in arr){printf "%s", $(genesymbol_index); for(i=4;i<(4+sample_nums);i++){if($(i) == "-/-"){printf "\t0"}else{printf "\t1"}};printf "\n"}}' $top30_degs_genelist - | sort -k1,1 | awk 'BEGIN{FS=OFS="\t";sym=""}{if(sym != $1){if(sym != ""){printf "%s", sym; for(i=2;i<=NF;i++){printf "\t%d", arr[i]}; printf "\n"};sym=$1; for(i=2;i<=NF;i++){arr[i]=$(i)};}else{ for(i=2;i<=NF;i++){arr[i]=arr[i]+$(i)}}}END{printf "%s", sym; for(i=2;i<=NF;i++){printf "\t%d", arr[i]}; printf "\n"}' > $oncoprint_input_count
		
		awk -v class_string=`my_join "," ${type_kind[@]}` -v sample_nums=`my_join "," ${type_len[@]}` 'BEGIN{FS=OFS="\t"; class_num=split(class_string, class_arr, ","); sample_num=split(sample_nums, sample_arr, ",")}{sum=0; for(i=2;i<=NF;i++){sum=sum+$(i)}; if(sum==0){next;}; printf "%s", $1; my_index=1; for(i=1;i<=class_num;i++){for(j=1;j<=sample_arr[i];j++){if($(my_index+j) > 0){printf "\t%s", class_arr[i]}else{printf "\t"}}; my_index = my_index + sample_arr[i]};printf "\n"}' $oncoprint_input_count > $oncoprint_input
	
		# make oncoprint R code
		class_list=`my_join ";" ${type_kind[@]}`

		color_temp_list=()
		for((i=0; i<${#type_kind[@]};i++)); do
			color_temp_list+=(${TEN_COLOR_CODE[$i]})
		done

		color_list=`my_join ";" ${color_temp_list[@]}`
		bash $bin_dir/gen_oncoprint_R.ver3.sh $class_list $color_list > $viz_result_dir/$work_file".oncoprint.R"
		
		# draw oncoprint plot
		$NEW_R_DIR/Rscript $viz_result_dir/$work_file".oncoprint.R" $oncoprint_input $class_list `my_join ";" ${header_string[@]}` "ALL_class SNP" $viz_result_dir/$work_file".onco.png" $bin_dir/oncoprint_1.6_bugfix.R

# move onco result to final directory
cp $viz_result_dir/*.onco.png $final_result_dir_inte


###############################################################
## Final Result JSON Create ##########################
###############################################################
function json_elem() {
  echo -e "\t\t\t\t\t{\"name\" : \"$1\"," ;
  echo -e "\t\t\t\t\t\"grid\" : $2," ;
  echo -e "\t\t\t\t\t\"type\" : \"$3\"" ;
  if [ "$4" -eq "1" ]; then
    echo -e "\t\t\t\t\t},"
  else
    echo -e "\t\t\t\t\t}"
  fi
}

function json_tag_elem() {
  echo -e "\t\t\t\t\t{\"name\" : \"$1\"," ;
  echo -e "\t\t\t\t\t\"grid\" : $2," ;
  echo -e "\t\t\t\t\t\"tag\" : \"$3\"," ;
  echo -e "\t\t\t\t\t\"type\" : \"$4\"" ;
  if [ "$5" -eq "1" ]; then
    echo -e "\t\t\t\t\t},"
  else
    echo -e "\t\t\t\t\t}"
  fi
}


function json_dummy() {
  echo -e "\t\t\t\t\t{\"grid\" : $1" ;
  if [ "$2" -eq "1" ]; then
    echo -e "\t\t\t\t\t},"
  else
    echo -e "\t\t\t\t\t}"
  fi
}

final_json=$final_result_root_dir/"final_result.json"

{
	echo -e "{" ;
	echo -e "\t\"result\" : [" ;
##############################
	# for gene_expression
################################
	echo -e "\t\t{\"gene_expression\" : [" ;
	
	for((k=0; k<${#type_kind[@]}; k++)); do
		for((l=$k+1; l<${#type_kind[@]}; l++)); do
			work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}
			
			#echo -e "\t\t\t{\"${type_kind[$k]}"_vs_"${type_kind[$l]}\" : [" ;
			echo -e "\t\t\t{\"$work_file\" : [" ;
	
			## title	
			json_elem "All sample status" 12 title 1 ;

			# MA plot
			json_elem $work_file".MA_plot.png" 6 image 1 ;
		
			# box plot
			json_elem $work_file".norm.boxplot.png" 6 image 1 ;
			
			# NGS plot for RNA-seq
			
			if [ "$ge_data_type" -eq "0" ]; then
				json_elem $work_file".NGS.genebody.png" 6 image 1 ;
				json_elem $work_file".NGS.cpgi.png" 6 image 1 ;
			fi
	
			## title	
			json_tag_elem "Differentially expressed genes (DEG)s" 12 DEG title 1 ;

			# DEG Top 100 table
			json_elem $work_file".DEG.Top100.table.txt" 12 table 1 ;
			
			#dummy
			json_elem "dummy" 6 "dummy" 1 ;

			# Total table
			json_elem $work_file".Total.table.txt" 3 link 1 ;

			# DEG table
			json_elem $work_file".DEG.table.txt" 3 link 1 ;

			## title	
			json_elem "Plots of DEGs" 12 title 1 ;

			# volcano plot
			json_elem $work_file".volcano.png" 4 image 1 ;

			# DEG Top 100 heatmap
			json_elem $work_file".DEG.MGD.heatmap.png" 4 image 1 ;

			# DEG PCA
			json_elem $work_file".DEG.PCA.png" 4 image 1 ;

			## title	
			json_tag_elem "Gene set analyses (GSA)s" 12 GENE_SET title 1 ;

			# read map list to array
			IFS=$'\r\n' map_list=($(cat $final_result_dir_ge/$work_file".deg.for_KEGG_fig_maplist.txt"))

			# DAVID. print final"," depends on final kegg existance
			if [ ${#map_list[@]} -gt 0 ]; then
				json_elem $work_file".DEG.list.DAVID_summary.txt" 12 link 1 ;
	
				## title	
				json_tag_elem "Pathway analysis" 12 PATHWAY title 1 ;

				# KEGG : list up all the genated images
				for (( i=0; i<${#map_list[@]}; i++ )); do
					if [ $i -eq $((${#map_list[@]} - 1)) ]; then
						json_elem ${map_list[$i]}"."$work_file".png" 2 image 0 ;
					else
						json_elem ${map_list[$i]}"."$work_file".png" 2 image 1 ;
					fi
				done

			else
				json_elem $work_file".DEG.list.DAVID_summary.txt" 12 link 0 ;
			fi

		#	for kegg_map in `cat \$workfile".deg.for_KEGG_fig_maplist.txt"`;do
				#hsa00983.TEST_BaA_vs_TEST_BaB.png
				
		#		json_elem $kegg_map"."$work_file".png" 2 image 1 ;
		#	done


		
			echo -e "\t\t\t\t]" ;
		
			if [ ${#type_kind[@]} -gt 2 ]; then
				echo -e "\t\t\t}, " ;
			else
				if [ $k -eq $((${#type_kind[@]} -2 )) ] && [ $l -eq $((${#type_kind[@]} -1 )) ]; then
					echo -e "\t\t\t}" ;
				else
					echo -e "\t\t\t}, " ;	
				fi
			fi
		done
	done

	if [ ${#type_kind[@]} -gt 2 ]; then
		# for all-class
		 echo -e "\t\t\t{\"ALL_class\" : [" ;	
		 work_file="ALL"

			## title	
			json_elem "Group DEGs" 12 title 1 ;

			# GDEG Top 100 table
			json_elem $work_file".DEG.Top100.table.txt" 12 table 1 ;
			
			#dummy
      json_elem "dummy" 6 dummy 1 ;

      # Total table
      json_elem $work_file".Total.table.txt" 3 link 1 ;

      # GDEG table
      json_elem $work_file".DEG.table.txt" 3 link 1 ;

			## title	
			json_elem "Plots of Group DEGs" 12 title 1 ;

			# GDEG Top 100 heatmap
      json_elem $work_file".DEG.MGD.heatmap.png" 6 image 1 ;
		
			# GDEG PCA
      json_elem $work_file".DEG.PCA.png" 6 image 1 ;
		
			## title	
			json_tag_elem "Gene set analyses (GSA)s" 12 GENE_SET title 1 ;

			# read map list to array
			IFS=$'\r\n' map_list=($(cat $final_result_dir_ge/$work_file".deg.for_KEGG_fig_maplist.txt"))

			# DAVID
			if [ ${#map_list[@]} -gt 0 ]; then
      	json_elem $work_file".DEG.list.DAVID_summary.txt" 12 link 1 ;

      # KEGG
				for (( i=0; i<${#map_list[@]}; i++ )); do
					if [ $i -eq $((${#map_list[@]} - 1)) ]; then
						json_elem ${map_list[$i]}"."$work_file".png" 2 image 0 ;
					else
						json_elem ${map_list[$i]}"."$work_file".png" 2 image 1 ;
					fi
				done
			else
      	json_elem $work_file".DEG.list.DAVID_summary.txt" 12 link 0 ;
			fi



      #json_elem $work_file".DEG.KEGG.txt" 6 link 0 ;
		
			echo -e "\t\t\t\t]" ;
			echo -e "\t\t\t} " ;	
	fi
##########################
	# for methylation
##########################
	echo -e "\t\t], \"methylation\" : [" ;
		echo -e "\t\t\t{\"Per_Sample\" : [" ;
		
		## title	
		json_elem "All sample methylation coverage and stats" 12 title 1 ;

		for((i=0; i<${#me_list[@]} ; i++)); do
			work_file=`basename \${me_list[\$i]}`
			file_extention=$(echo $work_file |awk -F . '{if (NF>1) {print $NF}}')
			filename_wo_ext=`basename $work_file "."$file_extention`

			if [ "$me_data_type" -eq 3 ]; then # BS-seq
				
				json_elem $filename_wo_ext"_bismark_bt2_pe.sorted.sam.CoverageStats.png" 2 image 1 ;
				if [ $i -eq $((${#me_list[@]} - 1)) ]; then
					json_elem $filename_wo_ext"_bismark_bt2_pe.sorted.sam.MethylStats.png" 2 image 0 ;
				else
					json_elem $filename_wo_ext"_bismark_bt2_pe.sorted.sam.MethylStats.png" 2 image 1 ;
				fi
			elif [ "$me_data_type" -eq 0 ]; then # MBD-seq
				if [ $i -eq $((${#me_list[@]} - 1)) ]; then
					json_elem $work_file".CoverageStats.png" 2 image 0 ;
				else
					json_elem $work_file".CoverageStats.png" 2 image 1 ;
				fi
			fi
		done
	
		echo -e "\t\t\t\t]" ;
		echo -e "\t\t\t}, " ;

	for((k=0;k<${#type_kind[@]};k++)); do
		for((l=$k+1;l<${#type_kind[@]};l++)); do
			work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}
			
			#echo -e "\t\t\t{\"${type_kind[$k]}"_vs_"${type_kind[$l]}\" : [" ;
			echo -e "\t\t\t{\"$work_file\" : [" ;
		

			# NGS plot for MBD-seq
			if [ "$me_data_type" -eq "0" ]; then
				## title	
				json_elem "Region specific average methylation" 12 title 1 ;

				json_elem $work_file".NGS.genebody.png" 6 image 1 ;
				json_elem $work_file".NGS.cpgi.png" 6 image 1 ;
			elif [ "$me_data_type" -eq "3" ]; then
				## title	
				json_elem "Pair-wise correlation" 12 title 1 ;

				json_elem $work_file".corr.png" 12 image 1 ;
			fi

			## title	
			json_elem "Pair-wise methylation density & Hierachcal clusturing with headmap" 12 title 1 ;

			# Density plot
			json_elem $work_file".DMR.density.png" 6 image 1 ;

			# promoter heatmap
			json_elem $work_file".DMR_pro.heatmap.png" 6 image 1 ;
		
			## title	
			json_elem "Regional DMR bin stats & PCA" 12 title 1 ;

			# DMR barplot
			json_elem $work_file".DMR.barplot.png" 6 image 1 ;

			# DMR PCA
			json_elem $work_file".DMR.PCA.png" 6 image 0 ;

			echo -e "\t\t\t\t]" ;
		
			if [ ${#type_kind[@]} -gt 2 ]; then
				echo -e "\t\t\t}, " ;
			else
				if [ $k -eq $((${#type_kind[@]} -2 )) ] && [ $l -eq $((${#type_kind[@]} -1 )) ]; then
					echo -e "\t\t\t}" ;
				else
					echo -e "\t\t\t}, " ;	
				fi

			fi


		done
	done

	if [ ${#type_kind[@]} -gt 2 ]; then
		# for all-class
		 echo -e "\t\t\t{\"ALL_class\" : [" ;	
		 work_file="ALL"

			## title	
			json_elem "All clase-wise Hierachcal clusturing with headmap" 12 title 1 ;

			# GDMR promoter heatamp
			json_elem $work_file".DMR_pro.heatmap.png" 12 image 1 ;
			
			## title	
			json_elem "Regional DMR bin stats & PCA" 12 title 1 ;

			# GDMR barplot
      json_elem $work_file".DMR.barplot.png" 6 image 1 ;

			# GDMR PCA
			json_elem $work_file".DMR.PCA.png" 6 image 0 ;
	
			echo -e "\t\t\t\t]" ;
		
			echo -e "\t\t\t} " ;	
	fi
##########################
	# for mutation
##########################
	echo -e "\t\t], \"mutation\" : [" ;

	for((k=0; k<${#type_kind[@]}; k++)); do
		for((l=$k+1; l<${#type_kind[@]}; l++)); do
			work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}
			
			#echo -e "\t\t\t{\"${type_kind[$k]}"_vs_"${type_kind[$l]}\" : [" ;
			echo -e "\t\t\t{\"$work_file\" : [" ;
			
			## title	
			json_elem "SNP and Indel status" 12 title 1 ;

			# total count
			json_elem $work_file".snp.tot.count.png" 6 image 1 ;
			json_elem $work_file".indel.tot.count.png" 6 image 1 ;
			
			# overlap AF
			json_elem $work_file".snp.overlap.AF.png" 6 image 1 ;
			json_elem $work_file".indel.overlap.AF.png" 6 image 1 ;

			## title	
			json_elem "TS/TV QUAL" 12 title 1 ;

			# ts_tv qual
			json_elem $work_file".ts_tv.qual.classA.png" 4 image 1 ;
			json_elem $work_file".ts_tv.qual.classB.png" 4 image 1 ;
			json_elem $work_file".ts_tv.qual.classAB.png" 4 image 1 ;

			## title	
			json_elem "Indel distribution" 12 title 1 ;

			# indel dist
			json_elem $work_file".indel_dist.classA.png" 4 image 1 ;
			json_elem $work_file".indel_dist.classB.png" 4 image 1 ;
			json_elem $work_file".indel_dist.classAB.png" 4 image 1 ;

			## title	
			json_elem "Depth disribution" 12 title 1 ;

			# depth dist
			json_elem $work_file".depth_dist.classA.png" 4 image 1 ;
			json_elem $work_file".depth_dist.classB.png" 4 image 1 ;
			json_elem $work_file".depth_dist.classAB.png" 4 image 1 ;

			## title	
			json_elem "Substitution types" 12 title 1 ;

			# subst type
			json_elem $work_file".subst_type.classA.png" 4 image 1 ;
			json_elem $work_file".subst_type.classB.png" 4 image 1 ;
			json_elem $work_file".subst_type.classAB.png" 4 image 0 ;

			echo -e "\t\t\t\t]" ;

			if [ $k -eq $((${#type_kind[@]} -2 )) ] && [ $l -eq $((${#type_kind[@]} -1 )) ]; then
				echo -e "\t\t\t}" ;
			else
				echo -e "\t\t\t}, " ;	
			fi

		done
	done

##########################
	# for integrated
##########################
	echo -e "\t\t], \"integrated\" : [" ;

	if [ "${#type_kind[@]}" -gt "2" ]; then
		for((k=0;k<${#type_kind[@]};k++)); do
			for((l=$k+1;l<${#type_kind[@]};l++)); do
				work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}
			
				echo -e "\t\t\t{\"$work_file\" : [" ;
				#echo -e "\t\t\t{\"${type_kind[$k]}"_vs_"${type_kind[$l]}\" : [" ;
				
				## title	
				json_elem "Genome wide methylation landscape with DEGs" 12 title 1 ;

				# circos
				json_elem $work_file".circos.png" 12 image 1 ;

				## title	
				json_elem "Mutation occurrence on DEGs" 12 title 1 ;

				# oncoprint
				json_elem $work_file".onco.png" 12 image 1 ;

				# dummy
				json_elem "dummy" 2 "dummy" 1 ;
				
				## title	
				json_elem "Visualization all analyzed data on UCSC genome browser" 12 title 1 ;

				# UCSC link
				json_elem $work_file".UCSC_hub_link.txt" 8 link 1 ;
			
				# dummy
				json_elem "dummy" 2 "dummy" 0 ;

				echo -e "\t\t\t\t]" ;

				echo -e "\t\t\t}, " ;	
			done
		done

		# for all-class
		echo -e "\t\t\t{\"ALL_class\" : [" ;	
		 	work_file="ALL"

			## title	
			json_elem "Genome wide methylation landscape with DEGs" 12 title 1 ;

			# circos
			json_elem $work_file".circos.png" 12 image 1 ;

			## title	
			json_elem "Mutation occurrence on DEGs" 12 title 1 ;

			# oncoprint
			json_elem $work_file".onco.png" 12 image 1 ;

			# mu_rate
			json_elem mu_rate.png 12 image 1 ;

      # dummy
      json_elem "dummy" 2 "dummy" 1 ;

			## title	
			json_tag_elem "Network analysis" 12 NETWORK title 1 ;

			# cyto 
			json_elem cyto.txt 8 link 1 ;
      
			# dummy
      json_elem "dummy" 2 "dummy" 1 ;
			 
			# dummy
      json_elem "dummy" 2 "dummy" 1 ;

			## title	
			json_elem "Visualization all analyzed data on UCSC genome browser" 12 title 1 ;

			# UCSC link
			json_elem $work_file".UCSC_hub_link.txt" 8 link 1 ;
 
			# dummy
      json_elem "dummy" 2 "dummy" 0 ;
		
			echo -e "\t\t\t\t]" ;
			echo -e "\t\t\t}" ;
	
	else
		 	work_file=${type_kind[0]}"_vs_"${type_kind[1]}

			echo -e "\t\t\t{\"${type_kind[$k]}"_vs_"${type_kind[$l]}\" : [" ;	

			## title	
			json_elem "Genome wide methylation landscape with DEGs" 12 title 1 ;

			# circos
			json_elem $work_file".circos.png" 12 image 1 ;

			## title	
			json_elem "Mutation occurrence on DEGs" 12 title 1 ;

			# oncoprint
			json_elem $work_file".onco.png" 12 image 1 ;

			## title	
			json_elem "Mutation rate over methylation changes" 12 title 1 ;

			# mu_rate
			json_elem mu_rate.png 12 image 1 ;

      # dummy
      json_elem "dummy" 2 "dummy" 1 ;

			## title	
			json_tag_elem "Network analysis" 12 NETWORK title 1 ;

			# cyto 
			json_elem cyto.txt 8 link 1 ;
      
			# dummy
      json_elem "dummy" 2 "dummy" 1 ;
			 
			# dummy
      json_elem "dummy" 2 "dummy" 1 ;

			## title	
			json_elem "Visualization all analyzed data on UCSC genome browser" 12 title 1 ;

			# UCSC link
			json_elem $work_file".UCSC_hub_link.txt" 8 link 1 ;
 
			# dummy
      json_elem "dummy" 2 "dummy" 0 ;
		
			echo -e "\t\t\t\t]" ;
			echo -e "\t\t\t}" ;
	fi	

	echo -e "\t\t]" ;
	echo -e "\t}]" ;

	## hallmark gene set
	#cat $HALLMARK_GENE_SET_JASON;

	echo -e "}" ;

} > $final_json

#formating the jason
python $bin_dir/formatting_data.py $result_dir

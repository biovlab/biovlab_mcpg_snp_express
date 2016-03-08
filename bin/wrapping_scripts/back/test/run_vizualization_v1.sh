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

# variables
type_kind=("lu" "baa" "bab")
####################################################################################
# 0. set parameters
####################################################################################

# TODO:set user-specific result dir to copy the result to user-specific web accessible dir


# set sample lists
mbd_list="..."

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
for (( i=0; i<${#cel_list[@]}; i++ )); do 
	# copy to web accessible directory
	# TODO :  this may be not necessary cuz already we are in the web accessible dir
	cp $cel_result_dir/`basename \${cel_list[$i]}`.bw $WEB_ACCESSIBLE_DIR
	# print UCSC url
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=test_SNP%20description=test2SNP%20visibility=dense%20bigDataUrl=$WEB_ACCESSIBLE_LOC/${cel_list[$i]}.bw"
done



############
# 6.1.2 ME
############

# TODO : set data kinds (MBD
# UCSC links for each samples : BigWig
count=0
for (( i=0; i<${#mbd_list[@]}; i++ )); do
	{
		met_level_file=$medips_from_mbd_result_dir/`basename \${mbd_list[$i]}`".met_level";

		# create bigWig file for UCSC
		awk 'NR==1{next;}{print}' $met_level_file | cut -f1,2,3,5 | bedtools intersect -a - -b $REF_HUMAN_GENOME_100_BIN_BED -sorted > $viz_result_dir/`basename \${mbd_list[$i]}`".met_level.bed" ;
		bedGraphToBigWig $viz_result_dir/`basename \${mbd_list[$i]}`".met_level.bed" $REF_HUMAN_CHR_SIZE $viz_result_dir/`basename \${mbd_list[$i]}`".met_level.bw";

		cp $viz_result_dir/`basename \${mbd_list[$i]}`".met_level.bw" $WEB_ACCESSIBLE_DIR;
		# print UCSC url
		echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=bigWig%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl=$WEB_ACCESSIBLE_LOC/`basename \${mbd_list[$i]}`.met_level.bw";
	} &
	let count+=1
	[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
done; wait

## for bisulfite samples
ls $methylkit_from_bs_result_dir/*.sorted.bw | xargs -I {} cp {} $WEB_ACCESSIBLE_DIR

# NGS plot
for (( k=0; k<${#type_kind[@]}; k++ )); do 
	echo -n "" > $viz_result_dir/ngsplot_config_ME.txt
	echo -n "" > $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt


	for (( i=0; i<${#mbd_list[@]}; i++ )); do 
		temp_filename_only=`basename \${mbd_list[$i]}`

		# create configure file for all multiple file plotting
		echo -e "$mbd_result_dir/$temp_filename_only.sorted.bam\t-1\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_all.txt
		#echo -e "$mbd_result_dir/$temp_filename_only.sorted.bam\tintron.bed\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_all.txt

		# create configure file for subtype multiple file plotting
		if [ ${MBD_TYPE_LIST[$i]} == ${type_kind[$k]} ]; then
			echo -e "$mbd_result_dir/$temp_filename_only.sorted.bam\t-1\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt
			#echo -e "$mbd_result_dir/$temp_filename_only.sorted.bam\tintron.bed\t$temp_filename_only" >> $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt
		fi
	done

	# NGS plot itself is parallelized, so do not need to be parallelzed
	#for region in "bed" ; do #genebody" "enhancer" "dhs" "cgi" "exon"; do 
	for region in "cgi" ; do #genebody" "enhancer" "dhs" "cgi" "exon"; do 

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
		/packages/test2/ngsplot-develop/bin/ngs.plot.r -SC $SC -YAS $YAS -G hg19 -R $region -C $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt -O $viz_result_dir/ME_${type_kind[$k]}_$region 

		# for manual intron region
		#Rscript /packages/test2/ngsplot-develop/bin/ngs.plot.r -G hg19 -R bed -E intron.bed -C $viz_result_dir/ngsplot_config_ME_${type_kind[$k]}.txt -O $viz_result_dir/ME_${type_kind[$k]}_$region 
	done
done

# for all sample
for region in "bed" ;do #"genebody" "enhancer" "dhs" "cgi"; do
	/packages/test2/ngsplot-develop/bin/ngs.plot.r -SC $SC -G hg19 -R $region -C $viz_result_dir/ngsplot_config_ME_all.txt -O $viz_result_dir/ME_all_$region -E /data/project/mcpg/lib/region_info/intron.bed
done






############
# 6.1.3 SNP
############
for (( i=0; i<${#mbd_list[@]}; i++ )); do
	cp $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.gz" $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.gz.tbi" $WEB_ACCESSIBLE_DIR
	# print UCSC url
	echo "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-1000000&hgct_customText=track%20type=vcfTabix%20name=ME_level%20description=ME_level%20visibility=dense%20bigDataUrl=$WEB_ACCESSIBLE_LOC/`basename \${mbd_list[$i]}`.vcf.gz"
done

i


















# 2. Region stat
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 2.2 SNP : dot plot 
#NUM_CPUS=5
#count=0 
#for (( i=0; i<${#mbd_list[@]}; i++ )); do 
#	bedtools intersect -wa -wb -a $snp_from_mbd_result_dir"/"`basename \${mbd_list[$i]}`".vcf.bed.all.cov" -b $human_refseq_various_ranges > $viz_result_dir/`basename \${mbd_list[$i]}`".SNP_bin_ranges" &
#	let count+=1
# 	[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
#done

# compute region based AVG values divide by 100
# example : chr10   95100   95200   0       chr10   95178   96678   -       NM_177987       TUBB8   promoter
# devide each range by relative 100 bins 
#awk '{interval=int(($3-$2)/100);if($4=="+"){for (i=0;i<100;i++){end=$2+interval*(i+1);start=$2+interval*i; if (end<$3){print $1"\t"start"\t"end"\t"i"\t"$7}else{print $1"\t"start"\t"$3"\t"i"\t"$7} }}else{for (i=99;i>=0;i--){end=$2+interval*(i+1);start=$2+interval*i; if (end<$3){print $1"\t"start"\t"end"\t"i"\t"$7}else{print $1"\t"start"\t"$3"\t"i"\t"$7} }} }' $human_refseq_various_ranges > $viz_result_dir/human_refseq_various_ranges.div_100

# create ID for stat matrix
echo "relative_range" > $viz_result_dir/prefix_percentage_ids
echo "5p" >> $viz_result_dir/prefix_percentage_ids
for ((i=1; i<99; i++)); do 
	echo $i"%" >> $viz_result_dir/prefix_percentage_ids
done
echo "3p" >> $viz_result_dir/prefix_percentage_ids

for range in "genebody" "cpgIsland" "exon"; do 
	# get range bins(1-100)
	grep $range $viz_result_dir/human_refseq_various_ranges.div_100 > $viz_result_dir/temp_range_div_100

	paste_list_all=()
	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		paste_list_sub_type=()
		count=0
		for (( i=0; i<${#mbd_list[@]}; i++ )); do 

???줄을 잃어버림
???줄을 잃어버림
???줄을 잃어버림
???줄을 잃어버림
???줄을 잃어버림
	python $bin_dir/profile_ME_GE.py $viz_result_dir/GDMR_GDEG_MGD_ME_GE_CORR_all.txt $viz_result_dir"/GDMR_GDEG_MGD_ME_GE_stat.txt" $sample_num $corr_threshold $viz_result_dir/GDMR_GDEG_ME_GE_MAT_THRESH".txt" $viz_result_dir/GDMR_GDEG_ME_GE_MAT_PNC".txt"

	# draw plot
	echo "[INFO] Draw bar plots"
	Rscript $bin_dir/stat_bar.r $viz_result_dir/GDMR_GDEG_ME_GE_MAT_PNC".txt" $viz_result_dir/BAR_GDMR_GDEG_ME_GE_MAT_PNC".jpg" "Among GDMR GDEG interseced bins over threshold, PC/NC ratio with GE in all subtypes" 2
	Rscript $bin_dir/stat_bar.r $viz_result_dir/GDMR_GDEG_ME_GE_MAT_THRESH".txt" $viz_result_dir/BAR_GDMR_GDEG_ME_GE_MAT_THRESH".jpg" "Ratio of GDMR GDEG intersected bins having CORR with GE over threshold in all subtypes" 1

	# copy to WEB_ACCESSIBLE_DIR
	cp $viz_result_dir/GDMR_GDEG_MGD_ME_GE_stat.txt $viz_result_dir/BAR_GDMR_GDEG_ME_GE_MAT_PNC".jpg" $viz_result_dir/BAR_GDMR_GDEG_ME_GE_MAT_THRESH".jpg" $WEB_ACCESSIBLE_DIR

# 3.2.10 DSNP & DEG , GDSNP & GDEG : Bar plot
	# DSNP & DEG
	for (( i=0; i<${#type_kind[@]}; i++ )); do 
		for (( j=$i+1; j<${#type_kind[@]}; j++ )); do 
			TYPE_PAIR=${type_kind[$i]}_${type_kind[$j]}
			TMP_MGD_LEN=$(( ${TYPE_LEN[$i]} + ${TYPE_LEN[$j]} ))

			DSNP_DEG_CORR=$PROF_DIR/DSNP_DEG_CORR_$TYPE_PAIR".txt"

			echo "[INFO] Compute % of correlation over threshold"
			python $bin_dir/profile_AF_GE.py $DSNP_DEG_CORR $viz_result_dir"/DSNP_DEG_MGD_AF_GE_stat_"$TYPE_PAIR".txt" $TMP_MGD_LEN $corr_threshold $viz_result_dir/DSNP_DEG_AF_GE_MAT_THRESH_$TYPE_PAIR".txt" $viz_result_dir/DSNP_DEG_AF_GE_MAT_PNC_$TYPE_PAIR".txt"

			# draw plot
			echo "[INFO] Draw bar plots"
			Rscript $bin_dir/stat_bar.r $viz_result_dir/DSNP_DEG_AF_GE_MAT_PNC_$TYPE_PAIR".txt" $viz_result_dir/BAR_DSNP_DEG_AF_GE_MAT_PNC_$TYPE_PAIR".jpg" "Among DSNP DEG interseced bins over threshold, PC/NC ratio with GE in $TYPE_PAIR" 2
			Rscript $bin_dir/stat_bar.r $viz_result_dir/DSNP_DEG_AF_GE_MAT_THRESH_$TYPE_PAIR".txt" $viz_result_dir/BAR_DSNP_DEG_AF_GE_MAT_THRESH_$TYPE_PAIR".jpg" "Ratio of DSNP DEG intersected bins having CORR between AF and GE over threshold in $TYPE_PAIR" 1

			# copy to WEB_ACCESSIBLE_DIR
			cp $viz_result_dir/DSNP_DEG_MGD_AF_GE_stat_$TYPE_PAIR".txt" $viz_result_dir/BAR_DSNP_DEG_AF_GE_MAT_THRESH_$TYPE_PAIR".jpg" $viz_result_dir/BAR_DSNP_DEG_AF_GE_MAT_PNC_$TYPE_PAIR".jpg" $WEB_ACCESSIBLE_DIR
		done
	done

	# GDSNP & GDEG 
	cut -f5- $PROF_DIR/GDEG_GDSNP.txt > $viz_result_dir/GDEG_GDSNP_CORR.txt
	GDEG_GDSNP_CORR=$viz_result_dir/GDEG_GDSNP_CORR.txt

	echo "[INFO] Compute % of correlation over threshold"
	python $bin_dir/profile_AF_GE.py $GDEG_GDSNP_CORR $viz_result_dir"/GDEG_GDSNP_MGD_AF_GE_stat.txt" $sample_num $corr_threshold $viz_result_dir/GDEG_GDSNP_AF_GE_MAT_THRESH.txt $viz_result_dir/GDEG_GDSNP_AF_GE_MAT_PNC.txt

	# draw plot
	echo "[INFO] Draw bar plots"
	Rscript $bin_dir/stat_bar.r $viz_result_dir/GDEG_GDSNP_AF_GE_MAT_PNC.txt $viz_result_dir/BAR_GDEG_GDSNP_AF_GE_MAT_PNC.jpg "Among GDSNP GDEG interseced bins over threshold, PC/NC ratio with GE in all subtypes" 2
	Rscript $bin_dir/stat_bar.r $viz_result_dir/GDEG_GDSNP_AF_GE_MAT_THRESH.txt $viz_result_dir/BAR_GDEG_GDSNP_AF_GE_MAT_THRESH.jpg "Ratio of GDSNP GDEG intersected bins having CORR between AF and GE over threshold in all subtypes" 1

	# copy to WEB_ACCESSIBLE_DIR
	cp $viz_result_dir/GDEG_GDSNP_MGD_AF_GE_stat.txt $viz_result_dir/BAR_GDEG_GDSNP_AF_GE_MAT_PNC.jpg $viz_result_dir/BAR_GDEG_GDSNP_AF_GE_MAT_THRESH.jpg $WEB_ACCESSIBLE_DIR
COMMENT

# create index page
sh $bin_dir/gnu-mirror-index-creator.sh $WEB_ACCESSIBLE_DIR /mcpg/ "147.46.15.115"
#################################
# end
#################################
<

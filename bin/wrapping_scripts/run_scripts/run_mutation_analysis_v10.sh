#!/bin/bash
source `dirname $0`/../../env.sh

# parse experiment information
source `dirname $0`/parse_exp_info.sh


NUM_CPUS=$mu_cpu

# directories
bin_dir="$WORK_DIR/bin"
lib_dir="$WORK_DIR/lib"
<<'COMMENT'
result_dir="$WORK_DIR/result"
COMMENT
rnaseq_result_dir=$result_dir"/mu/rna_seq"
dnaseq_result_dir=$result_dir"/mu/dna_seq"
array_result_dir=$result_dir"/mu/array"
genomeDir=$rnaseq_result_dir"/STAR_hg19_index/"
genomeDir_2pass=$rnaseq_result_dir"/2pass_genome_index"
runDir_2pass=$rnaseq_result_dir"/2pass"

final_result_dir=$final_result_root_dir/mutation
mkdir -p $final_result_dir

# variables
#genome_hg19=$genomeDir"/hg19.fa" -> REF_HUMAN
#hg19_gtf=$WORK_DIR"/lib/human_hg19_ucsc.gtf" -> REF_HUMAN_GTF
#genome_hg19_dict=$genomeDir"/hg19.dict" -> REF_HUMAN_DICT
merged_SJ=$rnaseq_result_dir"/merged_SJ.out.tab"
window_size=100

####################################################################################
# 0. parse input
####################################################################################

# type kind
mu_type_kind=($(echo ${mu_type_list[@]} | tr ' ' '\n' | sort | uniq | tr '\n' ' '))


# assign parsed values to each variable
cel_list_temp=(${mu_array_list[@]})
type_kind=(${mu_type_kind[@]})
type_list=(${mu_type_list[@]})
sample_fq1_list_temp=(${mu_sample_fq1_list[@]})
sample_fq2_list_temp=(${mu_sample_fq2_list[@]})
is_array=$mu_is_array
single_pair=$mu_single_pair
data_type=$mu_sample_format					# 0:DNAseq, 1:RNAseq, 2:HumanHap


# preprocess input / unzip
unzip=1
source `dirname $0`/preprocess_mu_input.sh

# read length
read_length_minus_1=100 #In case of reads of varying length, the ideal value is max(ReadLength)-1. In most cases, a generic value of 100 will work as well as the ideal value.

# directories
vcf_result_dir=""


# functions
function my_join { local IFS="$1"; shift; echo "$*"; }


file_extension=$(echo `basename \${sample_fq1_list[0]}` | awk -F . '{if (NF>1) {print $NF}}')



<<'COMMENT_TEMP'
#######################################################################################
# 1. vcf analysis
########################################################################################

# generate visual file
for(( i=0; i<${#vcf_list[@]}; i++ )); do
	vcf_filename_only=`basename \${vcf_list[\$i]}`
	output_prefix=`basename \$vcf_filename_only "."$file_extension`"_"
	
	# get snp count per window

	sh $bin_dir"/get_snp_count_per_window.sh" ${vcf_list[$i]} $vcf_result_dir $lib_dir/genome_wide_bed_window

	# compression
  $tabix_dir/bgzip -c ${vcf_list[$i]} > $vcf_result_dir"/"$vcf_filename_only".gz"

  # indexing
  $tabix_dir/tabix -p vcf $vcf_result_dir"/"$vcf_filename_only".gz"

	redirect_SNP_html=$vcf_result_dir"/"$output_prefix"SNP.html"

  # create redirect html
  echo -ne '<html><head><meta http-equiv="refresh" content="0; url='> $redirect_SNP_html
  echo -ne "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-10000000&hgct_customText=track%20type=vcfTabix%20name="$output_prefix"SNP%20description="$output_prefix"SNP%20visibility=dense%20bigDataUrl="$WEB_ACCESSIBLE_LOC""$vcf_filename_only".gz" >> $redirect_SNP_html
  echo -e '" /></head><body></body></html>' >> $redirect_SNP_html

  # copy to web_accessible_url
  cp $vcf_result_dir"/"$vcf_filename_only".gz.tbi" $vcf_result_dir"/"$vcf_filename_only".gz" $redirect_SNP_html $WEB_ACCESSIBLE_DIR
done	


# MERGE VCF FILES
for (( k=0; k<${#type_kind[@]}; k++)); do
	first_sample_vcf_list=()
	
	for (( i=0; i<${#vcf_list[@]}; i++)); do
		vcf_filename_only=`basename \${vcf_list[\$i]}`
		if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
			first_sample_vcf_list+=($vcf_result_dir"/"$vcf_filename_only".gz")
		fi
	done
	
	# merge vcf files by subtype

	first_sample_vcf_string=`my_join " " ${first_sample_vcf_list[@]}`
		
	$bcftools_dir/bcftools merge $first_sample_vcf_string -O z > $vcf_result_dir"/"${type_kind[$k]}".vcf.MGD.gz"
done

# SUBTYPE COMPARISON

for (( k=0; k<${#type_kind[@]}; k++ )); do
	for (( l=$k+1; l<${#type_kind[@]}; l++ )); do
		echo "[INFO] Currently processing ${type_kind[$k]} and ${type_kind[$l]}"
		work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}

		first_vcf_file=$vcf_result_dir"/"${type_kind[$k]}".vcf.MGD.gz"
		second_vcf_file=$vcf_result_dir"/"${type_kind[$l]}".vcf.MGD.gz"
		
    # compare two vcf files. *.log contains stat which is good for drawing venn diagram
    $vcftools_dir/vcftools --vcf $first_vcf_file --diff $second_vcf_file --out $vcf_result_dir"/"$work_file".vcf_diff"
    # or by vcf-compare
      # $vcftools_dir/vcf-compare $first_vcf_file".gz" $second_vcf_file".gz" > $runDir_2pass"/"$pair_prefix".vcf_compare"

    # stat by bcftools
    $bcftools_dir/bcftools stats $first_vcf_file $second_vcf_file >  $vcf_result_dir"/"$work_file".bcf_stats"

    # plot
    $bcftools_dir/plot-vcfstats $vcf_result_dir"/"$work_file".bcf_stats" -p $vcf_result_dir"/bcftools_compare_plots/"$work_file

    # cp
    cp $vcf_result_dir"/bcftools_compare_plots/"$work_file"-summary.pdf" $WEB_ACCESSIBLE_DIR
	done
done



#############################################################################################
#############################################################################################

COMMENT_TEMP



# check the input is array(human hap map), RNAseq, or DNA seq
#data_type=1

####################################################################################
# 3. Mutation analysis
####################################################################################


#####################
# 3.1 array : Hapmap
#####################
if [ "$data_type" -eq "2" ] ; then
	mkdir -p $array_result_dir
# TODO



#####################
# 3.2 RNAseq
#####################
elif [ "$data_type" -eq "1" ] ; then
	mkdir -p $rnaseq_result_dir
	clipped_dir=$rnaseq_result_dir/clipped/
	mkdir -p $clipped_dir

# clipping adapter and fastqc
	count=0
	NUM_CPUS=$(( $SYS_NUM_CPUS / 3 ))
	for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do
		if [ $single_pair == "single" ]; then
  		trim_galore --fastqc -o $clipped_dir ${sample_fq1_list[$i]} &
		else
    	trim_galore --paired --fastqc -o $clipped_dir ${sample_fq1_list[$i]} ${sample_fq2_list[$i]} &
		fi
		let count+=1
		[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
  done; wait
	NUM_CPUS=$mu_cpu

<<'COMMENT'
# generate STAR genome index
mkdir $genomeDir
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeDir/hg19.fa  --runThreadN $NUM_CPUS -sjdbGTFfile $REF_HUMAN_GTF

# genome prepare for gatk
$java_dir/java -jar $picard_dir/CreateSequenceDictionary.jar R= $REF_HUMAN O= $REF_HUMAN_DICT
COMMENT


# first pass alignment
echo -n "" > $merged_SJ


cd $rnaseq_result_dir
if [ $single_pair == "single" ]; then
  for (( i=0; i<${#sample_fq1_list[@]}; i++)); do
    sample_fq_filename_only=`basename \${sample_fq1_list[\$i]}`
    output_prefix=$rnaseq_result_dir"/"`basename $sample_fq_filename_only "."$file_extension`"_"
		
		STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $clipped_dir/`basename $sample_fq_filename_only "."$file_extension`"_trimmed.fq" --runThreadN $NUM_CPUS --outFileNamePrefix $output_prefix --alignEndsType EndToEnd
		
		# merge SJ
  	awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' $output_prefix"SJ.out.tab" >> $merged_SJ
		
	done
else
	
	for (( i=0; i<${#sample_fq1_list[@]}; i++)); do
		sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
    sample_fq2_filename_only=`basename \${sample_fq2_list[\$i]}`
    output_prefix=$rnaseq_result_dir"/"`basename $sample_fq1_filename_only "."$file_extension`"_"

    STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $clipped_dir/`basename $sample_fq1_filename_only "."$file_extension`"_val_1.fq" $clipped_dir/`basename $sample_fq2_filename_only "."$file_extension`"_val_2.fq" --runThreadN $NUM_CPUS --outFileNamePrefix $output_prefix --alignEndsType EndToEnd
 		
		# merge SJ
 	 	awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' $output_prefix"SJ.out.tab" >> $merged_SJ
		
	done
fi

cd $WORK_DIR

<<'COMMENT'
for (( i=0; i<${#clipped_fq_list[@]}; i++ )); do

  # variables
  output_prefix=${sample_prefix_list[$i]}"_"
  sample_fq=${clipped_fq_list[$i]}

  # first pass
  cd $rnaseq_result_dir
  STAR --genomeDir $genomeDir --readFilesIn $sample_fq --runThreadN $NUM_CPUS --outFileNamePrefix $output_prefix  --alignEndsType EndToEnd

  # merge SJ
  awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' $output_prefix"SJ.out.tab" >> $merged_SJ
  cd $WORK_DIR
done
COMMENT

# generate 2pass genome
mkdir -p $genomeDir_2pass
STAR --runMode genomeGenerate --genomeDir $genomeDir_2pass --genomeFastaFiles $REF_HUMAN --sjdbFileChrStartEnd $merged_SJ --sjdbOverhang $read_length_minus_1 --runThreadN $NUM_CPUS

echo "[INFO] Star 2nd align"

mkdir -p $runDir_2pass

cd $runDir_2pass

if [ $single_pair == "single" ]; then
  for (( i=0; i<${#sample_fq1_list[@]}; i++)); do
    sample_fq_filename_only=`basename \${sample_fq1_list[\$i]}`
    output_prefix=$runDir_2pass"/"`basename $sample_fq_filename_only "."$file_extension`"_"
		
		STAR --genomeDir $genomeDir_2pass --readFilesIn $clipped_dir/`basename $sample_fq_filename_only "."$file_extension`"_trimmed.fq" --runThreadN $NUM_CPUS --outFileNamePrefix $output_prefix --alignEndsType EndToEnd
		
	done
else
	for (( i=0; i<${#sample_fq1_list[@]}; i++)); do
		sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
    sample_fq2_filename_only=`basename \${sample_fq2_list[\$i]}`
    output_prefix=$runDir_2pass"/"`basename $sample_fq1_filename_only "."$file_extension`"_"

    STAR --genomeDir $genomeDir_2pass --readFilesIn $clipped_dir/`basename $sample_fq1_filename_only "."$file_extension`"_val_1.fq" $clipped_dir/`basename $sample_fq2_filename_only "."$file_extension`"_val_2.fq" --runThreadN $NUM_CPUS --outFileNamePrefix $output_prefix --alignEndsType EndToEnd
		
	done
fi

cd $WORK_DIR

mut_mapping_dir=$runDir_2pass

<<'COMMENT'
for (( i=0; i<${#clipped_fq_list[@]}; i++ )); do
  # variables
  output_prefix=${sample_prefix_list[$i]}"_"
  sample_fq=${clipped_fq_list[$i]}

  # second align
  mkdir $runDir_2pass
  cd $runDir_2pass
  STAR --genomeDir $genomeDir_2pass --readFilesIn $sample_fq --runThreadN $NUM_CPUS --outFileNamePrefix $output_prefix --alignEndsType EndToEnd
  cd $WORK_DIR
done
COMMENT



#####################
# 3.3 DNAseq
#####################
elif [ "$data_type" -eq "0" ] ; then
	mkdir -p $dnaseq_result_dir
	echo "[INFO] Start DNAseq analysis"
	
	#count=0
	for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do
		sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
		output_prefix=`basename $sample_fq1_filename_only "."$file_extension`"_"
		
		output_sam=$dnaseq_result_dir"/"$output_prefix"Aligned.out.sam"
	
		bash $bin_dir/bwa.sh ${sample_fq1_list[$i]} ${sample_fq2_list[$i]} $single_pair $dnaseq_result_dir $output_sam $NUM_CPUS
	
		#let count+=1
		#[[ $((count%$NUM_CPUS)) -eq 0 ]] && wait
	done; wait
	
	mut_mapping_dir=$dnaseq_result_dir
fi

########################
# Genearate VCF from mapping results
########################
num_cpus_for_gatk=$NUM_CPUS
mem_for_gatk_float=$(( $SYS_MAX_MEM / $num_cpus_for_gatk ))
mem_for_gatk=${mem_for_gatk_float%.*}
count=0
#echo "[INFO] mem_for_gatk="$mem_for_gatk
#mem_for_gatk="4"

# temporary test
#num_cpus_for_gatk=6
#mem_for_gatk="10"

if [ "$data_type" -eq "0" ] || [ "$data_type" -eq "1" ] ; then	

	for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do
	{	
		sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
	  output_prefix=`basename $sample_fq1_filename_only "."$file_extension`"_"
	
	 	second_star_out_sam=$mut_mapping_dir"/"$output_prefix"Aligned.out.sam"
	  second_star_out_bam=$mut_mapping_dir"/"$output_prefix"Aligned.out.bam"
		second_star_out_bam_sorted=$mut_mapping_dir"/"$output_prefix"Aligned.out.sorted.bam"
	
		split_bam=$mut_mapping_dir"/"$output_prefix"Aligned.out.split.bam"
	  second_star_out_dedupped_bam=$mut_mapping_dir"/"$output_prefix"Aligned.out.dedupped.bam"
	
		vcf_out=$mut_mapping_dir"/"$output_prefix"Aligned.out.split.vcf"

	  filtered_vcf_out=$mut_mapping_dir"/"$output_prefix"Aligned.out.split.filtered.vcf"
	  filtered_vcf_out_filename_only=$output_prefix"Aligned.out.split.filtered.vcf"
	  redirect_SNP_html=$mut_mapping_dir/$output_prefix"SNP.html"

	
		echo "[INFO] Add read groups, sort, mark duplicates, and create index" # Add read groups, sort, mark duplicates, and create index
	  java -jar -Xmx"$mem_for_gatk"g $picard_dir/AddOrReplaceReadGroups.jar I=$second_star_out_sam O=$second_star_out_bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

	#	echo "[INFO] Sort bam"
	#	java -Xmx20g -jar $picard_dir/SortSam.jar I=$second_star_out_bam O=$second_star_out_bam_sorted SO=coordinate

	#	echo "[INFO] Indexing"
	#	java -Xmx20g -jar $picard_dir/BuildBamIndex.jar I=$second_star_out_bam_sorted

	  echo "[INFO] Mark duplicates"
 	 java -jar -Xmx"$mem_for_gatk"g $picard_dir/MarkDuplicates.jar TMP_DIR=/tmp I=$second_star_out_bam O=$second_star_out_dedupped_bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics
	
	  echo "[INFO] Split'N'Trim and reassign mapping qualities" # Split'N'Trim and reassign mapping qualities
	  java -jar -Xmx"$mem_for_gatk"g $gatk_dir/GenomeAnalysisTK.jar -T SplitNCigarReads -R $REF_HUMAN -I $second_star_out_dedupped_bam -o $split_bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS # --fix_misencoded_quality_scores

	  echo "[INFO] Variant calling" #HaplotypeCaller
	  java -jar -Xmx"$mem_for_gatk"g $gatk_dir/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 4 -R $REF_HUMAN -I $split_bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $vcf_out

	  echo "[INFO] Variant filtering" # variant filtering
  # this filter causes errors on vcf comparing step (additional GD FS filter header lines in vcf file cuase the problem)
    #$java_dir/java -jar $gatk_dir/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_HUMAN -V $vcf_out -window 35 -cluster 3 -filterName FS -filter \""FS > 30.0\"" -filterName QD -filter \""QD < 2.0\"" -o $filtered_vcf_out
 	 java -jar -Xmx"$mem_for_gatk"g $gatk_dir/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_HUMAN -V $vcf_out -window 35 -cluster 3 -o $filtered_vcf_out

#---------------------------------------------------------------------------------------------------
<<'COMMENT'
		# get snp count per window
	
		bash $bin_dir"/get_snp_count_per_window.sh" $filtered_vcf_out $mut_mapping_dir $lib_dir/genome_wide_bed_window
COMMENT


		# convert vcf to bed
		awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+",$8}}' $filtered_vcf_out > $filtered_vcf_out".bed"

	  # compression
	  $tabix_dir/bgzip -c $filtered_vcf_out > $filtered_vcf_out".gz"

	  # indexing
	  $tabix_dir/tabix -p vcf $filtered_vcf_out".gz"

	  # create redirect html
	  echo -ne '<html><head><meta http-equiv="refresh" content="0; url='> $redirect_SNP_html
	  echo -ne "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:1-10000000&hgct_customText=track%20type=vcfTabix%20name="$output_prefix"SNP%20description="$output_prefix"SNP%20visibility=dense%20bigDataUrl="$final_result_url""$filtered_vcf_out_filename_only".gz" >> $redirect_SNP_html
	  echo -e '" /></head><body></body></html>' >> $redirect_SNP_html
	  # copy to web_accessible_url
	  cp $filtered_vcf_out".gz.tbi" $filtered_vcf_out".gz" $redirect_SNP_html $final_result_dir

	#---------------------------------------------------------------------------------------------------
	} &
	#TODO Check number of concurrent process
		let count+=1
		[[ $((count%$num_cpus_for_gatk)) -eq 0 ]] && wait
	done; wait

	echo "[INFO] Merge VCF file according to class"
	# MERGE VCF FILES
	for (( k=0; k<${#type_kind[@]}; k++)); do
		first_sample_vcf_list=()
	
	#	for (( i=0; i<${#vcf_list[@]}; i++)); do
		for (( i=0; i<${#sample_fq1_list[@]}; i++)); do
			sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
			output_prefix=`basename $sample_fq1_filename_only "."$file_extension`"_"
		
			vcf_file=$mut_mapping_dir"/"$output_prefix"Aligned.out.split.filtered.vcf.gz"
	
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				first_sample_vcf_list+=($vcf_file)
			fi
		done
	
		# merge vcf files by subtype

		first_sample_vcf_string=`my_join " " ${first_sample_vcf_list[@]}`
	
		$bcftools_dir/bcftools merge $first_sample_vcf_string -O z --force-samples > $mut_mapping_dir"/"${type_kind[$k]}".vcf.gz"
	
		# indexing
		$tabix_dir/tabix -p vcf $mut_mapping_dir"/"${type_kind[$k]}".vcf.gz"
	done


	# SUBTYPE COMPARISON

	for (( k=0; k<${#type_kind[@]}; k++ )); do
		for (( l=$k+1; l<${#type_kind[@]}; l++ )); do
			echo "[INFO] Currently processing ${type_kind[$k]} and ${type_kind[$l]}"
			work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}

			first_vcf_file=$mut_mapping_dir"/"${type_kind[$k]}".vcf.gz"
			second_vcf_file=$mut_mapping_dir"/"${type_kind[$l]}".vcf.gz"
	

<<'COMMENT' 
			echo "[INFO] Compare two vcf files"	
	    # compare two vcf files. *.log contains stat which is good for drawing venn diagram
  		# no diff options. Need to specific options  
			$vcftools_dir/vcftools --gzvcf $first_vcf_file --gzdiff $second_vcf_file --out $mut_mapping_dir"/"$work_file".vcf_diff"
   	 # or by vcf-compare
      # $vcftools_dir/vcf-compare $first_vcf_file".gz" $second_vcf_file".gz" > $runDir_2pass"/"$pair_prefix".vcf_compare"
COMMENT

			echo "[INFO] vcf stats"
   	 # stat by bcftools
  	  $bcftools_dir/bcftools stats $first_vcf_file $second_vcf_file >  $mut_mapping_dir"/"$work_file".bcf_stats"

			echo "[INFO] Drawing plot"
    # plot
  	  $bcftools_dir/plot-vcfstats $mut_mapping_dir"/"$work_file".bcf_stats" -p $mut_mapping_dir"/bcftools_compare_plots/"$work_file

    # cp
#  	  cp $mut_mapping_dir"/bcftools_compare_plots/"$work_file"*summary.pdf" $final_result_dir
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-venn_bars.snps.png" $final_result_dir/$work_file".snp.tot.count.png"
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-venn_bars.indels.png" $final_result_dir/$work_file".indel.tot.count.png"

			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-overlap_by_af.snps.png" $final_result_dir/$work_file".snp.overlap.AF.png"
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-overlap_by_af.indels.png"	$final_result_dir/$work_file".indel.overlap.AF.png"

			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-tstv_by_qual.0.png" $final_result_dir/$work_file".ts_tv.qual.classA.png"
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-tstv_by_qual.1.png" $final_result_dir/$work_file".ts_tv.qual.classB.png"
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-tstv_by_qual.2.png" $final_result_dir/$work_file".ts_tv.qual.classAB.png"

			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-indels.0.png" $final_result_dir/$work_file".indel_dist.classA.png"
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-indels.1.png" $final_result_dir/$work_file".indel_dist.classB.png"
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-indels.2.png" $final_result_dir/$work_file".indel_dist.classAB.png"
		
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-depth.0.png" $final_result_dir/$work_file".depth_dist.classA.png"
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-depth.1.png" $final_result_dir/$work_file".depth_dist.classB.png"
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-depth.2.png" $final_result_dir/$work_file".depth_dist.classAB.png"
			
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-substitutions.0.png" $final_result_dir/$work_file".subst_type.classA.png"
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-substitutions.1.png" $final_result_dir/$work_file".subst_type.classB.png"
			mv $mut_mapping_dir"/bcftools_compare_plots/"$work_file"-substitutions.2.png" $final_result_dir/$work_file".subst_type.classAB.png"
		done
	done
fi


<<'COMMENT'
# comparing

# comparing two vcf files.
echo "[INFO] compairng vcf files"
for (( i=0; i<${#clipped_fq_list[@]}-1; i++ )); do
  for (( j=i+1; j<${#clipped_fq_list[@]}; j++ )); do
    # variables
    first_output_prefix=${sample_prefix_list[$i]}"_"
    second_output_prefix=${sample_prefix_list[$j]}"_"
    pair_prefix=${sample_prefix_list[$i]}"_"${sample_prefix_list[$j]}
    first_vcf_file=$runDir_2pass"/"$first_output_prefix"Aligned.out.split.filtered.vcf"
    second_vcf_file=$runDir_2pass"/"$second_output_prefix"Aligned.out.split.filtered.vcf"

    # compare two vcf files. *.log contains stat which is good for drawing venn diagram
    $vcftools_dir/vcftools --vcf $first_vcf_file --diff $second_vcf_file --out $runDir_2pass"/"$pair_prefix".vcf_diff"
    # or by vcf-compare
      # $vcftools_dir/vcf-compare $first_vcf_file".gz" $second_vcf_file".gz" > $runDir_2pass"/"$pair_prefix".vcf_compare"

    # stat by bcftools
    $bcftools_dir/bcftools stats $first_vcf_file".gz" $second_vcf_file".gz" >  $runDir_2pass"/"$pair_prefix".bcf_stats"

    # plot
    $bcftools_dir/plot-vcfstats $runDir_2pass"/"$pair_prefix".bcf_stats" -p $runDir_2pass"/bcftools_compare_plots/"$pair_prefix

    # cp
    cp $runDir_2pass"/bcftools_compare_plots/"$pair_prefix"-summary.pdf" $WEB_ACCESSIBLE_DIR

  done
done
COMMENT


<<'COMMENT'
#####################
# 3.3 DNAseq
#####################
elif [ "$data_type" -eq "0" ] ; then
	mkdir -p $dnaseq_result_dir
	echo "[INFO] Start DNAseq analysis"
	
	for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do
		bash $bin_dir/mbd_ver2.sh ${sample_fq1_list[$i]} ${sample_fq2_list[$i]} $single_pair $dnaseq_result_dir		# align by STAR, so no need to be parallelization.
	done
	
	# run on multi cpu
	# NOTE : Require 1. bam file list, 2. index of bam file "samtools index bamfile"
	# NOTE : !!!! IF bam file list is used in mpileup, merged result is different at the result flag columns compare with one with run seperately. Only use for replicate. !!!

	# get SNP from MBDseq
	for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do 
		# get snp from mbd
		echo "[INFO] Get SNP from DNA-seq : "${sample_fq1_list[$i]}
		temp_filename_only=`basename \${sample_fq1_list[$i]}`
		filename_wo_ext=$(echo ${sample_fq1_list[$i]} | cut -d '.'  --complement -f2-)

		bash $bin_dir/snp_from_mbd.sh $dnaseq_result_dir"/"$temp_filename_only".sorted.bam" $dnaseq_result_dir $temp_filename_only

		# visualize by UCSC genome browser
		bgzip -cf $dnaseq_result_dir"/"$temp_filename_only".vcf" > $dnaseq_result_dir"/"$temp_filename_only".vcf.gz"
		tabix -p vcf $dnaseq_result_dir"/"$temp_filename_only".vcf.gz"

		# move to .vcf.gz .vcf.gz.tbi to web accessible loc
		cp $dnaseq_result_dir"/"$temp_filename_only".vcf.gz" $dnaseq_result_dir"/"$temp_filename_only".vcf.gz.tbi" $WEB_ACCESSIBLE_DIR
	
		# get snp count per genome wide window
		if [ $i -eq 0 ] ; then
			bash $bin_dir/create_per_chr_genome_window.sh $dnaseq_result_dir $window_size $dnaseq_result_dir/genome_wide_bed_window
		fi

		echo "[INFO] Get count per window  : "${sample_fq1_list[$i]}
		bash $bin_dir/get_snp_count_per_window.sh $dnaseq_result_dir"/"$temp_filename_only".vcf" $dnaseq_result_dir $window_size $dnaseq_result_dir/genome_wide_bed_window 

		# convert bed to bigbed to view on UCSC genome browser
		#sort -k1,1 -k2,2n $snp_from_mbd_result_dir"/"$temp_filename_only".vcf.bed.all.cov" > $snp_from_mbd_result_dir"/"$temp_filename_only".vcf.sorted.bed.all"
		#bedToBigBed $snp_from_mbd_result_dir"/"$temp_filename_only".vcf.sorted.bed.all" $REF_HUMAN_CHR_SIZE $snp_from_mbd_result_dir"/"$temp_filename_only".vcf.all.cov.bb"
	done
fi

COMMENT

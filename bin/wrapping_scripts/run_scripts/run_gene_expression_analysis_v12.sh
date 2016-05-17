#!/bin/bash
source `dirname $0`/../../env.sh

# parse experiment information
source `dirname $0`/parse_exp_info.sh

NUM_CPUS=$ge_cpu

# variables
bin_dir="$WORK_DIR/bin"
lib_dir="$WORK_DIR/lib"
#result_dir="$WORK_DIR/result"
#result_dir="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/rna_paired_test_data/"

# KEGG color code based on log2 foldchange
red4="#ff6666" # 5<x : red4
red3="#ff8888" # 5>x>2 : red3
red2="#ffaaaa" # 2>x>1 : red2
red1="#ffcccc" # 1>x>0 : red1
white="#ffffff" # x=0 : white
blue4="#6666ff" # -5>x : blue4
blue3="#8888ff" # -5<x<-2 : blue3
blue2="#aaaaff" # -2<x<-1 : blue2

# directories
rnaseq_result_dir=$result_dir"/gene_exp/rna_seq"
clipped_dir=$rnaseq_result_dir"/cliiped/"
deg_analysis_dir=$rnaseq_result_dir"/deg/"
cel_result_dir=$result_dir"/gene_exp/cel"
final_result_dir=$final_result_root_dir/gene_expression

mkdir -p $final_result_dir
#genomeDir=$WORK_DIR"/STAR_hg19_index/"

# functions
function my_join { local IFS="$1"; shift; echo "$*"; }

####################################################################################
# 0. parse input
####################################################################################

#echo ${exp_info[$i,$SAMPLE_TYPE]}

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

# preprcess input data / unzip
unzip=1
source `dirname $0`/preprocess_ge_input.sh

type_len=($(echo ${ge_type_list[@]} | tr ' ' '\n' | sort | uniq -c | awk '{print $1}'| tr '\n' ' '))


#######################
# set input, class, type lists
#######################


<<'COMMENT'
# input list
cell_test_dir="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/icbp/gene_expression_array/"
cel_list=($cell_test_dir"100730_s_1_export.txt.CEL" $cell_test_dir"100730_s_2_export.txt.CEL" $cell_test_dir"100730_s_4_export.txt.CEL" $cell_test_dir"100730_s_5_export.txt.CEL")


#sample_fq1_list
#sample_fq2_list
#rnsseq_test_dir="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/rna_single_test_data/"
#sample_fq1_list=($rnsseq_test_dir"MOM-CR_head_100000.fq" $rnsseq_test_dir"MOM-DX_head_100000.fq" $rnsseq_test_dir"SON-CR_head_100000.fq" $rnsseq_test_dir"SON-DX_head_100000.fq")
#sample_fq2_list=("NA" "NA" "NA" "NA")
rnaseq_test_dir="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/rna_paired_test_data/"
sample_fq1_list=($rnaseq_test_dir"C2_head_100000_1.fq" $rnaseq_test_dir"C3_head_100000_1.fq" $rnaseq_test_dir"Li1_head_100000_1.fq" $rnaseq_test_dir"Li5_head_100000_1.fq" $rnaseq_test_dir"Th1_head_100000_1.fq" $rnaseq_test_dir"Th5_head_100000_1.fq")
sample_fq2_list=($rnaseq_test_dir"C2_head_100000_2.fq" $rnaseq_test_dir"C3_head_100000_2.fq" $rnaseq_test_dir"Li1_head_100000_2.fq" $rnaseq_test_dir"Li5_head_100000_2.fq" $rnaseq_test_dir"Th1_head_100000_2.fq" $rnaseq_test_dir"Th5_head_100000_2.fq")



#adapter_seq_list=("AGATCGGAAGAGC" "AGATCGGAAGAGC" "AGATCGGAAGAGC" "AGATCGGAAGAGC" "AGATCGGAAGAGC" "AGATCGGAAGAGC")


# type list
#type_kind=("MOM_CR" "MOM_DX" "SON_CR" "SON_DX")
#type_list=("MOM_CR" "MOM_DX" "SON_CR" "SON_DX")
type_kind=("C" "Li" "Th")
type_list=("C" "C" "Li" "Li" "Th" "Th")
COMMENT

# fold change
fold_change=2

####################################################################################
# 2. Gene Expression
####################################################################################

if [ "$is_array" -eq "1" ] ; then
# 2.1 Microarray (.CELL)
	# limma phenotype group comparison. limma doesn't support 2 sample comparison unless if there are biological replicates -> just extract expression value and use fold changeAA
	# in case there is only one sample, limma.R only extract expression levels from input samples

	# create result directory
	mkdir -p $cel_result_dir
	
#	if [ $replicate_exist -eq 1 ]; then	
	if [ ${#type_kind[@]} -lt ${#type_list[@]} ]; then
	# FOR WITH REPLICATES
		for (( k=0; k<${#type_kind[@]}; k++ )); do 
			work_file=""
			for (( l=$k+1; l<${#type_kind[@]}; l++ )); do
				echo "[INFO] Currently processing ${type_kind[$k]} and ${type_kind[$l]}"

				work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}

				# init list
				first_sample_class_name=()
				second_sample_class_name=()
				first_sample_class=()
				second_sample_class=()

				# create sample list per subtype
				for (( i=0; i<${#cel_list[@]}; i++ )); do
					if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
						first_sample_class_name+=(${cel_list[$i]})
						first_sample_class+=(${type_list[$i]})
					elif [ ${type_list[$i]} == ${type_kind[$l]} ]; then
						second_sample_class_name+=(${cel_list[$i]})
						second_sample_class+=(${type_list[$i]})
					fi 
				done
				sample_class_list=( ${first_sample_class[@]} ${second_sample_class[@]} ) 

				echo "[INFO] Input1 sample lists are : ${first_sample_class_name[@]}"
				echo "[INFO] Input2 sample lists are : ${second_sample_class_name[@]}"
				echo "[INFO] sample_class_list : ${sample_class_list[@]}"
				$NEW_R_DIR/Rscript $bin_dir"/"limma.R -a `my_join ";" ${first_sample_class_name[@]}` -b `my_join ";" ${second_sample_class_name[@]}` -r $cel_result_dir -l `my_join ";" ${sample_class_list[@]}` -p $work_file --pvalue $P_VALUE_CUT
				# get deg
<<'COMMENT'
				awk -v fc=$fold_change '{if (($3<(-1*log(fc)/log(2))) || ($3>(log(fc)/log(2)))) print}' $cel_result_dir/$work_file".limma.txt" > $cel_result_dir/$work_file".deg"
COMMENT
				awk -v pval=$P_VALUE_CUT 'NR==1{print $0}{if ($7 < pval){ print $0}}' $cel_result_dir/$work_file".limma.txt" > $cel_result_dir/$work_file".deg"
			
				tail -n+2 $cel_result_dir/$work_file".deg" | cut -f2 | sort | uniq > $cel_result_dir/$work_file".DEG.list"

				# create DVAID links based on DEG lists
				bash $bin_dir/gsea3.sh $cel_result_dir/$work_file".DEG.list" $cel_result_dir "GENE_SYMBOL" $work_file".DEG.list"


				# create KEGG input gene_symbol, color for up & down regulation
				# This will be depricated since we will generate KEGG figure locally below
				awk -v red4=$red4 -v red3=$red3 -v red2=$red2 -v red1=$red1 -v blue4=$blue4 -v blue3=$blue3 -v blue2=$blue2 -v blue1=$blue1 -v white=$white \
				'NR==1{next;}{OFS="\t"; if ($3 > 5) print $2, red4;
					else if (($3 <= 5) && ($3 > 2)) print $2, red3;
					else if (($3 <= 2) && (1 < $3)) print $2, red2;
					else if (($3 <= 1) && (0 < $3)) print $2, red1;
					else if ($3 < -5) print $2, blue4;
					else if (($3 >= -5)  && (-2 > $3)) print $2, blue3;
					else if (($3 >= -2) && (-1 > $3)) print $2, blue2;
					else if (($3 >-1) && (0 > $3)) print $2, blue1;
					else print $2, white;}'  $cel_result_dir/$work_file".deg" > $cel_result_dir/$work_file".DEG.KEGG.txt"

				# generate KEGG figures based on DEG and their fold change or pvalue
				cut -f2,3 $cel_result_dir/$work_file".deg" | grep -v -w "NA" > $cel_result_dir/$work_file".deg.for_KEGG_fig.txt"
				#$NEW_R_DIR/Rscript $bin_dir/create_kegg_fig.R -i $cel_result_dir/$work_file".deg.for_KEGG_fig.txt" -s $work_file
				#$NEW_R_DIR/Rscript $bin_dir/create_kegg_fig.R $cel_result_dir/$work_file".deg.for_KEGG_fig.txt" $work_file $cel_result_dir $cel_result_dir/$work_file".deg.for_KEGG_fig_maplist.txt"

				# create json for generate figure


				#cp $cel_result_dir/$work_file".deg.for_KEGG_fig.txt" $final_result_dir

				# copy all final results
#				cp $cel_result_dir/*.DEG.* $cel_result_dir/*.KEGG $cel_result_dir/*.UP_DOWN $final_result_dir
			done
		done

	 #make DEG exp table & heatmap
		for (( k=0; k<${#type_kind[@]}; k++ )); do 
			work_file=""
			for (( l=$k+1; l<${#type_kind[@]}; l++ )); do
				echo "[INFO] Currently processing ${type_kind[$k]} and ${type_kind[$l]}"
	
					work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}
	
					# init list
					first_sample_class_file=()
					second_sample_class_file=()
					first_sample_class_header=()
					second_sample_class_header=()

					first_sample_class=()
					second_sample_class=()
					# create sample list per subtype
					for (( i=0; i<${#cel_list[@]}; i++ )); do
						exp_file_filename_only=`basename \${cel_list[$i]}`".exp"

						if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
							first_sample_class_file+=($cel_result_dir"/"$exp_file_filename_only)
							first_sample_class_header+=($exp_file_filename_only)
							first_sample_class+=(${type_kind[$k]})
						elif [ ${type_list[$i]} == ${type_kind[$l]} ]; then
							second_sample_class_file+=($cel_result_dir"/"$exp_file_filename_only)
							second_sample_class_header+=($exp_file_filename_only)
							second_sample_class+=(${type_kind[$l]})
						fi 
					done
		
					deg_limma_file=$cel_result_dir"/"$work_file".limma.txt"

#					class1_header_string=$(IFS=$'\t'; echo "${first_sample_class_header[*]}")
					class1_header_string=`my_join "/" ${first_sample_class_header[@]}`
#					class2_header_string=$(IFS=$'\t'; echo "${second_sample_class_header[*]}")
					class2_header_string=`my_join "/" ${second_sample_class_header[@]}`
	
					work_file_exp_MGD=$cel_result_dir"/"$work_file".exp.MGD"

					# Total table
					paste ${first_sample_class_file[@]} ${second_sample_class_file[@]} | awk 'BEGIN{FS=OFS="\t"}{printf "%s\t%s", $1,$2; for(i=3;i<=NF;i=i+3){printf "\t%.4g", $(i)};printf "\n"}' >  $work_file_exp_MGD
					
					{
						echo -e "ID\tGeneSymbol\t"$class1_header_string"\t"$class2_header_string"\tP.val\tadj.P.val" | tr '/' '\t' ;
					awk 'BEGIN{FS=OFS="\t"}FILENAME==ARGV[1]{arr[$1]=$0} FILENAME==ARGV[2]{if(NR>1){if($1 in arr){printf "%s\t%.4g\t%.4g\n", arr[$1],$6,$7}}}' $work_file_exp_MGD $deg_limma_file ;
					} > $cel_result_dir"/"$work_file".Total.table.txt"
				
					# DEG table
					awk -v threshold=$P_VALUE_CUT 'BEGIN{FS=OFS="\t"}NR==1{print $0}{if($(NF) < threshold){print $0}}' $cel_result_dir"/"$work_file".Total.table.txt" > $cel_result_dir"/"$work_file".DEG.table.txt"
					
					# calculate top 30 degs for oncoprint
					awk 'BEGIN{FS=OFS="\t";i=0}NR>1{if(i<30){if(!($2 in arr) && ($2 != "NA")){i=i+1;arr[$2]=""}}else{exit}}END{for(x in arr){print x}}' $cel_result_dir"/"$work_file".DEG.table.txt" > $cel_result_dir"/"$work_file".DEG.Top30.genelist"
					# DEG Top 100
					head -n 101 $cel_result_dir"/"$work_file".DEG.table.txt" > $cel_result_dir"/"$work_file".DEG.Top100.table.txt"

					# draw heatmap
					DEG_top100_heatmap_input=$cel_result_dir"/"$work_file".DEG.Top100.table.txt.heatmap.input"
					awk -v first_num=${type_len[$k]} -v second_num=${type_len[$l]} 'BEGIN{FS=OFS="\t"}NR==1{print "ID","class1", "class2";next}{printf "%s", $1"|"$2; sum=0; for(i=3;i<(3+first_num);i++){sum=sum+$(i);};printf "\t%f", sum/first_num; sum=0; for(i=(3+first_num);i<(3+first_num+second_num);i++){sum=sum+$(i)};printf "\t%f\n", sum/second_num}' $cel_result_dir"/"$work_file".DEG.Top100.table.txt" > $DEG_top100_heatmap_input

					$R_DIR/Rscript $bin_dir/make_heatmap3_ver2_ge_v2.r $DEG_top100_heatmap_input ${type_kind[$k]}","${type_kind[$l]} "Class" ${type_kind[$k]}","${type_kind[$l]} $cel_result_dir"/"$work_file".DEG.MGD.heatmap.png" "Top100_DEG"
				
				 # draw PCA
				 DEG_PCA_input=$cel_result_dir"/"$work_file".PCA.input"
				 cut -f1,3-$((${type_len[$k]} + ${type_len[$l]}+2)) $cel_result_dir"/"$work_file".DEG.table.txt" > $DEG_PCA_input

				 $NEW_R_DIR/Rscript $bin_dir/pca.r $DEG_PCA_input `my_join "," ${first_sample_class_header[@]}`","`my_join "," ${second_sample_class_header[@]}` `my_join "," ${first_sample_class[@]}`","`my_join "," ${second_sample_class[@]}` $cel_result_dir"/"$work_file".DEG.PCA.stat.txt" $cel_result_dir"/"$work_file".DEG_PCA.var.png" $cel_result_dir"/"$work_file".DEG.PCA.png"  
			done
		done
		
		# move to final result directory
		cp $cel_result_dir"/"*.deg.for_KEGG_fig.txt $cel_result_dir"/"*.MA_plot.png $cel_result_dir"/"*.norm.boxplot.png $cel_result_dir"/"*DEG.Top100.table.txt $cel_result_dir"/"*.DEG.table.txt $cel_result_dir"/"*.Total.table.txt $cel_result_dir"/"*.volcano.png $cel_result_dir"/"*.DEG.MGD.heatmap.png $cel_result_dir"/"*.DEG.PCA.png $cel_result_dir"/"*.DEG.list.DAVID_*.txt $cel_result_dir"/"*.DEG.KEGG.txt $final_result_dir

	else

	# FOR NO REPLICATES 
	# only extract expression level since no replicates
		for (( i=0; i<${#cel_list[@]}; i=i+2 )); do
			sample_num_list=$i";"$(($i+1))
			$NEW_R_DIR/Rscript $bin_dir"/"limma.R -a ${cel_list[$i]} -b ${cel_list[$i+1]} -r $cel_result_dir -l $sample_num_list --pvalue $P_VALUE_CUT
		done
		
		if [ `echo ${#cel_list[@]}" % 2" | bc` -eq 1 ]; then
			sample_num_list="0;"$((${#cel_list[@]}-1))
			$NEW_R_DIR/Rscript $bin_dir"/"limma.R -a ${cel_list[0]} -b ${cel_list[$((${#cel_list[@]}-1))]} -r $cel_result_dir -l $sample_num_list
		fi
		
		# compute DEG by foldchange
		for (( i=0; i<${#cel_list[@]}; i++ )); do 
			for (( j=i+1; j<${#cel_list[@]}; j++ )); do 
				exp_file1_filename_only=`basename \${cel_list[$i]}`".exp"
				exp_file2_filename_only=`basename \${cel_list[$j]}`".exp"

				# get gene symbol avrage gene expression level file
				cut -f2,3 $cel_result_dir/$exp_file1_filename_only | sort | awk '{sum_exp[$1]=((sum_exp[$1]*count[$1] + $2)/(count[$1]+1));count[$1] = count[$1] + 1}END{for(key in sum_exp) print key"\t"sum_exp[key]}' > $cel_result_dir/$exp_file1_filename_only".geneSymbol_avg_exp.txt"
				cut -f2,3 $cel_result_dir/$exp_file2_filename_only | sort | awk '{sum_exp[$1]=((sum_exp[$1]*count[$1] + $2)/(count[$1]+1));count[$1] = count[$1] + 1}END{for(key in sum_exp) print key"\t"sum_exp[key]}' > $cel_result_dir/$exp_file2_filename_only".geneSymbol_avg_exp.txt"

				paste $cel_result_dir"/"$exp_file1_filename_only".geneSymbol_avg_exp.txt" $cel_result_dir"/"$exp_file2_filename_only".geneSymbol_avg_exp.txt" | awk '{print $1"\t"$2"\t"$4"\t"log($4/$2)/log(2)}' > $cel_result_dir"/"$exp_file1_filename_only"_vs_"$exp_file2_filename_only".fc"
				awk -v fc=$fold_chage '{if ($4 > (log(fc)/log(2)) || $4 < (-1 * (log(fc)/log(2)))) print $0}' $cel_result_dir"/"$exp_file1_filename_only"_vs_"$exp_file2_filename_only".fc" > $cel_result_dir"/"${type_list[$i]}"_vs_"${type_list[$j]}".deg"

				# Gene_symbol	expression1	expression2	log2_fold_change
				# AXL     4.57813287184114        10.586776865637 1.20943

				cut -f1  $cel_result_dir"/"${type_list[$i]}"_vs_"${type_list[$j]}".deg" | sort | uniq > $cel_result_dir"/"${type_list[$i]}"_vs_"${type_list[$j]}".DEG.list"

			done
		done
	fi 

	# get expression level bigWig to visualize 
	for (( i=0; i<${#cel_list[@]}; i=i+1 )); do
		exp_file1_filename_only=`basename \${cel_list[$i]}`".exp"

		if [ ${#type_kind[@]} -lt ${#type_list[@]} ] ; then
			cut -f2,3 $cel_result_dir/$exp_file1_filename_only | sort | awk '{sum_exp[$1]=((sum_exp[$1]*count[$1] + $2)/(count[$1]+1));count[$1] = count[$1] + 1}END{for(key in sum_exp) print key"\t"sum_exp[key]}' > $cel_result_dir/$exp_file1_filename_only".geneSymbol_avg_exp.txt"
		fi
		# convert expression values to genomic coordinate to visualize expression level for each sample 
#		cp $cel_result_dir/$exp_file1_filename_only".geneSymbol_avg_exp.txt" $cel_result_dir/$exp_file1_filename_only".gene_symbol_exp.txt"
		awk -v input_file=$cel_result_dir"/"$exp_file1_filename_only".geneSymbol_avg_exp.txt" 'BEGIN { while ((getline < input_file ) > 0) data[$1] = $2 }{if (data[$5]) print $1"\t"$3"\t"$4"\t"$5"\t"data[$5] }' $lib_dir"/"gene_symbol_chr_start_end.txt > $cel_result_dir"/"$exp_file1_filename_only".gene_symbol_exp.bed"
		# convert bed to bigWig to view on UCSC genome browser
		sort -k1,1 -k2,2n $cel_result_dir"/"$exp_file1_filename_only".gene_symbol_exp.bed" | bedtools merge -c 4,5 -o distinct,mean -i - | cut -f1,2,3,5 > $cel_result_dir"/"$exp_file1_filename_only".gene_symbol_exp.sorted.cut.bed"
		$bin_dir"/"bedGraphToBigWig $cel_result_dir"/"$exp_file1_filename_only".gene_symbol_exp.sorted.cut.bed" $REF_HUMAN_CHR_SIZE $cel_result_dir"/"`basename \${cel_list[$i]}`".bw"
	done

<<'COMMENT'
#TODO
	# make DEG bw
	for ((k=0; k<${#type_kind[@]}; k++)); do
		for ((l=$k+1; l<${#type_kind[@]}; l++)); do
			work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}
			deg_file=$cel_result_dir/$work_file".deg"
COMMENT

else

# 2.2 RNA-seq (.fastq) using STAR 
<<'COMMENT'
	# variable based on sample number
	for (( i=0; i<${#clipped_fq_list[@]}; i++ )); do
		# variables
		output_prefix=${sample_prefix_list[$i]}
		htseq_out=$deg_analysis_dir/$output_prefix".htseq"
		count_file_list+=($htseq_out".geneSymbol")
	done
COMMENT

mkdir -p $rnaseq_result_dir
mkdir -p $clipped_dir

if [ "$single_pair" == "single" ]; then
	# clipping adapter and fastqc
	for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do
		# variables
		sample_fq=${sample_fq1_list[$i]}

		# clip & QC
		trim_galore --fastqc -o $clipped_dir $sample_fq &
	done; wait
else
	# clipping adapter and fastqc
	for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do
		# variables
		sample_fq1=${sample_fq1_list[$i]}
		sample_fq2=${sample_fq2_list[$i]}
		# clip & QC
		trim_galore --paired --fastqc -o $clipped_dir $sample_fq1 $sample_fq2&
	done; wait
fi
	######################
	# generate STAR genome index
	######################
<<'COMMENT'
	mkdir $genomeDir
	STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeDir/hg19.fa  --runThreadN $NUM_CPUS -sjdbGTFfile $hg19_gtf

	# genome prepare for gatk (using picard)
	$java_dir/java -jar CreateSequenceDictionary.jar R= $genome_hg19 O= $genome_hg19_dict

COMMENT

file_extension=$(echo `basename \${sample_fq1_list[0]}` | awk -F . '{if (NF>1) {print $NF}}')
	
# first pass alignment
#	echo -n "" > $merged_SJ
cd $rnaseq_result_dir

if [ $single_pair == "single" ]; then
	for (( i=0; i<${#sample_fq1_list[@]}; i++)); do
		sample_fq_filename_only=`basename \${sample_fq1_list[\$i]}`
		output_prefix=$rnaseq_result_dir"/"`basename $sample_fq_filename_only "."$file_extension`"_"

		STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $clipped_dir/`basename $sample_fq_filename_only "."$file_extension`"_trimmed.fq" --runThreadN $NUM_CPUS --outFileNamePrefix $output_prefix --alignEndsType EndToEnd

	done
else
	for (( i=0; i<${#sample_fq1_list[@]}; i++)); do
		sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
		sample_fq2_filename_only=`basename \${sample_fq2_list[\$i]}`
		output_prefix=$rnaseq_result_dir"/"`basename $sample_fq1_filename_only "."$file_extension`"_"
		
		STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $clipped_dir/`basename $sample_fq1_filename_only "."$file_extension`"_val_1.fq" $clipped_dir/`basename $sample_fq2_filename_only "."$file_extension`"_val_2.fq" --runThreadN $NUM_CPUS --outFileNamePrefix $output_prefix --alignEndsType EndToEnd
	done
fi

cd $WORK_DIR

<<'COMMENT'
	for (( i=0; i<${#clipped_fq_list[@]}; i++ )); do

		# variables
		output_prefix=${sample_prefix_list[$i]}"_"
		sample_fq=${clipped_fq_list[$i]}

		# first pass
		mkdir -p $rnaseq_result_dir
		cd $rnaseq_result_dir
#		STAR --genomeDir $genomeDir --readFilesIn $sample_fq --runThreadN $NUM_CPUS --outFileNamePrefix $output_prefix  --alignEndsType EndToEnd
		
		STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $sample_fq --runThreadN $NUM_CPUS --outFileNamePrefix $output_prefix  --alignEndsType EndToEnd
		cd $WORK_DIR
	done
COMMENT
	
######################
	# ngs plot config file
	######################
	echo "[DEBUG] NOT_RUN"

	ngs_all_config=$rnaseq_result_dir/ngsplot_config_all.txt
	echo -ne "" > $ngs_all_config

	for (( k=0; k<${#type_kind[@]}; k++ )); do 
		for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do
			ngs_subtype_config=$rnaseq_result_dir/ngsplot_config_GE_${type_kind[$k]}.txt

			sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
			output_prefix=`basename $sample_fq1_filename_only "."$file_extension`"_"
	#		output_prefix=${sample_prefix_list[$i]}"_"
			second_star_out_sam=$rnaseq_result_dir"/"$output_prefix"Aligned.out.sam"
			second_star_out_bam_uniq=$rnaseq_result_dir"/"$output_prefix"Aligned.out.uniq.bam"

			config_line="$second_star_out_bam_uniq\t-1\t$output_prefix"

			echo -e $config_line >> $ngs_all_config
			if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
				# filter only unique
				samtools view -Sh $second_star_out_sam | awk '{if (substr($1, 1, 1)=="@" || $5 == 255) print  }' | samtools view -bS - | samtools sort - $rnaseq_result_dir"/"$output_prefix"Aligned.out.uniq"  &

				echo -e $config_line >> $ngs_subtype_config
			fi
		done; wait

		for region in "genebody" "cgi"; do
			ngs.plot.r -G hg19 -R $region -C $rnaseq_result_dir/ngsplot_config_GE_${type_kind[$k]}.txt -O $rnaseq_result_dir/GE_${type_kind[$k]}_$region

			# convert pdf to png
			convert -density 150 $rnaseq_result_dir/GE_${type_kind[$k]}_$region".avgprof.pdf" -quality 90 $rnaseq_result_dir/GE_${type_kind[$k]}_$region".avgprof.png"
		done
	done
	#for region in "genebody" "cgi"; do
		#ngs.plot.r -G hg19 -R $region -C $rnaseq_result_dir/ngsplot_config_all.txt -O $rnaseq_result_dir/GE_all_$region
	#done
	cp $rnaseq_result_dir/*.png $final_result_dir


	######################
	# DEG analysis
	######################
	echo "[INFO] DEG extraction"

	# get expression levels using HTseq-count
	mkdir -p $deg_analysis_dir
	for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do
	{
		# variables
		sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
		output_prefix=`basename $sample_fq1_filename_only "."$file_extension`
		htseq_out=$deg_analysis_dir/$output_prefix".htseq"
		count_file_list+=($htseq_out".geneSymbol")
		htseq-count -s no $rnaseq_result_dir"/"$output_prefix"_Aligned.out.sam" $REF_HUMAN_GTF | sort -k1,1 > $htseq_out

		# convert refseq to geneSymbol and sum based on same gene symbol
		join $htseq_out $REF_HUMAN_GENESYMBOL | awk '{print $3"\t"$2}' | sort | awk '{data[$1]=data[$1]+$2}END{for( j in data){print j"\t"data[j]}}' | sort -k1,1 > $htseq_out".geneSymbol" &
	} &
	done; wait

	# visualize bw
	for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do
		sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
		output_prefix=`basename $sample_fq1_filename_only "."$file_extension`
		htseq_out=$deg_analysis_dir/$output_prefix".htseq"
			
		awk -v input_file=$htseq_out".geneSymbol" 'BEGIN { while ((getline < input_file ) > 0) data[$1] = $2 }{if (data[$5]) print $1"\t"$3"\t"$4"\t"$5"\t"data[$5] }' $lib_dir"/"gene_symbol_chr_start_end.txt > $htseq_out".gene_symbol_exp.bed"
		# convert bed to bigWig to view on UCSC genome browser
		sort -k1,1 -k2,2n $htseq_out".gene_symbol_exp.bed" | bedtools merge -c 4,5 -o distinct,mean -i - | cut -f1,2,3,5 > $htseq_out".gene_symbol_exp.sorted.cut.bed"
		$bin_dir"/"bedGraphToBigWig $htseq_out".gene_symbol_exp.sorted.cut.bed" $REF_HUMAN_CHR_SIZE $deg_analysis_dir"/"$sample_fq1_filename_only".bw"
	done
	
	# SUBTYPE COMPARISON

	for (( k=0; k<${#type_kind[@]}; k++)); do
		for (( l=$k+1; l<${#type_kind[@]}; l++)); do
			
			# init file list
			first_sample_class_list=()
      first_sample_class_name=()
			first_sample_class_num=()
			first_sample_class_header=()
			first_sample_class_wo_ext=()

      second_sample_class_list=()
      second_sample_class_name=()	
			second_sample_class_num=()
			second_sample_class_header=()
			second_sample_class_wo_ext=()	
		
			for (( i=0; i<${#sample_fq1_list[@]}; i++)); do
				if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
					sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
					output_prefix=`basename $sample_fq1_filename_only "."$file_extension`

					first_sample_class_list+=($deg_analysis_dir/$output_prefix".htseq.geneSymbol")
					first_sample_class_name+=($sample_fq1_filename_only)
					first_sample_class_num+=(${type_kind[$k]})
					first_sample_class_header+=($sample_fq1_filename_only"_exp_level")
					first_sample_class_wo_ext+=($output_prefix)
		
				elif [ ${type_list[$i]} == ${type_kind[$l]} ]; then
					sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
					output_prefix=`basename $sample_fq1_filename_only "."$file_extension`

					second_sample_class_list+=($deg_analysis_dir/$output_prefix".htseq.geneSymbol")
					second_sample_class_name+=($sample_fq1_filename_only)
					second_sample_class_num+=(${type_kind[$l]})
					second_sample_class_header+=($sample_fq1_filename_only"_exp_level")
					second_sample_class_wo_ext+=($output_prefix)
				fi
			done

			work_file=${type_kind[$k]}"_vs_"${type_kind[$l]}	
			
			sample_list=`my_join ";" ${first_sample_class_name[@]}`";"`my_join ";" ${second_sample_class_name[@]}`
			sample_num_list=`my_join ";" ${first_sample_class_num[@]}`";"`my_join ";" ${second_sample_class_num[@]}`
			

			merged_name_with_path=$deg_analysis_dir/$work_file
			count_table=$merged_name_with_path".cnt"

			deseq_result=$merged_name_with_path".deseq"
			deseq_w_count=$merged_name_with_path".deseq.w_cnt"

			fold_change_cut_file=$deseq_w_count".FC"$fold_change
			DEG_file_name=$merged_name_with_path".FC"$fold_change".DEGs"
			DEG_w_UP_DOWN=$DEG_file_name".UP_DOWN"

			
			# merge count value to create count table for deseq
			cut_column=()
			temp_column=0
			for ((i=0; i<$((${#first_sample_class_list[@]}+${#second_sample_class_list[@]})); i++)); do
				temp_column=$(($temp_column+2))
				cut_column+=($temp_column)
			done
			
			paste -d "\t" ${first_sample_class_list[@]} ${second_sample_class_list[@]} | cut -f1,`my_join "," ${cut_column[@]}` | sort > $count_table
			

			#run deseq2			
			$R_DIR/Rscript $bin_dir/run_DESeq.R -s $sample_list -r $deg_analysis_dir -o $deseq_result -l $sample_num_list -c $count_table --output_prefix $work_file --pval_cut $P_VALUE_CUT

#			class_1_header=$(IFS=$'\t'; echo "${first_sample_class_header[*]}")
#			class_2_header=$(IFS=$'\t'; echo "${second_sample_class_header[*]}")
			class_1_header=`my_join "/" ${first_sample_class_header[@]}`
			class_2_header=`my_join "/" ${second_sample_class_header[@]}`

			deseq_w_count_header="Gene_Symbol\t"$class_1_header"\t"$class_2_header"\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tUp/Down"
			
			sample_total_num=$((${#first_sample_class_list[@]}+${#second_sample_class_list[@]}))
	
			# merge deseq & get expression level, cut by foldchange, DEG file with UP/DOWN, DEG
			(echo -e $deseq_w_count_header | tr '/' '\t' && awk -v input_file=$deseq_result 'BEGIN { while ((getline < input_file ) > 0) {data[$1] = $0}; OFS="\t"; }{if (data[$1]) print $0, data[$1]}' $count_table | awk -v sample_cnt=$sample_total_num '{OFS="\t";
									if ($(4+sample_cnt)=="NA") {print $0, "NA";} else if ($(4+sample_cnt) > 0){ print $0,  "UP";} else if ($(4+sample_cnt) <0){ print $0, "Down";} else {print $0, "-";}}' | cut --complement -f$((2+$sample_total_num)) )|tee  $deseq_w_count |
													awk -v sample_cnt=$sample_total_num -v fc=$fold_change '{if (($(3+sample_cnt)>(log(fc)/log(2)) || $(3+sample_cnt)<-(log(fc)/log(2)))&& $(3+sample_cnt)!="NA") print }' |tee $fold_change_cut_file |
										awk '{print $1"\t"$NF}' | sort | uniq |tee $DEG_w_UP_DOWN | cut -f1 > $DEG_file_name
			sort $DEG_file_name | uniq > $merged_name_with_path".DEG.list"

			# DEG table 
			pvalue_column=$((${type_len[$k]} + ${type_len[$l]}+6))
			adj_pvalue_column=$(($pvalue_column+1))

			total_table=$deg_analysis_dir/$work_file".Total.table.txt"
			cut -f1-$((${type_len[$k]} + ${type_len[$l]}+1)),$pvalue_column,$adj_pvalue_column $deseq_w_count | sort -k$((${type_len[$k]} + ${type_len[$l]}+3)),$((${type_len[$k]} + ${type_len[$l]}+3))g | awk 'BEGIN{FS=OFS="\t"}NR==1{print $0;next}{printf "%s", $1; for(i=2;i<=NF;i++){printf "\t%.4g", $(i)}}' >  $total_table

			deg_table=$deg_analysis_dir/$work_file".DEG.table.txt"
			awk -v pval=$P_VALUE_CUT 'BEGIN{FS=OFS="\t"}NR==1{print $0;next}{if(($(NF) < pval) && ($(NF) != "NA")){print $0}' > $deg_table
			
			deg_top100_table=$deg_analysis_dir/$work_file".DEG.Top100.table.txt"
			head -n 101 $deg_table >  $deg_top100_table

			# calculate top 30 degs for oncoprint
			awk 'BEGIN{FS=OFS="\t";i=0}NR>1{if(i<30){if(!($1 in arr) && ($1 != "NA")){i=i+1;arr[$1]=""}}else{exit}}END{for(x in arr){print x}}' $total_table > $deg_analysis_dir"/"$work_file".DEG.Top30.genelist"

			# draw heatmap
			DEG_top100_heatmap_input=$deg_analysis_dir"/"$work_file".DEG.Top100.table.txt.heatmap.input"
		  awk -v first_num=${type_len[$k]} -v second_num=${type_len[$l]} 'BEGIN{FS=OFS="\t"}NR==1{print "ID","class1", "class2";next}{printf "%s", $1; sum=0; for(i=2;i<(2+first_num);i++){sum=sum+$(i);};printf "\t%f", sum/first_num; sum=0; for(i=(2+first_num);i<(2+first_num+second_num);i++){sum=sum+$(i)};printf "\t%f\n", sum/second_num}' $deg_analysis_dir"/"$work_file".DEG.Top100.table.txt" > $DEG_top100_heatmap_input

			$R_DIR/Rscript $bin_dir/make_heatmap3_ver2_ge_v2.r $DEG_top100_heatmap_input ${type_kind[$k]}","${type_kind[$l]} "Class" ${type_kind[$k]}","${type_kind[$l]} $deg_analysis_dir"/"$work_file".DEG.MGD.heatmap.png" "Top100_DEG"
				
			# draw PCA
			DEG_PCA_input=$deg_analysis_dir"/"$work_file".PCA.input"
			cut -f1,2-$((${type_len[$k]} + ${type_len[$l]}+1)) $deg_analysis_dir"/"$work_file".DEG.table.txt" > $DEG_PCA_input

		 $NEW_R_DIR/Rscript $bin_dir/pca.r $DEG_PCA_input `my_join "," ${first_sample_class_wo_ext[@]}`","`my_join "," ${second_sample_class_wo_ext[@]}` `my_join "," ${first_sample_class_num[@]}`","`my_join "," ${second_sample_class_num[@]}` $deg_analysis_dir"/"$work_file".DEG.PCA.stat.txt" $deg_analysis_dir"/"$work_file".DEG_PCA.var.png" $deg_analysis_dir"/"$work_file".DEG.PCA.png"  


			# create DVAID links based on DEG lists
			sh $bin_dir/gsea3.sh $merged_name_with_path".DEG.list" $deg_analysis_dir "GENE_SYMBOL" $work_file".DEG.list"

			# create KEGG input gene_symbol, color for up & down regulation
			# This will be depricated since we will generate KEGG figure locally below
			awk -v pval=$P_VALUE_CUT 'BEGIN{FS=OFS="\t"}NR>1{if($(NF) < pval) print $0}' $deseq_result | awk -v red4=$red4 -v red3=$red3 -v red2=$red2 -v red1=$red1 -v blue4=$blue4 -v blue3=$blue3 -v blue2=$blue2 -v blue1=$blue1 -v white=$white \
				'NR>1{OFS="\t"; if ($3 > 5) print $1, red4;
					else if (($3 <= 5) && ($3 > 2)) print $1, red3;
					else if (($3 <= 2) && (1 < $3)) print $1, red2;
					else if (($3 <= 1) && (0 < $3)) print $1, red1;
					else if ($3 < -5) print $1, blue4;
					else if (($3 >= -5)  && (-2 > $3)) print $1, blue3;
					else if (($3 >= -2) && (-1 > $3)) print $1, blue2;
					else if (($3 >-1) && (0 > $3)) print $1, blue1;
					else print $1, white;}' > $deg_analysis_dir/$work_file".DEG.KEGG.txt"

			# generate KEGG figures based on DEG and their fold change or pvalue
				cut -f1,2 $deseq_result | grep -v -w "NA" > $deg_analysis_dir/$work_file".deg.for_KEGG_fig.txt"
				#$NEW_R_DIR/Rscript $bin_dir/create_kegg_fig.R -i $deg_analysis_dir/$work_file".deg.for_KEGG_fig.txt" -s $work_file
				# moved to viz
				#$NEW_R_DIR/Rscript $bin_dir/create_kegg_fig.R $deg_analysis_dir/$work_file".deg.for_KEGG_fig.txt" $work_file $deg_analysis_dir $deg_analysis_dir/$work_file".deg.for_KEGG_fig_maplist.txt"

		done
	done
	# move to final result directory
		cp $deg_analysis_dir"/"*.MA_plot.png $deg_analysis_dir"/"*.norm.boxplot.png $deg_analysis_dir"/"*DEG.Top100.table.txt $deg_analysis_dir"/"*.DEG.table.txt $deg_analysis_dir"/"*.Total.table.txt $deg_analysis_dir"/"*.volcano.png $deg_analysis_dir"/"*.DEG.MGD.heatmap.png $deg_analysis_dir"/"*.DEG.PCA.png $deg_analysis_dir"/"*.DEG.list.DAVID_*.txt $deg_analysis_dir"/"*.DEG.KEGG.txt $deg_analysis_dir"/"*.deg.for_KEGG_fig.txt $final_result_dir

fi

<<'COMMENT'

	sample_list=""
	for (( i=0; i<$sample_num; i++ )); do
		for (( j=i+1; j<$sample_num; j++ )); do

			# variables
			output_prefix=${sample_prefix_list[$i]}; output_prefix2=${sample_prefix_list[$j]};
			htseq_w_geneSymbol=$deg_analysis_dir/$output_prefix".htseq.geneSymbol"; htseq_w_geneSymbol2=$deg_analysis_dir/$output_prefix2".htseq.geneSymbol"

			pair_prefix=${sample_prefix_list[$i]}"_"${sample_prefix_list[$j]}
			sample_pair=${sample_prefix_list[$i]}";"${sample_prefix_list[$j]}
			samplenum_1=$(( $i + 1 )); samplenum_2=$(( $j + 1 ))
			sample_num_list="$samplenum_1;$samplenum_2"

			merged_name_with_path=$deg_analysis_dir/$pair_prefix
			count_table=$merged_name_with_path".cnt"

			deseq_result=$merged_name_with_path".deseq"
			deseq_w_count=$merged_name_with_path".deseq.w_cnt"
			deseq_w_count_header="Gene_Symbol\t${sample_prefix_list[$i]}_exp_level\t${sample_prefix_list[$j]}_exp_level\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tUp/Down"

			fold_change_cut_file=$deseq_w_count".FC"$fold_change
			DEG_file_name=$merged_name_with_path".FC"$fold_change".DEGs"
			DEG_w_UP_DOWN=$DEG_file_name".UP_DOWN"


			# merge count value to create count table for deseq
			paste $htseq_w_geneSymbol $htseq_w_geneSymbol2 | cut -f1,2,4 | head -n -5 | sort > $count_table

			# run deseq
			$R_DIR/Rscript $bin_dir/run_DESeq.R -s $sample_pair -r $deg_analysis_dir -o $deseq_result -l $sample_num_list -c $count_table

			# merge deseq & get expression level, cut by foldchange, DEG file with UP/DOWN, DEG
			(echo -e $deseq_w_count_header && awk -v input_file=$deseq_result 'BEGIN { while ((getline < input_file ) > 0) {data[$1] = $0}; OFS="\t"; }{if (data[$1]) print $0, data[$1]}' $count_table | awk '{OFS="\t";
									if ($6=="NA") {print $0, "NA";} else if ($6 > 0){ print $0,  "UP";} else if ($6 <0){ print $0, "Down";} else {print $0, "-";}}' | cut -f1,2,3,5- )|tee  $deseq_w_count |
													awk -v fc=$fold_change '{if (($5>(log(fc)/log(2)) || $5<-(log(fc)/log(2)))&& $5!="NA") print }' |tee $fold_change_cut_file |
										awk '{print $1"\t"$NF}' | sort | uniq |tee $DEG_w_UP_DOWN | cut -f1 > $DEG_file_name

			# create deg file list string
			deg_file_list_array+=($DEG_file_name)

			# create DVAID links based on DEG lists
			sh $bin_dir/gsea3.sh $DEG_file_name $deg_analysis_dir "GENE_SYMBOL" $pair_prefix".FC"$fold_change".DEGs"

			# create annotation link for each deg
			sh $bin_dir/create_annotation_link_html_based_on_gene_list.sh $DEG_file_name > $DEG_file_name".annot.html"

			# create KEGG input gene_symbol, color for up & down regulation
			awk -v red4=$red4 -v red3=$red3 -v red2=$red2 -v red1=$red1 -v blue4=$blue4 -v blue3=$blue3 -v blue2=$blue2 -v blue1=$blue1 -v white=$white \
				'{OFS="\t"; if ($5 > 5) print $1, red4;
					else if (($5 <= 5) && ($5 > 2)) print $1, red3;
					else if (($5 <= 2) && (1 < $5)) print $1, red2;
					else if (($5 <= 1) && (0 < $5)) print $1, red1;
					else if ($5 < -5) print $2, blue4;
					else if (($5 >= -5)  && (-2 > $5)) print $1, blue3;
					else if (($5 >= -2) && (-1 > $5)) print $1, blue2;
					else if (($5 >-1) && (0 > $5)) print $1, blue1;
					else print $1, white;}'  $fold_change_cut_file > $fold_change_cut_file".KEGG"

			# NEED TO TF TG INFORMATION
			# create TG TF KEGG
			# red3 = TF, blue3 = TG
			#refseq_tf_list="lib/refseq_mouse_grn_tf.txt"
			#awk -v red3=$red3 -v blue3=$blue3 -v input_file=$refseq_tf_list 'BEGIN { while ((getline < input_file ) > 0) data[$1] = "tf" }{if (data[$1]) {print $2"\t"red3;} else {print $2"\t"blue3;} }' $fold_change_cut_file > $fold_change_cut_file".KEGG_TG_TF"

			# create tf,tf file for cytoscape from deg file
			#sh bin/create_cytoscape_input.sh $DEG_file_name $fold_change_cut_file
		done
	done
COMMENT

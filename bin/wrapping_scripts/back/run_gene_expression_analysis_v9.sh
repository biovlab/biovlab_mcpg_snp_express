#!/bin/bash
source `dirname $0`/../../env.sh

# parse experiment information
source `dirname $0`/parse_exp_info.sh

# variables
bin_dir="$WORK_DIR/bin"
lib_dir="$WORK_DIR/lib"
result_dir="$WORK_DIR/result"
#result_dir="/usr/local/apache2/htdocs/biovlab_mcpg_snp_express/test_data/rna_paired_test_data/"

# directories
rnaseq_result_dir=$result_dir"/gene_exp/rna_seq"
clipped_dir=$rnaseq_result_dir"/cliiped/"
deg_analysis_dir=$rnaseq_result_dir"/deg/"
cel_result_dir=$result_dir"/gene_exp/cel"
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

data_dir=""

if [ "$is_array" -eq "1" ]; then
	data_dir=`dirname \${cel_list_temp[0]}`
else
	data_dir=`dirname \${sample_fq1_list_temp[0]}`
fi

cel_list=()
sample_fq1_list=()
sample_fq2_list=()

if [ "$is_array" -eq "1" ]; then
	for ((i=0;i<${#cel_list_temp[@]};i++)); do
		file_extension1=$(echo `basename \${cel_list_temp[\$i]}` | awk -F . '{if(NF>1){print $NF}}' | tr '[:upper:]' '[:lower:]')
	
		if [ $file_extension1 == "cel" ]; then
			cel_list+=(${cel_list_temp[$i]})
		elif [ $file_extension1 == "gz" ]; then
			filename1=$data_dir"/"`basename ${cel_list_temp[$i]} ".gz"`
			
			cel_list+=($filename1)
			
			gunzip -d -c ${cel_list_temp[$i]} > $filename1
		else 
			echo "Invalid file format : "$file_extension1" (Allowed : CEL, gz)"
			exit
		fi
	done

else
	for ((i=0;i<${#sample_fq1_list_temp[@]};i++)); do
		file_extension1=$(echo `basename \${sample_fq1_list_temp[\$i]}` | awk -F . '{if(NF>1){print $NF}}' | tr '[:upper:]' '[:lower:]')
	
		if [ $file_extension1 == "fq" ]; then
			sample_fq1_list+=(${sample_fq1_list_temp[$i]})
		elif [ $file_extension1 == "gz" ]; then
			filename1=$data_dir"/"`basename ${sample_fq1_list_temp[$i]} ".gz"`
		
			sample_fq1_list+=($filename1)

			gunzip -d -c ${sample_fq1_list_temp[$i]} > $filename1
		else
			echo "Invalid file format : "$file_extension1" (Allowed : fq, gz)"
			exit
		fi

		if [ $single_pair == "pair" ]; then
			file_extension2=$(echo `basename \${sample_fq2_list_temp[\$i]}` | awk -F . '{if(NF>1){print $NF}}' | tr '[:upper:]' '[:lower:]')
			if [ $file_extension2 == "fq" ]; then
				sample_fq2_list+=(${sample_fq2_list_temp[$i]})
			elif [ $file_extension2 == "gz" ]; then
				filename2=$data_dir"/"`basename ${sample_fq2_list_temp[$i]} ".gz"`
		
				sample_fq2_list+=($filename2)

				gunzip -d -c ${sample_fq2_list_temp[$i]} > $filename2
			else
				echo "Invalid file format : "$file_extension2" (Allowed : fq, gz)"
				exit
			fi
		else
			sample_fq2_list+=("NA")
		fi
	done
fi

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

if [ $is_array -eq 1 ] ; then
# 2.1 Microarray (.CELL)
	# limma phenotype group comparison. limma doesn't support 2 sample comparison unless if there are biological replicates -> just extract expression value and use fold changeAA
	# in case there is only one sample, limma.R only extract expression levels from input samples

	# create result directory
	mkdir -p $cel_result_dir
	
	if [ $replicate_exist -eq 1 ]; then	
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
				$R_DIR/Rscript $bin_dir"/"limma.R -a `my_join ";" ${first_sample_class_name[@]}` -b `my_join ";" ${second_sample_class_name[@]}` -r $cel_result_dir -l `my_join ";" ${sample_class_list[@]}` -p $work_file

				# get deg
				awk -v fc=$fold_change '{if (($3<(-1*log(fc)/log(2))) || ($3>(log(fc)/log(2)))) print}' $cel_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".limma.txt" > $cel_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".deg"
				
				tail -n+2 $cel_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".deg" | cut -f2 | sort | uniq > $cel_result_dir/${type_kind[$k]}"_vs_"${type_kind[$l]}".DEG.list"
			done
		done
	else

	# FOR NO REPLICATES 
	# only extract expression level since no replicates
		for (( i=0; i<${#cel_list[@]}; i=i+2 )); do
			sample_num_list=$i";"$(($i+1))
			$R_DIR/Rscript $bin_dir"/"limma.R -a ${cel_list[$i]} -b ${cel_list[$i+1]} -r $cel_result_dir -l $sample_num_list
		done
		
		if [ `echo ${#cel_list[@]}" % 2" | bc` -eq 1 ]; then
			sample_num_list="0;"$((${#cel_list[@]}-1))
			$R_DIR/Rscript $bin_dir"/"limma.R -a ${cel_list[0]} -b ${cel_list[$((${#cel_list[@]}-1))]} -r $cel_result_dir -l $sample_num_list
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

		if [ "$replicate_exist" -eq 1 ] ; then
			cut -f2,3 $cel_result_dir/$exp_file1_filename_only | sort | awk '{sum_exp[$1]=((sum_exp[$1]*count[$1] + $2)/(count[$1]+1));count[$1] = count[$1] + 1}END{for(key in sum_exp) print key"\t"sum_exp[key]}' > $cel_result_dir/$exp_file1_filename_only".geneSymbol_avg_exp.txt"
		fi
		# convert expression values to genomic coordinate to visualize expression level for each sample 
		cp $cel_result_dir/$exp_file1_filename_only".geneSymbol_avg_exp.txt" $cel_result_dir/$exp_file1_filename_only".gene_symbol_exp.txt"
		awk -v input_file=$cel_result_dir"/"$exp_file1_filename_only".gene_symbol_exp.txt" 'BEGIN { while ((getline < input_file ) > 0) data[$1] = $2 }{if (data[$5]) print $1"\t"$3"\t"$4"\t"$5"\t"data[$5] }' $lib_dir"/"gene_symbol_chr_start_end.txt > $cel_result_dir"/"$exp_file1_filename_only".gene_symbol_exp.bed"
		# convert bed to bigWig to view on UCSC genome browser
		sort -k1,1 -k2,2n $cel_result_dir"/"$exp_file1_filename_only".gene_symbol_exp.bed" | bedtools merge -c 4,5 -o distinct,mean -i - | cut -f1,2,3,5 > $cel_result_dir"/"$exp_file1_filename_only".gene_symbol_exp.sorted.cut.bed"
		$bin_dir"/"bedGraphToBigWig $cel_result_dir"/"$exp_file1_filename_only".gene_symbol_exp.sorted.cut.bed" $REF_HUMAN_CHR_SIZE $cel_result_dir"/"`basename \${cel_list[$i]}`".bw"
	done

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

if [ $single_pair == "single" ]; then
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
	echo -ne "" > $rnaseq_result_dir/ngsplot_config_all.txt

	for (( i=0; i<${#sample_fq1_list[@]}; i++ )); do
		sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
		output_prefix=`basename $sample_fq1_filename_only "."$file_extension`"_"
#		output_prefix=${sample_prefix_list[$i]}"_"
		second_star_out_sam=$rnaseq_result_dir"/"$output_prefix"Aligned.out.sam"
		second_star_out_bam_uniq=$rnaseq_result_dir"/"$output_prefix"Aligned.out.uniq.bam"

		# filter only unique
		samtools view -Sh $second_star_out_sam | awk '{if (substr($1, 1, 1)=="@" || $5 == 255) print  }' | samtools view -bS - | samtools sort - $rnaseq_result_dir"/"$output_prefix"Aligned.out.uniq"  &

		echo -e "$second_star_out_bam_uniq\t-1\t$output_prefix" >> $rnaseq_result_dir/ngsplot_config_all.txt
	done; wait

	#for region in "genebody" "cgi"; do
		#ngs.plot.r -G hg19 -R $region -C $rnaseq_result_dir/ngsplot_config_all.txt -O $rnaseq_result_dir/GE_all_$region
	#done


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

      second_sample_class_list=()
      second_sample_class_name=()	
			second_sample_class_num=()
			second_sample_class_header=()	
		
			for (( i=0; i<${#sample_fq1_list[@]}; i++)); do
				if [ ${type_list[$i]} == ${type_kind[$k]} ]; then
					sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
					output_prefix=`basename $sample_fq1_filename_only "."$file_extension`

					first_sample_class_list+=($deg_analysis_dir/$output_prefix".htseq.geneSymbol")
					first_sample_class_name+=($sample_fq1_filename_only)
					first_sample_class_num+=(${type_kind[$k]})
					first_sample_class_header+=($sample_fq1_filename_only"_exp_level")
				
				elif [ ${type_list[$i]} == ${type_kind[$l]} ]; then
					sample_fq1_filename_only=`basename \${sample_fq1_list[\$i]}`
					output_prefix=`basename $sample_fq1_filename_only "."$file_extension`

					second_sample_class_list+=($deg_analysis_dir/$output_prefix".htseq.geneSymbol")
					second_sample_class_name+=($sample_fq1_filename_only)
					second_sample_class_num+=(${type_kind[$l]})
					second_sample_class_header+=($sample_fq1_filename_only"_exp_level")
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
			
			paste -d "\t" ${first_sample_class_list[@]} ${second_sample_class_list[@]} | cut -f1,`my_join "," ${cut_column[@]}` | head -n -5 | sort > $count_table
			

			#run deseq2			
			$R_DIR/Rscript $bin_dir/run_DESeq.R -s $sample_list -r $deg_analysis_dir -o $deseq_result -l $sample_num_list -c $count_table


			deseq_w_count_header="Gene_Symbol\t"`my_join "\t" ${first_sample_class_header[@]}`"\t"`my_join "\t" ${second_sample_class_header[@]}`"\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tUp/Down"
			
			sample_total_num=$((${#first_sample_class_list[@]}+${#second_sample_class_list[@]}))
	
			# merge deseq & get expression level, cut by foldchange, DEG file with UP/DOWN, DEG
			(echo -e $deseq_w_count_header && awk -v input_file=$deseq_result 'BEGIN { while ((getline < input_file ) > 0) {data[$1] = $0}; OFS="\t"; }{if (data[$1]) print $0, data[$1]}' $count_table | awk -v sample_cnt=$sample_total_num '{OFS="\t";
									if ($(4+sample_cnt)=="NA") {print $0, "NA";} else if ($(4+sample_cnt) > 0){ print $0,  "UP";} else if ($(4+sample_cnt) <0){ print $0, "Down";} else {print $0, "-";}}' | cut --complement -f$((2+$sample_total_num)) )|tee  $deseq_w_count |
													awk -v sample_cnt=$sample_total_num -v fc=$fold_change '{if (($(3+sample_cnt)>(log(fc)/log(2)) || $(3+sample_cnt)<-(log(fc)/log(2)))&& $(3+sample_cnt)!="NA") print }' |tee $fold_change_cut_file |
										awk '{print $1"\t"$NF}' | sort | uniq |tee $DEG_w_UP_DOWN | cut -f1 > $DEG_file_name
			sort $DEG_file_name | uniq > $merged_name_with_path".DEG.list"
		done
	done
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

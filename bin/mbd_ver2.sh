#!bin/bash
source `dirname $0`/../env.sh

input_fastq1=$1
input_fastq2=$2
single_pair=$3
result_dir=$4
window_size=$5

# variables
input_dir=`dirname $input_fastq1`
input_fq1_filename_only=`basename \$input_fastq1`
input_fq2_filename_only=`basename \$input_fastq2`
file_extention=$(echo $input_fq1_filename_only |awk -F . '{if (NF>1) {print $NF}}')
filename1_wo_ext=`basename $input_fq1_filename_only "."$file_extention`
filename2_wo_ext=`basename $input_fq2_filename_only "."$file_extention`

################################################
# Script for MBD-seq data analysis
################################################

# multi line comment
<<'COMMENT'

echo "[INFO] Align MBDseq using Bowtie2"
if [ "$single_pair" == "pair" ]; then
	echo "[INFO] This is paired-end data"
	#ln -s $input_fastq2 $result_dir/temp/lane1_read2.fastq
	bowtie2 -q -p $NUM_CPUS -x $BOWTIE2_INDEX_HUMAN -1 $input_fastq1 -2 $input_fastq2 -S $result_dir/$input_fq1_filename_only".sam"
else
	echo "[INFO] This is single-end data"
	bowtie2 -q -p $NUM_CPUS -x $BOWTIE2_INDEX_HUMAN -U $input_fastq1 -S $result_dir/$input_fq1_filename_only".sam"
fi

# covert sam to sorted bam
echo "[INFO] Convert SAM to Sorted BAM"
samtools view -ubS $result_dir/$input_fq1_filename_only".sam" | samtools sort -m $MAX_MEM_IN_BYTE - $result_dir/$input_fq1_filename_only".sorted"

echo "[INFO] Index sorted BAM"
samtools index $result_dir/$input_fq1_filename_only".sorted.bam"



COMMENT
echo "[INFO] Align MBDseq using STAR"
if [ $single_pair == "single" ]; then
	trim_galore -o $result_dir $input_fastq1	

	cd $result_dir
	STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $result_dir/$filename1_wo_ext"_trimmed.fq" --runThreadN $SYS_NUM_CPUS --alignIntronMax 1 --outFileNamePrefix $result_dir/$input_fq1_filename_only"." --alignEndsType EndToEnd --seedSearchStartLmax 20 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 # --outFilterMatchNmin 15 # --alignEndsType EndToEnd --seedSearchStartLmax 30 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 options are for increasing sensitivity for less than 50bp reads. MBDseq ICBP 36bp length 
	cd $WORK_DIR
else
	trim_galore --paired -o $result_dir $input_fastq1 $input_fastq2
		
	cd $result_dir
	STAR --genomeDir $REF_HUMAN_DIR --readFilesIn $result_dir/$filename1_wo_ext"_val_1.fq" $result_dir/$filename2_wo_ext"_val_2.fq" --runThreadN $SYS_NUM_CPUS --alignIntronMax 1 --outFileNamePrefix $result_dir/$input_fq1_filename_only"." --alignEndsType EndToEnd --seedSearchStartLmax 20 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 # --outFilterMatchNmin 15 # --alignEndsType EndToEnd --seedSearchStartLmax 30 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 options are for increasing sensitivity for less than 50bp reads. MBDseq ICBP 36bp length 
	cd $WORK_DIR
fi


# covert sam to sorted bam
echo "[INFO] Convert SAM to Sorted BAM"
samtools view -ubS -@ $SYS_NUM_CPUS $result_dir/$input_fq1_filename_only".Aligned.out.sam" | samtools sort -@ $SYS_NUM_CPUS -m 3G - $result_dir/$input_fq1_filename_only".sorted"

echo "[INFO] Index sorted BAM"
samtools index $result_dir/$input_fq1_filename_only".sorted.bam"

<<'COMMENT'
# backup

# remove fastq additional IDs cause by GEO data remove all except '+' 
#awk '{if (NR%4==3){a=substr($1,1,1); if (a=="+") {print "+";} else {print $0;}} else {print $0;}}' BrCa-03.fastq > BrCa-03_modid.fastq

# remove previous links
#unlink $result_dir/temp/lane1_read1.fastq
#unlink $result_dir/temp/lane1_read2.fastq

# create temp dir
#mkdir -p $result_dir/temp

# create symbolic link
#ln -s $input_fastq1 $result_dir/temp/lane1_read1.fastq

	#bedtools nuc -C -pattern CG -fi /data/project/mcpg/ref_data/human/ucsc/human_hg19_ucsc_all.fa -bed ./human_hg19_ucsc_all.w100.bed > ./refgenome_profile.txt

#2 10m
#bedtools coverage -abam ../BrCa-02_head10000000.fastq.bam -b human_hg19_ucsc_all.w100.bed | bedtools sort -i - > ./BrCa-02_head10000000.fastq.cov




# Align MBDseq reads with bowtie2
#isaac-align -t $result_dir/temp -r $ref_isaac_index -o $result_dir/temp -b $result_dir/temp -m $MAX_MEM --base-calls-format fastq

# move bam to result dir 
#mv $result_dir/temp/Project/default/default/sorted.bam $result_dir/$input_fq1_filename_only.bam
#mv $result_dir/temp/Project/default/default/sorted.bai $result_dir/$input_fq1_filename_only.bai

# Analyze bam
#isaac-align -r $ref_isaac_index -b ./sorted.bam -m 40 --base-calls-format bam


#head -n 1 $result_dir/refgenome_profile.chr1.frq | awk "{print \$0\"\t\"\"methylation_level\"}" > $result_dir/temp_header.txt
#cat $result_dir/chr_list.txt  | xargs -I {} -P 13 sh -c '(head -n 1 "$2"/refgenome_profile."$1".frq | awk "{print \$0\"\t\"\"methylation_level\"}" && tail -n +2 "$2"/refgenome_profile."$1".frq | paste "$3"."$1".cov - | awk "{if (\$14==0 || \$4==0){me_level=0}else{me_level=\$4/\$14}}; print \"\$1\t\$2\t\$3\t\$4\t\$14\t\"me_level}") > "$2"/final_result.txt ' -- {} $result_dir $result_dir/$input_fq1_filename_only".sorted.bam"

# merge coverage file
#echo "[INFO] Merge coverage results"
#ls $result_dir/*.cov | xargs -I {} cat {} > $result_dir/$input_fq1_filename_only".sorted.bam.cov.sorted.temp"

# sort coverage file
#echo "[INFO] Sort coverage result"
#bedtools sort -i $result_dir/$input_fq1_filename_only".sorted.bam.cov.temp" > $result_dir/$input_fq1_filename_only".sorted.bam.cov.sorted.temp

# add header
echo "[INFO] Add header to coverage result file"
echo -e "chr\tstart\tend\tread_count\t#of_non-zero_bases\twindow_size\t%of_non-zero_cov" > $result_dir/temp_header.txt
cat $result_dir/temp_header.txt $result_dir/$input_fq1_filename_only".sorted.bam.cov.sorted.temp" > $result_dir/$input_fq1_filename_only".sorted.bam.cov.sorted.all"

# merge cov + genome frequency
echo "[INFO] Merge coverate + genome profile result"
paste $result_dir/$input_fq1_filename_only".sorted.bam.cov.sorted.all" $result_dir/refgenome_profile.sorted.all > $result_dir/$input_fq1_filename_only".sorted.bam.cov.frq.all"

	# tt contains all info (cov + genome frequency)

# extract data and compute methylation level
echo "[INFO] compute genome-wide methylation level per window"
awk '{if(NR==1){me_level="methylation_level"}else{if ($11==0 || $4==0){me_level=0}else{me_level=$4/$11}};print $1"\t"$2"\t"$3"\t"$4"\t"$11"\t"me_level}' $result_dir/$input_fq1_filename_only".sorted.bam.cov.frq.all" > $result_dir/$input_fq1_filename_only".met_level"

# custom
# compute coverage in parallel 1m
echo "[INFO] Compute genome wide coverage"
cat $result_dir/chr_list.txt  | xargs -I {} -P 13 sh -c 'samtools view -ubh "$2" {} |  bedtools coverage -abam - -b $3/genome_w"$4".{}.bed | cut -f1,2,3,4 | bedtools sort -i - > "$2"."$1".cov' -- {} $result_dir/$input_fq1_filename_only".sorted.bam" $result_dir $window_size

# compute methylation level
echo "[INFO] Compute methylation level"
cat $result_dir/chr_list.txt  | xargs -I {} -P 13 sh -c 'tail -n +2 "$2"/refgenome_profile."$1".frq | paste "$3"."$1".cov - | awk "{if (\$11==0 || \$4==0){me_level=0}else{me_level=\$4/\$11}; print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$11\"\t\"me_level;}" > "$4"."$1".met_level ' -- {} $result_dir $result_dir/$input_fq1_filename_only".sorted.bam" $result_dir/$input_fq1_filename_only

ls $result_dir/$input_fq1_filename_only.*.met_level | xargs -I {} cat {} > $result_dir/$input_fq1_filename_only".met_level.merged"




COMMENT

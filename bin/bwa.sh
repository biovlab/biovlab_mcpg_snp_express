#!bin/bash
source `dirname $0`/../env.sh

input_fastq1=$1
input_fastq2=$2
single_pair=$3
result_dir=$4
output_name=$5
num_cpus=$6

# variables
input_dir=`dirname $input_fastq1`
input_fq1_filename_only=`basename \$input_fastq1`
input_fq2_filename_only=`basename \$input_fastq2`
file_extention=$(echo $input_fq1_filename_only |awk -F . '{if (NF>1) {print $NF}}')
filename1_wo_ext=`basename $input_fq1_filename_only "."$file_extention`
filename2_wo_ext=`basename $input_fq2_filename_only "."$file_extention`

################################################
# Script for DNA-seq data analysis
################################################




echo "[INFO] Align DNAseq using bwa"
if [ $single_pair == "single" ]; then
	trim_galore -o $result_dir $input_fastq1	

	bwa mem -t $num_cpus $REF_HUMAN $result_dir/$filename1_wo_ext"_trimmed.fq" > $output_name

else
	trim_galore --paired -o $result_dir $input_fastq1 $input_fastq2
	
	bwa mem -t $num_cpus $REF_HUMAN $result_dir/$filename1_wo_ext"_val_1.fq" $result_dir/$filename2_wo_ext"_val_2.fq" > $output_name	
fi

<<'COMMENT'
# covert sam to sorted bam
echo "[INFO] Convert SAM to Sorted BAM"
samtools view -ubS -@ $NUM_CPUS $result_dir/$input_fq1_filename_only".Aligned.out.sam" | samtools sort -@ $NUM_CPUS -m 3G - $result_dir/$input_fq1_filename_only".sorted"

echo "[INFO] Index sorted BAM"
samtools index $result_dir/$input_fq1_filename_only".sorted.bam"
COMMENT

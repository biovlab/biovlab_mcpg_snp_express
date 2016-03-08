#!bin/bash
source `dirname $0`/../env.sh

sorted_bam=$1
result_dir=$2
filename_only=$3

	vcfutils.pl splitchr -l 10000000000 $REF_HUMAN".fai" | xargs --verbose -I {} -P $NUM_CPUS bash -c 'samtools mpileup -uf "$2" -r "$1" "$3" | bcftools view -bcvg - > "$4"/"$5".part-"$1".bcf' -- {} $REF_HUMAN $sorted_bam $result_dir $filename_only
	ls $result_dir/$filename_only"."part-*.bcf |sort | xargs bcftools cat > $result_dir/$filename_only".all.bcf.merged"
	bcftools view $result_dir/$filename_only".all.bcf.merged" >  $result_dir/$filename_only".vcf"
	# remove temp part bcf files
	rm -rf $result_dir/$filename_only"."part-chr*.bcf
	rm -rf $result_dir/$filename_only"."all.bcf.merged


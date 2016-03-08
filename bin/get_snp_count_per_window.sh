source `dirname $0`/../env.sh 


vcf_file=$1
result_dir=$2
window_size=$3
genome_wide_window_dir=$4

bed_file=$vcf_file".bed"

# convert vcf to bed
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+",$8}}' $vcf_file > $bed_file

# compute coverage
cat $result_dir/chr_list.txt  | xargs -I {} -P 13 sh -c 'grep -w {} $2 | bedtools coverage -a - -b $3/genome_w"$4".{}.bed | cut -f1,2,3,4 | bedtools sort -i - > "$2"."$1".cov' -- {} $bed_file $genome_wide_window_dir $window_size 

# merge coverage file
ls $bed_file"."chr*.cov | xargs -I {} cat {} > $bed_file"."all.cov

# remove chr cov files
rm -rf $bed_file"."chr*.cov

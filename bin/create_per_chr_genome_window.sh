source `dirname $0`/../env.sh 


result_dir=$1
window_size=$2
genome_wide_window_dir=$3

# create genome wide bed window file
mkdir -p $genome_wide_window_dir

bedtools makewindows -g "$REF_HUMAN".fai -w $window_size | bedtools sort -i - > $genome_wide_window_dir/genome_w"$window_size".bed

# create chr list
cut -f1 $genome_wide_window_dir/genome_w"$window_size".bed | uniq > $result_dir/chr_list.txt

# create per chr genome window file
cat $result_dir/chr_list.txt | xargs -I {} -P 13 sh -c 'grep -w {} "$2" > "$4"/genome_w"$3"."$1".bed' -- {} $genome_wide_window_dir/genome_w"$window_size".bed $window_size $genome_wide_window_dir

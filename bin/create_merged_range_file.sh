#!/bin/bash
source `dirname $0`/../env.sh

region_info_dir="$WORK_DIR/lib/region_info/"
human_refseq_various_ranges="$region_info_dir/human_refseq_various_ranges.txt"

	echo -n "" > $human_refseq_various_ranges
	for range in "3utr" "5utr" "exon" "genebody" "intron" "promoter"; do 
		# chr1	10373	11873	+	NR_046018	DDX11L1	TSS_methyl_range
		awk -v range_kind=$range '{print $1"\t"$2"\t"$3"\t"$6"\t"$5"\t"$4"\t"range_kind}' $region_info_dir/$range".bed" >> $human_refseq_various_ranges
	done


	for range in "cpgIsland" "cpgShelf" "cpgShore"; do 
		# cpg island		
		# chr1	135124	135563	CpG:30	chr1	134772	140566	LOC729737	NR_039983	-
		# chr1	762416	763445	CpG:115	chr1	762970	794826	+	NR_047522	LOC643837	0
		bedtools closest -a $region_info_dir/$range".bed" -b $WORK_DIR/lib/region_info/genebody.bed | awk -v range_kind=$range '{print $1"\t"$2"\t"$3"\t"$10"\t"$9"\t"$8"\t"range_kind}'  >> $human_refseq_various_ranges
	done


	awk -v TSS_range=$TSS_range_for_snp -v range_kind="SNP_TSS_TSE_flanking_range+-250kb" '{if ($2-TSS_range <0) {start=0;} else {start=$2-TSS_range;}; print $1"\t"start"\t"$3+TSS_range"\t"$6"\t"$5"\t"$4"\t"range_kind;} ' $WORK_DIR/lib/region_info/genebody.bed >> $human_refseq_various_ranges

	# sort for later bedtools intersect & groupby
	sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 -k6,6 -k7,7 $human_refseq_various_ranges > $region_info_dir/temp.txt && mv $region_info_dir/temp.txt $human_refseq_various_ranges


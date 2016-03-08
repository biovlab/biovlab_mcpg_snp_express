#!/bin/bash
source `dirname $0`/../../env.sh

# directories
bin_dir="$WORK_DIR/bin"
result_dir="$WORK_DIR/result"
profile_result_dir="$WORK_DIR/profile"

# variables
type_kind=("lu" "baa" "bab")

####################################################################################
# 0. parse input
####################################################################################
# set parameters here




####################################################################################
# 5. Gene set analysis : GDEGs, GDEGs related with GDMR, GDEGs related with DSNPs
####################################################################################

# 5.1 Functional annotation & Pathway analysis by DAVID
# 5.1.1 Subtype pairwise
id_type="GENE_SYMBOL"
	for (( i=0; i<${#type_kind[@]}; i++ )); do 
		for (( j=$i+1; j<${#type_kind[@]}; j++ )); do 
			# deg file 
			type_pair_prefix=${type_kind[$i]}"_vs_"${type_kind[$j]}
			deg_file=$cel_result_dir/$type_pair_prefix".deg"
			awk 'NR==1{next;}{print $2}' $deg_file | sort | uniq > $profile_result_dir/$type_pair_prefix".deg_GSym"

			sh $bin_dir"/gsea3.sh" $profile_result_dir/$type_pair_prefix".deg_GSym" $profile_result_dir $id_type $type_pair_prefix".deg"

		done
	done

# 5.1.2 all together
cut -f4 $profile_result_dir/GDEG_w_GSym.txt | sort | uniq > $profile_result_dir/GDEG_only.txt
echo "**** GDEG ****"
sh $bin_didr"/gsea3.sh" $profile_result_dir/GDEG_only.txt $profile_result_dir $id_type "GDEG"
echo ""

echo "**** GDEG that have GDMR in genebody ****"
sh $bin_dir"/gsea3.sh" $profile_result_dir/GDEGs_that_have_GDMR_in_genebody.txt $profile_result_dir $id_type "GDEGs_that_have_GDMR_in_genebody"
echo ""

echo "**** GDEG that have GDMR in promoter ****"
sh $bin_dir"/gsea3.sh" $profile_result_dir/GDEGs_that_have_GDMR_in_promoter.txt $profile_result_dir $id_type "GDEGs_that_have_GDMR_in_promoter"



# copy results to web directory
cp $profile_result_dir/*.DAVID.clustering.html $profile_result_dir/*.DAVID_table.html $profile_result_dir/*.DAVID_chart.html $WEB_ACCESSIBLE_DIR



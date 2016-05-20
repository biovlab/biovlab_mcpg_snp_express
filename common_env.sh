#!/bin/sh
[ -z "$WORK_DIR" ] && WORK_DIR=$(readlink -f $(dirname ${BASH_SOURCE[0]}))

# This is common (local, cloud) environmental setup for BioVLAB-MCPG-SNP_EXPRESS

# reference genomes
BOWTIE2_INDEX_HUMAN="$WORK_DIR/ref_data/human/ucsc/numeric_order/hg19"
REF_HUMAN="$WORK_DIR/ref_data/human/ucsc/numeric_order/hg19.fa"
REF_HUMAN_DICT="$WORK_DIR/ref_data/human/ucsc/numeric_order/hg19.dict"
REF_HUMAN_CHR_SIZE="$WORK_DIR/ref_data/human/ucsc/hg19.chrom.sizes"
REF_HUMAN_GENE_INFO="$WORK_DIR/lib/human_refseq_info.txt"
REF_HUMAN_GENE_RANGE_INFO="$WORK_DIR/lib/region_info/human_refseq_various_ranges.txt.only_NM.txt"
REF_HUMAN_GENE_RANGE_INFO_GROUP_BY_REFSEQ="$WORK_DIR/lib/region_info/human_refseq_various_ranges.txt.only_NM_groupby_refseq.txt"
REF_HUMAN_GENE_RANGE_INFO_DIV_100="$WORK_DIR/lib/region_info/human_refseq_various_ranges.txt.only_NM.div_100.txt"
REF_HUMAN_GENE_GENEBODY_INFO="$WORK_DIR/ref_data/human/ucsc/region_info/genebody.bed.GeneSymbol.txt.only_NM.txt"
REF_HUMAN_GENE_GENEBODY_INFO_SORTED="$WORK_DIR/ref_data/human/ucsc/region_info/genebody.bed.GeneSymbol.txt.only_NM.txt.sorted"
REF_HUMAN_GENE_TSS_TSE_INFO="$WORK_DIR/lib/human_refseq_TSS_TSE_list_fix.txt"
REF_HUMAN_GENOME_100_BIN_BED="$WORK_DIR/ref_data/human/ucsc/hg19_100bp.bed"
REF_HUMAN_GENOME_10M_BIN_BED="$WORK_DIR/ref_data/human/ucsc/hg19.chrom.sizes_10000000.bed"
REF_HUMAN_GENOME_1M_BIN_BED="$WORK_DIR/ref_data/human/ucsc/hg19.chrom.sizes_1000000.bed"
HUMAN_GENE_SYMBOL_LIST="$WORK_DIR/lib/gene_symbol_list.txt"
REF_HUMAN_DIR="$WORK_DIR/ref_data/human/ucsc/numeric_order/"
BISMARK_REF_DIR_HUMAN="$WORK_DIR/ref_data/human/ucsc/numeric_order/bisulfite_index/"
#BISMARK_REF_DIR_HUMAN="$WORK_DIR/ref_data/human/ucsc/test_ref/bisulfite_index/"
METHYLKIT_BED_REF_HUMAN="$WORK_DIR/ref_data/human/ucsc/human_hg19_ucsc.bed.txt"
METHYLKIT_BED_REF_HUMAN_CPGI="$WORK_DIR/ref_data/human/ucsc/human_hg19_CpGi_ucsc.bed.txt"
DHS_CLUSTER_UCSC="$WORK_DIR/lib/region_info/DnaseCluster.bed"
REF_HUMAN_GTF="$WORK_DIR/ref_data/human/ucsc/human_hg19_ucsc.gtf"
REF_HUMAN_GENESYMBOL="$WORK_DIR/ref_data/human/ucsc/refseq2geneSymbol.txt.sorted"
HUMAN_HG19_CPGSITE="$WORK_DIR/ref_data/human/ucsc/hg19_cpg_site.bed"
ENCODE_TF_GENE_INFO="$WORK_DIR/lib/region_info/hg19_promoter_only_NM_refseq_Encode_TF.bed"
GENOME_VERSION="hg19"
REF_HUMAN_PROMOTER="$WORK_DIR/lib/region_info/promoter_NM_only.bed.valid_chr.bed"

TEN_COLOR_CODE=('#fb8072' '#80b1d3' '#fdb462' '#b3de69' '#fccde5' '#bc80bd' '#8dd3c7' '#ffffb3' '#bebada' '#ccebc5')
GENESYMBOL_GENEBODY_INFO="$WORK_DIR/ref_data/human/ucsc/region_info/GeneSymbol_genebody.bed"
HALLMARK_GENE_SET_JASON="/data/project/MSG/lib/hallmark_genes_msigdb.json"

# MYSQL DB setting
MYSQL_DB_NAME="msg"
MSQL_USER_ID="airavata"
MYSQL_PASSWD="airavata"
MYSQL_QUERY="SELECT c.class_name, p.patient_name, s.sample_type, s.sample_format, s.is_pair, s.replicates, f.file_path, f.url, f.file_name, f.pair_number, f.replicate_number, e.* FROM experiments e INNER JOIN experiment_files ef ON e.id = ef.experiment_id INNER JOIN files f ON f.id = ef.file_id INNER JOIN sample_groups s ON s.id = f.sample_group_id INNER JOIN patients p ON p.id = s.patient_id INNER JOIN classes c ON c.id = p.class_id WHERE e.uid="
MYSQL_HOSTNAME="bhi2.snu.ac.kr"
MYSQL_HOST_PORT="3306"
# NOTE : ADD airavata remote account and grand all with set passwd
# after root login
# grant all on *.* to airavata@'%';
# SET PASSWORD FOR airavata@'%' = PASSWORD('airavata');

# MYSQL field name matching
CLASS_NAME="class_name"
PATIENT_NAME="patient_name"
SAMPLE_TYPE="sample_type"
SAMPLE_FORMAT="sample_format"
IS_PAIR="is_pair"
FILE_PATH="file_path"
FILE_NAME="file_name"
PAIR_NUMBER="pair_number"
REPLICATE_NUMBER="replicate_number"
REPLICATES="replicates"
URL="url"
ID="id"
PGA_UID="uid"
EXPRESSION_TYPE="expression_type"
METHYLATION_TYPE="methylation_type"
MUTATION_TYPE="mutation_type"
P_VALUE_CUT=0.15


# BioVLAB
BIOVLAB_EMAIL=snu.biovlab@gmail.com
export WEB_ACCESSIBLE_LOC="http://bhi2.snu.ac.kr:3000/"

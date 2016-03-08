# load MEDIPS library
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(optparse)

# make options
option_list <- list( 
    make_option(c("-a", "--sample_list1"), action="store"),
    make_option(c("-r", "--result_dir"), action="store")
    )

opt <- parse_args(OptionParser(option_list=option_list))
# assumes input1;input1_rep1;input1_rep2,input2;input2_rep1,input3;input3_rep1;input3_rep2

# parameters
genome="BSgenome.Hsapiens.UCSC.hg19"
uniq=TRUE
extend=300
shift=0
window_size=100

sample_list1 <- unlist(strsplit(opt$sample_list1, ";"))
result_dir <- opt$result_dir

# create MEDIP set1 
for (replicate in sample_list1) {

	temp_MEDIPS1=list() 
	# create MEDIP set using all replicates for single phenotype

	# saturation analysis per replicate
	sr = MEDIPS.saturation(file=replicate, BSgenome=genome, uniq=uniq, extend=extend, shift=shift, window_size=window_size, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE)

		# draw saturation plot
		result_fig <- paste(replicate,"_sat.jpg", collapse = "", sep="") 
		jpeg(result_fig)
		MEDIPS.plotSaturation(sr)
		dev.off()


	# coverage analysis
	cr = MEDIPS.seqCoverage(file=replicate, BSgenome=genome, pattern = "CG", uniq=uniq, extend=extend, shift=shift)
		# draw plot
		result_fig <- paste(replicate,"_cov.jpg", collapse = "", sep="") 
		jpeg(result_fig)
		MEDIPS.plotSeqCoverage(seqCoverageObj = cr, type = "pie", cov.level = c(0,1,2,3,4,5))
		dev.off()

	# compute CpG enrichment
	er = MEDIPS.CpGenrich(file=replicate, BSgenome=genome, extend=extend, shift=shift, uniq=uniq)
	cat("[INFO] er\n")

	# save 
	result_file <- paste(replicate,"_cpg_er.txt", collapse = "", sep="") 
	con <- file(result_file)
	sink(con, append=FALSE, split=TRUE)
	print(er)
	# Restore output to console
	sink() 

}

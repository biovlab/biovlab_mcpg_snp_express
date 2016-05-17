# load MEDIPS library
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(optparse)

# make options
option_list <- list( 
    make_option(c("-a", "--sample_list1"), action="store"),
    make_option(c("-b", "--sample_list2"), action="store"),
    make_option(c("-r", "--result_dir"), action="store"),
    make_option(c("-p", "--result_prefix"), action="store"),
		make_option(c("--pvalue"), action="store", default=0.05)
#    make_option(c("-l", "--sample_num_list"), action="store"),
    #make_option(c("-c", "--count_table"), action="store")
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
sample_list2 <- unlist(strsplit(opt$sample_list2, ";"))
result_dir <- opt$result_dir
result_prefix <- opt$result_prefix
pvalue_cut <- as.numeric(opt$pvalue)

# functions
to_file <- function(filename, content){

	write.table(content, filename, sep="\t", quote=FALSE, row.names=F)
	#save(content,file=filename)

	#con <- file(filename)
	#sink(con, append=FALSE, split=FALSE)

	# max print 
	#options(max.print=1.0E8)
	#print(content)

	# Restore output to console
	#sink() 

	return("")
}
if(FALSE){
# Sequence Pattern Coverage
cat("[INFO] Draw Sequence pattern coverage plots")
for (replicate in sample_list1) {
	cr = MEDIPS.seqCoverage(file = replicate, pattern = "CG", BSgenome = genome, extend = extend, shift = shift,uniq = uniq)
	
	result_file <- paste(result_dir,"/", result_prefix, ".CoverageStats.png", collapse = "", sep="")
	
	png(result_file)
	MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", cov.level = c(0,1, 2, 3, 4, 5))
	dev.off()
}
for (replicate in sample_list2) {
	cr = MEDIPS.seqCoverage(file = replicate, pattern = "CG", BSgenome = genome, extend = extend, shift = shift,uniq = uniq)
	
	result_file <- paste(result_dir,"/", result_prefix, ".CoverageStats.png", collapse = "", sep="")
	
	png(result_file)
	MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", cov.level = c(0,1, 2, 3, 4, 5))
	dev.off()
}
}
# create MEDIP set1 
cat("[INFO] Create input1 MEDIP dataset")
temp_MEDIPS1=list() 
for (replicate in sample_list1) {

	# create MEDIP set using all replicates for single phenotype

	# merge to MEDIP set
	temp_MEDIPS1=c(temp_MEDIPS1,MEDIPS.createSet(file=replicate, BSgenome=genome, uniq=uniq, extend=extend, shift=shift, window_size=window_size))

#	cat(temp_MEDIPS1)
}

# For CpG density dependent normalization of MeDIP-seq data, we need to generate a coupling set
cat("[INFO] Create coupling set for CpG densigy dependent normalization")
CS = MEDIPS.couplingVector(pattern = "CG", refObj = temp_MEDIPS1[[1]])

# create MEDIP set2 
cat("[INFO] Create input2 MEDIP dataset")
temp_MEDIPS2=list() 
for (replicate in sample_list2) {

	# create MEDIP set using all replicates for single phenotype
	temp_MEDIPS2=c(temp_MEDIPS2,MEDIPS.createSet(file=replicate, BSgenome=genome, uniq=uniq, extend=extend, shift=shift, window_size=window_size))

#	cat(temp_MEDIPS2)
}

# coverage methylation profilies and differential coverage
mr.edgeR = MEDIPS.meth(MSet1 = temp_MEDIPS1, MSet2 = temp_MEDIPS2,  CSet = CS, p.adj = "bonferroni",  diff.method = "edgeR", prob.method = "poisson", MeDIP = T, CNV = F, type = "rpkm", minRowSum = 1)
	# print result	
	cat("[INFO] mr.edgeR\n")

	# FIXME !!! name change required
	result_file <- paste(result_dir,"/", result_prefix, ".mr.edgeR.txt", collapse = "", sep="") 

	# store to file
	to_file(result_file, mr.edgeR)
	cat("[INFO] mr.edgeR done\n")

# select significant window
mr.edgeR.s = MEDIPS.selectSig(results = mr.edgeR, p.value = pvalue_cut,adj = T, ratio = NULL, bg.counts = NULL, CNV = F)
	cat("[INFO] mr.edgeR.s\n")

	# FIXME !!! name change required
	result_file <- paste(result_dir,"/", result_prefix, ".mr.edgeR.s.txt", collapse = "", sep="") 

	# store to file
	to_file(result_file, mr.edgeR.s)
	cat("[INFO] mr.edgeR.s done\n")

# merge frame 
mr.edgeR.s.gain = mr.edgeR.s[which(mr.edgeR.s[, grep("logFC", colnames(mr.edgeR.s))] > 0), ]
	cat("[INFO] mr.edgeR.gain\n")

	# FIXME !!! name change required
	result_file <- paste(result_dir,"/", result_prefix, ".mr.edgeR.s.gain.txt", collapse = "", sep="") 

	# store to file
	to_file(result_file, mr.edgeR.s.gain)
	cat("[INFO] mr.edgeR.gain done\n")

mr.edgeR.s.gain.m = MEDIPS.mergeFrames(frames = mr.edgeR.s.gain, distance = 1)
	cat("[INFO] mr.edgeR.gain.m\n")

	# FIXME !!! name change required
	result_file <- paste(result_dir,"/", result_prefix, ".mr.edgeR.s.gain.m.txt", collapse = "", sep="") 

	# store to file
	to_file(result_file, mr.edgeR.s.gain.m)
	cat("[INFO] mr.edgeR.gain.m done\n")

# WiG format export
#result_file <- paste(result_dir,"/CS.wig", collapse = "", sep="") 
#MEDIPS.exportWIG(Set = CS[[1]], file=result_file, format="pdensity", descr="")
if(FALSE){
result_file <- paste(result_dir,"/", result_prefix, ".file1.wig", collapse = "", sep="") 
MEDIPS.exportWIG(Set = temp_MEDIPS1[[1]], file=result_file, format="rpkm", descr="")

result_file <- paste(result_dir,"/", result_prefix, ".file2.wig", collapse = "", sep="") 
MEDIPS.exportWIG(Set = temp_MEDIPS2[[1]], file=result_file, format="rpkm", descr="")
}

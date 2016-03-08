# load MEDIPS library
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(optparse)

# make options
option_list <- list( 
    make_option(c("-a", "--sample_list1"), action="store"),
    make_option(c("-b", "--sample_list2"), action="store"),
    make_option(c("-r", "--result_dir"), action="store"),
    #make_option(c("-o", "--deseq_result"), action="store"),
    make_option(c("-l", "--sample_num_list"), action="store")
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

# parsing arguments
sample_list1 <- unlist(strsplit(opt$sample_list1, ";"))
sample_list2 <- unlist(strsplit(opt$sample_list2, ";"))

# create MEDIP set1 
for (replicate in sample_list1) {

	temp_MEDIPS1=list() 
	# create MEDIP set using all replicates for single phenotype

	# saturation analysis per replicate
	sr = MEDIPS.saturation(file=replicate, BSgenome=genome, uniq=uniq, extend=extend, shift=shift, window_size=window_size, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE)

		# draw saturation plot
		result_fig <- paste(result_dir,"/",replicate,"_sat.jpg", collapse = "", sep="") 
		jpeg(result_fig)
		MEDIPS.plotSaturation(sr)
		dev.off()


	# coverage analysis
	cr = MEDIPS.seqCoverage(file=replicate, BSgenome=genome, pattern = "CG", uniq=uniq, extend=extend, shift=shift)
		# draw plot
		result_fig <- paste(result_dir,"/",replicate,"_cov.jpg", collapse = "", sep="") 
		jpeg(result_fig)
		MEDIPS.plotSeqCoverage(seqCoverageObj = cr, type = "pie", cov.level = c(0,1,2,3,4,5))
		dev.off()

	# compute CpG enrichment
	er = MEDIPS.CpGenrich(file = replicate, BSgenome=genome, extend=extend, shift=shift, uniq=uniq)
	cat "[INFO] er\n"
	er

	# merge to MEDIP set
	temp_MEDIPS1=c(temp_MEDIPS1,MEDIPS.createSet(file=replicate, BSgenome=genome, uniq=uniq, extend=extend, shift=shift, window_size=window_size))
}

# For CpG density dependent normalization of MeDIP-seq data, we need to generate a coupling set
CS = MEDIPS.couplingVector(pattern = "CG", refObj = temp_MEDIPS1[[1]])

# create MEDIP set2 
for (replicate in sample_list2) {

	temp_MEDIPS2=list() 

	# saturation analysis per replicate
	sr = MEDIPS.saturation(file=replicate, BSgenome=genome, uniq=uniq, extend=extend, shift=shift, window_size=window_size, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE)
		result_fig <- paste(result_dir,"/",replicate,"_sat.jpg", collapse = "", sep="") 
		jpeg(result_fig)
		# draw saturation plot
		MEDIPS.plotSaturation(sr)
		dev.off()

	# coverage analysis
	cr = MEDIPS.seqCoverage(file=replicate, BSgenome=genome, pattern = "CG", uniq=uniq, extend=extend, shift=shift)
		# draw plot
		result_fig <- paste(result_dir,"/", replicate,"_cov.jpg", collapse = "", sep="") 
		jpeg(result_fig)
		MEDIPS.plotSeqCoverage(seqCoverageObj = cr, type = "pie", cov.level = c(0,1,2,3,4,5))
		dev.off()

	# compute CpG enrichment
	er = MEDIPS.CpGenrich(file = replicate, BSgenome=genome, extend=extend, shift=shift, uniq=uniq)
	cat "[INFO] er\n"
	er

	# create MEDIP set using all replicates for single phenotype
	temp_MEDIPS2=c(temp_MEDIPS2,MEDIPS.createSet(file=replicate, BSgenome=genome, uniq=uniq, extend=extend, shift=shift, window_size=window_size))
}

# coverage methylation profilies and differential coverage
mr.edgeR = MEDIPS.meth(MSet1 = temp_MEDIPS1, MSet2 = temp_MEDIPS2,  CSet = CS, ISet1 = temp_MEDIPS1, ISet2 = temp_MEDIPS2, p.adj = "bonferroni",  diff.method = "edgeR", prob.method = "poisson", MeDIP = T, CNV = F, type = "rpkm", minRowSum = 1)
	# print result	
	cat "[INFO] mr.edgeR\n"
	mr.edgeR



# select significant window
mr.edgeR.s = MEDIPS.selectSig(results = mr.edgeR, p.value = 0.1,adj = T, ratio = NULL, bg.counts = NULL, CNV = F)
	cat "[INFO] mr.edgeR.s\n"
	mr.edgeR.s

# merge frame 
mr.edgeR.s.gain = mr.edgeR.s[which(mr.edgeR.s[, grep("logFC", colnames(mr.edgeR.s))] > 0), ]
	cat "[INFO] mr.edgeR.gain\n"
	mr.edgeR.s.gain


mr.edgeR.s.gain.m = MEDIPS.mergeFrames(frames = mr.edgeR.s.gain, distance = 1)
	cat "[INFO] mr.edgeR.gain.m\n"
	mr.edgeR.s.gain.m

















# create data set
#input1=MEDIPS.createSet(file="/data/project/mcpg/result/BrCa-02.fastq.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", uniq=TRUE, extend=300, shift=0, window_size=100)
	# add list
	#  CS = MEDIPS.couplingVector(pattern = "CG", refObj = hESCs_MeDIP[[1]])

	# in case there are replicates
	#input1_rep1=MEDIPS.createSet(file="/data/project/mcpg/result/BrCa-02_rep1.fastq.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", uniq=TRUE, extend=300, shift=0, window_size=100)

	# create MEDIPS SET using replicates
	#input1_MEDIP = c(input1, input1_rep1)

#input2=MEDIPS.createSet(file="/data/project/mcpg/result/result/BrCa-03.fastq.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", uniq=TRUE, extend=300, shift=0, window_size=100)
#Reading bam alignment BrCa-03.fastq.bam
	# add list
	#  CS = MEDIPS.couplingVector(pattern = "CG", refObj = hESCs_MeDIP[[1]])

	#input1_rep2=MEDIPS.createSet(file="/data/project/mcpg/result/BrCa-03_rep1.fastq.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", uniq=TRUE, extend=300, shift=0, window_size=100)

	# create MEDIPS SET using replicates
	#input2_MEDIP = c(input2, input2_rep1)

# saturation analysis
#sr = MEDIPS.saturation(file="/data/project/mcpg/result/BrCa-02.fastq.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", uniq=TRUE, extend=300, shift=0, window_size=100, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE)

	# draw saturation plot
	#MEDIPS.plotSaturation(sr)

# coverage analysis
#cr = MEDIPS.seqCoverage(file="/data/project/mcpg/result/BrCa-02.fastq.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", pattern = "CG", uniq=TRUE, extend=300, shift=0)

	# draw plot
	#MEDIPS.plotSeqCoverage(seqCoverageObj = cr, type = "pie", cov.level = c(0,1,2,3,4,5))

# compute CpG enrichment
#er = MEDIPS.CpGenrich(file = "/data/project/mcpg/result/BrCa-02.fastq.bam", BSgenome="BSgenome.Hsapiens.UCSC.hg19", extend=300, shift=0, uniq=TRUE)



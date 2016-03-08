# load MEDIPS library
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(optparse)
require(foreach)
require(doMC)
# make options
option_list <- list( 
    make_option(c("-a", "--sample_list1"), action="store"),
    make_option(c("-b", "--sample_list2"), action="store"),
    make_option(c("-r", "--result_dir"), action="store"),
    make_option(c("-p", "--result_prefix"), action="store"),
    make_option(c("-P", "--num_process"), action="store"),
		make_option(c("--pvalue"), action="store"),
		make_option(c("-w", "--window_size"), action="store")
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
num_process <- as.numeric(opt$num_process)
window_size <- as.numeric(opt$window_size)
registerDoMC(cores=num_process)
pval <- as.numeric(opt$pvalue)

# functions
to_file <- function(filename, content){

	write.table(content, filename, sep="\t", quote=FALSE, row.names=F)
	return("")
}

MEDIPS_meth_parallel <- function (MSet1 = NULL, MSet2 = NULL, CSet = NULL, ISet1 = NULL, 
    ISet2 = NULL, chr = NULL, p.adj = "bonferroni", diff.method = "ttest", 
    prob.method = "poisson", CNV = FALSE, MeDIP = TRUE, type = "rpkm", 
    minRowSum = 1) 
{
    nMSets1 = length(MSet1)
    nMSets2 = length(MSet2)
    nISets1 = length(ISet1)
    nISets2 = length(ISet2)
    if (is.list(CSet)) 
        CSet = CSet[[1]]
    if (!is.list(MSet1)) 
        MSet1 = c(MSet1)
    if (!is.list(MSet2)) 
        MSet2 = c(MSet2)
    if (!is.list(ISet1)) 
        ISet1 = c(ISet1)
    if (!is.list(ISet2)) 
        ISet2 = c(ISet2)
    if (MeDIP) {
        if (class(CSet) != "COUPLINGset") {
            stop("You have to state a COUPLINGset object!")
        }
    }
    controlSet = MSet1[[1]]
    for (i in 1:nMSets1) {
        if (class(MSet1[[i]]) != "MEDIPSset" & class(MSet1[[i]]) != 
            "MEDIPSroiSet") {
            stop("You have to state a MEDIPSset or MEDIPSroiSet object!")
        }
        if (length(genome_count(MSet1[[i]])) != length(genome_count(controlSet))) 
            stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
    }
    if (!is.null(MSet2)) {
        for (i in 1:nMSets2) {
            if (class(MSet2[[i]]) != "MEDIPSset" & class(MSet2[[i]]) != 
                "MEDIPSroiSet") {
                stop("You have to state a MEDIPSset or MEDIPSroiSet object!")
            }
            if (length(genome_count(MSet2[[i]])) != length(genome_count(controlSet))) 
                stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
        }
    }
    if (!is.null(ISet1)) {
        for (i in 1:nISets1) {
            if (class(ISet1[[i]]) != "MEDIPSset" & class(ISet1[[i]]) != 
                "MEDIPSroiSet") {
                stop("You have to state a MEDIPSset or MEDIPSroiSet object!")
            }
            if (length(genome_count(ISet1[[i]])) != length(genome_count(controlSet))) 
                stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
        }
    }
    if (!is.null(ISet2)) {
        for (i in 1:nISets2) {
            if (class(ISet2[[i]]) != "MEDIPSset" & class(ISet2[[i]]) != 
                "MEDIPSroiSet") {
                stop("You have to state a MEDIPSset or MEDIPSroiSet object!")
            }
            if (length(genome_count(ISet2[[i]])) != length(genome_count(controlSet))) 
                stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
        }
    }
    if (class(controlSet) == "MEDIPSset") {
        window_size = window_size(controlSet)
        no_chr_windows = ceiling(chr_lengths(controlSet)/window_size(controlSet))
        supersize_chr = cumsum(no_chr_windows)
        GRanges.genome = MEDIPS.GenomicCoordinates(supersize_chr, 
            no_chr_windows, chr_names(controlSet), chr_lengths(controlSet), 
            window_size(controlSet))
    }
    else if (class(controlSet) == "MEDIPSroiSet") {
        GRanges.genome = rois(controlSet)
        window_size = as.numeric(width(GRanges.genome))
    }
    if (MeDIP) {
        base = data.frame(chr = as.vector(seqnames(GRanges.genome)), 
            start = start(GRanges.genome), stop = end(GRanges.genome), 
            CF = genome_CF(CSet), stringsAsFactors = F)
    }
    else {
        base = data.frame(chr = as.vector(seqnames(GRanges.genome)), 
            start = start(GRanges.genome), stop = end(GRanges.genome), 
            stringsAsFactors = F)
    }
    rm(controlSet)
    gc()
    counts.medip = NULL
    rpkm.medip = NULL
    rms = NULL
    prob = NULL
    counts.input = NULL
    rpkm.input = NULL
    if (!is.null(MSet1)) {
        cat("Preprocessing MSet1 in parallel\n")
	tryCatch({
	counts.medip = as.matrix(foreach(i=1:nMSets1, .combine='cbind', .inorder=TRUE) %dopar% {
			       genome_count(MSet1[[i]])})
	colnames(counts.medip) = rep("MSet1", nMSets1)
	}, error=function(err){
		cat("[Error] Fail to count.medip\n")
		print(nMSets1)
		stop(err)
	}, finally={
	})
	
	tryCatch({
	rpkm.medip = as.matrix(foreach(i=1:nMSets1, .combine='cbind', .inorder=TRUE) %dopar%{
			     ((genome_count(MSet1[[i]]) * 10^9)/(window_size * number_regions(MSet1[[i]])))})
	}, error=function(err){
		cat("[Error] Fail to rpkm.medip\n")
		print(nMSets1)
		stop(error)
	}, finally={
	})

        if (MeDIP) {
		registerDoMC(cores=3)
		tryCatch({
		rms = as.matrix(foreach(i=1:nMSets1, .combine='cbind', .inorder=TRUE) %dopar%{
            		cat(paste("Processing MSet1 SET rms", i, " ...\n")) 
                	MEDIPS.rms(MSet1[[i]], CSet, ccObj = MEDIPS.calibrationCurve(MSet = MSet1[[i]], CSet = CSet, input = F))
		})
		}, error=function(err){
			cat("[Error] Fail to rms\n")
			print(nMSets1)
			stop(err)
		}, finally={
		})

                if (prob.method == "poisson") {
			tryCatch({
			prob = as.matrix(foreach(i=1:nMSets1, .combine='cbind', .inorder=TRUE) %dopar%{
            			cat(paste("Processing MSet1 SET prob", i, " ...\n")) 
				MEDIPS.pois(MSet1[[i]], CSet, ccObj = MEDIPS.calibrationCurve(MSet = MSet1[[i]], CSet = CSet, input = F))
			})
			}, error=function(err){
				cat("[Error] Fail to prob\n")
				print(nMSets1)
				stop(err)
			}, finally={
			})
		}
                else if (prob.method == "negBinomial") {
			tryCatch({
			prob = as.matrix(foreach(i=1:nMSets1, .combine='cbind', .inorder=TRUE) %dopar%{
                		ccObj = MEDIPS.calibrationCurve(MSet = MSet1[[i]], CSet = CSet, input = F)
				MEDIPS.negBin(MSet1[[i]], CSet, ccObj = ccObj)
			})
			}, error=function(err){
				cat("[Error] Fail to prob 2\n")
				print(nMSets1)
				stop(err)
			}, finally={
			})
		}
                else { stop(paste("Method ", prob.method, " not supported.", sep = "")) }
	}
	registerDoMC(cores=num_process)
    }
    if (!is.null(MSet2)) {
        cat("Preprocessing MSet2 in parallel\n")
	counts.medip = cbind(counts.medip,foreach(i=1:nMSets2, .combine='cbind', .inorder=TRUE) %dopar% {
			       genome_count(MSet2[[i]])})

	colnames(counts.medip) = c(rep("MSet1", nMSets1),rep("MSet2", nMSets2))

	rpkm.medip = cbind(rpkm.medip, foreach(i=1:nMSets2, .combine='cbind', .inorder=TRUE) %dopar%{
			     ((genome_count(MSet2[[i]]) * 10^9)/(window_size * number_regions(MSet2[[i]])))})


        if (MeDIP) {
		registerDoMC(cores=3)
		rms = cbind(rms, foreach(i=1:nMSets2, .combine='cbind', .inorder=TRUE) %dopar%{
            		cat(paste("Processing MSet2 SET rms", i, " ...\n")) 
                	MEDIPS.rms(MSet2[[i]], CSet, ccObj = MEDIPS.calibrationCurve(MSet = MSet2[[i]], CSet = CSet, input = F))
		})
                if (prob.method == "poisson") {
			prob = cbind(prob, foreach(i=1:nMSets2, .combine='cbind', .inorder=TRUE) %dopar%{
            			cat(paste("Processing MSet2 SET prob", i, " ...\n")) 
				MEDIPS.pois(MSet2[[i]], CSet, ccObj = MEDIPS.calibrationCurve(MSet = MSet2[[i]], CSet = CSet, input = F))
			})
		}
                else if (prob.method == "negBinomial") {
			prob = cbind(prob,foreach(i=1:nMSets2, .combine='cbind', .inorder=TRUE) %dopar%{
                		ccObj = MEDIPS.calibrationCurve(MSet = MSet2[[i]], CSet = CSet, input = F)
				MEDIPS.negBin(MSet2[[i]], CSet, ccObj = ccObj)
			})
		}
                else { stop(paste("Method ", prob.method, " not supported.", sep = "")) }

	}
	registerDoMC(cores=num_process)
    }
    if (!is.null(ISet1)) {
        for (i in 1:nISets1) {
            cat(paste("Preprocessing INPUT SET ", i, " in ISet1...\n", 
                sep = ""))
            counts.input = cbind(counts.input, ISet1 = genome_count(ISet1[[i]]))
            rpkm.input = cbind(rpkm.input, ((genome_count(ISet1[[i]]) * 
                10^9)/(window_size * number_regions(ISet1[[i]]))))
        }
    }
    if (!is.null(ISet2)) {
        for (i in 1:nISets2) {
            cat(paste("Preprocessing INPUT SET ", i, " in ISet2...\n", 
                sep = ""))
            counts.input = cbind(counts.input, ISet2 = genome_count(ISet2[[i]]))
            rpkm.input = cbind(rpkm.input, ((genome_count(ISet2[[i]]) * 
                10^9)/(window_size * number_regions(ISet2[[i]]))))
        }
    }
    if (!is.null(chr)) {
        fi = base[, 1] %in% chr
        cat("Extracting data for", chr, "...\n", sep = " ")
        if (length(fi) == 0) {
            stop("Stated chromosome does not exist in the COUPLING SET.")
        }
        if (!is.null(counts.medip)) {
            counts.medip = counts.medip[fi, ]
            rpkm.medip = rpkm.medip[fi, ]
            rms = rms[fi, ]
            prob = prob[fi, ]
        }
        if (!is.null(counts.input)) {
            counts.input = counts.input[fi, ]
            rpkm.input = rpkm.input[fi, ]
        }
        base = base[fi, ]
        cat(nrow(base), "windows on", chr, "\n", sep = " ")
    }
    col.names.count = NULL
    col.names.rpkm = NULL
    col.names.rms = NULL
    col.names.prob = NULL
    if (nMSets1 != 0) {
        for (i in 1:nMSets1) {
            col.names.count = c(col.names.count, paste(sample_name(MSet1[[i]]), 
                ".counts", sep = ""))
            col.names.rpkm = c(col.names.rpkm, paste(sample_name(MSet1[[i]]), 
                ".rpkm", sep = ""))
            if (MeDIP) {
                col.names.rms = c(col.names.rms, paste(sample_name(MSet1[[i]]), 
                  ".rms", sep = ""))
                col.names.prob = c(col.names.prob, paste(sample_name(MSet1[[i]]), 
                  ".prob", sep = ""))
            }
        }
    }
    if (nMSets2 != 0) {
        for (i in 1:nMSets2) {
            col.names.count = c(col.names.count, paste(sample_name(MSet2[[i]]), 
                ".counts", sep = ""))
            col.names.rpkm = c(col.names.rpkm, paste(sample_name(MSet2[[i]]), 
                ".rpkm", sep = ""))
            if (MeDIP) {
                col.names.rms = c(col.names.rms, paste(sample_name(MSet2[[i]]), 
                  ".rms", sep = ""))
                col.names.prob = c(col.names.prob, paste(sample_name(MSet2[[i]]), 
                  ".prob", sep = ""))
            }
        }
    }
    if (nMSets1 != 0 | nMSets2 != 0) {
        counts.medip = data.frame(counts.medip)
        colnames(counts.medip) = col.names.count
        rpkm.medip = data.frame(rpkm.medip)
        colnames(rpkm.medip) = col.names.rpkm
        if (MeDIP) {
            rms = data.frame(rms)
            colnames(rms) = col.names.rms
            prob = data.frame(prob)
            colnames(prob) = col.names.prob
        }
    }
    col.names.count.input = NULL
    col.names.rpkm.input = NULL
    if (nISets1 != 0) {
        for (i in 1:nISets1) {
            col.names.count.input = c(col.names.count.input, 
                paste(sample_name(ISet1[[i]]), ".counts", sep = ""))
            col.names.rpkm.input = c(col.names.rpkm.input, paste(sample_name(ISet1[[i]]), 
                ".rpkm", sep = ""))
        }
    }
    if (nISets2 != 0) {
        for (i in 1:nISets2) {
            col.names.count.input = c(col.names.count.input, 
                paste(sample_name(ISet2[[i]]), ".counts", sep = ""))
            col.names.rpkm.input = c(col.names.rpkm.input, paste(sample_name(ISet2[[i]]), 
                ".rpkm", sep = ""))
        }
    }
    if (nISets1 != 0 | nISets2 != 0) {
        counts.input = data.frame(counts.input)
        colnames(counts.input) = col.names.count.input
        rpkm.input = data.frame(rpkm.input)
        colnames(rpkm.input) = col.names.rpkm.input
    }
    if (!is.null(MSet1) & !is.null(MSet2)) {
        cat(paste("Differential coverage analysis...\n", sep = " "))
        if ((nMSets1 < 3 | nMSets2 < 3) & diff.method == "ttest") {
            stop("Method 'ttest' is not valid for less than 3 replicates per group. Method 'edgeR' can be applied in this case.")
        }
        if (diff.method == "edgeR") {
            if (!is.null(MSet1)) {
                n.r.M1 = NULL
                for (i in 1:nMSets1) {
                  n.r.M1 = c(n.r.M1, number_regions(MSet1[[i]]))
                }
            }
            if (!is.null(MSet2)) {
                n.r.M2 = NULL
                for (i in 1:nMSets2) {
                  n.r.M2 = c(n.r.M2, number_regions(MSet2[[i]]))
                }
            }
            diff.results.list = MEDIPS.diffMeth(base = base, 
                values = counts.medip, diff.method = "edgeR", 
                nMSets1 = nMSets1, nMSets2 = nMSets2, p.adj = p.adj, 
                n.r.M1 = n.r.M1, n.r.M2 = n.r.M2, MeDIP = MeDIP, 
                minRowSum = minRowSum)
        }
        else if (diff.method == "ttest") {
            if (type == "rpkm") {
                diff.results.list = MEDIPS.diffMeth(base = base, 
                  values = rpkm.medip, diff.method = "ttest", 
                  nMSets1 = nMSets1, nMSets2 = nMSets2, p.adj = p.adj, 
                  MeDIP = MeDIP, minRowSum = minRowSum)
            }
            else if (type == "rms") {
                if (MeDIP) {
                  diff.results.list = MEDIPS.diffMeth(base = base, 
                    values = rms, diff.method = "ttest", nMSets1 = nMSets1, 
                    nMSets2 = nMSets2, p.adj = p.adj, MeDIP = MeDIP, 
                    minRowSum = minRowSum)
                }
                else {
                  stop("Invalid specification for parameter type because parameter MeDIP is FALSE (no rms values have been calculated).")
                }
            }
            else {
                stop("Unknown specification for parameter type.")
            }
        }
        else {
            stop("Selected method for calculating differential coverage not supported")
        }
        cat("Please note, log2 ratios are reported as log2(MSet1/MSet2).\n")
        diff.results = diff.results.list$diff.results
        diff.index = diff.results.list$diff.index
        rm(diff.results.list)
        gc()
    }
    else {
        cat("No differential coverage will be calculated- only one group of MEDIPS SETs given.\n")
    }
    if (CNV) {
        if (!is.null(ISet1) & !is.null(ISet2)) {
            cat(paste("CNV analysis...\n", sep = " "))
            cnv.combined = MEDIPS.cnv(base = base, rpkm.input = rpkm.input, 
                nISets1 = nISets1, nISets2 = nISets2)
        }
        else {
            cat("Cannot perform CNV analysis- please specify two groups of INPUT SETs!\n")
        }
    }
    cat(paste("Creating results table...\n", sep = " "))
    if (!is.null(counts.medip)) {
        if (MeDIP) {
            results = data.frame(base, counts.medip, rpkm.medip, 
                rms, prob, stringsAsFactors = F)
        }
        else {
            results = data.frame(base, counts.medip, rpkm.medip, 
                stringsAsFactors = F)
        }
    }
    if (!is.null(counts.input)) {
        if (!is.null(counts.medip)) {
            results = data.frame(results, counts.input, rpkm.input, 
                stringsAsFactors = F)
        }
        else {
            results = data.frame(base, counts.input, rpkm.input, 
                stringsAsFactors = F)
        }
    }
    set1idx = 1:(nMSets1)
    counts.mean.C = numeric(dim(counts.medip)[1])
    rpkm.mean.C = numeric(dim(rpkm.medip)[1])
    for (i in set1idx) {
        counts.mean.C = counts.mean.C + counts.medip[, i]
        rpkm.mean.C = rpkm.mean.C + rpkm.medip[, i]
    }
    counts.mean.C = counts.mean.C/nMSets1
    rpkm.mean.C = rpkm.mean.C/nMSets1
    if (MeDIP) {
        rms.mean.C = numeric(dim(rms)[1])
        prob.mean.C = numeric(dim(prob)[1])
        for (i in set1idx) {
            rms.mean.C = rms.mean.C + rms[, i]
            prob.mean.C = prob.mean.C + prob[, i]
        }
        rms.mean.C = rms.mean.C/nMSets1
        prob.mean.C = prob.mean.C/nMSets1
        results = data.frame(results, MSets1.counts.mean = counts.mean.C, 
            MSets1.rpkm.mean = rpkm.mean.C, MSets1.rms.mean = rms.mean.C, 
            MSets1.prob.mean = prob.mean.C, stringsAsFactors = F)
        rm(counts.mean.C, rpkm.mean.C, set1idx, rms.mean.C, prob.mean.C)
    }
    else {
        results = data.frame(results, MSets1.counts.mean = counts.mean.C, 
            MSets1.rpkm.mean = rpkm.mean.C, stringsAsFactors = F)
        rm(counts.mean.C, rpkm.mean.C, set1idx)
    }
    if (nMSets2 > 0) {
        set2idx = (nMSets1 + 1):(nMSets1 + nMSets2)
        counts.mean.T = numeric(dim(counts.medip)[1])
        rpkm.mean.T = numeric(dim(rpkm.medip)[1])
        for (i in set2idx) {
            counts.mean.T = counts.mean.T + counts.medip[, i]
            rpkm.mean.T = rpkm.mean.T + rpkm.medip[, i]
        }
        counts.mean.T = counts.mean.T/nMSets2
        rpkm.mean.T = rpkm.mean.T/nMSets2
        if (MeDIP) {
            rms.mean.T = numeric(dim(rms)[1])
            prob.mean.T = numeric(dim(prob)[1])
            for (i in set2idx) {
                rms.mean.T = rms.mean.T + rms[, i]
                prob.mean.T = prob.mean.T + prob[, i]
            }
            rms.mean.T = rms.mean.T/nMSets2
            prob.mean.T = prob.mean.T/nMSets2
            results = data.frame(results, MSets2.counts.mean = counts.mean.T, 
                MSets2.rpkm.mean = rpkm.mean.T, MSets2.rms.mean = rms.mean.T, 
                MSets2.prob.mean = prob.mean.T, stringsAsFactors = F)
            rm(counts.mean.T, rpkm.mean.T, set2idx, rms.mean.T, 
                prob.mean.T)
        }
        else {
            results = data.frame(results, MSets2.counts.mean = counts.mean.T, 
                MSets2.rpkm.mean = rpkm.mean.T, stringsAsFactors = F)
            rm(counts.mean.T, rpkm.mean.T, set2idx)
        }
    }
    if (nISets1 > 1) {
        setI1idx = 1:(nISets1)
        counts.input.mean.C = counts.input[, setI1idx[1]]
        rpkm.input.mean.C = rpkm.input[, setI1idx[1]]
        for (i in setI1idx[-1]) {
            counts.input.mean.C = counts.input.mean.C + counts.input[, 
                i]
            rpkm.input.mean.C = rpkm.input.mean.C + rpkm.input[, 
                i]
        }
        counts.input.mean.C = counts.input.mean.C/nISets1
        rpkm.input.mean.C = rpkm.input.mean.C/nISets1
        results = data.frame(results, ISets1.counts.mean = counts.input.mean.C, 
            ISets1.rpkm.mean = rpkm.input.mean.C, stringsAsFactors = F)
        rm(counts.input.mean.C, rpkm.input.mean.C, setI1idx)
    }
    if (nISets2 > 1) {
        setI2idx = (nISets1 + 1):(nISets1 + nISets2)
        counts.input.mean.T = counts.input[, setI2idx[1]]
        rpkm.input.mean.T = rpkm.input[, setI2idx[1]]
        for (i in setI2idx[-1]) {
            counts.input.mean.T = counts.input.mean.T + counts.input[, 
                i]
            rpkm.input.mean.T = rpkm.input.mean.T + rpkm.input[, 
                i]
        }
        counts.input.mean.T = counts.input.mean.T/nISets2
        rpkm.input.mean.T = rpkm.input.mean.T/nISets2
        results = data.frame(results, ISets2.counts.mean = counts.input.mean.T, 
            ISets2.rpkm.mean = rpkm.input.mean.T, stringsAsFactors = F)
        rm(counts.input.mean.T, rpkm.input.mean.T, setI2idx)
    }
    if (MeDIP) {
        rm(base, counts.medip, rpkm.medip, rms, prob, counts.input, 
            rpkm.input)
    }
    else {
        rm(base, counts.medip, rpkm.medip, counts.input, rpkm.input)
    }
    gc()
    if (nMSets1 != 0 & nMSets2 != 0) {
        cat(paste("Adding differential coverage results...\n", 
            sep = " "))
        dummy.results = matrix(ncol = ncol(diff.results), nrow = nrow(results))
        if (diff.method == "edgeR") {
            c.names = colnames(diff.results)
            diff.results <- matrix(unlist(diff.results), ncol = ncol(diff.results), 
                byrow = FALSE)
            colnames(diff.results) = c.names
            rm(c.names)
        }
        dummy.results[diff.index, ] = diff.results
        colnames(dummy.results) = colnames(diff.results)
        results = data.frame(results, dummy.results, stringsAsFactors = F)
        rm(diff.results, dummy.results, diff.index)
        gc()
    }
    if (!is.null(ISet1) & !is.null(ISet2)) {
        if (CNV) {
            cat(paste("Adding CNV results...\n", sep = " "))
            dummy.results = matrix(ncol = 1, nrow = (nrow(results)))
            for (i in 1:nrow(cnv.combined)) {
                dummy.results[cnv.combined[i, 1]:cnv.combined[i, 
                  2]] = cnv.combined[i, 3]
            }
            colnames(dummy.results) = "CNV.log2.ratio"
            results = data.frame(results, dummy.results, stringsAsFactors = F)
            rm(dummy.results)
            gc()
        }
    }
    rownames(results) = seq(1, nrow(results))
    gc()
    return(results)
}

# create MEDIP set1 
cat("[INFO] Create input1 MEDIP dataset")
temp_MEDIPS1=foreach(replicate=sample_list1, .combine='c', .inorder=TRUE) %dopar% {
	MEDIPS.createSet(file=replicate, BSgenome=genome, uniq=uniq, extend=extend, shift=shift, window_size=window_size)}
temp_MEDIPS1=as.list(c(temp_MEDIPS1))

# For CpG density dependent normalization of MeDIP-seq data, we need to generate a coupling set
cat("[INFO] Create coupling set for CpG densigy dependent normalization")
CS = MEDIPS.couplingVector(pattern = "CG", refObj = temp_MEDIPS1[[1]])

# create MEDIP set2 
cat("[INFO] Create input2 MEDIP dataset")
temp_MEDIPS2=foreach(replicate=sample_list2, .combine='c', .inorder=TRUE) %dopar% {
	MEDIPS.createSet(file=replicate, BSgenome=genome, uniq=uniq, extend=extend, shift=shift, window_size=window_size)}
temp_MEDIPS2=as.list(c(temp_MEDIPS2))

# coverage methylation profilies and differential coverage

	chr_length=length(chr_names(temp_MEDIPS1[[1]]))
	temp_interval<-seq(from = 1, to=chr_length-1, by = 2)
	#cat("[INFO] Start MEDIPS.meth\n")
	#mr.edgeR=foreach(chr=temp_interval, .combine='rbind', .inorder=TRUE) %dopar% {
	#	next_chr=chr+1
	#	MEDIPS.meth(MSet1 = temp_MEDIPS1, MSet2 = temp_MEDIPS2,  CSet = CS, ISet1 = temp_MEDIPS1, ISet2 = temp_MEDIPS2, p.adj = "bonferroni",  diff.method = "edgeR", prob.method = "poisson", MeDIP = T, CNV = F, type = "rpkm", minRowSum = 1, chr_names(temp_MEDIPS1[[1]])[chr:next_chr])
	#}
	#cat("[INFO] End MEDIPS.meth\n")

	cat("[INFO] Start MEDIPS.meth(1/2)\n")
mr.edgeR = MEDIPS_meth_parallel(MSet1 = temp_MEDIPS1, MSet2 = temp_MEDIPS2,  CSet = CS, ISet1 = temp_MEDIPS1, ISet2 = temp_MEDIPS2, p.adj = "bonferroni",  diff.method = "edgeR", prob.method = "poisson", MeDIP = T, CNV = F, type = "rpkm", minRowSum = 1, chr_names(temp_MEDIPS1[[1]])[1:11])

	cat("[INFO] Start MEDIPS.meth(2/2)\n")
mr.edgeR2 = MEDIPS_meth_parallel(MSet1 = temp_MEDIPS1, MSet2 = temp_MEDIPS2,  CSet = CS, ISet1 = temp_MEDIPS1, ISet2 = temp_MEDIPS2, p.adj = "bonferroni",  diff.method = "edgeR", prob.method = "poisson", MeDIP = T, CNV = F, type = "rpkm", minRowSum = 1, chr_names(temp_MEDIPS1[[1]])[12:chr_length])

	
	mr.edgeR = rbind(mr.edgeR, mr.edgeR2)
	rm(mr.edgeR2)


	# print result	
	cat("[INFO] mr.edgeR data storing start\n")
	# FIXME !!! name change required
	result_file <- paste(result_dir,"/", result_prefix, ".mr.edgeR.txt", collapse = "", sep="") 

	# store to file
	write.table(mr.edgeR, result_file, sep="\t", quote=FALSE, row.names=F)
	#to_file(result_file, mr.edgeR)
	cat("[INFO] mr.edgeR done\n")
	cat("[INFO] Start MEDIPS.selectSig\n")
# select significant window
mr.edgeR.s = MEDIPS.selectSig(results = mr.edgeR, p.value = pval,adj = T, ratio = NULL, bg.counts = NULL, CNV = F)
	cat("[INFO] done MEDIPS.selectSig\n")

	# FIXME !!! name change required
	result_file <- paste(result_dir,"/", result_prefix, ".mr.edgeR.s.txt", collapse = "", sep="") 

	cat("[INFO] mr.edgeR.s data storing start\n")
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

result_file <- paste(result_dir,"/", result_prefix, ".file1.wig", collapse = "", sep="") 
MEDIPS.exportWIG(Set = temp_MEDIPS1[[1]], file=result_file, format="rpkm", descr="")

result_file <- paste(result_dir,"/", result_prefix, ".file2.wig", collapse = "", sep="") 
MEDIPS.exportWIG(Set = temp_MEDIPS2[[1]], file=result_file, format="rpkm", descr="")

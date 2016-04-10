# load DESeq library
library(DESeq2)
library(optparse)

resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}

# make options
option_list <- list( 
    make_option(c("-s", "--sample_list"), action="store"),
    make_option(c("-r", "--result_dir"), action="store"),
    make_option(c("-o", "--deseq_result"), action="store"),
    make_option(c("-l", "--sample_num_list"), action="store"),
    make_option(c("-c", "--count_table"), action="store"),
		make_option(c("--output_prefix"), action="store")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# parsing arguments
result_dir <- opt$result_dir
sample_list <- unlist(strsplit(opt$sample_list, ";"))
sample_num_list <- unlist(strsplit(opt$sample_num_list, ";"))
count_table_file <- opt$count_table
final_result_file <- opt$deseq_result

#output_prefix <- paste(unique(sample_num_list), collapse="_")
output_prefix <- opt$output_prefix

print("**** arguments ****")
result_dir
sample_list
count_table_file
final_result_file
sample_num_list
output_prefix
print("**** arguments ****")

# read couting table
mycounts <- read.table(count_table_file, header=FALSE, row.names=1)

# set sample information
# temp <- c(test[3:length(test)]) handle input argument like this to make below
		# for DESeq condition list	

#samples <- data.frame(row.names=c("Nip_2D_0_1_head100000.fastq","Nip_2D_6_1_head100000.fastq"), condition=as.factor(c(rep("Nip_2D_0_1_head100000.fastq",1),rep("Nip_2D_6_1_head100000.fastq"))))
samples <- data.frame(row.names=sample_list, condition=as.factor(sample_num_list))

# create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = mycounts, colData=samples, design=~condition)

# run DESeq analysis
dds <- DESeq(dds)

# store result
#res <- results(dds, contrasts=c("conditions", sample_num_list))
res <- results(dds)

#res <- results(dds)
res <- res[order(res$padj),]


res

# create MA plot
result_dir
result_fig <- paste(result_dir, "/",output_prefix,".MA_plot.png", collapse = "", sep="") 
result_fig
png(result_fig)
plotMA(dds,ylim=c(-2,2),main=output_prefix)
dev.off()

#cooks distance boxplot
boxplot_fig <- paste(result_dir, "/",output_prefix,".norm.boxplot.png", collapse = "", sep="")

par(mar=c(8,5,2,2))
png(boxplot_fig)
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

par(resetPar())

# write result to file
write.table(res,file=final_result_file, quote= FALSE , sep="\t")

mcols(res, use.names=TRUE)

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

#Heatmaps showing the expression data of the 30 most highly expressed genes
#raw counts 
result_fig <- paste(result_dir, "/",output_prefix,"_heatmap1.jpg", collapse = "", sep="") 
png(result_fig)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,8), cexCol=1)
dev.off()

#Heatmaps showing the expression data of the 30 most highly expressed genes
#from regularized log transformation
result_fig <- paste(result_dir, "/",output_prefix,"_heatmap2.jpg", collapse = "", sep="") 
png(result_fig)
heatmap.2(assay(rld)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 8),cexCol=1)
dev.off()

#Heatmaps showing the expression data of the 30 most highly expressed genes
#from regularized log transformation
result_fig <- paste(result_dir, "/",output_prefix,"_heatmap3.jpg", collapse = "", sep="") 
png(result_fig)
heatmap.2(assay(vsd)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 8),cexCol=1)
dev.off()

if (length(sample_list) >2){
	result_fig <- paste(result_dir, "/",output_prefix,"_heatmap4.jpg", collapse = "", sep="") 
	png(result_fig)
	distsRL <- dist(t(assay(rld)))
	mat <- as.matrix(distsRL)
	rownames(mat) <- colnames(mat) <- with(colData(dds), condition)
	heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
	dev.off()

	result_fig <- paste(result_dir, "/",output_prefix,"_heatmap5.jpg", collapse = "", sep="") 
	png(result_fig)
	print(plotPCA(rld, intgroup=c("condition")))
	dev.off()

}
result_fig <- paste(result_dir, "/",output_prefix,"_heatmap6.jpg", collapse = "", sep="") 
png(result_fig)
plotDispEsts(dds)
dev.off()

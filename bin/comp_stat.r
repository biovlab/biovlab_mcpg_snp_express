library(optparse)
library(methylKit)
library(graphics)

#make options
option_list = list(
		make_option(c("-s", "--sample_list"), action="store"),
		make_option(c("-i", "--sample.id_list"), action="store"),
		make_option(c("-c", "--class_list"), action="store"),
		make_option(c("-a", "--genome_assembly"), action="store",default="hg19"),
		make_option(c("-m", "--methy_context"), action="store", default="CpG"),
		make_option("--reference", action="store", default="/data/project/mcpg/ref_data/human/ucsc/human_hg19_ucsc.bed.txt"),
		make_option("--cpgi_reference", action="store",default="/data/project/mcpg/ref_data/human/ucsc/human_hg19_CpGi_ucsc.bed.txt"),
		make_option(c("-r", "--result_dir"), action="store"),
		make_option(c("-o", "--commeRbund_result"), action="store"),
		make_option(c("-w", "--window_size"), action="store", default=100),
		make_option(c("-n", "--normalize"), action="store", default=FALSE),
		make_option(c("--pval"), action="store", default=0.01),
		make_option("--out_prefix", action="store"),
		make_option(c("-l", "--sample_num_list"), action="store")
		)

opt = parse_args(OptionParser(option_list=option_list))

#parsing arguments
sample_list = lapply(unlist(strsplit(opt$sample_list, ",")),function(x) x)
sample.id_list = lapply(unlist(strsplit(opt$sample.id_list, ",")), function(x) x)
class_list = unlist(lapply(unlist(strsplit(opt$class_list, ",")),as.numeric))
genome_assembly = opt$genome_assembly
m_context = opt$methy_context
result_dir = opt$result_dir
ref = opt$reference
cpgi_ref = opt$cpgi_reference
window = as.numeric(opt$window_size)
normal_flag = opt$normalize
pval_cut = as.numeric(opt$pval)

comp_sample = opt$out_prefix

#file read

myobj = read(
		location=sample_list,
		sample.id=sample.id_list,
		assembly=genome_assembly,
		pipeline='bismark',
		context=m_context,
		treatment=class_list
		)

#sample statistics

basename_list = lapply(sample_list,basename)
#stat_out_list = lapply(basename_list, function(x) paste(result_dir,x,".stat.txt",sep=""))

comp_out_file = paste(result_dir,"/", comp_sample, ".stat.txt",sep="")

##Normalizing steps
if(normal_flag)
{
	#filtering
	filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)

	#normalize
	normalized.myobj = normalizeCoverage(filtered.myobj,method="median")

	myobj = normalized.myobj
}

#tiling - w/o normalize
tiles = tileMethylCounts(myobj, win.size=window, step.size=window)

#file union
meth = unite(myobj, destrand=FALSE)
methr = unite(tiles, destrand=FALSE)

#mr_out_file = paste(result_dir,"/",comp_sample,".MR_list.bed",sep="")
#write.table(getData(methr), mr_out_file, quote=FALSE, row.name=FALSE, col.name=FALSE, sep="\t")



sink(comp_out_file)

#sample correlation
cat("[Correlation Analysis - base]\n")
#getCorrelation(meth, plot=FALSE)
png(paste(result_dir,"/",comp_sample, ".Correlation.png",sep=""))
getCorrelation(meth, plot=TRUE)
dev.off()
cat("\n")

#cat("[Correlation Analysis - region]\n")
#getCorrelation(methr, plot=FALSE)
#getCorrelation(methr, plot=TRUE)
#cat("\n")



#clustering samples
cat("[Clustering Analysis - base]\n")
png(paste(result_dir,"/",comp_sample, ".Cluster.png",sep=""))
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
cat("\n")

#cat("[Clustering Analysis - region]\n")
#clusterSamples(methr, dist="correlation", method="ward", plot=TRUE)
#cat("\n")


#PCA analysis
#cat("[PCA Analysis - base]\n")
png(paste(result_dir,"/",comp_sample, ".PCA_screeplot.png",sep=""))
PCASamples(meth, screeplot=TRUE)
dev.off()
png(paste(result_dir,"/",comp_sample, ".PCA_plot.png",sep=""))
PCASamples(meth) #PCA plot with PC1 and PC2 axis and a scatter plot
dev.off()
cat("\n")

#cat("[PCA Analysis - region]\n")
#PCASamples(methr, screeplot=TRUE)
#PCASamples(methr) #PCA plot with PC1 and PC2 axis and a scatter plot
#cat("\n")


##DMR analysis
cat("[DMR Analysis]\n")
#q-value calculation
myDiff = calculateDiffMeth(methr)

#q-value<0.01 and percent methylation difference larger than 25%.
#get hyper methylated bases
myDiff25p.hyper=get.methylDiff(myDiff, difference=25,qvalue=pval_cut,type="hyper")
#get hypo methylated bases
myDiff25p.hypo=get.methylDiff(myDiff, difference=25,qvalue=pval_cut,type="hypo")
#get all differentially methylated bases
myDiff25p = get.methylDiff(myDiff, difference=25, qvalue=pval_cut)

cat("[Hyper Methylated Region (versus control)]\n")
myDiff25p.hyper
cat("\n")

cat("[Hypo Methylated Region (versus control)]\n")
myDiff25p.hypo
cat("\n")

cat("[All Differentially Methylated Region (versus control)]\n")
myDiff25p
cat("\n")

#distribution of hypo/hyper-methylated bases/regions per chromosome

diffMethPerChr(myDiff, plot=FALSE, qvalue.cutoff=0.01, meth.cutoff=25)
png(paste(result_dir,"/",comp_sample, ".diffMethPerChr.png",sep=""))
diffMethPerChr(myDiff, plot=TRUE, qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()

sink()

#make DMR list
dmr_out_file = paste(result_dir,"/",comp_sample,".DMR_list.txt",sep="")
hyper_dmr_out_file = paste(result_dir,"/",comp_sample,".hyper_DMR_list.txt",sep="")
hypo_dmr_out_file = paste(result_dir,"/",comp_sample,".hypo_DMR_list.txt",sep="")
dmr_merge_out_file = paste(result_dir,"/",comp_sample,".DMR_merge_list.txt",sep="")


dmr_data = getData(myDiff25p)

write.table(dmr_data, dmr_out_file, quote=FALSE, row.names=FALSE, col.name=FALSE, sep="\t")



write.table(getData(myDiff25p.hyper), hyper_dmr_out_file, quote=FALSE, row.names=FALSE, col.name=FALSE, sep="\t")



write.table(getData(myDiff25p.hypo), hypo_dmr_out_file, quote=FALSE, row.names=FALSE, col.name=FALSE, sep="\t")

sink(dmr_merge_out_file)

n = dim(dmr_data)[1]

end = -1
chr = ""

#cat("chr\tstart\tend\n")

for (i in 1:n)
{
	if (dmr_data$start[i]-1 != end | as.character(dmr_data$chr[i]) != chr)	
	{
		if(end != (-1))
		{
			cat(end)
			cat("\n")
		}
		cat(as.character(dmr_data$chr[i]))
		cat("\t")
		cat(dmr_data$start[i])
		cat("\t")
	}

	chr = as.character(dmr_data$chr[i])
	end = dmr_data$end[i]
}
cat(end)

sink()

sink(comp_out_file, append=TRUE)

## not running ??

if (FALSE){

#annotating
gene.obj = read.transcript.features(ref)

head(gene.obj)

diffAnn=annotate.WithGenicParts(myDiff25p,gene.obj)

cat("[Summary of Ref annotation]\n")
diffAnn
cat("\n")

#cpgi_annotating
cpg.obj=read.feature.flank(cpgi_ref, feature.flank.name=c("CpGi","shores"))

diffCpGann=annotate.WithFeature.Flank(myDiff25p,cpg.obj$CpGi,cpg.obj$shores,
		feature.name="CpGi",flank.name="shores")

cat("[Summary of CpGi annotation]\n")
diffCpGann
cat("\n")

#regional analysis
cat("[Get distance to nearest TSS and gene id from annotationByGenicPart]\n")
cat("#target.row is the row number in myDiff25p\n")
gaTTS = getAssociationWithTSS(diffAnn)
head(gaTTS)
cat("\n")

cat("[get percentage/number of diffrentially methylated regions that overlap with intron/exon/promoters]\n")
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
png(paste(result_dir,"/",comp_sample, ".DM_annotation.png",sep=""))
plotTargetAnnotation(diffAnn,precedence=TRUE,main="differential methylation annotation")
dev.off()
cat("[percentage of differentially methylated bases are on CpG islands,CpG island shores and other regions]\n")
png(paste(result_dir,"/",comp_sample, ".CpGi_DM_annotation.png",sep=""))
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
		main="differential methylation annotation")
dev.off()
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)

sink()

#overlap promoters with All MR(not DMR)

promoters=regionCounts(methr, gene.obj$promoter)

dmr_pro_file =  paste(result_dir,"/",comp_sample,".MR_overlap_promoters.txt",sep="")

write.table(getData(promoters), dmr_pro_file, quote=FALSE, row.name=FALSE, col.name=FALSE, sep="\t")

#overlap CpGi with All MR(not DMR)

cpg_islands=regionCounts(methr, cpg.obj$CpGi)
dmr_cpg_file = paste(result_dir,"/",comp_sample,".MR_overlap_CpGisland.txt",sep="")
write.table(getData(cpg_islands), dmr_cpg_file, quote=FALSE, row.name=FALSE, col.name=FALSE, sep="\t")

}

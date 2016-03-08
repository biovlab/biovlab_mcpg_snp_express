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
		make_option("--reference", action="store"),
		make_option("--cpgi_reference", action="store"),
		make_option(c("-r", "--result_dir"), action="store"),
		make_option(c("-w", "--window"), action="store", default=100),
		make_option(c("-o", "--commeRbund_result"), action="store"),
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
window = as.numeric(opt$window)

#file read
myobj = read.bismark(
		location=sample_list,
		sample.id=sample.id_list,
		assembly=genome_assembly,
		save.folder=result_dir,
		save.context=NULL,
		read.context=m_context,
		treatment=class_list
		)


#sample statistics

basename_list = lapply(sample_list,basename)
stat_out_list = lapply(basename_list, function(x) paste(result_dir,"/",x,".stat.txt",sep=""))

for(i in 1:length(sample.id_list)){
#	cat(paste("<", sample.id_list[[i]],">","\n", sep=" "))

	sink(stat_out_list[[i]])

#	getMethylationStats(myobj[[i]], plot=FALSE, both.strands=FALSE)
	png(paste(result_dir, "/", basename_list[[i]], ".MethylationStats.png",sep=""))
	getMethylationStats(myobj[[i]], plot=TRUE, both.strands=FALSE)
	dev.off()

#	getCoverageStats(myobj[[i]],plot=FALSE,both.strands=FALSE)
	png(paste(result_dir, "/", basename_list[[i]], ".CoverageStats.png",sep=""))
	getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
	dev.off()

	sink()
}

output_list = lapply(basename_list, function(x) paste(result_dir,"/",x,".bed", sep=""))

for(i in 1:length(sample.id_list)){
	bedgraph(myobj[[i]],file.name=output_list[[i]],'perc.meth')
}

#single methylated region stastics
tiles = tileMethylCounts(myobj, win.size=window, step.size=window)

output_MR_list = lapply(basename_list, function(x) paste(result_dir, "/", x, "_MR.bed", sep=""))
output_MR_coverage_list = lapply(basename_list, function(x) paste(result_dir, "/", x, "_MR_coverage.bed", sep=""))

for(i in 1:length(sample.id_list)){
	bedgraph(tiles[[i]], file.name=output_MR_list[[i]], 'perc.meth')
	bedgraph(tiles[[i]], file.name=output_MR_coverage_list[[i]], 'coverage')
}


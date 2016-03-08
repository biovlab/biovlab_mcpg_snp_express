library(optparse)
library(parallel)

#make options
option_list = list(
		make_option(c("-i", "--sample_file"), action="store"),
		make_option(c("-o","--output"),action="store")
		)

opt = parse_args(OptionParser(option_list=option_list))

#parsing arguments

sample_file = opt$sample_file
output = opt$output

#read sample data
sample_data <- read.table(opt$sample_file,header=F,sep="\t")

#spearman
cat(cor(as.numeric(sample_data[[1]]),as.numeric(sample_data[[2]]),method="spearman"))

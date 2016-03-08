library(optparse)
library(parallel)

#make options
option_list = list(
		make_option(c("-i", "--sample_file"), action="store"),
		make_option(c("-n", "--sample_number"), action="store"),
		make_option(c("--s1"), action="store"),
		make_option(c("--s2"), action="store"),
		make_option(c("-o","--output"),action="store"),
		make_option(c("-p", "--cpu_num"),action="store")
		)

opt = parse_args(OptionParser(option_list=option_list))

#parsing arguments

sample_file = opt$sample_file
sample_num = as.numeric(opt$sample_number)
s1_start = as.numeric(opt$s1)
s2_start = as.numeric(opt$s2)
cpu_num = as.numeric(opt$cpu_num)
output = opt$output

#read sample data
#sample_data <- read.table(opt$sample_file,header=F,sep="\t")
sample_data <- read.table(opt$sample_file,header=T,sep="\t")

#spearman
calc_spearman <- function(x,n, start1, start2){
		return(cor(as.numeric(x[start1:(start1+n-1)]),as.numeric(x[start2:(start2+n-1)]),method="spearman"))
}

#parallel
numWorkers <- cpu_num
cl <- makeCluster(numWorkers, type="PSOCK")

#get spearman's rho
spearman_result <- parRapply(cl, sample_data, calc_spearman, sample_num, s1_start, s2_start)

#make output result
new_sample_data <- cbind(sample_data, spearman_result)

#write output
write.table(new_sample_data, file=output, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

stopCluster(cl)


# input example
#stat	cpg_island_shelf	cpg_island 
#Corr_>_threshold	13.5	12.1
#Corr_<_threshold	66.5	87.9

library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)

input <- args[1]
output <- args[2]
title <- args[3]
flag <- args[4]

data <- read.table(input, sep="\t", header=T)
#re_data <- melt(data, colnames(data)[1])
if (flag == "1"){
	color=c("#CCCCFF", "#FF6666")
}else{
	color=c("#FF6666", "#CCCCFF")
}
jpeg(output)
ggplot(data=data, aes(x=variable, y=value)) + geom_bar(aes(fill=stat), stat="identity") + theme(legend.title=element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1)) + ylab("Percentage (%)") + ggtitle(title) + scale_fill_manual(values=color)
dev.off()

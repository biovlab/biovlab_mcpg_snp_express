library(reshape2)
library(ggplot2)

args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]
ylabel <- args[3]

data <- read.table(input, sep="\t",header=T)

#re_data <- melt(data, colnames(data)[1])
#re_data[,1] <- factor(re_data[,1], levels=re_data[1:100,1], ordered=TRUE)

pdf(output)
ggplot(data=data, aes(x=factor(1), y=counts, fill=factor(Region), stat="identify" ),) + geom_bar(stat="identity") + facet_grid(paste(".~",colnames(data)[3],sep="")) + ylab(ylabel) + xlab("")
dev.off()

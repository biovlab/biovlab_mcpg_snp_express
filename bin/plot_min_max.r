library(reshape2)
library(ggplot2)

args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]

ylim_min <- args[3]
ylim_max <- args[4]

data <- read.table(input, sep="\t",header=T)

re_data <- melt(data, colnames(data)[1])
re_data[,1] <- factor(re_data[,1], levels=re_data[1:100,1], ordered=TRUE)

jpeg(output, width=7, height=5, units="in", res=500)
ggplot(data=re_data, aes(x=relative_range,y=value,colour=variable)) + geom_line(aes(colour=variable,group=variable)) + scale_x_discrete(breaks=c("5p","10%","20%","30%","40%","50%","60%","70%","80%","90%","3p")) + ylim(as.numeric(ylim_min), as.numeric(ylim_max))
dev.off()

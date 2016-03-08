library(reshape2)
library(ggplot2)
library(cowplot)

args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]
ylabel <- args[3]

data <- read.table(input, sep="\t",header=T,quote="\"")

#re_data <- melt(data, colnames(data)[1])
#re_data[,1] <- factor(re_data[,1], levels=re_data[1:100,1], ordered=TRUE)

pdf(output)
ggplot(data=data, aes(x=factor(1), y=counts, fill=factor(Region), stat="identify" ),) + geom_bar(stat="identity") + facet_grid(paste(".~",colnames(data)[3],sep="")) + ylab(ylabel) + xlab("")



#################################################################
########### Require cowplot package ############################
####### Must be generalized for numbers ########################

ggplot(data=data, aes(x=subtype_pair, y=counts, fill=factor(Region), stat="identify" ),) + geom_bar(stat="identity") + coord_flip() + ylab("Number of differentially methylated bins") + xlab("") + theme_classic() -> b
b + scale_fill_manual("Region",values=c("#0F73B7","#D46128","#1D9F77","#CC78A9","#E79E24")) -> b2
b2  + theme(legend.position="left") -> b3
ggdraw(switch_axis_position(b3,'xy')) -> b4
b4

dev.off()

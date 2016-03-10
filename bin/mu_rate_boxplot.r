library(ggplot2)

args <- commandArgs(TRUE)

data <- read.table(args[1], sep="\t",header=T)

data$DNA_methylation_level <- factor(data$DNA_methylation_level,levels=data$DNA_methylation_level)
ggplot(data=data, aes(x=factor(DNA_methylation_level),y=Mutation_rate,fill=factor(Subtype))) + geom_boxplot(outlier.shape=NA) -> m2

split(data, data[,c('Subtype', 'DNA_methylation_level')]) -> split_list

lapply(split_list,function(x) boxplot.stats(x$Mutation_rate)$stats) -> qu


lapply(qu,function(x) x[2]) -> t

min(unlist(t)) -> lo

max(unlist(qu))-> up

m2 + coord_cartesian(ylim=c(lo/2,up*1.05)) -> m3

m3 + theme_classic() + xlab("DNA methylation level (RMS)") + ylab("Mutation rate") -> m4

m4 + guides(fill=guide_legend(title="Class")) -> m5

#m4 + scale_fill_manual("Subtype",values=c("#0F73B7","#E79E24","#D46128")) -> m5


#m5 + theme(legend.position=c(0,1), legend.justification=c(-0.5,1.5), legend.text=element_text(size=12), legend.key.size=unit(1.0,"cm")) -> m6

m5 + theme(axis.text=element_text(size=13),axis.title=element_text(size=15),axis.title.y = element_text(vjust=1)) -> m6

m6 + theme(legend.title=element_text(size=13),legend.text=element_text(size=13), legend.key.size=unit(1.3,"cm")) -> m7 # , legend.position=c(0,1), legend.justification=c(-0.8,1.15)) -> m8

ggsave(args[2], m7, width=10, height=8)

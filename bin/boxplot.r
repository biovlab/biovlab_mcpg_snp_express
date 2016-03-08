library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)

input_list <- args[1]
sample_class_list <- args[2]
sample_class_kind <- args[3]
region_output <- args[4]
cgi_output <- args[5]

data_list <- unlist(strsplit(input_list, ";"))
class_list <- unlist(lapply(unlist(strsplit(sample_class_list, ";")), as.numeric))
class_kind <- unlist(strsplit(sample_class_kind, ";"))

grp <- c()

for(i in 1:length(class_list))
{
	grp <- c(grp, rep(i, class_list[i]))
}

# column header
names <- c("grp","3'UTR","5'UTR","CpG Island","CpGI Shelf", "CpGI Shore", "Exon", "Intron", "Promoter")

# data frame
all_df <- data.frame(grp=grp)

# read all files
for(i in 1:length(data_list))
{
	data <- read.table(data_list[i],sep="\t",header=F)
	all_df <- cbind(all_df, data[,1])
}
colnames(all_df) <- names


all_df.m <- melt(all_df, id="grp")
all_df.m[which(all_df.m$variable=="CpG Island" | all_df.m$variable =="CpGI Shelf" | all_df.m$variable=="CpGI Shore"),] -> all_df.m.cpg
all_df.m[-which(all_df.m$variable=="CpG Island" | all_df.m$variable =="CpGI Shelf" | all_df.m$variable=="CpGI Shore"),] -> all_df.m.genomic

#all_df.m.cpg$Group="CpGI related"

#all_df.m.genomic$Group="Genomic Region related"

all_df.m.cpg$variable <- factor(all_df.m.cpg$variable, levels=c("CpG Island","CpGI Shore","CpGI Shelf"))


# plot for gemoic region
g <- ggplot(data=all_df.m.genomic, aes(x=variable,y=value,fill=factor(grp))) + geom_boxplot(outlier.shape=NA) + scale_fill_discrete(name="Class",labels=class_kind) + theme_classic()
gg <- g + xlab("")+ylab("Number of subtype specific mutation") + theme(axis.text=element_text(size=30),axis.title.y=element_text(size=30,margin=margin(0,20,0,0))) + theme(legend.text=element_text(size=30),legend.title=element_text(size=30), legend.key.size=unit(3.0,"cm")) 
ggsave(region_output ,gg,width=20,height=20)


# plot for gemoic region
g <- ggplot(data=all_df.m.cpg, aes(x=variable,y=value,fill=factor(grp))) + geom_boxplot(outlier.shape=NA) + scale_fill_discrete(name="Class",labels=class_kind) + theme_classic()
gg <- g + xlab("")+ylab("Number of subtype specific mutation") + theme(axis.text=element_text(size=30),axis.title.y=element_text(size=30,margin=margin(0,20,0,0))) + theme(legend.text=element_text(size=30),legend.title=element_text(size=30), legend.key.size=unit(3.0,"cm")) 
ggsave(cgi_output ,gg,width=20,height=20)




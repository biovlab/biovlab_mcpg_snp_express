library(gplots)
library(RColorBrewer)

args <- commandArgs(TRUE)
cat(args)
input <- args[1]
column <- args[2]
range <- eval(parse(text=args[3]))
sample_break <- args[4]
color <- eval(parse(text=args[5]))
out <- args[6]
title <- args[7]


sep_break <- cumsum(as.vector(unlist(lapply(unlist(strsplit(sample_break,",")),as.numeric))))
sep_break <- sep_break[-length(sep_break)]

data <- as.matrix(read.table(input, sep="\t",header=F)[,2:(range+1)])

row_blank <- rep("", nrow(data))
cnames <- as.vector(unlist(strsplit(column,",")))

rownames(data) <- row_blank
colnames(data) <- cnames

#my_p = colorRampPalette(c("blue", "red"))(n=color)

png(filename=out, width = 1480, height = 5000, units="px") 
par(cex.main=2.0)
heatmap.2(data, lwid=c(0.5,4), lhei=c(0.1,4), Rowv=NULL, Colv=NULL, main=title, density.info="none",trace="none", dendrogram="none", scale="row",col=bluered,keysize=1,colsep=sep_break, sepcolor="white",sepwidth=c(0.05,0.05),na.color="black",labRow=NULL,srtRow=NULL,cexRow=3.0, cexCol=3.0)
dev.off()


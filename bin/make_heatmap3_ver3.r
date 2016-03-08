library(heatmap3)
library(RColorBrewer)

args <- commandArgs(TRUE)
input <- args[1]
#column <- args[2]
#range <- eval(parse(text=args[3]))
#sample_break <- args[4]
#color <- eval(parse(text=args[5]))
add_info_name <- args[2]
add_info_name <- "Subtype"

add_info_column <- args[3]
out <- args[4]
title <- args[5]

#ini_start <- eval(parse(text=args[8]))

#sep_break <- cumsum(as.vector(unlist(lapply(unlist(strsplit(sample_break,",")),as.numeric))))
#sep_break <- sep_break[-length(sep_break)]

#data <- as.matrix(read.table(input, sep="\t",header=F)[,(ini_start):(range+ini_start-1)])
#data <- as.matrix(read.table(input, sep="\t", header=T))
data_table <- read.table(input,sep="\t",header=T)

#data <- log(as.matrix(data_table[,-c(1:3)]))
data <- as.matrix(data_table[,-c(1)])

#row_blank <- rep("", nrow(data))
#cnames <- as.vector(unlist(strsplit(column,",")))
#cnames <- colnames(data_table)[-1]
col_blank <- rep("", ncol(data))

#data <- na.omit(data)

#####rainbow(n)

rownames(data) <- data_table[,1]
colnames(data) <- col_blank


add_info_table <- read.table(add_info_column, sep="\t",header=T)

add_info_table_factor_num <- data.frame(as.double(factor(add_info_table[,2])))

factor_count <- nlevels(factor(add_info_table[,2]))

pallet <- topo.colors(factor_count)

make_pallet <- function(x){
	sapply(x, function(y) pallet[y])
}

ColSideColors <- do.call(cbind,lapply(add_info_table_factor_num, make_pallet))
colnames(ColSideColors) <- add_info_name


insert_blank <- function(x,n,elem,foo){
		result <- c()
		temp_result <- c()
				
		for(i in 1:length(x)){
			if((i!=1)&&(i%%n)==1){
				result <- c(result, unique(temp_result), elem)
				temp_result <- c(foo(x[i]))
			}
			else{
				temp_result <- c(temp_result, foo(x[i]))
			}
		}
		result <- c(result,unique(temp_result))
		result
}

colside_label <- insert_blank(add_info_table[,2], length(col_blank), "", function(x) x)
colside_color <- insert_blank(add_info_table_factor_num, length(col_blank), "#FFFFFFFF", function(x) pallet[x])


pdf(out)
pdf("test.pdf")
par(cex.main=1)
heatmap3(data,margin=c(6,12),Rowv=TRUE,Colv=TRUE,scale="row",ColSideColors=ColSideColors,showRowDendro=FALSE)##,main=title)
legend("topright", legend=colside_label, fill=colside_color, border=FALSE, bty="n", y.intersp=0.7, cex=0.7)
dev.off()

#my_p = colorRampPalette(c("blue","white" ,"red"))(n=color)

#png(filename=out, width = 1480, height = 5000, units="px")
#pdf(out, width = 1480, height = 5000)
#par(cex.main=2.0)
#heatmap.2(data, lwid=c(0.5,4), lhei=c(0.1,4), Rowv=TRUE, Colv=TRUE, main=title, density.info="none",trace="none", dendrogram="column", scale="row",col=bluered,keysize=1,labRow=NULL,srtRow=NULL,cexRow=3.0, cexCol=3.0)
#dev.off()


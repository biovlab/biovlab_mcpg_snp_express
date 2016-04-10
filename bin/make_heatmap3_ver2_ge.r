library(heatmap3)
library(RColorBrewer)

args <- commandArgs(TRUE)
input <- args[1]
column <- args[2]
#range <- eval(parse(text=args[3]))
#sample_break <- args[4]
#color <- eval(parse(text=args[5]))
add_info_name <- args[3]
add_info_column <- args[4]
out <- args[5]
title <- args[6]

#ini_start <- eval(parse(text=args[8]))

#sep_break <- cumsum(as.vector(unlist(lapply(unlist(strsplit(sample_break,",")),as.numeric))))
#sep_break <- sep_break[-length(sep_break)]

#data <- as.matrix(read.table(input, sep="\t",header=F)[,(ini_start):(range+ini_start-1)])
#data <- as.matrix(read.table(input, sep="\t", header=T))
data_table <- read.table(input,sep="\t",header=T)

row_blank <- data_table[,1]
#data <- log(as.matrix(data_table[,-c(1:3)]))
data <- as.matrix(data_table[,-c(1)])

#row_blank <- rep("", nrow(data))
cnames <- as.vector(unlist(strsplit(column,",")))
#cnames <- colnames(data_table)[-1]

#data <- na.omit(data)

#####rainbow(n)

rownames(data) <- row_blank
colnames(data) <- cnames

add_info_name <- unlist(strsplit(add_info_name,','))
add_info_column <- unlist(strsplit(add_info_column,','))
add_info_column_num <- as.double(factor(add_info_column))
add_info_column_num_df <- as.data.frame(matrix(add_info_column_num, nrow=length(cnames)))

factor_count <- nlevels(factor(add_info_column))

pallet <- topo.colors(factor_count)

make_pallet <- function(x){
	sapply(x, function(y) pallet[y])
}

ColSideColors <- do.call(cbind,lapply(add_info_column_num_df, make_pallet))
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

colside_label <- insert_blank(add_info_column, length(cnames), "", function(x) x)
colside_color <- insert_blank(add_info_column_num, length(cnames), "#FFFFFFFF", function(x) pallet[x])


max_len <- max(sapply(cnames, nchar))
bottom_margin <- max_len * 5 / 9

if (bottom_margin < 5){
  bottom_margin <- 5
}

png(out)
#par(cex.main=1,mar = c(30,4,4,2) + 0.1)
par(cex.main=1)
heatmap3(data,margins=c(bottom_margin,5),Rowv=TRUE,ColSideColors=ColSideColors,showRowDendro=FALSE,main=title)
legend("topright", legend=colside_label, fill=colside_color, border=FALSE, bty="n", y.intersp=0.7, cex=0.7)
dev.off()


#my_p = colorRampPalette(c("blue","white" ,"red"))(n=color)

#png(filename=out, width = 1480, height = 5000, units="px")
#pdf(out, width = 1480, height = 5000)
#par(cex.main=2.0)
#heatmap.2(data, lwid=c(0.5,4), lhei=c(0.1,4), Rowv=NULL, Colv=TRUE, main=title, density.info="none",trace="none", dendrogram="column", scale="row",col=bluered,keysize=1,labRow=NULL,srtRow=NULL,cexRow=3.0, cexCol=3.0)
#dev.off()


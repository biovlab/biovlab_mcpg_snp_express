library(devtools)
library(ggbiplot)
library(ggplot2)
#library(rgl)

args <- commandArgs(TRUE)
input <- args[1]
sample_list <- unlist(strsplit(args[2], ','))
class_list <- factor(unlist(strsplit(args[3], ',')))
stat <- args[4]
pc_plot <- args[5]
ggbiplot_file <- args[6]
#pc3d_plot <- args[7]

data <-read.table(input, sep="\t", header=T)
id_list <- data[,1]

data <- data[,-1]
#rownames(data) <- id_list

data <- t(data)

#log.data <-log(data)
data <- data[,which(colSums(data)>0)]

data.pca <- prcomp(data, center=TRUE, scale.=TRUE, na.action=na.omit)

sink(stat)
print(data.pca)

cat("\n")

sink()

png(pc_plot)
plot(data.pca, type="l")
dev.off()

sink(stat, append=TRUE)
summary(data.pca)
sink()

png(ggbiplot_file)
g <- ggbiplot(data.pca, obs.scale=0.5, var.scale=1, groups=class_list, circle=TRUE,var.axes=FALSE,labels=sample_list)
g <- g + scale_color_discrete(name='')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)
dev.off()




#theta <- seq(0,2*pi,length.out = 100)
#circle <- data.frame(x = cos(theta), y = sin(theta))

#p <- ggplot(circle,aes(x,y)) + geom_path()

#loadings <- data.frame(data.pca$rotation, .names = row.names(data.pca$rotation))
#p + geom_text(data=loadings, mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) + coord_fixed(ratio=1) + labs(x = "PC1", y = "PC2")
#pdf(loading_plot)
#p
#dev.off()

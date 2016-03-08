library(reshape2)
library(ggplot2)

# execute the following code to create a theme_publish object 
# for example, ggplot(mtcars,aes(mpg,hp,size=wt))+geom_point()+theme_publish()

# note: for discrete colored data, use gray scale: 
# > palette(gray(0:3 / 3))
# > greys <- palette()




theme_publish <- function(base_size = 12) {
  structure(list(
    axis.line =         theme_blank(),
    axis.text.x =       theme_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 1),
    axis.text.y =       theme_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
    axis.ticks =        theme_segment(colour = "black", size = 0.2),
    axis.title.x =      theme_text(size = base_size, vjust = 1),
    axis.title.y =      theme_text(size = base_size, angle = 90, vjust = 0.5),
    axis.ticks.length = unit(0.3, "lines"),
    axis.ticks.margin = unit(0.5, "lines"),
    
    legend.background = theme_rect(colour=NA), 
    legend.key =        theme_rect(colour = "grey80"),
    legend.key.size =   unit(1.2, "lines"),
    legend.text =       theme_text(size = base_size * 0.8,family = "sans"),
    legend.title =      theme_text(size = base_size * 0.8, face = "bold", hjust = 0),
    legend.position =   "right",
    
    panel.background =  theme_rect(fill = "white", colour = NA), 
    panel.border =      theme_rect(fill = NA, colour="grey50"), 
    panel.grid.major =  theme_blank(),
    panel.grid.minor =  theme_blank(),
    panel.margin =      unit(0.25, "lines"),
    
    strip.background =  theme_rect(fill = "grey80", colour = "grey50"), 
    strip.text.x =      theme_text(size = base_size * 0.8),
    strip.text.y =      theme_text(size = base_size * 0.8, angle = -90),
    
    plot.background =   theme_rect(colour = NA),
    plot.title =        theme_text(size = base_size * 1.2),
    plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")
  ), class = "options")
}
args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]
ylabel <- args[3]

data <- read.table(input, sep="\t",header=T)

#re_data <- melt(data, colnames(data)[1])
#re_data[,1] <- factor(re_data[,1], levels=re_data[1:100,1], ordered=TRUE)

pdf(output)
#ggplot(data=data, aes(x=Region, y=counts, fill=factor(Subtype), stat="identify" ),) + geom_bar(position="dodge",stat="identity") + facet_grid(paste(".~",colnames(data)[1],sep="")) + ylab(ylabel) + theme(legend.title=element_blank, axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))
ggplot(data=data, aes(x=Region, y=counts)) + geom_bar(position="dodge", aes(fill=Subtype), stat="identity") + theme(legend.title=element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(size=10,face="bold",angle = 60, hjust = 1)) + ylab(ylabel) 
dev.off()

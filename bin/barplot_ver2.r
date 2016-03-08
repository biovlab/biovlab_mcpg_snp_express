library(ggplot2)
library(cowplot)

args <- commandArgs(TRUE)

input <- args[1]
output <- args[2]
type <- args[3]

if(type == "region"){
	data <- read.table(input, sep="\t", header=TRUE, quote="\"")
	
	ggplot(data=data, aes(x=subtype_pair, y=counts, fill=factor(Region), stat="identify" ),) + geom_bar(stat="identity") + coord_flip() + ylab("Number of differentially methylated bins") + xlab("") + theme_classic() -> fig
	
	fig + scale_fill_manual("Region",values=c("#0F73B7","#D46128","#1D9F77","#CC78A9","#E79E24")) -> fig2
	
	fig2  + theme(legend.position="left",legend.title = element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=15), axis.title.x=element_text(size=15)) -> fig3
	
	ggdraw(switch_axis_position(fig3,'xy')) -> fig4
	
	ggsave(output,fig4, width=10, height=5)
}
else{
	c_data <- read.table(input, sep="\t", header=TRUE, quote="\"")
	
	ggplot(data=c_data, aes(x=subtype_pair, y=counts, fill=factor(Region), stat="identify" ),) + geom_bar(stat="identity") + coord_flip() + ylab("Number of differentially methylated bins") + xlab("") + theme_classic() -> c_fig
	
	c_fig + scale_fill_manual("Region",values=c("#0F73B7","#E79E24","#D46128")) -> c_fig2

	c_fig2 + theme(legend.position="left",legend.title = element_text(size=15),legend.text=element_text(size=15),axis.text = element_text(size=15), axis.title.x=element_text(size=15)) -> c_fig3

	ggdraw(switch_axis_position(c_fig3,'xy')) -> c_fig4

	ggsave(output,c_fig4, width=10, height=5)
}

#!/bin/sh

# variables
old_IFS=$IFS

IFS=';' read -r -a class_list <<< "$1"

IFS=";" read -r -a color_list <<< "$2"

IFS=$old_IFS

echo '''library(ComplexHeatmap)

args <- commandArgs(TRUE)

data_input <- args[1]
class_list <- unlist(strsplit(args[2], ";"))
gene_list <- unlist(strsplit(args[3], ";"))
output <- args[4]

data <- read.table(data_input, header=T, sep="\t")
data[is.na(data)] = ""
rownames(data) <- gene_list
data <- as.matrix(data)

alter_fun_list = list(
	background = function(x, y, w, h) {
		grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
	},
'''

for((i=0;i<${#class_list[@]}-1;i++)); do

echo	"	"${class_list[$i]}''' = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "'''${color_list[$i]}'''", col = NA))
  },
'''
done

echo  " "${class_list[$((${#class_list[@]}-1))]}''' = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "'''${color_list[$((${#class_list[@]}-1))]}'''", col = NA))
  }
)
'''

echo "col = c("

for((i=0;i<${#class_list[@]}-1;i++)); do
		echo '"'${class_list[$i]}'" = "'${color_list[$i]}'",'
done

echo '"'${class_list[$((${#class_list[@]}-1))]}'" = "'${color_list[$((${#class_list[@]}-1))]}'"'
echo ")"


echo '''sample_order=colnames(data)


png(output)
oncoPrint(data, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun_list = alter_fun_list, col = col, 
          row_order=NULL,
          column_order = sample_order, show_column_names =TRUE,
          column_title = "OncoPrint for Test",
          heatmap_legend_param = list(title = "Class", at = class_list, 
                                      labels = class_list))
dev.off()
'''



library(pathview)
library(gage)
library(gageData)
library(gageData)
#library(dplyr)

data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

# variables
# filter kegg pathway by qvalue
#qval_threshold=0.1
qval_threshold=1
top_num=10

# function
convertSymbol2Entrez <- function(x){
  ids <- AnnotationDbi::select(org.Hs.eg.db, x[1], c("ENTREZID","GENENAME"), "ALIAS")$ENTREZID
  fcs <- rep(as.numeric(x[2]), length(ids))
  names(fcs) <- ids
  fcs
}

args <- commandArgs(TRUE)
input_file  <- args[1]
output_suffix <- args[2]
output_dir <- args[3]
output_map_list_file <- args[4]

# input should contains "geneid	foldchange" with header
# read data
data<-read.table(input_file, header=TRUE)

foldchanges <- unlist(apply(data, 1, convertSymbol2Entrez))

# get kegg pathway
keg_maps = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)


print("[INFO] Extract kegg map by genes and foldchange")
print(summary(keg_maps))
print(summary(keg_maps$greater))
print(head(keg_maps$greater))
print(head(keg_maps$greater[,4]))

print("[INFO] Fingering by qvalue")
filtered_kegg_maps <- keg_maps$greater[which(keg_maps$greater[,4] < qval_threshold),]
filtered_kegg_maps <- head(keg_maps$greater) # for temp testing

print("[INFO] Fingered kegg maps")
print(head(filtered_kegg_maps))

filtered_kegg_maps_sorted <- filtered_kegg_maps[order(as.numeric(filtered_kegg_maps[,4])),]
print("[INFO] Sorted kegg maps")
print((filtered_kegg_maps_sorted))

# top N filter
#filtered_kegg_maps_sorted_topN = rownames(data.frame(id=rownames(filtered_kegg_maps_sorted), filtered_kegg_maps_sorted))
filtered_kegg_maps_sorted_topN <- rownames(head(filtered_kegg_maps_sorted, top_num))

print("[INFO] Top N maps")
print(filtered_kegg_maps_sorted_topN)

# Get the IDs.
kegg_maps_ids = substr(filtered_kegg_maps_sorted_topN, start=1, stop=8)

# for testing
print("kegg_maps_ids")
print(kegg_maps_ids)
write(kegg_maps_ids, output_map_list_file)

setwd(output_dir)
# plot pathways
#plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
kegg_template_download_dir="."
pv.out.list = sapply(kegg_maps_ids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", out.suffix = output_suffix)) # NOTE : This will generate mapid.suffix.png files

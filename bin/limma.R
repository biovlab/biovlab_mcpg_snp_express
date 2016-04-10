
library(affy)
library(limma)
library(annotate)
library(optparse)

# make options
option_list <- list( 
    make_option(c("-a", "--sample_list1"), action="store"),
    make_option(c("-b", "--sample_list2"), action="store"),
    make_option(c("-r", "--result_dir"), action="store"),
    make_option(c("-o", "--limma_result"), action="store"),
    make_option(c("-p", "--output_prefix"), action="store"),
		make_option(c("--pvalue", action="store")),
    make_option(c("-l", "--sample_num_list"), action="store")
    #make_option(c("-c", "--count_table"), action="store")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# parsing arguments
sample_list1 <- unlist(strsplit(opt$sample_list1, ";"))
sample_list2 <- unlist(strsplit(opt$sample_list2, ";"))
sample_num_list <- unlist(strsplit(opt$sample_num_list, ";"))
limma_result <- opt$limma_result
output_prefix <- opt$output_prefix
result_dir <- opt$result_dir
pvalue_cut <- as.numeric(opt$pvalue)

# variables
#file1_1="/home/heechae/projects/wd_back/user_data/gene_expression/100730_s_1_export.txt.CEL"
#file1_2="/home/heechae/projects/wd_back/user_data/gene_expression/100730_s_2_export.txt.CEL"
#file2_1="/home/heechae/projects/wd_back/user_data/gene_expression/100730_s_4_export.txt.CEL"
#file2_2="/home/heechae/projects/wd_back/user_data/gene_expression/100730_s_5_export.txt.CEL"
replicate=1
if (length(sample_list1)==1){
	replicate=0
}



cell_file_path=""

#file1
#file2
filenames=c(sample_list1, sample_list2)
#filenames
unlist(filenames)
#c(file1_1,file1_2, file2_1, file2_2)
#data<-ReadAffy(filenames=c(file1_1,file1_2, file2_1, file2_2))
data<-ReadAffy(filenames=c(filenames))
#a<-ReadAffy( filenames="/home/heechae/projects/wd_back/user_data/gene_expression/100730_s_6_export.txt.CEL","/home/heechae/projects/wd_back/user_data/gene_expression/100730_s_7_export.txt.CEL")

cat("[INFO] Data set\n")
data

eset<-rma(data)
cat("[INFO] eset\n")
eset

cat("[INFO] draw box plot\n")
box_fig <-  paste(result_dir, "/", output_prefix, ".norm.boxplot.png", collapse = "", sep="")
png(box_fig)
boxplot(exprs(eset))
dev.off()


# sample names
cat("[INFO] eset sample names\n")
sampleNames(eset)

# load platform annotation
platform=annotation(eset)
mydb=paste(platform,".db", sep="")
cat("[INFO] mydb\n")
mydb
library(mydb,character.only = TRUE)
ID<-featureNames(eset)
gene_symbol <- getSYMBOL(ID,mydb)

# write expression values 
for(i in 1:ncol(eset)){	
	mtemp<-matrix(0,nrow=nrow(eset),ncol=3)
	mtemp[,1]<-featureNames(eset)
	mtemp[,2]<-gene_symbol
	mtemp[,3]<-assayData(eset)[["exprs"]][,i]
	#write.table(mtemp,file=paste(result_dir,sampleNames(eset)[i],".txt",sep=""),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(mtemp,file=paste(result_dir,"/",sampleNames(eset)[i],".exp",sep=""),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
}


# merge gene_symbol to data set
	## Additional informations
	#Name <- as.character(lookUp(ID, mydb, "GENENAME"))
	#Ensembl <- as.character(lookUp(ID, mydb, "ENSEMBL"))
	#tmp <- data.frame(ID=ID, Symbol=gene_symbol, Name=Name, Ensembl=Ensembl, stringsAsFactors=F)
	#tmp <- data.frame(ID=ID, Symbol=gene_symbol, Ensembl=Ensembl, stringsAsFactors=F)
tmp <- data.frame(ID=ID, Symbol=gene_symbol, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA 
fData(eset) <- tmp

cat("[INFO] Chip Platform: ")
platform

cat("[INFO] design\n")
design <- model.matrix(~ 0+factor(sample_num_list))
design
colnames(design) <- c("group1", "group2" )
design

fit <- lmFit(eset, design)
cat("[INFO] fit\n")
fit
#fit
#fit$coefficients
#rownames(fit$coefficients)
#fit$genes$Symbol <- getSYMBOL(fit$genes$ProbeName,annotation)

contrast.matrix <- makeContrasts(group2 - group1, levels=design)
#cat("***contrast.matrix***\n")
#contrast.matrix

fit2 <- contrasts.fit(fit, contrast.matrix)

#fit2

#cat("***log2 foldchange***\n")
#log2fold<-fit2$coefficients

#t=data.frame(log2fold, gene_symbol)
#t
if (replicate==0){
	q()
}

# statistical significance test based on replicates
fit2 <- eBayes(fit2)


cat("[INFO] Store result to files\n")
limma_result1=topTable(fit2, n=Inf, coef=1, adjust="BH")
write.table(limma_result1, file=paste(result_dir,"/",output_prefix,".limma.txt", sep=""), row.names =FALSE, col.names = TRUE, sep ="\t", quote=FALSE)

results <- decideTests(fit2)
results

#draw MA plot
cat("MA plot\n")

# get P.adj with orignal order
limma_result2=topTable(fit2, n=Inf, coef=1, adjust="BH",sort.by="none")

# get deg numbers
degs_num <- length(which(limma_result1$adj.P.Val < pvalue_cut))

degs <- order(limma_result2$adj.P.Val)[1:degs_num]


MA_plot <- paste(result_dir, "/", output_prefix, ".MA_plot.png", collapse = "", sep="")
png(MA_plot)
plotMD(fit2, main=output_prefix)
abline(0,0,col="red")
points(fit2$Amean[degs], fit2$coef[degs], cex=0.8, pch=16, col="red")
dev.off()


# draw volcano plot
volcano <- paste(result_dir, "/", output_prefix, ".volcano.png", collapse = "", sep="")
png(volcano)
volcanoplot(fit2, main=output_prefix)
points(fit2$coef[degs], fit2$lods[degs], cex=0.8, pch=16, col="red")
dev.off()

# print out results
#tops <- topTable(fit2, n=Inf)

# Draw venndiagram
cat("[INFO] Create vennDiagram\n")
result_fig <- paste(result_dir, "/", output_prefix, ".limma_vennDiagram.jpg", collapse = "", sep="") 
png(result_fig)
vennDiagram(results)
dev.off()

# print out results
#cat("[INFO] Store result2 to files\n")
#limma_result2=topTableF(fit2, n=Inf, paste(result_dir,"/limma_result2.txt"))
#write.table(limma_result2, file=paste(result_dir,"/limma_result2.txt", sep=""), row.names =FALSE, col.names = TRUE, sep ="\t", quote=FALSE)


# tops [which(tops$logFC > 0), ] [1:25,] # up reg top 25
#fit2$coefficients
#fit


#contrast.matrix
#
#
#fit2
#
#fit2 <- eBayes(fit2)
#
#fit2

#topTable(fit2, coef=1, adjust="BH")

#results <- decideTests(fit2)

#vennDiagram(results)




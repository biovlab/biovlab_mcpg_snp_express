library(affy)
library(limma)
library(annotate)
library("gplots")
#sample_list1<-unlist(strsplit("/home/heechae/projects/wd_back/user_data/gene_expression/100730_s_1_export.txt.CEL;/home/heechae/projects/wd_back/user_data/gene_expression/100730_s_2_export.txt.CEL", ";"))
#sample_list1
#sample_list2<-unlist(strsplit("/home/heechae/projects/wd_back/user_data/gene_expression/100730_s_4_export.txt.CEL;/home/heechae/projects/wd_back/user_data/gene_expression/100730_s_5_export.txt.CEL", ";"))
#sample_list2
#files=c(sample_list1, sample_list2)
#files
#unlist(files)
#data<-ReadAffy(filenames=c(files))
#data
#eset<-rma(data)
#eset
#sampleNames(eset)
#platform=annotation(eset)
#platform
#ID<-featureNames(eset)
#ID
#mydb=paste(platform,".db", sep="")
#mydb
#gene_symbol <- getSYMBOL(ID,mydb)
#library(mydb,character.only = TRUE)
#eset
#library(mydb,character.only = TRUE)
#mydb
#eset
#fData(eset)
#sample_num_list <- unlist(strsplit("1;1;2;2;", ";"))
#sample_num_list
#design <- model.matrix(~ 0+factor(sample_num_list))
#design
#colnames(design) <- c("group1", "group2" )
#design
#fit <- lmFit(eset, design)
#fit
#topTable(fit)
#head(fit)
#topTable(fit)
#contrast.matrix <- makeContrasts(group2 - group1, levels=design)
#fit2 <- contrasts.fit(fit, contrast.matrix)
#fit2
#fit2 <- eBayes(fit2)
#topTable(fit2)
#topTable(fit2, coef=2)
#topTable(fit2)
#selected  <- p.adjust(fit2$p.value[, 2]) <0.05
#selected<-p.adjust(fit2$p.value[, 2]) <0.05
#fit2$p.value
#fit2$p.value[,2]
#selected  <- p.adjust(fit2$p.value) <0.05
#selected
#esetSel <- eset [selected, ]
#esetSel
#topTable(esetSel)
#heatmap(exprs(esetSel))
#eset
#heatmap(eset)
#heatmap
#topTable(fit2)
#exprs(eset)
#head(exprs(eset))
#heatmap(exprs(eset))
#heatmap(exprs(eset[1:100,]))
#q()
#heatmap(exprs(eset[1:100,]))
#library(affy)
#library(limma)
#library(annotate)
#heatmap(exprs(eset[1:100,]))
#library("gplots")
#heatmap.2(exprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors,
#heatmap.2(exprs(eset[1:100,]), col=redgreen(75), scale="row", ColSideColors=patientcolors,
#key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#heatmap.2(exprs(eset[1:100,]), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#label(eset)
#feature(eset)
#colnames(eset)
#eset
#eset
#patientcolors <- c("#FF0000","#FF0000","#0000FF","#0000FF")
#heatmap.2(exprs(eset[1:100,]), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#heatmap.2(exprs(eset[1:1000,]), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#fit2
#topTable(fit2)
#fit2$adj.P.Val
#colnames(fit2)
#pvalue(fit2)
#p-value(fit2)
#fit2$pvalue
#fit2$p-value
#fit2$p.value
#head(fit2$p.value)
#topTable(fit2)
#head(fit2$q.value)
#head(fit2$p.value<0.05)
#rownames(fit2)
#attr(fit2)
#attrib(fit2)
#fit2
#selected<-head(fit2$p.value<0.05)
#esetSel[eset[selected,]]
#esetSel<-[eset[selected,]]
#esetSel<-eset[selected,]
#esetSel
#heatmap.2(exprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#heatmap.2(exprs(xprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#heatmap.2(exprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#esetSel
#colnames(esetSel)
#rownames(esetSel)
#head(esetSel)
#heatmap.2(exprs(esetSel[1:100,]), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#heatmap.2(exprs(esetSel[1:100,]), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#fit2
#topTable(eset)
#topTable(fit2)
#topTable(fit2,coef=1, adjust="BH")
# results <- decideTests(fit2)
#vennDiagram(results)
#results
#cel_lu1=/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_1_export.txt.CEL
cel_lu1="/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_1_export.txt.CEL"
cel_lu2="/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_7_export.txt.CEL"
cel_lu3="/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_8_export.txt.CEL"
cel_lu4="/data/project/mcpg/test_data/icbp/gene_expression_array/100803_s_5_export.txt.CEL"
cel_lu5="/data/project/mcpg/test_data/icbp/gene_expression_array/100803_s_6_export.txt.CEL"
cel_lu6="/data/project/mcpg/test_data/icbp/gene_expression_array/100812_s_3_export.txt.CEL"
cel_lu7="/data/project/mcpg/test_data/icbp/gene_expression_array/100812_s_4_export.txt.CEL"
cel_lu8="/data/project/mcpg/test_data/icbp/gene_expression_array/100824_s_2_export.txt.CEL"
cel_lu9="/data/project/mcpg/test_data/icbp/gene_expression_array/100831_s_5_export.txt.CEL"
cel_lu10="/data/project/mcpg/test_data/icbp/gene_expression_array/100902_s_7_export.txt.CEL"
cel_lu11="/data/project/mcpg/test_data/icbp/gene_expression_array/100908_s_3_export.txt.CEL"
cel_lu12="/data/project/mcpg/test_data/icbp/gene_expression_array/100908_s_4_export.txt.CEL"
cel_lu13="/data/project/mcpg/test_data/icbp/gene_expression_array/100908_s_8_export.txt.CEL"
cel_baa1="/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_5_export.txt.CEL"
cel_baa2="/data/project/mcpg/test_data/icbp/gene_expression_array/100831_s_6_export.txt.CEL"
cel_baa3="/data/project/mcpg/test_data/icbp/gene_expression_array/100831_s_8_export.txt.CEL"
cel_baa4="/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_4_export.txt.CEL"
cel_baa5="/data/project/mcpg/test_data/icbp/gene_expression_array/100812_s_2_export.txt.CEL"
cel_baa6="/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_6_export.txt.CEL"
cel_baa7="/data/project/mcpg/test_data/icbp/gene_expression_array/100908_s_5_export.txt.CEL"
cel_bab1="/data/project/mcpg/test_data/icbp/gene_expression_array/100730_s_2_export.txt.CEL"
cel_bab2="/data/project/mcpg/test_data/icbp/gene_expression_array/100812_s_1_export.txt.CEL"
cel_bab3="/data/project/mcpg/test_data/icbp/gene_expression_array/100803_s_3_export.txt.CEL"
cel_bab4="/data/project/mcpg/test_data/icbp/gene_expression_array/100803_s_7_export.txt.CEL"
cel_bab5="/data/project/mcpg/test_data/icbp/gene_expression_array/100812_s_5_export.txt.CEL"
cel_bab6="/data/project/mcpg/test_data/icbp/gene_expression_array/100824_s_4_export.txt.CEL"
cel_bab7="/data/project/mcpg/test_data/icbp/gene_expression_array/100831_s_4_export.txt.CEL"
cel_bab8="/data/project/mcpg/test_data/icbp/gene_expression_array/100902_s_6_export.txt.CEL"
cel_bab9="/data/project/mcpg/test_data/icbp/gene_expression_array/100908_s_6_export.txt.CEL"
cel_bab10="/data/project/mcpg/test_data/icbp/gene_expression_array/100910_s_4_export.txt.CEL"

filenames=c(cel_lu1, cel_lu2, cel_lu3, cel_lu4, cel_lu5, cel_lu6, cel_lu7, cel_lu8, cel_lu9, cel_lu10, cel_lu11, cel_lu12, cel_lu13, cel_baa1, cel_baa2, cel_baa3, cel_baa4, cel_baa5, cel_baa6, cel_baa7, cel_bab1, cel_bab2, cel_bab3, cel_bab4, cel_bab5, cel_bab6, cel_bab7, cel_bab8, cel_bab9, cel_bab10)
data<-ReadAffy(filenames=c(filenames))
eset<-rma(data)
sample_num_list<-c("lu", "lu", "lu", "lu", "lu", "lu", "lu", "lu", "lu", "lu", "lu", "lu", "lu", "baa", "baa", "baa", "baa", "baa", "baa", "baa", "bab", "bab", "bab", "bab", "bab", "bab", "bab", "bab", "bab", "bab")
design <- model.matrix(~ 0+factor(sample_num_list))

design
colnames(design) <- c("lu", "baa", "bab" )
fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(baa-lu, bab-baa, bab-lu, levels=design)
contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)
#vennDiagram(results)

#selected<-fit2$p.value<0.05
#selected
#selected
#selected
#selected[,2]
#selected[2,]
selected <-apply(fit2$p.value[,1:3], 1, function(x) any(x<0.05))
selected
esetSel<-eset[selected,]
#selected  <- p.adjust(fit2$p.value[, 2]) <0.05
selected
head(fit2)
#heatmap.2(exprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#patientcolors
patientcolors<-c("#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#00FF00", "#00FF00", "#00FF00", "#00FF00", "#00FF00", "#00FF00", "#00FF00", "#00FF00", "#00FF00", "#00FF00")
heatmap.2(exprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)

#heatmap.2(exprs(eset[1:100,]), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#heatmap.2(exprs(eset[1:1000,]), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#q()
#heatmap.2(exprs(eset[1:1000,]), col=redgreen(75), scale="row", ColSideColors=patientcolors,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
#savehistory(file = ".Rhistory")
#q()

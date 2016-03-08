library(limma)

mat <- matrix(c(5.098, 5.076, 5.072, 4.677,8.906,8.903,6.700,6.653,7.409,7.498,5.392,5.432), nrow=3, byrow=TRUE)

Group <- factor(c("p1","p1","p2","p2"))

design <- model.matrix(~0+Group)

colnames(design) <- gsub("Group","",colnames(design))

fit <- lmFit(mat,design)

contrast.matrix<-makeContrasts(p1-p2,levels=design)

fit2 <- contrasts.fit(fit,contrast.matrix)

fit2<-eBayes(fit2)

tab <- topTable(fit2)

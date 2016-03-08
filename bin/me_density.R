library(hexbin)
library(RColorBrewer)
library(gplots)
library(ggplot2)

# color setting
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

args <- commandArgs(TRUE)

input1 <- args[1]
input2 <- args[2]
class1 <- args[3]
class2 <- args[4]
max_val <- as.numeric(args[5])
step_val <- as.numeric(args[6])
output <- args[7]

# max setting
setting_max = function(x){ if (x >= max_val) max_val else x}


data1 <- read.table(input1, header=T, sep="\t")
data2 <- read.table(input2, header=T, sep="\t")

merged_data <- cbind(data1[,1],data2[,1])
colnames(merged_data) <- c(class1, class2)

merged_data <- merged_data[complete.cases(merged_data),]

merged_data_modified <- structure(vapply(merged_data, setting_max, numeric(1)),dim=dim(merged_data))
colnames(merged_data_modified) <- c(class1, class2)

merged_data_modified <- as.data.frame(merged_data_modified)

nbins <- 25
x.bin <- seq(min(merged_data_modified[,1]), max(merged_data_modified[,1]), length=nbins)
y.bin <- seq(min(merged_data_modified[,2]), max(merged_data_modified[,2]), length=nbins)

freq <-  as.data.frame(table(findInterval(merged_data_modified[,1], x.bin),findInterval(merged_data_modified[,2], y.bin)))
freq[,1] <- as.numeric(freq[,1])
freq[,2] <- as.numeric(freq[,2])

freq2D <- diag(nbins)*0
freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]

df <- expand.grid(x=x.bin, y=y.bin)
df$counts <- freq[,3]
colnames(df) <- c(class1, class2,"counts")


g <- ggplot(df, aes_string(class1, class2)) +
  geom_raster(aes(fill =counts), interpolate = TRUE) + theme_classic() + scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_distiller(palette = "Spectral",trans="log10", na.value="#3386BC") +
  geom_vline(xintercept=seq(0,max_val,by=step_val), linetype="longdash", colour="white") +
  geom_hline(yintercept=seq(0,max_val,by=step_val), linetype="longdash", colour="white") +
  geom_abline(intercept=0,slope=1,colour="white")

ggsave(output, g, width=10, height=8)

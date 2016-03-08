library(optparse)
library(parallel)

#make options
option_list = list(
		make_option(c("-i", "--sample_file"), action="store"),
		make_option(c("-c", "--class_list"), action="store"),
		make_option(c("-m", "--mode"), action="store"),
		make_option(c("--o1"), action="store"),
		make_option(c("--o2"), action="store"),
		make_option(c("-n", "--cpu_num"),action="store"),
		make_option(c("-a", "--alpha"), action="store", default=0.05)
		)

opt = parse_args(OptionParser(option_list=option_list))

#parsing arguments

sample_file = opt$sample_file
sample_file_name = basename(sample_file)
group <- factor(unlist(lapply(unlist(strsplit(opt$class_list, ",")), function(x) x)))
mode = opt$mode
alpha <- as.numeric(opt$alpha)
cpu_num <- as.numeric(opt$cpu_num)
#output_dir <- opt$output_dir

#output <- paste(output_dir, sample_file_name, ".kw_test.txt",sep="")
#output2 <- paste(output_dir, sample_file_name, ".kw_test.signif.cut", alpha,".txt",sep="")

output <- opt$o1
output2 <- opt$o2

#group
#(group2<-factor(c(rep(1,10),rep(2,10),rep(3,10))))


#read sample data
sample_data <- read.table(opt$sample_file,header=F,sep="\t")
#sample_data <- cbind(sample_data, c(rep(opt$class_list, nrow(sample_data))))

#make kruskal-wallis test  // maked for all difference
kw_test <- function(x,group, alpha,mode){
	#value=unlist(x)
	value <- x
#	value <- x[-length(x)]
#	group <- factor(unlist(lapply(unlist(strsplit(x[length(x)], ",")), function(x) x)))
	data <- data.frame(group,value)
	kw <- kruskal.test(value ~ group, data=data)
	
	result <- ""

	if(is.nan(kw$p.value) | kw$p.value >= alpha){
		if(mode == 0){
			result <- paste(kw$p.value,0)
		}else{
			result <- paste(kw$p.value,0,".")
		}
	}else{
		post_hoc <- pairwise.wilcox.test(value, group, p.adj="bonferroni", exact=F)

		flag <- TRUE

	 	if(mode ==0){
			func <- function(x) return (is.nan(x) | (x >= alpha))
			flag <- Reduce("|", sapply(post_hoc$p.value[row(post_hoc$p.value) >= col(post_hoc$p.value)], func))
			
			if(flag)
			{
				result <- paste(kw$p.value, 0)
			}else{
				result <- paste(kw$p.value, 1)
			}
		}else{
			save_flag <- FALSE
			class <- ""
			for(k in 1:(length(levels(group)))){
				i <- 1
				j <- (k-1)
				l <- i
				temp_flag <- FALSE
				if(l <= j){
					for(n in l:j){
						if(is.nan(post_hoc$p.value[j,n]) | (post_hoc$p.value[j,n] >= alpha)){
							temp_flag <- TRUE
							break
						}
						l <- n
					}
					l <- l+1
				}
				if(temp_flag){
					next
				}
				j <- j+1
				
				if(j <= (length(levels(group))-1)){
					for(m in j:(length(levels(group))-1)){
						if(is.nan(post_hoc$p.value[m,l]) | (post_hoc$p.value[m,l] >= alpha)){
						temp_flag <- TRUE
						break
						}
					}
				}

				if(!temp_flag)
				{
					if(save_flag){
						flag <- TRUE
						break
					}else{
						save_flag <- TRUE
						class <- levels(group)[k]
						flag <- FALSE
					}
				}

			}

			if(flag){
				result <- paste(kw$p.value, 0, ".")
			}else{
				result <- paste(kw$p.value, 1, class)
			}
		}	
	}
	return (result)
}

r_split <- function(x,mode){
	result <- matrix()
	if(mode==0){
		result <- matrix(unlist(strsplit(x," ")),ncol=2,byrow=T)
	}else{
		result <- matrix(unlist(strsplit(x," ")),ncol=3,byrow=T)
	}

	return (result)
}
		
#data <- data.frame(group,value=unlist(sample_data[1,-1]))
numWorkers <- cpu_num
cl <- makeCluster(numWorkers, type="PSOCK")
#print(data)
#if(mode == 0) # find all different
#{

#	kw_test_result <- as.data.frame(t(parRapply(cl,sample_data[,-1],kw_test,group,alpha)))
	kw_test_result <- parRapply(cl,sample_data[,-1],kw_test,group,alpha,mode)
#	kw_test_result <- as.data.frame(t(apply(sample_data[,-1],1,kw_test)))
	kw_test_result <- r_split(kw_test_result,mode)
#	kw_test_result <- parRapply(cl,kw_test_result,r_split)	
	if(mode==0){
		colnames(kw_test_result) <- c("p-value","flag")
	}else{
		colnames(kw_test_result) <- c("p-value","flag","class")
	}
	new_sample_data <- cbind(sample_data, kw_test_result)
	write.table(new_sample_data, file=output,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
#	subset_data <- new_sample_data[which(new_sample_data$flag == 1)]
#	subset_data <- subset(m, m[,4]>=2)
	subset_data <- subset(new_sample_data, new_sample_data$flag == 1)

	write.table(subset_data, file=output2, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
#	if(mode==0){
#	}else{
#		write.table(subset_data[,-(ncol(subset_data-1))],file=output2,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
#	}
#}else{
#}

stopCluster(cl)

#/bin/sh

####################################################################################
############# 1. methylation level range is ambiguous   #############################
# 2. Need to refine heatmap image
#####################################################################################


#$1 : input file (without header, only data, first column : chr_start_end , second~ : data)
#$2 : class label list sep by "," (ex "br1,br2,br3,co1,co2,co3")
#$3 : mode (0 or 1) (ex 0 : all diff, 1 : 1 class diff)
#$4 : alpha (p-value threshold)
#$5 : cpu number 
#$6 : output file 1 for all data (path+name)
#$7 : output file 2 for only significant result in mode (path+name)


#$8 : column id file by 1 column
#(Ex. col_id.txt)
# br1
# br2
# br3
# co1
# co2

#$9 : # of samples
#$10 : # of classes (ex "13,7,10" or "5,5,5,5,5")
#####################################$11 : methylation level range (ex. if level 0~100, then 100)  ######################################
#$12 : heatmap png output file (path+name)
#$13 : heatamp title


#run kruskal wallis test
Rscript p_kruskal.r -i $1 -c $2 -m $3 -a $4 -n $5 --o1 $6 --o2 $7

#make heatmap
Rscript make_heatmap_flex.r $6 $8 $9 $10 $11 $12 $13

	

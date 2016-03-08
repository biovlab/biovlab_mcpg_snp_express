#!/bin/sh

#$1 : input file
#$2 : sample numbers sep by "," (ex, 13,7,10 / 5,5,5,5,5,5)
#$3 : pseudo count
#$4 : output file

awk -v str="$2" 'BEGIN{class_num=split(str,numbers,",")} {total=0; arr[class_num]; for(i=1; i<=class_num; i++){arr[i] = 0}; c=1; margin=1+numbers[c]; for(i=2;i<=NF;i++){ arr[c]=arr[c]+$(i)+ps; if(i==margin && i!=NF){c=c+1;margin=margin+numbers[c]}} ; for(i=1;i<=class_num;i++){arr[i]=arr[i]/numbers[i]; total=total+arr[i]}; ent = 0; for(i=1;i<=class_num;i++){ if(arr[i] != 0) ent=ent-((arr[i]/total)*(log(arr[i]/total)/log(2)))}; print $0"\t"(ent/(log(class_num)/log(2)))}' ps="$3" $1 > $4

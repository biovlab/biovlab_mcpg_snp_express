#!/usr/bin/awk -f

BEGIN { while ((getline < input_file ) > 0) data[$1] = $2 }
	#NR==1{print header; next;}{ 
	NR==1{next;}{ 
		#print temp_iter;
		#print NF;
		gene_symbol_list=$(NF-2);
		split(gene_symbol_list, tokens, ","); 

		# remove last ","
		if(tokens[length(tokens)]=="") {
			tokens_length=length(tokens)-1;
		} 
		else {
			tokens_length=length(tokens);
		} 

		#for (j=1; j<=tokens_length; j++){ print tokens[j]}

		# for each gene symbol store to array
		for (j=1; j<=tokens_length; j++){
			if (j==1) {
				if (data[tokens[j]]) {printf "%s", data[tokens[j]];} else {printf "%s", "n/a";};
			}
			else{
				if (data[tokens[j]]) {printf "%s", ","data[tokens[j]];} else {printf "%s", ",n/a";};
			}
		}
		print ""
	}  

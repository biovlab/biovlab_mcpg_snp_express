#!/usr/bin/awk -f

## Calculates Pearson Correlation Coefficient in one-pass
## pseudocode from http://en.wikipedia.org/wiki/Talk:Correlation_and_dependence

## Input: 2 cols from stdin
{
	sum_sq_x = 0
	sum_sq_y = 0
	sum_coproduct = 0
	N = 0
	mean_x = 0
	mean_y = 0
	
	counter=0
	for (i=1; i<=sample_num; i++){

		if (i==1) counter=0
		first=$(first_value_column+i-1)
		#print first
		second=$(NF - sample_num + i)
		# now only consider first expression value for multiple gene expresion value e.g : value,value 
		split($(NF - sample_num + i), tokens, ",")
		if (length(tokens) > 1) second=tokens[1]

		#print second
		# if value is not available, skip
		if (first == "-") continue
		if (second == "n/a") continue

		counter++
		N++

		if(mean_x==0){
			mean_x = first
			mean_y = second
		}else{
			sweep = (counter - 1.0) /counter 
			delta_x = first - mean_x
			delta_y = second - mean_y
			sum_sq_x += delta_x * delta_x * sweep
			sum_sq_y += delta_y * delta_y * sweep
			sum_coproduct += delta_x * delta_y * sweep
			mean_x += delta_x / counter
			mean_y += delta_y / counter
		}
	}


        if(N > 0){
		#print sum_sq_x
		#print sum_sq_y
		#print "here"
		pop_sd_x = sqrt( sum_sq_x/N )
		#print pop_sd_x
		pop_sd_y = sqrt( sum_sq_y/N )
		cov_x_y = sum_coproduct/N
		#print pop_sd_y

	}
        if(pop_sd_x * pop_sd_y != 0){
		correlation = cov_x_y / (pop_sd_x * pop_sd_y)
		print correlation
	}else{
	    	print "-"
	}
}

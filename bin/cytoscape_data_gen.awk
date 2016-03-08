function round(x,   ival, aval, fraction)
{
   ival = int(x)    # integer part, int() truncates

   # see if fractional part
   if (ival == x)   # no fraction
      return ival   # ensure no decimals

   if (x < 0) {
      aval = -x     # absolute value
      ival = int(aval)
      fraction = aval - ival
      if (fraction >= .5)
         return int(x) - 1   # -2.5 --> -3
      else
         return int(x)       # -2.3 --> -2
   } else {
      fraction = x - ival
      if (fraction >= .5)
         return ival + 1
      else
         return ival
   }
}

function value_ratio(start, N, sum,arr,   i,temp_sum,s)
{
  for(i=1;i<N;i++)
  {
    arr[i]=round($(start+i-1)*100/sum);
    temp_sum+=arr[i];
  };

  arr[N]=100-temp_sum;
	s=arr[1];
	for(i=2;i<=N;i++)
	{	
		s=s","arr[i];
	}
	return s;
}


BEGIN{
	FS=OFS="\t";
	ii=0;
}
{ 
	methyl_sum = 0.0;
	for(i=2;i<sample_num+2;i++)
	{
		methyl_sum += $(i);
	};
		
	exp_sum = 0.0;
	for(i=sample_num+2; i<2*sample_num+2; i++)
	{
		exp_sum += $(i);
	};

	if (!($(2*sample_num +2) in tf_arr))
	{
		tf_sum = 0.0;
		for(i=2*sample_num+3;i<3*sample_num+3;i++)
		{
			tf_sum += $(i);
		};
		
		temp_string=value_ratio(2*sample_num+3, sample_num, tf_sum, temp_arr);
		tf_arr[$(2*sample_num +2)] = temp_string;
	};

	temp_methyl_string=value_ratio(2,sample_num, methyl_sum, methyl_temp_arr);
	temp_exp_string=value_ratio(sample_num+2,sample_num, exp_sum, exp_temp_arr);

	deg_arr[$1]=temp_exp_string;
	dmr_arr[$(1)"_Promoter"]=temp_methyl_string;

	temp_edge_arr= $(1)"_Promoter,"$1;

	edge_arr[ii]= temp_edge_arr;
	ii++;

	temp_edge_arr = $(2*sample_num+2)","$(1)"_Promoter";
	
	edge_arr[ii]= temp_edge_arr;
	ii++;
}
END{
	print "var elements = {";
	print "\tnodes: [";
	for(elem in deg_arr)
	{
		printf "\t\t{ data: { id: \"%s\", name: \"%s\", type: \"deg\", ", elem, elem;
		split(deg_arr[elem], temp_arr,",");
		for(i=1;i<sample_num;i++)
		{
			printf "c%d: %d, ", i, temp_arr[i];
		};
		printf "c%d: %d } },", sample_num, temp_arr[sample_num];
		printf "\n";
	};
	printf "\n";
				
	for(elem in dmr_arr)
	{
		printf "\t\t{ data: { id: \"%s\", name: \"%s\", type: \"dmr\", ", elem, elem;
		split(dmr_arr[elem], temp_arr, ",");
		for(i=1;i<sample_num;i++)
		{
			printf "c%d: %d, ", i, temp_arr[i];
		};
		printf "c%d: %d } },", sample_num, temp_arr[sample_num];
		printf "\n";
	};
	printf "\n";
				
	for(elem in tf_arr)
	{
		printf "\t\t{ data: { id: \"%s\", name: \"%s\", type: \"tf\", ", elem, elem;
		split(tf_arr[elem], temp_arr, ",");
		for(i=1;i<sample_num;i++)
		{
			printf "c%d: %d, ", i, temp_arr[i];
		};
		printf "c%d: %d } },", sample_num, temp_arr[sample_num];
		printf "\n";
	};				
	print "\t],";
	print "edges: [";

	for(elem in edge_arr)
	{
		split(edge_arr[elem], temp_arr, ",");
		printf "\t\t{ data: { id: \"%s-%s\", source: \"%s\", target: \"%s\" }},", temp_arr[1], temp_arr[2], temp_arr[1], temp_arr[2];
		printf "\n";
	}
	print "\t]"
	print "};"
}


	
	

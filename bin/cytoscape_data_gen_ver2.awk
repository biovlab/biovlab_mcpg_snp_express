BEGIN{
	FS=OFS="\t";
	ii=0;
	tf_index=1;
	deg_index=(2+class_num);
	promoter_index=(3+2*class_num);
}
{ 
	if(!($(tf_index) in tf_arr))
	{
		temp_string=$(tf_index+1);
		for(i=2;i<=class_num;i++)
		{
			temp_string=temp_string","$(tf_index+i)
		};
		tf_arr[$(tf_index)]=temp_string
	};
	if(!($(deg_index) in deg_arr))
	{
		temp_string=$(deg_index+1);
		for(i=2;i<=class_num;i++)
		{
			temp_string=temp_string","$(deg_index+i)
		};
		deg_arr[$(deg_index)]=temp_string;
		
		temp_string=$(promoter_index+1);
		for(i=2;i<=class_num;i++)
		{
			temp_string=temp_string","$(promoter_index+i)
		};
		promoter_arr[$(promoter_index)]=temp_string;
	
	
		temp_edge_arr= $(promoter_index)","$(deg_index);

		edge_arr[ii]= temp_edge_arr;
		ii++;

	};
	
	temp_edge_arr = $(tf_index)","$(promoter_index);
	
	edge_arr[ii]= temp_edge_arr;
	ii++;
}
END{
	print "{\n\t\"elements\" : {";
	print "\t\t\"nodes\" : [";
	for(elem in deg_arr)
	{
		printf "\t\t\t{ \"data\" : { \"id\" : \"%s\", \"name\" : \"%s\", \"type\" : \"deg\", ", elem, elem;
		split(deg_arr[elem], temp_arr,",");
		for(i=1;i<class_num;i++)
		{
			printf "\"c%d\" : %d, ", i, temp_arr[i];
		};
		printf "\"c%d\" : %d } },", class_num, temp_arr[class_num];
		printf "\n";
	};
	printf "\n";
				
	for(elem in promoter_arr)
	{
		printf "\t\t\t{ \"data\" : { \"id\" : \"%s\", \"name\" : \"%s\", \"type\" : \"dmr\", ", elem, elem;
		split(promoter_arr[elem], temp_arr, ",");
		for(i=1;i<class_num;i++)
		{
			printf "\"c%d\": %d, ", i, temp_arr[i];
		};
		printf "\"c%d\": %d } },", class_num, temp_arr[class_num];
		printf "\n";
	};
	printf "\n";
	tf_len=length(tf_arr);
	my_index=1;				
	for(elem in tf_arr)
	{
		printf "\t\t\t{ \"data\" : { \"id\" : \"%s\", \"name\" : \"%s\", \"type\" : \"tf\", ", elem, elem;
		split(tf_arr[elem], temp_arr, ",");
		for(i=1;i<class_num;i++)
		{
			printf "\"c%d\": %d, ", i, temp_arr[i];
		};
		if(my_index < tf_len){
			printf "\"c%d\": %d } },", class_num, temp_arr[class_num];
		}else{
			printf "\"c%d\": %d } }", class_num, temp_arr[class_num];
		};
		printf "\n";
		my_index=my_index+1;
	};				
	print "\t\t],";
	print "\t\t\"edges\": [";
	edge_len=length(edge_arr);
	my_index=1;
	for(elem in edge_arr)
	{
		if(my_index < edge_len){
			split(edge_arr[elem], temp_arr, ",");
			printf "\t\t\t{ \"data\": { \"id\": \"%s-%s\", \"source\" : \"%s\", \"target\" : \"%s\" }},", temp_arr[1], temp_arr[2], temp_arr[1], temp_arr[2];
			printf "\n";
		} else{
			split(edge_arr[elem], temp_arr, ",");
			printf "\t\t\t{ \"data\": { \"id\": \"%s-%s\", \"source\" : \"%s\", \"target\" : \"%s\" }}", temp_arr[1], temp_arr[2], temp_arr[1], temp_arr[2];
      printf "\n";
		};
		my_index=my_index+1;
	}
	print "\t\t]";
	print "\t}";
	print "}"
}


	
	

BEGIN{
	FS=OFS="\t"
}
{
	printf "{ data: { id: \"%s\", name : \"%s\", type: \"%s\", ";




print "var elements = {";
	print "\tnodes: [";
	for(elem in deg_arr)
	{
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
	

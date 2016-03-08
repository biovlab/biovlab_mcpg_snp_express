#!bin/bash

gene_list_file=$1


name_only=`basename $1`


# create html file
echo '
<html>
<head>
</head>
<body>
	<table border="1" align="center">
	<tr><th colspan="2">'$name_only'</th></tr>
	<tr><th width="100">Gene</th><th width="200">Annotaion link</th></tr>'
	for gene in `sort "$gene_list_file" `
	do 
		echo "<tr><td align='center'>$gene</td><td align='center'><a href='http://www.ncbi.nlm.nih.gov/nuccore/$gene'>View</a></td></tr>"
	done
	echo '</table>
</body>
</html>'




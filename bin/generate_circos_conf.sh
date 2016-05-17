#!/bin/bash
source `dirname $0`/../env.sh

# variables
class_num=$1
data_file_list=$2
diff_file=$3
gene_file=$4
circos_config_loc=$5
min=$6
max=$7
final_result_url=$8

# circos parameters
#min="0"
#max="10"
top="0.98"

#diff_min="-0.05"
diff_min=$8
#diff_max="0.05"
diff_max=$9

#interval=`echo - | awk -v cl_num=$class_num '{print 0.4/cl_num}'`
#interval=`echo - | awk -v cl_num=$class_num '{temp=0.4/cl_num; printf("%.2f\n", temp)}'`
interval=`echo - | awk -v cl_num=$class_num '{printf("%.2f\n", 0.4/cl_num)}'`

#echo $interval

# directories
bin_dir="$WORK_DIR/bin"




# set default parameters
echo "
karyotype = data/karyotype/karyotype.human.txt
chromosomes_display_default = yes
#chromosomes                 = hs1
chromosomes_units = 1000000
<links>
<link>
</link>
</links>

#<<include $circos_config_loc/ideogram.conf>>
<ideogram>

<spacing>
default = 0.005r
</spacing>

radius           = 0.90r
thickness        = 20p
fill             = yes
stroke_color     = dgrey
stroke_thickness = 2p

show_label       = yes

label_font       = default
label_radius     = dims(image,radius) - 60p
label_size       = 30
label_parallel   = yes
ideogram_url = $final_result_dir_inte/circos_[chr].html
#ideogram_url = http://147.46.15.115/temp/circos.html

</ideogram>

<<include $circos_config_loc/ticks.conf>>

<image>
<<include etc/image.conf>>                
image_map_use      = yes
image_map_name     = circosmap
</image>
<plots>
"
# set each plots

IFS=';' read -ra data_files <<< "$data_file_list"
for data_file in "${data_files[@]}"; do
	echo "
		<plot>
			<<include $circos_config_loc/colors.ucsc.conf>>
			type = histogram
			fill             = yes
			stroke_type = outline
			stroke_thickness = 1p
			thickness   = 3
			color=lgrey
			extend_bin  = no
			#fill_color = blue
			file = "$data_file

# compute radius of each band for each data 
r1=$top"r"
r0=`echo - | awk -v temp=$top -v intv=$interval '{print temp-intv+0.02"r"}'`
top=`echo - | awk -v temp=$top -v intv=$interval '{print temp-intv}'`

	echo -n "
			r1   = $r1
			r0   = $r0
			min=$min
			max=$max
			<rules>
		# This is filter rule
		#	<rule>
		#	condition = var(value) > 0.10
		#	color = red
		#	</rule>
			<rule>
			condition  = 1
			fill_color = eval(qw(vvlblue vlblue lblue blue dblue vdblue)[remap_int(var(value),$min,$max,0,5)])
			</rule>
			</rules>

			<<include $circos_config_loc/axes.conf>>
			#<<include $circos_config_loc/backgrounds.conf>>
		</plot>
	"
done
echo "
<plot>
	type    = heatmap
	 r1 = 0.58r
	 r0 = 0.55r
	file = "$diff_file
echo "
#	color  = rdylbu-9-div-rev
	color  = reds-7-seq
	stroke_thickness = 0
	min              = "$diff_min
echo "
	max              = "$diff_max
echo "
	<rules>
	<rule>
	condition  = 1
	fill_color = eval(qw(vvlblue vlblue lblue blue dblue vdblue)[remap_int(var(value),"$diff_min","$diff_max",0,5)])
	#condition = var(value) > 0.10
	#color = red
	</rule>
	</rules>


</plot>

<plot>
type             = text
color            = black
file             = "$gene_file
echo "
# inside circle
 # r0 = 0.4r
 # r1 = 0.8r

 # on tick scale
 r0 = 0.30r
 r1 = 0.48r

 show_links     = yes
 link_dims      = 0p,0p,50p,0p,10p
 #link_dims      = 4p,4p,8p,4p,4p
 link_thickness = 2p
 link_color     = red

 label_size   = 24p
 label_font   = condensed

 padding  = 0p
 rpadding = 0p

 </plot>



</plots>



<<include etc/colors_fonts_patterns.conf>> 
<<include etc/housekeeping.conf>>

"




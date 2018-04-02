#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
@result_files = ();

#$ENV{rsat_echo}=1;

### Read the CGI query
$query = new CGI;

## Open result page
&RSA_header("feature-map result", "results");
&ListParameters() if ($ENV{rsat_echo} >= 2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();


### intialization
$feature_map_command = "$SCRIPTS/feature-map ";
$prefix = "feature-map";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); $tmp_file_name = &ShortFileName($tmp_file_path);
$tmp_color_path = &RSAT::util::make_temp_file("","color_file", 1); $tmp_color_file = &ShortFileName($tmp_color_path);



$features_from_swissprot_cmd = $SCRIPTS."/features-from-swissprot";
$features_from_msf_cmd = $SCRIPTS."/features-from-msf";
$features_from_gibbs_cmd = $SCRIPTS."/features-from-gibbs";
$features_from_fugue_cmd = $SCRIPTS."/features-from-fugue";
$features_from_dssp_cmd = $SCRIPTS."/features-from-dssp";
$features_from_matins_cmd = $SCRIPTS."/features-from-matins";
$features_from_sigscan_cmd = $SCRIPTS."/features-from-sigscan";
$features_from_dnapat_cmd = $SCRIPTS."/features-from-dnapat";
$features_from_tffact_cmd = $SCRIPTS."/features-from-tffact";
$features_from_tfgene_cmd = $SCRIPTS."/features-from-tfgene";
$features_from_patser_cmd = $SCRIPTS."/features-from-patser";

$title = "feature-map result";

#### read parameters ####
$parameters = "";

### dynamic map
if (lc($query->param('htmap')) eq "on") {
    $parameters .= " -htmap ";
}

### feature thickness
### proportional to score
if (lc($query->param('scorethick')) eq "on") {
    $parameters .= " -scorethick ";
}

### max feature thickness
if (&IsNatural($query->param('maxfthick'))) {
    $parameters .= " -xcmaxfthick ".$query->param('maxfthick');
}

### min feature thickness
if (&IsNatural($query->param('minfthick'))) {
    $parameters .= " -minfthick ".$query->param('minfthick');
}


### legend ###
if (lc($query->param('legend')) eq "on") {
    $parameters .= " -legend ";
}

### scale bar ###
if ($query->param('scalebar') eq "on") {
    $parameters .= " -scalebar ";
}

### sequence names
unless ($query->param('seq_names') eq "on") {
    $parameters .= " -no_name ";
}

### horizontal format ###
if ($query->param('orientation') =~ /vertic/i) {
    $parameters .= " -vertical ";
} else {
    $parameters .= " -horizontal ";
}


### scale bar step ###
if ($query->param('scalestep') =~ /^\d+$/) {
  $parameters .= " -scalestep ".$query->param('scalestep');
}

### title ###
if ($query->param('title') ne "") {
    $title = $query->param('title');
    $title =~ s/\'//g;
    $query->param('title',$title);
    $parameters .= " -title \'".$query->param('title'). "\'";
}

my $title_position = $query->param('group_name');
my $view_scalebar = $query->param('view_scalebar');
my $view_score = $query->param('view_score');
my $view_legends = $query->param('view_legends');
my $legends_position = $query->param('legends');
my $export_file = $query->param('print_file');




print( $group );

### display limits ###
if ($query->param('from') ne "") {
    $parameters .= " -from ".$query->param('from');
}
if ($query->param('to') ne "") {
    $parameters .= " -to ".$query->param('to');
}
if ($query->param('origin') ne "") {
    $parameters .= " -origin ".$query->param('origin');
}

### map size ###
if ($query->param('mlen') ne "") {
    $parameters .= " -mlen ".$query->param('mlen');
}
if ($query->param('mapthick') =~ /\d+/) {
    $parameters .= " -mapthick ".$query->param('mapthick');
}
if ($query->param('mspacing') =~ /\d+/) {
    $parameters .= " -mspacing ".$query->param('mspacing');
}

### handle ###
if (lc($query->param('handle')) =~ /dot/) {
    $parameters .= " -dot ";
} elsif (lc($query->param('handle')) =~ /symbol/) {
    $parameters .= " -symbol ";
}

### image format
if ($query->param('img_format')) {
    $image_format = $query->param('img_format');
} else {
    $image_format = $ENV{rsat_img_format} || "png";
}
$parameters .= " -format ".$image_format;

### palette
if (lc($query->param('palette')) =~ /mono/i) {
    $parameters .= " -mono ";
}

### bg color
if (lc($query->param('bgcolor'))) {
    $parameters .= " -bgcolor ".$query->param('bgcolor');
}


### label keys ###
$label = "";
if (lc($query->param('label_strand')) eq "on") {
    $label .= "strand,";
}
if (lc($query->param('label_pos')) eq "on") {
    $label .= "pos,";
}
if (lc($query->param('label_id')) eq "on") {
    $label .= "id,";
}
if (lc($query->param('label_descr')) eq "on") {
    $label .= "descr,";
}
if (lc($query->param('label_score')) eq "on") {
    $label .= "score,";
}
$label =~ s/,$//;
$parameters .= " -label $label " unless ($label eq "");

### id selection ###
@selected_ids = $query->param('id_selection');
if ($#selected_ids >= 0) {
    for $i (0..$#selected_ids) {
	$selected{$selected_ids[$i]} = 1;
	$selected_ids[$i] = "'".$selected_ids[$i]."'";
    }
    $id_selection = join ",", @selected_ids;
    unless ($selected{'*all*'}) {
	$parameters .= " -select $id_selection ";
    }
}

## color file ##
if ($query->param('color_file')) {
  $color_file = $tmp_file_path."_colors.tab";
  push @result_files, "Colors", $color_file;
  open COLOR, ">$color_file";
  $upload_color_file = $query->param('color_file');
  while (<$upload_color_file>) {
    print COLOR;
  }
  close COLOR;
  $parameters .= " -colors $color_file ";
}

### data file ####
if ($query->param('feature_file') =~ /\S/) {

  ### file on the server
  $feature_file = $query->param('feature_file');

} elsif (($query->param('uploaded_file')) ||
	 ($query->param('data') =~ /\S/)) {

  $feature_format = $query->param('format');
  $feature_file = $tmp_file_path.".".$feature_format;
  $feature_file =~ s/\s+/_/;

  ### convert data towards feature-map format
  if ($feature_format =~ /swiss/i) {
    open DATA, "| $features_from_swissprot_cmd -o $feature_file";
  } elsif ($feature_format =~ /transfac factor/i) {
    open DATA, "| $features_from_tffact_cmd -o $feature_file";
  } elsif ($feature_format =~ /transfac gene/i) {
    open DATA, "| $features_from_tfgene_cmd -o $feature_file";
  } elsif ($feature_format =~ /msf/i) {
    open DATA, "| $features_from_msf_cmd -o $feature_file";
    $parameters .= " -aacolors ";
    $parameters .= " -horiz ";
  } elsif ($feature_format =~ /dssp/i) {
    open DATA, "| $features_from_dssp_cmd -o $feature_file";
  } elsif ($feature_format =~ /matins/i) {
    open DATA, "| $features_from_matins_cmd -o $feature_file";
  } elsif ($feature_format =~ /dna\-pattern/i) {
    open DATA, "| $features_from_dnapat_cmd -o $feature_file";
  } elsif ($feature_format =~ /patser/i) {
    open DATA, "| $features_from_patser_cmd -o $feature_file";
  } elsif ($feature_format =~ /signal scan/i) {
    open DATA, "| $features_from_sigscan_cmd -o $feature_file";
  } elsif ($feature_format =~ /gibbs/i) {
    open DATA, "| $features_from_gibbs_cmd -o $feature_file";
  } elsif ($feature_format =~ /fugue/i) {
    open DATA, "| $features_from_fugue_cmd -o $feature_file";
  } else {
    open DATA, ">$feature_file";
  }

  if ($query->param('uploaded_file')) {
    $upload_feature_file = $query->param('uploaded_file');
    $type = $query->uploadInfo($upload_feature_file)->{'Content-Type'};
    #	&Info($feature_file, "\n", $upload_feature_file, "\n", $type);
    while (<$upload_feature_file>) {
      #	    print $_;
      print DATA;
    }
    close DATA;

    ### upload file from the client
    #	$fh = $query->param('uploaded_file');
    #	while (<$fh>) {
    #	    print DATA;
    #	}
  } else {
    ### data from the textarea
    print DATA $query->param('data');
    #	print "<PRE>\$feature_file\t$feature_file</PRE>";
  }
  close DATA;

} else {
    &cgiError("The feature list should not be empty.");
}

push @result_files, "Input features (.ft)", $feature_file;
$parameters .= " -i $feature_file ";

### map file ###
$graph_file = $tmp_file_path.".".$image_format;
#print( $graph_file . "\n" );
push @result_files, "Map file ($image_format)", $graph_file;
$htmap_file = $tmp_file_path.".html";
my $htmap_file_new = $tmp_file_path."_new.html";
push @result_file, "Html report (html)", $htmap_file;
$parameters .= " -o $graph_file > $htmap_file";

$feature_map_command .= " ".$parameters;

## Report the command
&ReportWebCommand($feature_map_command);

### execute the command
system($feature_map_command);
&DelayedRemoval($feature_file);
&DelayedRemoval($graph_file);
&DelayedRemoval($htmap_file);

my $json_file = $tmp_file_path . ".json";
###execute command to make hierarchical JSON
my $json_cmd = "./script_tsv_to_JSON.pl -a $feature_file > $json_file";
#print( $json_cmd );
my $make_json = system( $json_cmd );
if( $make_json ){
	print( "$make_json: $0\n" );
}

my $json_URL = $ENV{rsat_www}."/tmp/"; $json_URL .= &RSAT::util::RelativePath(&RSAT::util::get_pub_temp(), $json_file );
### display the result ###
my $graph_URL = $ENV{rsat_www}."/tmp/"; $graph_URL .= &RSAT::util::RelativePath(&RSAT::util::get_pub_temp(), $graph_file);
#print( "PATH ÄTH " . $graph_URL . "\n" );
my $html_URL = $ENV{rsat_www}."/tmp/"; $html_URL .= &RSAT::util::RelativePath(&RSAT::util::get_pub_temp(), $htmap_file);
my $short_graph_file = &ShortFileName($graph_file);
my $short_feature_file = &ShortFileName($feature_file);


##$htmap_file_new = '<html><body><script type="text/javascript">' . join( "", @html1 ) . "d3.json( \"$tmp_file_path.json\", function(data){" . join( "", @html2 ) . '</script></body></html>';

&print_html( $title, $json_URL , $file_name, $title_position, $view_scalebar, $legends_position, $export_file);

if (($image_format ne 'ps')
    && (lc($query->param('htmap')) eq "on")) {
  my $htmap_content = `cat $htmap_file`;
  $htmap_content =~ s/<\/*html>//gi;
  $htmap_content =~ s/<\/*body>//gi;
  $htmap_content =~ s/<\/*head>//gi;
  $htmap_content =~ s/<title>.*<\/title>//gi;
  $htmap_content =~ s/${short_graph_file}/${graph_URL}/g;
  $htmap_content =~ s/${short_feature_file}/${data_URL}/g;
#  $htmap_content =~ s/</&lt;/g;
#  $htmap_content =~ s/>/&gt;/g;
# print $htmap_content;
  # print $htmap_file_new;
} else {
  #print "<center><a href='".$html_URL."'><img src='".$html_URL."'></a></center><P>\n";
}

&PrintURLTable(@result_files);


print "<hr>";

# die join "\n",
#   "graph_file = ".$graph_file,
#   "graph_URL = ".$graph_URL,
#   "short_graph_file = ".$short_graph_file;

exit(0);


sub print_html{
	my ( $title, $json_URL, $file_name, $title_position, $view_scalebar, $legends_position, $export_to_svg, $export_to_png, $export_to_jpg) = @_;
	#print( "Title " . $title . ", ". "path " . $json_URL . " $title_position" . "$legends_position" ."$view_scalebar" . "$export_file");
   #print("$json_file");

	my $html = qq(
      <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
      <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js"></script>
      <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
      <script src="https://d3js.org/d3.v4.min.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/randomcolor/0.5.2/randomColor.js" ></script>
      <script src="https://cdn.rawgit.com/eligrey/canvas-toBlob.js/f1a01896135ab378aa5c0118eadd81da55e698d8/canvas-toBlob.js"></script>
      <script src="https://cdn.rawgit.com/eligrey/FileSaver.js/e9d941381475b5df8b7d7691013401e171014e89/FileSaver.min.js"></script>

      <!-- Defines the limits of the graph according to the size of the grid in which the SVG is created -->
      <div class="container">
        <div class="row">
          <div class="col-md-0"></div>
          <div class="col-md-12">
            <div id="feature_map_div"></div>
          </div>
          <div class="col-md-0"></div>
        </div>
      </div>

      <canvas id="myCanvas" width="0" height="0" style="border:1px solid #d3d3d3;">
      	Your browser does not support the HTML5 canvas tag.
      </canvas>

      <script type="text/javascript">
        drawFeatureMap();

        function drawFeatureMap(){
          // JSON file generated by RSAT - script_tsv_to_JSON.pl
          var jsonFile = jsonFileName("$json_URL");

          // GRAPH TITLE OPTIONS
          // Title position variables - default positions
          var titlePositionLeft = "false";
          var titlePositionCenter = "true";
          var titlePositionRight = "false";

          // If to re-define title position variables according to user input
          if( "$title_position" == "left"){
            titlePositionLeft = "true";
            titlePositionCenter = "false";
          }
          else if( "$title_position" == "right"){
            titlePositionRight = "true";
            titlePositionCenter = "false";
          }

          //EXPORT OPTIONS
          // Export file variables - default positions
          var svg_export = "false";
          var png_export = "false";
          var jpg_export = "false";

          // If to re-define export variables according to user input
          if("$export_file" == "svg"){
            svg_export = "true";
          }
          else if("$export_file" == "png"){
            png_export = "true";
          }
          else if("$export_file" == "jpg"){
            jpg_export = "true";
          }

          //ADVANCED OPTIONS VARIABLES
          // View legend and legend position options
          var view_legend = "true";
          var legends_top = "true";
          var legends_bottom = "false";
          var legends_left = "false";
          var legends_right = "false";

          // If to re-define legends options variables
          if("$legends_position" == "left"){
            legends_top = "false";
            legends_left = "true";
          }

          // View scale options
          var viewScale;

          // If to re-define view scalebar options variables
          if("$view_scalebar" == "true"){
            viewScale = "true";
          }
          else{
            viewScale = "false";
          }

          // View score options
          var viewScore = "true";

          //Other options
          var defaultFont = "Arial";
          var titleFontSize = 20;
          var userTitleFontSize = 18;
          var transparency = "true";
          var extraMargin = 100;

          //Open json file
          d3.json( "$json_URL", function(data){

            // Some graph variables
            var p_title_x = 25;
            var p_title_y = 20;
            var p_title_width = 0;
            var p_title_height = 0;
            var p_title_margin_left = 0;
            var p_title_margin_right = 40;
            var p_title_space = 26;
            var p_legends_left_x = p_title_x;

            //Calculate max size of the name of genes

            var geneNameList = [];
            for (i = 0; i < data.length; i++){
              geneNameList.push(data[i].Gene_name);
            }

            var maxSizeGeneName = calculateMaxSizeOfTextList(geneNameList, defaultFont, 16);


            // Relative width of the svg according to the grid in which the feature_map_div resides
            var svgWidth = document.getElementById("feature_map_div").offsetWidth;

            if(maxSizeGeneName >= (.5 * svgWidth)){
              svgWidth = svgWidth + maxSizeGeneName + 50; 
            }

            // Make function to calculate
            var gene_num = 0;
            var dedicated_space = 60;
            var title_space = 100;

            // To know how many genes are there to graph
            for(i=0; i<data.length; i++){
               gene_num = gene_num +1;
            }

            // FUNCTION TO CALCULATE GENE Number
            var geneNumber = calculateGeneNumber(data);
            console.log("EL NUMERO DE GENES ES " + geneNumber);

            function calculateGeneNumber(data){
              var geneNameArray = [];

              for(i = 0; i < data.length; i++){
                geneNameArray.push(data[i].Gene_name);
              }

              var geneNameArraySorted = geneNameArray.sort();

              var geneListCount = 1;

              for(i = 0; i < geneNameArraySorted.length - 1; i++){
                if(geneNameArraySorted[i+1] != geneNameArraySorted[i]){
                  geneListCount++;
                }
              }
              return geneListCount;
            }



            //Calculate height according to file and user preferences
            var total_svg_height = ((geneNumber * dedicated_space) + title_space + 200);
            //var total_svg_height = 4000;
		        console.log("Total_svg_height = "+ total_svg_height);

            var width = svgWidth, height = total_svg_height;

            //Creation of the svg
            var graph_svg = d3.select("#feature_map_div")
               .append("svg")
               .attr("x",0)
               .attr("y",0)
               .attr("width",svgWidth)
               .attr("height",total_svg_height);

               // TITLE BAR CODE
               var g_title_x = 20;
               var g_title_y = 40;
               var g_title_anchor = "start";

               //Print of the application title
               graph_svg.append("text")
                  .attr("x",function(){
                     if(titlePositionLeft == "true"){
                        return g_title_x;
                     }
                     else if(titlePositionCenter == "true"){
                        return (svgWidth/2);
                     }
                   else if(titlePositionRight == "true"){
                        return (svgWidth - p_title_margin_right);
                     }
                  })
                  .attr("y",g_title_y)
                  .attr("font-size", titleFontSize)
                  .attr("text-anchor",function(){
                     if(titlePositionLeft == "true"){
                        return "start";
                     }
                     else if(titlePositionCenter == "true"){
                        return "middle";
                     }
                     else if(titlePositionRight == "true"){
                        return "end";
                     }
                  })
                  .attr("font-family", defaultFont)
                  .attr("fill","rgb(51,51,51)")
                  .text("");

                  //Print of the User title
               graph_svg.append("text")
                  .attr("x",function(){
                     if(titlePositionLeft == "true"){
                        return g_title_x;
                     }
                     else if(titlePositionCenter == "true"){
                        return (svgWidth/2);
                     }
                     else if(titlePositionRight == "true"){
                        return (svgWidth - p_title_margin_right);
                     }
                  })
                  .attr("y",g_title_y + p_title_space)
                  .attr("font-size", userTitleFontSize)
                  .attr("text-anchor",function(){
                     if(titlePositionLeft == "true"){
                        return "start";
                     }
                     else if(titlePositionCenter == "true"){
                        return "middle";
                     }
                     else if(titlePositionRight == "true"){
                        return "end";
                     }
                  })
                  .attr("font-family", defaultFont)
                  .attr("fill","rgb(51,51,51)")
                  .text('$title');

                  //LEGENDS BAR
               //Get all the legends and extract discard the repeated ones
               var label_array_total = [];

               for(i=0; i < data.length; i++){
            		for(f=0; f< data[i].Features.length; f++){
                     label_array_total.push(data[i].Features[f].seq);
                     //console.log(data[i].Features[f].seq);
                  }
               }

               var sorted_label_array_total = label_array_total.sort();
               var selected_label_array_total = [];

               for(i=0; i<sorted_label_array_total.length; i++){

                  if(sorted_label_array_total[i] != sorted_label_array_total[i+1]){
                     selected_label_array_total.push(sorted_label_array_total[i]);
                  }
               }

               console.log("SIZE OF ARRAY SELECTED " + selected_label_array_total.length)

               //Color assignment DEFAULT
               var default_colors = ["#f6412d","#3e4eb7","#48af4b","#ff9800","#9c1cb1","#02BBD5","#FFEC19","#9D9D9D","#F6412D","#3E4EB7","#48AF4B","#FF9800","#9C1CB1","#02BBD5","#FFEC19","#9D9D9D"];
               var label_array_color = [];
               var label_color;
               var rand_color;
               var default_color;

               for(i=0; i < selected_label_array_total.length; i++){
               if(i < 8){
                  //Asigna de los de default
                  default_color = default_colors[i];
                  label_color = {label:selected_label_array_total[i], color:default_color};
                  label_array_color.push(label_color);
               }
               else{
                  //Asigna los random
                  rand_color = randomColor();
                  label_color = {label:selected_label_array_total[i], color:rand_color};
                  label_array_color.push(label_color);
               }
            }

            //Color assignment RANDOM
            // var label_array_color = [];// Array of labels and randonmly asigned colors.
            // var label_color;
            // var rand_color;
            // //
            // for(i=0; i < sorted_label_array_total.length; i++){
            //    if(sorted_label_array_total[i] != sorted_label_array_total[i+1]){
            //       rand_color = randomColor();
            //       label_color = {label:sorted_label_array_total[i], color:rand_color};
            //       label_array_color.push(label_color);
            //       console.log(i);
            //    }
            // }

            //Printing of legends

            var legends_x = 25;
            var legends_y = 110;

            //Legends title
               graph_svg.append("text")
                  .attr("x",legends_x + 20)
                  .attr("y",legends_y)
                  .attr("font-size", 16)
                  .attr("text-anchor","start")
                  .attr("font-family", defaultFont)
                  .attr("fill","rgb(51,51,51)")
                  .text("Legends");

                  //Code for printing the legends automatically

                        var max_size = 0;
                        var rect_size = 14;

                        //Get the size of the longest label
                        for(i=0; i< label_array_color.length; i++){
                           size_2 = calculate_text_width(label_array_color[i].label,defaultFont,14);

                           if (max_size < size_2){
                              max_size = size_2;
                           }
                        }
                        //console.log("Size " + max_size);

                        var tab_size = 40;
                        var auto_columns = Math.floor((svgWidth / (tab_size + (max_size + 20) + rect_size)));
         //Printing of the labels and their color squares
                        var labels_x = legends_x;
                        var labels_y = legends_y + 16;
                        var labels_left_y = legends_y + 16;
                        var labels_left_y_pass = labels_y;

                        var labels_x_jump = 10+ rect_size + tab_size + max_size;
                        var labels_y_jump = 40;

                        var cont = 0;


            if(legends_top == "true"){
                        for(i=0; i< label_array_color.length; i++){
                           graph_svg.append("rect")
                              .attr("x",labels_x + 20)
                              .attr("y",labels_y)
                              .attr("width",rect_size)
                              .attr("height",rect_size)
                              .attr("fill",label_array_color[i].color);

                           graph_svg.append("text")
                              .attr("x",labels_x + rect_size +4 + 20)
                              .attr("y",labels_y + rect_size -1)
                              .attr("font-size", 11)
                              .attr("text-anchor","start")
                              .attr("font-family", defaultFont)
                              .attr("fill","rgb(51,51,51)")
                              .text(label_array_color[i].label);

                           cont = cont +1;
                           labels_x = labels_x + labels_x_jump;

                           if(cont >= auto_columns){
                              labels_y = labels_y + 24;
                              labels_x = legends_x;
                              cont = 0;
                           }
                        }
                  }

               else if(legends_left == "true"){
   for(i=0; i< label_array_color.length; i++){
      graph_svg.append("rect")
         .attr("x",labels_x)
         .attr("y",labels_left_y)
         .attr("width",rect_size)
         .attr("height",rect_size)
         .attr("fill",label_array_color[i].color);

      graph_svg.append("text")
         .attr("x",labels_x + rect_size +4)
         .attr("y",labels_left_y + rect_size -1)
         .attr("font-size", 16)
         .attr("text-anchor","start")
         .attr("font-family", defaultFont)
         .attr("fill","rgb(51,51,51)")
         .text(label_array_color[i].label);

      cont = cont +1;
      labels_left_y = labels_left_y + labels_y_jump;

      // if(cont >= auto_columns){
      //    labels_y = labels_y + 24;
      //    labels_x = legends_x;
      //    cont = 0;
      // }
   }

   graph_svg.append("line")
   .attr("x1",labels_x + rect_size +40 + max_size)
   .attr("y1",labels_left_y_pass)
   .attr("x2",labels_x + rect_size +40 + max_size)
   .attr("y2", total_svg_height)
   .attr("stroke","#000000")
   .attr("stroke-width",.5);
}

//Print of the main scale bar
var scale_bar_x;

if(legends_top == "true"){
scale_bar_x = 120;

}
else if(legends_left == "true"){
scale_bar_x = max_size + 220;
}


var scale_bar_y;
if (legends_top == "true"){
scale_bar_y =labels_y + 80;
}
else if(legends_left == "true"){
scale_bar_y = labels_left_y_pass;
}
var scale_width = svgWidth - 20 - scale_bar_x;
//


// //Obtener el mayor rango de los genes para calcular el rango de la main scale_bar
var range = 0;
//
for(i=0; i<data.length; i++){
 if(data[i].Start > range){
    range = data[i].Start;
 }
}
// console.log("range "+ range);
//
// //Graficar ticks de la barra principal
// //Primero hacer calculo de la escala
var line_width = (svgWidth - scale_bar_x - 50);
var tick_space = 50; // Puede ser definida por el usuario
var tick_jump = ((tick_space * line_width) / range);
var value = range;



// CODIGO PRUEBA NUEVA ESCALA
var geneNameSpace = {
  x: (p_title_x + 20),
  y: scale_bar_y,
  width: (p_title_x + 20 + maxSizeGeneName)
};

var featuresSpace = geneNameSpace.width;

var scaleBarX = featuresSpace + 20;
var svgInnerRightMargin = 50;

var scaleBarEnd = svgWidth - (svgInnerRightMargin + geneNameSpace.x + maxSizeGeneName);

var scalePrueba = d3.scaleLinear()
  .domain([range,0])
  .range([0,scaleBarEnd]);

var xAxis = d3.axisBottom()
  .scale(scalePrueba);

var marginLeft = 0;
var marginTop = scale_bar_y;

graph_svg.append("g")
  .attr("transform", "translate(" + scaleBarX + "," + scale_bar_y + ")")
  .call(xAxis);

if(viewScale == "true"){
    //  graph_svg.append("line")
      //   .attr("x1",scale_bar_x)
        // .attr("y1",scale_bar_y)
         //.attr("x2",svgWidth -50)
         //.attr("y2",scale_bar_y)
         //.attr("stroke","#000000")
         //.attr("stroke-width",.5);

        // while(scale_bar_x < (svgWidth - 49)){
        //
          //  graph_svg.append("line")
            //   .attr("x1",scale_bar_x)
            //   .attr("y1",scale_bar_y-3)
            //   .attr("x2",scale_bar_x)
            //   .attr("y2",scale_bar_y+3)
            //   .attr("stroke","#000000")
            //   .attr("stroke-width",1);

          // graph_svg.append("text")
            //  .attr("x",scale_bar_x)
            //  .attr("y",scale_bar_y+12)
            //  .attr("font-size", 10)
            //  .attr("text-anchor","middle")
            //  .attr("font-family", defaultFont)
            //  .attr("fill","rgb(51,51,51)")
            //  .text(value);
        //
          //  scale_bar_x = scale_bar_x + tick_jump;
          //  value = value - tick_space;
        // }
   }

  //Print of each of the gene names
   var gene_name_x;
   var gene_name_y = scale_bar_y + 40;
   var gene_jump = 40;
   console.log(data);
   graph_svg.selectAll("text.genes")
      .data(data)
      .enter()
      .append("text")
      .attr("x",function(){
      if(legends_top == "true"){
         return p_title_x + 20;
      }
      else if(legends_left == "true"){
         return p_title_x + max_size + extraMargin;
      }
   })
   .attr("y", function(d,i){
      return gene_name_y + ((i+1)*gene_jump);
   })
   .attr("font-size", 16)
   .attr("text-anchor","start")
   .attr("font-family", defaultFont)
   .attr("fill","rgb(51,51,51)")
   .text(function(d){
      return d.Gene_name;
   });

   // //Print each of the gene horizontal lines
   var gene_line_x2 = svgWidth - 50;
   var gene_line_y ;
   var start_y = 95;
   //
      graph_svg.selectAll("line.bar2")
         .data(data)
         .enter()
         .append("line")
         .attr("x1", function(){
           return (scalePrueba(0) + scaleBarX);
           })
         .attr("y1",function(d,i){
            return (gene_name_y + ((i+1)*gene_jump))-10;
            })
         .attr("x2",function(d){
           console.log("D.START ES " + d.Start);
            return (scalePrueba(d.Start) + scaleBarX);
            })
          .attr("y2",function(d,i){
             return (gene_name_y + ((i+1)*gene_jump))-10;
             })
         .attr("stroke","#333333")
         .attr("stroke-width",1);

         graph_svg.selectAll("line.tickCero")
            .data(data)
            .enter()
            .append("line")
            .attr("x1", function(){
              return (scalePrueba(0) + scaleBarX);
              })
            .attr("y1",function(d,i){
               return ((gene_name_y + ((i+1)*gene_jump))-10) - 5 ;
               })
            .attr("x2",function(){
              return (scalePrueba(0) + scaleBarX);
               })
             .attr("y2",function(d,i){
                return ((gene_name_y + ((i+1)*gene_jump))-10) + 5;
                })
            .attr("stroke","#333333")
            .attr("stroke-width",1);

            graph_svg.selectAll("line.tickFinal")
               .data(data)
               .enter()
               .append("line")
               .attr("x1", function(d){
                 return (scalePrueba(d.Start) + scaleBarX);
                 })
               .attr("y1",function(d,i){
                  return ((gene_name_y + ((i+1)*gene_jump))-10) - 5 ;
                  })
               .attr("x2",function(d){
                 return (scalePrueba(d.Start) + scaleBarX);
                  })
                .attr("y2",function(d,i){
                   return ((gene_name_y + ((i+1)*gene_jump))-10) + 5;
                   })
               .attr("stroke","#333333")
               .attr("stroke-width",1);

      //var gene_line_ticks_y;
   //var gene_line_ticks_x = svgWidth - 50;

   //for(i=0; i < data.length; i ++){

    //  gene_line_ticks_y = gene_name_y + ((i+1)*gene_jump)-10;

      //Tick on cero
      // graph_svg.append("line")
        //  .attr("x1",svgWidth - 50)
        //  .attr("y1",gene_line_ticks_y - 3)
        //  .attr("x2",svgWidth - 50)
        //  .attr("y2",gene_line_ticks_y + 3)
        //  .attr("stroke","#333333")
        //  .attr("stroke-width",1);

      //Tick on Start
      //graph_svg.append("line")
      //  .attr("x1",svgWidth - 50 - ((data[i].Start * line_width)/range) + maxSizeGeneName)
      //  .attr("y1",gene_line_ticks_y - 3)
      //  .attr("x2",svgWidth - 50 - ((data[i].Start * line_width)/range) + maxSizeGeneName)
      //  .attr("y2",gene_line_ticks_y + 3)
      //  .attr("stroke","#333333")
      //  .attr("stroke-width",.5);

      //while(gene_line_ticks_x > (svgWidth - 50 - ((data[i].Start * line_width)/range) + maxSizeGeneName)){
      //   graph_svg.append("line")
      //       .attr("x1",gene_line_ticks_x)
      //       .attr("y1",gene_line_ticks_y - 3)
      //       .attr("x2",gene_line_ticks_x)
      //       .attr("y2",gene_line_ticks_y + 3)
      //       .attr("stroke","#333333")
      //       .attr("stroke-width",.5);

        //     gene_line_ticks_x = gene_line_ticks_x - ((50 * line_width)/range);
    //  }
    //  gene_line_ticks_x = svgWidth - 50;
   //}



   for(i=0; i < data.length; i++){

   for(f=0; f < data[i].Features.length; f++){
      graph_svg.append("rect")
         .attr("class","feature_rect")
         .attr("x", function(){
           return (scalePrueba(data[i].Features[f].start) + scaleBarX);
           })
         .attr("y",function(){
            if(data[i].Features[f].strand == "D"){
               return (gene_name_y + ((i+1)*gene_jump) - (data[i].Features[f].score * 3)-10);
            }
          else if(data[i].Features[f].strand == "R"){
               return (gene_name_y + ((i+1)*gene_jump)-10);
            }
          else if(data[i].Features[f].strand == "DR"){
               return (gene_name_y + ((i+1)*gene_jump) - (data[i].Features[f].score * 3)-10);
            }
         })
         .attr("width", function(){
          console.log("fStart " + data[i].Features[f].start + " fEnd " + data[i].Features[f].end);
          var vOne = scalePrueba(data[i].Features[f].start);
          var vTwo = scalePrueba(data[i].Features[f].end);
          console.log("VStart " + vOne + " VEnd " + vTwo);
          return (scalePrueba(data[i].Features[f].end) - scalePrueba(data[i].Features[f].start));
         })
         .attr("height", function(){
           if(data[i].Features[f].strand == "D"){
              return (data[i].Features[f].score * 3);
           }
         else if(data[i].Features[f].strand == "R"){
              return (data[i].Features[f].score * 3);
           }
         else if(data[i].Features[f].strand == "DR"){
              return ((data[i].Features[f].score * 3)*2);
           }
          })
         .attr("opacity",function(){
            if(transparency == "true"){
               return 0.7;
            }
         })
         .attr("fill", function(){
            for(c=0; c < label_array_color.length; c++){
               if(data[i].Features[f].seq == label_array_color[c].label){
                  //console.log("SEC: " + data[i].Features[f].seq + "LABEL: " + label_array_color[c].label + "Color " + label_array_color[c].color);
                  return label_array_color[c].color;
               }
            }
         })
         .on("mouseover", function(){

           console.log("DATA " + data[i].Features[f].seq);

			   });
   }
}



d3.select("#saveButton").on("click", function(){
   var svgString = getSVGString(graph_svg.node());
   svgString2Image( svgString, 2*width, 2*height, "png", save ); // passes Blob and filesize String to the callback

   function save( dataBlob, filesize ){
   saveAs( dataBlob, "D3 vis exported to PNG.png" ); // FileSaver.js function
   }
});

d3.select("#saveButton_jpeg").on("click", function(){
	var svgString = getSVGString(graph_svg.node());
	svgString2Image_jpeg( svgString, 2*width, 2*height, "jpeg", save ); // passes Blob and filesize String to the callback

	function save( dataBlob, filesize ){
		saveAs( dataBlob, "D3 vis exported to JPEG.jpeg" ); // FileSaver.js function
	}
});



//Export options
if (svg_export == "true"){
create_svg(jsonFile);
}
else if(png_export == "true"){

      var svgString = getSVGString(graph_svg.node());
      svgString2Image( svgString, 2*width, 2*height, "png", save ); // passes Blob and filesize String to the callback

      function save( dataBlob, filesize ){
      saveAs( dataBlob, jsonFile + ".png" ); // FileSaver.js function
      }

}
else if(jpg_export == "true"){
   var svgString = getSVGString(graph_svg.node());
   svgString2Image_jpeg( svgString, 2*width, 2*height, "jpeg", save ); // passes Blob and filesize String to the callback

   function save( dataBlob, filesize ){
      saveAs( dataBlob, jsonFile + ".jpeg" ); // FileSaver.js function
   }
}



         });

      }//Closing of function feature map

   //FUNCTIONS

   function save_to_png(){
      var svgString = getSVGString(graph_svg.node());
      svgString2Image( svgString, 2*width, 2*height, "png", save ); // passes Blob and filesize String to the callback

      function save( dataBlob, filesize ){
      saveAs( dataBlob, "D3 vis exported to PNG.png" ); // FileSaver.js function
      }
   }

   function calculate_text_width(text, font, size){
      var c = document.getElementById("myCanvas");
      var ctx = c.getContext("2d");
      ctx.font = size+"px "+font;
      var txt = text;
      var size = ctx.measureText(txt).width;
      return size;
   }

   function calculateSvgSize(data){
     return 4000;
   }

   //------------------------

   function getSVGString( svgNode ) {
	svgNode.setAttribute("xlink", "http://www.w3.org/1999/xlink");
	var cssStyleText = getCSSStyles( svgNode );
	appendCSS( cssStyleText, svgNode );

	var serializer = new XMLSerializer();
	var svgString = serializer.serializeToString(svgNode);
	svgString = svgString.replace(/(\w+)?:?xlink=/g, 'xmlns:xlink='); // Fix root xlink without namespace
	svgString = svgString.replace(/NS\d+:href/g, 'xlink:href'); // Safari NS namespace fix

	return svgString;

	function getCSSStyles( parentElement ) {
		var selectorTextArr = [];

		// Add Parent element Id and Classes to the list
		selectorTextArr.push( "#"+parentElement.id );
		for (var c = 0; c < parentElement.classList.length; c++)
				if ( !contains("."+parentElement.classList[c], selectorTextArr) )
					selectorTextArr.push( "."+parentElement.classList[c] );

		// Add Children element Ids and Classes to the list
		var nodes = parentElement.getElementsByTagName("*");
		for (var i = 0; i < nodes.length; i++) {
			var id = nodes[i].id;
			if ( !contains("#"+id, selectorTextArr) )
				selectorTextArr.push( "#"+id );

			var classes = nodes[i].classList;
			for (var c = 0; c < classes.length; c++)
				if ( !contains("."+classes[c], selectorTextArr) )
					selectorTextArr.push( "."+classes[c] );
		}

		// Extract CSS Rules
		var extractedCSSText = "";
		for (var i = 0; i < document.styleSheets.length; i++) {
			var s = document.styleSheets[i];

			try {
			    if(!s.cssRules) continue;
			} catch( e ) {
		    		if(e.name !== "SecurityError") throw e; // for Firefox
		    		continue;
		    	}

			var cssRules = s.cssRules;
			for (var r = 0; r < cssRules.length; r++) {
				if ( contains( cssRules[r].selectorText, selectorTextArr ) )
					extractedCSSText += cssRules[r].cssText;
			}
		}


		return extractedCSSText;

		function contains(str,arr) {
			return arr.indexOf( str ) === -1 ? false : true;
		}

	}

	function appendCSS( cssText, element ) {
		var styleElement = document.createElement("style");
		styleElement.setAttribute("type","text/css");
		styleElement.innerHTML = cssText;
		var refNode = element.hasChildNodes() ? element.children[0] : null;
		element.insertBefore( styleElement, refNode );
	}
}


//----------------

function svgString2Image( svgString, width, height, format, callback ) {
	var format = format ? format : "png";

	var imgsrc = "data:image/svg+xml;base64,"+ btoa( unescape( encodeURIComponent( svgString ) ) ); // Convert SVG string to data URL

	var canvas = document.createElement("canvas");
	var context = canvas.getContext("2d");

	canvas.width = width;
	canvas.height = height;

	var image = new Image();
	image.onload = function() {
		context.clearRect ( 0, 0, width, height );
		context.drawImage(image, 0, 0, width, height);

		canvas.toBlob( function(blob) {
			var filesize = Math.round( blob.length/1024 ) + ' KB';
			if ( callback ) callback( blob, filesize );
		});


	};

	image.src = imgsrc;
}

//--------------
function svgString2Image_jpeg( svgString, width, height, format, callback ) {
	var format = format ? format : "jpeg";

	var imgsrc = "data:image/svg+xml;base64,"+ btoa( unescape( encodeURIComponent( svgString ) ) ); // Convert SVG string to data URL

	var canvas = document.createElement("canvas");
	var context = canvas.getContext("2d");

	canvas.width = width;
	canvas.height = height;

	var image = new Image();
	image.onload = function() {
		context.clearRect ( 0, 0, width, height );
		context.drawImage(image, 0, 0, width, height);

		canvas.toBlob( function(blob) {
			var filesize = Math.round( blob.length/1024 ) + " KB";
			if ( callback ) callback( blob, filesize );
		});


	};

	image.src = imgsrc;
}

//----------------------


function jsonFileName(json_url){
  var data = json_url;
  var name = data.substring(data.indexOf("feature-map_")+12,data.indexOf(".json"));
  return name;
}

function calculateMaxSizeOfTextList (list, geneFont, geneFontSize){
  var maxSizeOfList = 0;
  for (i = 0; i < list.length; i++){
    if (calculate_text_width(list[i], geneFont, geneFontSize) > maxSizeOfList) {
      maxSizeOfList = calculate_text_width(list[i]);
    }
  }
  return maxSizeOfList;
}


function create_svg(jsonFile){
   console.log("CREATE SVG OUTPUT " + jsonFile );

   var doctype = '<?xml version="1.0" standalone="no"?><!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">';

   window.URL = (window.URL || window.webkitURL);

   var body = document.body;

   var prefix = {
     xmlns: "http://www.w3.org/2000/xmlns/",
     xlink: "http://www.w3.org/1999/xlink",
     svg: "http://www.w3.org/2000/svg"
   }

   initialize();

   function initialize() {
     var documents = [window.document],
         SVGSources = [];
         iframes = document.querySelectorAll("iframe"),
         objects = document.querySelectorAll("object");

     [].forEach.call(iframes, function(el) {
       try {
         if (el.contentDocument) {
           documents.push(el.contentDocument);
         }
       } catch(err) {
         console.log(err)
       }
     });

     [].forEach.call(objects, function(el) {
       try {
         if (el.contentDocument) {
           documents.push(el.contentDocument);
         }
       } catch(err) {
         console.log(err)
       }
     });

     documents.forEach(function(doc) {
       var styles = getStyles(doc);
       var newSources = getSources(doc, styles);
       // because of prototype on NYT pages
       for (var i = 0; i < newSources.length; i++) {
         SVGSources.push(newSources[i]);
       };
     })
     if (SVGSources.length > 1) {
       createPopover(SVGSources);
     } else if (SVGSources.length > 0) {
       download(SVGSources[0]);
     } else {
       alert("The Crowbar couldn’t find any SVG nodes.");
     }
   }

   function createPopover(sources) {
     cleanup();

     sources.forEach(function(s1) {
       sources.forEach(function(s2) {
         if (s1 !== s2) {
           if ((Math.abs(s1.top - s2.top) < 38) && (Math.abs(s1.left - s2.left) < 38)) {
             s2.top += 38;
             s2.left += 38;
           }
         }
       })
     });

     var buttonsContainer = document.createElement("div");
     body.appendChild(buttonsContainer);

     buttonsContainer.setAttribute("class", "svg-crowbar");
     buttonsContainer.style["z-index"] = 1e7;
     buttonsContainer.style["position"] = "absolute";
     buttonsContainer.style["top"] = 0;
     buttonsContainer.style["left"] = 0;



     var background = document.createElement("div");
     body.appendChild(background);

     background.setAttribute("class", "svg-crowbar");
     background.style["background"] = "rgba(255, 255, 255, 0.7)";
     background.style["position"] = "fixed";
     background.style["left"] = 0;
     background.style["top"] = 0;
     background.style["width"] = "100%";
     background.style["height"] = "100%";

     sources.forEach(function(d, i) {
       var buttonWrapper = document.createElement("div");
       buttonsContainer.appendChild(buttonWrapper);
       buttonWrapper.setAttribute("class", "svg-crowbar");
       buttonWrapper.style["position"] = "absolute";
       buttonWrapper.style["top"] = (d.top + document.body.scrollTop) + "px";
       buttonWrapper.style["left"] = (document.body.scrollLeft + d.left) + "px";
       buttonWrapper.style["padding"] = "4px";
       buttonWrapper.style["border-radius"] = "3px";
       buttonWrapper.style["color"] = "white";
       buttonWrapper.style["text-align"] = "center";
       buttonWrapper.style["font-family"] = "'Helvetica Neue'";
       buttonWrapper.style["background"] = "rgba(0, 0, 0, 0.8)";
       buttonWrapper.style["box-shadow"] = "0px 4px 18px rgba(0, 0, 0, 0.4)";
       buttonWrapper.style["cursor"] = "move";
       buttonWrapper.textContent =  "SVG #" + i + ": " + (d.id ? "#" + d.id : "") + (d.class ? "." + d.class : "");

       var button = document.createElement("button");
       buttonWrapper.appendChild(button);
       button.setAttribute("data-source-id", i)
       button.style["width"] = "150px";
       button.style["font-size"] = "12px";
       button.style["line-height"] = "1.4em";
       button.style["margin"] = "5px 0 0 0";
       button.textContent = "Download";

       button.onclick = function(el) {
         // console.log(el, d, i, sources)
         download(d);
       };

     });

   }

   function cleanup() {
     var crowbarElements = document.querySelectorAll(".svg-crowbar");

     [].forEach.call(crowbarElements, function(el) {
       el.parentNode.removeChild(el);
     });
   }


   function getSources(doc, styles) {
     var svgInfo = [],
         svgs = doc.querySelectorAll("svg");

     styles = (styles === undefined) ? "" : styles;

     [].forEach.call(svgs, function (svg) {

       svg.setAttribute("version", "1.1");

       var defsEl = document.createElement("defs");
       svg.insertBefore(defsEl, svg.firstChild); //TODO   .insert("defs", ":first-child")
       // defsEl.setAttribute("class", "svg-crowbar");

       var styleEl = document.createElement("style")
       defsEl.appendChild(styleEl);
       styleEl.setAttribute("type", "text/css");


       // removing attributes so they aren't doubled up
       svg.removeAttribute("xmlns");
       svg.removeAttribute("xlink");

       // These are needed for the svg
       if (!svg.hasAttributeNS(prefix.xmlns, "xmlns")) {
         svg.setAttributeNS(prefix.xmlns, "xmlns", prefix.svg);
       }

       if (!svg.hasAttributeNS(prefix.xmlns, "xmlns:xlink")) {
         svg.setAttributeNS(prefix.xmlns, "xmlns:xlink", prefix.xlink);
       }

       var source = (new XMLSerializer()).serializeToString(svg).replace('</style>', '<![CDATA[' + styles + ']]></style>');
       var rect = svg.getBoundingClientRect();
       svgInfo.push({
         top: rect.top,
         left: rect.left,
         width: rect.width,
         height: rect.height,
         class: svg.getAttribute("class"),
         id: svg.getAttribute("id"),
         childElementCount: svg.childElementCount,
         source: [doctype + source]
       });
     });
     return svgInfo;
   }

   function download(source) {
     var filename = jsonFile;
      console.log("Y AQUI "+ filename);
     if (source.id) {
       filename = jsonFile;
     } else if (source.class) {
       filename = jsonFile;
     } else if (window.document.title) {
       filename = jsonFile;
     }

     var url = window.URL.createObjectURL(new Blob(source.source, { "type" : "text\/xml" }));

     var a = document.createElement("a");
     body.appendChild(a);
     a.setAttribute("class", "svg-crowbar");
     a.setAttribute("download", filename + ".svg");
     a.setAttribute("href", url);
     a.style["display"] = "none";
     a.click();

     setTimeout(function() {
       window.URL.revokeObjectURL(url);
     }, 10);
   }

   function getStyles(doc) {
     var styles = "",
         styleSheets = doc.styleSheets;

     if (styleSheets) {
       for (var i = 0; i < styleSheets.length; i++) {
         processStyleSheet(styleSheets[i]);
       }
     }

     function processStyleSheet(ss) {
       if (ss.cssRules) {
         for (var i = 0; i < ss.cssRules.length; i++) {
           var rule = ss.cssRules[i];
           if (rule.type === 3) {
             // Import Rule
             processStyleSheet(rule.styleSheet);
           } else {
             // hack for illustrator crashing on descendent selectors
             if (rule.selectorText) {
               if (rule.selectorText.indexOf(">") === -1) {
                 styles += "\\n" + rule.cssText;
               }
             }
           }
         }
       }
     }
     return styles;
   }
}




   </script>

);

print( $html );

}

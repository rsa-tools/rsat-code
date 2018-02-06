#!/usr/bin/perl
#### this cgi script fills the HTML form for the feature map cgi program

if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
### Default values for filling the new feature map form
### read values to fill the form ###
$default{title} = "";
$default{data} = "";
$default{format} = 'feature map';
$default{img_format} = $ENV{rsat_img_format} || "jpg";
$default{from} = 'auto';
$default{to} = 'auto';
$default{handle} = 'none';
$default{origin} = '0';
$default{map_len} = 500;
$default{spacing} = 2;
$default{thickness} = 25;
$default{file_name} = '';

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
    if ($query->param($key)) {
	$default{$key} = $query->param($key);
    }
}

if (-e $query->param('feature_file')) {
    $file = $query->param('feature_file');
    $default{data} = `cat $file`;
} else {
    $default{data} = $query->param('data');
}
##REVISAR ESTA LINEA
#$default{data} =~ s/\"//g; #### remove quotes for security reasons (avoid imbedded command)

&ListParameters() if ($ENV{rsat_echo} >= 2);

################################################################
### print the form ###

################################################################
### header
&RSA_header_bootstrap("feature map 2", "form");

##RO-REVISAR Creo que es para linkear form y html cgiss
print $query->start_multipart_form(-action=>"feature-map2.cgi");
print '
   <!-- Form with bootstrap -->
   <div class="container">
	   <div class="row">
         <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 bhoechie-tab-container">
            <div class="col-lg-2 col-md-3 col-sm-3 col-xs-3 bhoechie-tab-menu">
               <div class="list-group">

                  <a href="#" class="list-group-item active text-center">
                     <h4 class="glyphicon"><i class="fa fa-info-circle fa-2x"></i></h4><br/>Feature map 2
                  </a>

                  <a href="#" class="list-group-item text-center">
                     <h4 class="glyphicon"><i class="fa fa-tag fa-2x"></i></h4><br/>Mandatory inputs
                  </a>

                  <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-tasks fa-2x"></i></h4><br/>Advanced options
                </a>

                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-play-circle fa-2x"></i></h4><br/>Run analysis
                </a>

               </div>
            </div>

            <div class="col-lg-9 col-md-9 col-sm-9 col-xs-9 bhoechie-tab">
               <!-- ################################################################ -->
               <!-- ### info ### -->

               <div class="bhoechie-tab-content active">
                  <h2> <img src="images/RSAT_logo.jpg" style="max-width:150px;max-height:60px;padding-bottom:10px" alt="RSAT server" border="0"></img>
                  feature map 2</h2>

                  <span class="fa-stack fa-lg">
  							<i class="fa fa-info-circle fa-stack-1x"></i>
					   </span>
   					Generates a graphical map of features localized on one or several sequences.<br>

                  <span class="fa-stack fa-lg">
 							<i class="fa fa-user fa-stack-1x"></i>
					   </span>
					   Raul Ossio, Daniela Robles-Espinoza, Alejandra Medina-Rivera<br>

					   <span class="fa-stack fa-lg">
  					      <i class="fa fa-folder-open fa-stack-1x"></i>
					   </span>
					   Sample output<br>

					   <span class="fa-stack fa-lg">
  						   <i class="fa fa-book fa-stack-1x"></i>
					   </span>
					   <a class="iframe" href="help.feature-map.html">User Manual</a><br>

					   <span class="fa-stack fa-lg">
  						   <i class="fa fa-twitter fa-stack-1x"></i>
					   </span>
					   <a href="https://twitter.com/rsatools" target="_blank">Ask a question to the RSAT team</a><br>
                  <span class="fa-stack fa-lg">
  						   <i class="fa fa-laptop fa-stack-1x"></i>
					   </span>
					   <a href="http://embnet.ccg.unam.mx/rsat/feature-map_form.cgi" target="_blank">Use previous version of feature map</a><br>
               </div>

               <!-- ################################################################ -->
               <!-- ### mandatory inputs ### -->
               <div class="bhoechie-tab-content">
                  <!-- title -->
                  <div class="panel panel-danger">
                     <div class="panel-heading">Graph Title
                     </div>
                     <div class="panel-body">
                        <div class="form-group" id="fm_form">';

                           print
                              $query->textfield(-id=>'html_title',
                              -name=>'title', -class=>'form-control',-placeholder=>'Provide a Title for this graph ', -required=>'true',
                              -default=>$default{title}) .'

                        </div>
			            </div>
		            </div>


                  <!-- Input options -->
                  <div class="panel panel-danger">
                     <div class="panel-heading">List of features to display
                       <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Input here the list of features (relative coordinates) to be displayed on the map, you can either paste the list in the text box or upload it from your computer" data-original-title=""></i>
                     </div>
		               <div class="panel-body">
                        <div class="form-inline">
                           <div class="form-group">';

  						print "<label for='format'>Select a format</label>\n";
  						print "<A class='badge badge-primary iframe' HREF='help.feature-map.html#formats'>Info</a>";

                              print
                                 $query->popup_menu(-name=>'format',
                                 -Values=>['feature map',
                                 'dna-pattern',
                                 'Patser',
                                 'Matinspector',
                                 'Signal scan',
                                 'Swissprot',
                                 'Transfac factor',
                                 'Transfac gene',
                                 'GCG msf',
                                 'DSSP',
                                 'Gibbs sampler',
                                 'Fugue'],
                                 -default=>$default{format});

print '</div></div>
<div class="form-group">';

  print $query->textarea(-id=>'data',
    -name=>'data', -class=>'form-control',-placeholder=>'Paste here feature list, or select a file to upload below ',
			 -default=>$default{data},
			 -rows=>4,
			 -columns=>60,
			 -wrap=>'soft');


  ## Option to upload the matrix file from the client machine
  print "<b>Or</b> select a file to upload <BR>\n";
  print $query->filefield(-name=>'uploaded_file',
			  -default=>'',
			  -size=>45);
    print "<BR></div>\n";

                              print
                           '</div>
                     </div>
                  </div>

				<!-- ### Advanced options ### -->

                  <div class="bhoechie-tab-content">
                     <div class="panel panel-warning">
                        <div class="panel-heading">Title position
                        </div>

                        <div class="panel-body" style="font-weight:normal;">

                           <label class="radio-inline" style="color:black; font-weight:normal;">';

                              print $query->radio_group(-name=>'group_name',
                              -values=>["left"],
                              -default=>'center'
                              );

                           print '
                           </label>
                           <label class="radio-inline">';

                              print $query->radio_group(-name=>'group_name',
                              -values=>["center"],
                              -default=>'center'
                              );

                              print '
                           </label>
                           <label class="radio-inline">';

                              print $query->radio_group(-name=>'group_name',
                              -values=>["right"],
                              -default=>'center',
                              );

                              print '
                           </label>
                        </div>
                     </div>

                     <div class="panel panel-warning">
                        <div class="panel-heading">Display options
                        </div>
                        <div class="panel-body">
                           <label class="checkbox-inline">';

                              print $query->checkbox(
                                   -name    => 'view_scalebar',
                                   -checked => 1,
                                   -value   => 'true',
                                   -label   => 'View scalebar',
                               );

                              print'
                           </label>&nbsp&nbsp
                           <label class="checkbox-inline">';

                              print $query->checkbox(
                                   -name    => 'view_score',
                                   -checked => 1,
                                   -value   => 'true',
                                   -label   => 'View score',
                               );

                           print '
                           </label>
                        </div>
                     </div>

                     <div class="panel panel-warning">
                        <div class="panel-heading">Legends
                        </div>

                        <div class="panel-body" style="font-weight:normal;">
                           <label class="checkbox-inline">';

                           print $query->checkbox(
                                -name    => 'view_legends',
                                -checked => 1,
                                -value   => 'true',
                                -label   => 'View legends',
                            );


                           print '
                              </label>&nbsp&nbsp
                              <label class="radio-inline" style="color:black; font-weight:normal;">
                           ';

                           print $query->radio_group(-name=>'legends',
                              -values=>["left"],
                              -default=>'left',
                           );

                           print '
                           </label>
                           <label class="radio-inline">';

                           print $query->radio_group(-name=>'legends',
                              -values=>["top"],
                              -default=>'top',
                           );


                           print '
                           </label>
                        </div>
                     </div>
                  </div>

				 <!--################################################################-->
				 <!--### output & run ###-->
                  <div class="bhoechie-tab-content">
                     <div class="panel panel-info">
                        <div class="panel-heading">Output options
                        </div>
                        <div class="panel-body">
                           <div class="form-group">
                           ';

                           print

                           # $query->textfield(-id=>'file_name',
                           # -name=>'title', -class=>'form-control',-placeholder=>'Provide a name for this file ',
                           # -default=>$default{file_name})

                           .'

                           </div>
 							<div class="form-group">
 							<label>File type </label>

                        <label class="radio-inline">';


                           print $query->radio_group(-name=>'print_file',
                           -values=>["none"],
                           -default=>'none'
                           );

                        print '
                        </label>
                        <label class="radio-inline" style="color:black; font-weight:normal;">';

                        print $query->radio_group(-name=>'print_file',
                        -values=>["svg"],
                        -default=>'none'
                        );

                        print '
                        </label>
                        <label class="radio-inline">';

                           print $query->radio_group(-name=>'print_file',
                           -values=>["png"],
                           -default=>'none'
                           );

                           print '
                        </label>
                        <label class="radio-inline">';

                           print $query->radio_group(-name=>'print_file',
                           -values=>["jpg"],
                           -default=>'none'
                           );

                           print '
                           </label>
                     </div>
                  </div>
                  </div>';

                            print $query->submit(-label=>"GO", -class=>"btn btn-success", id=>"go_button", -type=>"button");
                            print " ";
							print $query->reset(-id=>"reset",-class=>"btn btn-warning", -type=>"button");
print '
                        </div>
                  </div>
               </div>';

################################################################
## Demo area
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr' class='col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container'></div>";

                        print '

                        <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container">
                  <button type="button" class="btn btn-info" id="demo_button" onclick="demo_fill()">DEMO</button>
         </div>';


            my $demo_data = "";
            open ($fh, "demo_files/feature-map_demo_data.ft");
            while($row = <$fh>){
                chomp $row;
                $demo_data .= $row;
                $demo_data .= "\\n";
            }

            print'
            <script>

               function demo_fill(){
                  document.getElementById("html_title").value = "Motifs discovered in upstream sequences of 19 MET genes";
                  document.getElementById("data").value ="' .$demo_data .'";
                  document.getElementById("file_name").value = "feature_map";
                  document.getElementById("view_scalebar").checked = "true";
                  document.getElementById("view_legends").checked = "true";
                  document.getElementById("view_scores").checked = "true";
                  document.getElementById("title_center").checked = "true";
                  descr_demo = "<blockquote class =\'demo_3 blockquote text-justify small\'>\
     Check the panel <b>Mandatory inputs</b> and then <b>Run analysis</b></p>\n</blockquote>";
                  demo_descr.innerHTML = descr_demo;
                  document.getElementById("demo") = descr_demo;
               }

            </script>';

            print $query->end_html;
            exit(0);

#
################################################################
#################### SUBROUTINE DEFINITIONS  ###################
################################################################

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

################################################################
### print the form ###

################################################################
### header
&RSA_header_bootstrap("feature map", "form");

##RO-REVISAR Creo que es para linkear form y html cgiss
print $query->start_multipart_form(-action=>"feature-map2.cgi");
print '
   <!-- Form with bootstrap -->
   <div class="container">
	   <div class="row">
         <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 bhoechie-tab-container">
            <div class="col-lg-2 col-md-3 col-sm-3 col-xs-3 bhoechie-tab-menu">
               <div class="list-group">

                  <a href="#" class="list-group-item active text-center" onclick="hide_buttons()">
                     <h4 class="glyphicon"><i class="fa fa-info-circle fa-2x"></i></h4><br/>Feature map
                  </a>

                  <a href="#" class="list-group-item text-center" onclick="unhide_buttons()">
                     <h4 class="glyphicon"><i class="fa fa-tag fa-2x"></i></h4><br/>Input options<br>( mandatory )
                  </a>

                  <a href="#" class="list-group-item text-center onclick="unhide_buttons()">
                     <h4 class="glyphicon"><i class="fa fa-tags fa-2x"></i></h4><br/>Export options
                  </a>

                  <a href="#" class="list-group-item text-center onclick="unhide_buttons()">
                     <h4 class="glyphicon"><i class="fa fa-tags fa-2x"></i></h4><br/>Advanced options
                  </a>

               </div>
            </div>

            <div class="col-lg-9 col-md-9 col-sm-9 col-xs-9 bhoechie-tab">
               <!-- ################################################################ -->
               <!-- ### info ### -->

               <div class="bhoechie-tab-content active">
                  <h2> <img src="images/RSAT_logo.jpg" style="max-width:150px;max-height:60px;padding-bottom:10px" alt="RSAT server" border="0"></img>
                  feature map</h2>

                  <span class="fa-stack fa-lg">
  							<i class="fa fa-info-circle fa-stack-1x"></i>
					   </span>
   					Generates a graphical map of features localized on one or several sequences.<br>

                  <span class="fa-stack fa-lg">
 							<i class="fa fa-user fa-stack-1x"></i>
					   </span>
					   Raul Ossio, Daniela Robles-Espinoza, Morgane Thomas-Chollier<br>

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
                  <div class="panel panel-primary">
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
                  <div class="panel panel-primary">
                     <div class="panel-heading">Input options
                     </div>
		               <div class="panel-body">
                        <div class="form-inline">
                           <div class="form-group">';

                              print
                                 "<b>1. Select a file to upload</b> <BR>\n";

                              print
                                 $query->filefield(-name=>'uploaded_file',
                                 -default=>'',
                                 -size=>45);

                              print
                                 "</div><br>or<br>\n";


                              print
                                 "<b>2. Provide a feature list</b><BR>\n";


                              print
                                 '<label for="format">Select a format  &nbsp&nbsp</label>';

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

                              print
                                 "<BR>\n";

                              print
                                 $query->textarea(-id=>'data',
                                 -name=>'data', -class=>'form-control',-placeholder=>'Paste here your feature list',
             			           -default=>$default{data},
             			           -rows=>4,
             			           -columns=>70,
                                 -wrap=>'soft');


                              print
                           '</div>
                        </div>
                     </div>
                  </div>


                  <!-- ################################################################ -->
                  <!-- ### Export options ### -->
                  <div class="bhoechie-tab-content">
                     <div class="panel panel-primary">
                        <div class="panel-heading">File name
                        </div>
                        <div class="panel-body">
                           <div class="form-group">';

                           print
                           $query->textfield(-id=>'file_name',
                           -name=>'title', -class=>'form-control',-placeholder=>'Provide a name for this file ',
                           -default=>$default{file_name}) .'

                           </div>
                        </div>
                     </div>

                     <div class="panel panel-primary">
                        <div class="panel-heading">File type
                        </div>
                        <div class="panel-body">
                        <label class="radio-inline" style="color:black; font-weight:normal;">';

                           print $query->radio_group(-name=>'print_file',
                           -values=>["svg"],
                           -default=>'svg'
                           );

                        print '
                        </label>
                        <label class="radio-inline">';

                           print $query->radio_group(-name=>'print_file',
                           -values=>["png"],
                           -default=>'svg'
                           );

                           print '
                        </label>
                        <label class="radio-inline">';

                           print $query->radio_group(-name=>'print_file',
                           -values=>["jpg"],
                           -default=>'svg'
                           );

                           print '
                           </label>
                        </div>
                     </div>
                  </div>';

                  print '
                  <div class="bhoechie-tab-content">
                     <div class="panel panel-primary">
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

                     <div class="panel panel-primary">
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

                     <div class="panel panel-primary">
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

                  <div class="container">
                     <div class="row">
                        <div class="col-lg-6 col-md-5 text-center">';

                           print
                           $query->submit(-label=>"GO", -class=>"btn btn-success", id=>"go_button", -style=>"display:none");

                           print '
                              <button type="button" class="btn btn-info" id="demo_button" onclick="demo_fill()" style="display:none">DEMO</button>
                              <button type="button" class="btn btn-warning" id="reset_button" onclick="reset()" style="display:none" >RESET</button>
                        </div>

                        <div class="row"><BR></div>
                        <div class="row"><BR></div>';

                        print'
                     </div>
                  </div>
               </div>';

                        print '
            </div>
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


                  ';


               print'
               }

               function reset(){
                  document.getElementById("html_title").value = "";
                  document.getElementById("data").value = "";
                  document.getElementById("file_name").value = "";
                  document.getElementById("Legends_top").checked = "true";
                  document.getElementById("Legends_bottom").checked = "false";
                  document.getElementById("Legends_left").checked = "false";
                  document.getElementById("Legends_right").checked = "false";
                  document.getElementById("view_scalebar").checked = "true";
                  document.getElementById("view_legends").checked = "true";
                  document.getElementById("view_scores").checked = "true";
                  document.getElementById("title_left").checked = "false";
                  document.getElementById("title_center").checked = "true";
                  document.getElementById("title_right").checked = "false";

               }
               function unhide_buttons(){
                  var go = document.getElementById("go_button");
                  go.style.display = "inline";

                  var x = document.getElementById("demo_button");
                  x.style.display = "inline";

                  var y = document.getElementById("reset_button");
                  y.style.display = "inline";
               }

               function hide_buttons(){
                  document.getElementById("go_button").style.display = "none";
                  document.getElementById("demo_button").style.display = "none";
                  document.getElementById("reset_button").style.display = "none";
               }


            </script>';

            print $query->end_html;
            exit(0);


################################################################
#################### SUBROUTINE DEFINITIONS  ###################
################################################################

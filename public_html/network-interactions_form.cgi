#!/usr/bin/env perl
#### this cgi script fills the HTML form for the program network-interactions
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
    require "RSA.lib";
}
#if ($0 =~ /([^(\/)]+)$/) {
#    push (@INC, "$`lib/");
#}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";
use RSAT::matrix;
use RSAT::MatrixReader;

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

### Read the CGI query
$query = new CGI;

################################################################
### default values for filling the form
$default{html_title} = "";
$default{tf_selection} = "";
$default{cre_selection} = "";
$default{genome_v} = "";
$default{net_selection} = "";
$default{matrixscan_score} = 5;

### replace defaults by parameters from the cgi call, if defined
foreach $key (keys %default) {
  if ($query->param($key)) {
    $default{$key} = $query->param($key);
  }
  if ($query->param($key) =~ /checked/i) {
    $checked{$key} = "CHECKED";
  }
}

&ListParameters() if ($ENV{rsat_echo} >= 2);

################################################################
### print the form ###

################################################################
### header
&RSA_header_bootstrap("network-interactions", "form");

print $query->start_multipart_form(-action=>"network-interactions.cgi");

print '
 <!-- Form with bootstrap -->
<div class="container">
	<div class="row">
        <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 bhoechie-tab-container">
            <div class="col-lg-2 col-md-3 col-sm-3 col-xs-3 bhoechie-tab-menu">
              <div class="list-group">
                <a href="#" class="list-group-item active text-center">
                  <h4 class="glyphicon"><i class="fa fa-info-circle fa-2x"></i></h4><br/>Network Interactions
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-tag fa-2x"></i></h4><br/>Mandatory inputs
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-tasks fa-2x"></i></h4><br/>Mandatory options
                </a>
                <a href="#" class="list-group-item text-center">
                  <h4 class="glyphicon"><i class="fa fa-tags fa-2x"></i></h4><br/>Optional inputs
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
                      network-interactions</h2>
                    <span class="fa-stack fa-lg">
  							<i class="fa fa-info-circle fa-stack-1x"></i>
					</span>
          Constructs and compares Gene Regulatory Networks based on pattern-matching. <br>
                    <span class="fa-stack fa-lg">
 							 <i class="fa fa-user fa-stack-1x"></i>
					</span>

					Monica Padilla Galvez and <a target="_blank" href="http://liigh.unam.mx/amedina/index.html">Alejandra Medina Rivera</a>.<br>

					<span class="fa-stack fa-lg">
  							<i class="fa fa-book fa-stack-1x"></i>
					</span>
					<a class="iframe" href="help.network-interactions.html">User Manual</a><br>

					<span class="fa-stack fa-lg">
  							<i class="fa fa-graduation-cap fa-stack-1x"></i>
					</span>
					<a class="iframe" href="tutorials/tut_network-interactions.html">Tutorial</a><br>

					<span class="fa-stack fa-lg">
  							<i class="fa fa-twitter fa-stack-1x"></i>
					</span>
					<a href="https://twitter.com/rsatools" target="_blank">Ask a question to the RSAT team</a><br>

          <! publication pencil >

                </div>

 <!-- ################################################################ -->
 <!-- ### mandatory inputs ### -->
<div class="bhoechie-tab-content">

  <!-- Title -->
  <div class="panel panel-danger">
    <div class="panel-heading">Analysis Title
      <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Title that will be displayed at the top of the report page." data-original-title=""></i>
    </div>

    <div class="panel-body">
      <div class="form-group">';
        print $query->textfield(-id=>'html_title',-name=>'html_title', -class=>'form-control',-placeholder=>'Provide a Title for this analysis ', -required=>'true',
  			 -default=>$default{html_title}) .'
  	   </div>
  	</div>
  </div>

  <!-- TFs -->
  <div class="panel panel-danger">
    <div class="panel-heading">Transcription Factors
      <i class="fa fa-info-circle" rel="popover" data-container="body" data-trigger="hover" data-placement="right" data-content="Transcription Factors list needed to construct the Gene Regulatory Network, you can either write the TFs list or upload it from your computer."></i>
    </div>

    <div class="panel-body">
      <div class="form-group">';

        print "Specify TFs below\n";
        print $query->textarea(-id=>"tf_selection", -name=>"tf_selection", -default=>$default{tfs_selection}, -rows=>5, -columns=>60, class=>"form-control");

        ### option to upload a file with the gene list from the client machine
        print "<BR>Or upload TFs list from file<BR>\n";
        print $query->filefield(-id=>'uploaded_tf_file',-name=>'uploaded_tf_file',-default=>'',-size=>45,-maxlength=>200);

  print '</div></div></div>

  <!-- Regulatory Sequences BED File -->
  <div class="panel panel-danger">
    <div class="panel-heading"> Regulatory Sequences BED File
      <i class="fa fa-info-circle" rel="popover" data-container="body" data-trigger="hover" data-placement="right" data-content="BED file indicating the regulatory sequences coordinates for all genes (on fourth column) of interest, including TFs. You can either paste the BED in the text box or upload it from your computer."></i>
    </div>

    <div class="panel-body">
      <div class="form-group">';

        print " Specify BED file below\n";
        print $query->textarea(-id=>"cre_selection", -name=>"cre_selection", -default=>$default{cre_selection}, -rows=>5, -columns=>60, class=>"form-control");

        ### option to upload a file with the gene list from the client machine
        print "<BR>Or upload BED file<BR>\n";
        print $query->filefield(-id=>'uploaded_cre_file',-name=>'uploaded_cre_file',-default=>'',-size=>45,-maxlength=>200);

    print '</div></div></div>

</div>
 <!-- ################################################################ -->
 <!-- ### mandatory options ### -->
 <div class="bhoechie-tab-content">

 <!-- Genome version -->
 <div class="panel panel-danger">
   <div class="panel-heading">Genome
     <i class="fa fa-info-circle" rel="popover" data-container="body" data-trigger="hover" data-placement="right" data-content="Genome version according to coordinates in BED file."></i>
   </div>

   <div class="panel-body">
     Select a genome version<br>
     <div class="form-group">';

     print &UCSCGenomePopUpSelectable('genome_v', 'genome_v');
     print '

     </div>
   </div>
 </div>
  <!-- Directory Name -->
 <div class="panel panel-danger">
   <div class="panel-heading"> Directory Name
     <i class="fa fa-info-circle" data-container="body" data-toggle="tooltip" data-placement="top" title="Directory Name" data-original-title=""></i>
   </div>

   <div class="panel-body">
     <div class="form-group">';
       print $query->textfield(-id=>'dir_name',-name=>'dir_name', -class=>'form-control',-placeholder=>'Provide a directory name for analysis output ', -required=>'true',
        -default=>$default{dir_name}) .'
      </div>
   </div>
 </div>
</div>

 <!-- ################################################################-->
 <!-- ### optional inputs ###-->
 <div class="bhoechie-tab-content">

   <!-- Network -->
   <div class="panel panel-warning">
     <div class="panel-heading">Network
      <i class="fa fa-info-circle" rel="popover" data-container="body" data-trigger="hover" data-placement="right" data-content="Previously done network to allow its comparison to the one generated by the program."></i>
     </div>

     <div class="panel-body">
       <div class="form-group">';

         print " Specify the network below\n";
         print $query->textarea(-id=>"net_selection", -name=>"net_selection", -default=>$default{net_selection}, -rows=>5, -columns=>60, class=>"form-control");

         ### option to upload a file with the gene list from the client machine
         print "<BR>Or upload the network from file<BR>\n";
         print $query->filefield(-id=>'uploaded_net_file', -name=>'uploaded_net_file',-default=>'',-size=>45,-maxlength=>200);

         print '
      </div>
    </div>
   </div>
 </div>

 <!-- ################################################################-->
 <!-- ### advanced options ###-->

 <!-- ADVANCED OPTIONS -->
 <div class="bhoechie-tab-content">

 <div class="panel panel-warning">
 <div class="panel-heading"><em>matrix-scan</em> parameters</div>
 <div class="panel-body"> <br>';

print " Score lower threshold \n";
print '
 <div class="panel-body">
   <div class="form-group">';
     print $query->textfield(-id=>'matrixscan_score',-name=>'matrixscan_score', -class=>'form-control',-placeholder=>'5', -required=>'true',
      -default=>$default{matrixscan_score}) .'
    </div>
 </div>';

 print " P-value upper threshold \n";
 print '
  <div class="panel-body">
    <div class="form-group">';
      print $query->textfield(-id=>'matrixscan_pval',-name=>'matrixscan_pval', -class=>'form-control',-placeholder=>'0.0005', -required=>'true',
       -default=>$default{matrixscan_pval}) .'
     </div>
  </div>';

 print '</div>
 </div>

 <div class="panel panel-warning">
 <div class="panel-heading">Motif Database</div>
 <div class="panel-body"> <br>';

my %args = ();
#$args{space} = 1;

print &MotifSelection(%args);

 print '</div>
 </div>

 </div>

 <!--################################################################-->
 <!--### output & run ###-->

<div class="bhoechie-tab-content">

 <!-- ## Specific options for output files-->
  <div class="panel panel-info">
   <div class="panel-heading">Output options
   </div>

  <div class="panel-body">
    <div class="form-group">';


################################################################
## Send results by email only
print "<p>\n";
&SelectOutput("display");
print " </div>
  </div> </div>";


################################################################
## Action buttons
#print "<TABLE class='formbutton'>\n";
#print "<TR VALIGN=MIDDLE>\n";
#print "<TD>", $query->submit(-label=>"GO", -class=>"btn btn-success"), "</TD>\n";
print $query->submit(-label=>"GO", -class=>"btn btn-success", -type=>"button");
print " ";
print $query->reset(-id=>"reset",-class=>"btn btn-warning", -type=>"button");
print $query->end_form;

 print ' </div>
        </div>
  </div>
</div>
';

################################################################
## Demo area
print "<textarea id='demo' style='display:none'></textarea>";
print "<div id='demo_descr' class='col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container'></div>";



#
################################################################
## Demo 1 data

$tfs_file = "demo_files/TFs_short.txt";
$tfs = "";
open(my $fh1, $tfs_file);
while(my $row = <$fh1>){
    chomp $row;
    $tfs .= $row;
    $tfs .= "\\n";
}

$reg_file = "demo_files/enhancers_short.bed";
$bed = "";
open(my $fh2, $reg_file);
while(my $row = <$fh2>){
    chomp $row;
    $bed .= $row;
    $bed .= "\\n";
}

$net_file = "demo_files/net_short.tsv";
$net = "";
open(my $fh3, $net_file);
while(my $row = <$fh3>){
    chomp $row;
    $net .= $row;
    $net .= "\\n";
}

print '<script>
function setDemo1(tfs, bed, net){
    $("#reset").trigger("click");
    descr_1 = "<H4>Demonstration: Network Comparison</H4>\n \
    <blockquote class =\'blockquote text-justify small\'>\
    In this demo, we will reconstruct a Gene Regulatory Network \
    for a list of genes and transcription factors expressed during \
    chicken ear development and compare it with experimentally validated \
    TF-gene interactions. Data source: Anwar, M. et al., 2017. \
    Check the panel <b>Mandatory inputs</b>, \
    the panel <b>Mandatory options</b>, \
    the panel <b>Optional inputs</b>, \
    and then <b>Run analysis</b></blockquote>";

    demo_descr.innerHTML = descr_1;
    html_title.value = "\'Demo 01\'";
    matrixscan_pval.value = "0.0005";
    genome_v.value = "galGal4";
    dir_name.value = "Demo01";
    tf_selection.value = tfs;
    cre_selection.value = bed;
    net_selection.value = net;
    demo.value = descr_1;
}
</script>';

print ' <div class="col-lg-9 col-md-5 col-sm-8 col-xs-9 demo-buttons-container">

<button type="button" class="btn btn-info" onclick="setDemo1('. "'$tfs', '$bed', '$net'" .')">DEMO</button> ';


################################################################
print $query->end_html;

exit(0);



################################################################
#################### SUBROUTINE DEFINITIONS  ###################
################################################################

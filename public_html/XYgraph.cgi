#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "cgi-lib.pl";
require "RSA.lib.pl";
require "RSA_cgi_lib.pl";
$XYgraph_command = "$SCRIPTS/XYgraph";
$tmp_file_name = sprintf "XYgraph.%s", &AlphaDate;

MAIN:
{
  ### Read the content of the form 
  &ReadParse(*input);

  #### update log file ####
  &UpdateLogFile;

  #### read parameters ####
  $parameters = "-v ";
  
  ### general parameters ###
  if ($input{'title1'} ne "") {
    $parameters .= "-title1 \"$input{'title1'}\" ";
  }
  if ($input{'title2'} ne "") {
    $parameters .= "-title2 \"$input{'title2'}\" ";
  }
  if ($input{'pointsize'} =~ /\d+/) {
    $parameters .= "-pointsize \"$input{'pointsize'}\" ";
  }
  if ($input{'lines'} eq "yes") {
    $parameters .= "-lines ";
  }
  if ($input{'symbols'} eq "yes") {
    $parameters .= "-symbols ";
  }
  if ($input{'legend'} eq "yes") {
    $parameters .= "-legend ";
  }
  if (&IsNatural($input{'label_col'})) {
    $parameters .= "-lc $input{'label_col'} ";
  }

  $accepted_bg{'black'} = 1;
  $accepted_bg{'white'} = 1;
  $accepted_bg{'blue'} = 1;
  $accepted_bg{'gray'} = 1;
  if ($accepted_bg{$input{'bg'}}) {
    $parameters .= "-bg $input{'bg'} ";
  }


  
  ### X axis parameters ###
  if ($input{'xcol'} ne "") {
    $parameters .= "-xcol $input{'xcol'} ";
  }
  if ($input{'xlog'} eq "yes") {
    $parameters .= "-xlog ";
  }
  if ($input{'xleg1'} ne "") {
    $parameters .= "-xleg1 \"$input{'xleg1'}\" ";
  }
  if ($input{'xleg2'} ne "") {
    $parameters .= "-xleg2 \"$input{'xleg2'}\" ";
  }
  if ($input{'xmin'} ne "") {
    $parameters .= "-xmin $input{'xmin'} ";
  }
  if ($input{'xmax'} ne "") {
    $parameters .= "-xmax $input{'xmax'} ";
  }
  if ($input{'xgstep1'} ne "") {
    $parameters .= "-xgstep1 $input{'xgstep1'} ";
  }
  if ($input{'xgstep2'} ne "") {
    $parameters .= "-xgstep2 $input{'xgstep2'} ";
  }
  if ($input{'xsize'} ne "") {
    $parameters .= "-xsize $input{'xsize'} ";
  }

  ### Y axis parameters ###
  if ($input{'ycol'} ne "") {
    $parameters .= "-ycol $input{'ycol'} ";
  }
  if ($input{'ylog'} eq "yes") {
    $parameters .= "-ylog ";
  }
  if ($input{'yleg1'} ne "") {
    $parameters .= "-yleg1 \"$input{'yleg1'}\" ";
  }
  if ($input{'yleg2'} ne "") {
    $parameters .= "-yleg2 \"$input{'yleg2'}\" ";
  }
  if ($input{'ymin'} ne "") {
    $parameters .= "-ymin $input{'ymin'} ";
  }
  if ($input{'ymax'} ne "") {
    $parameters .= "-ymax $input{'ymax'} ";
  }
  if ($input{'ygstep1'} ne "") {
    $parameters .= "-ygstep1 $input{'ygstep1'} ";
  }
  if ($input{'ygstep2'} ne "") {
    $parameters .= "-ygstep2 $input{'ygstep2'} ";
  }
  if ($input{'ysize'} ne "") {
    $parameters .= "-ysize $input{'ysize'} ";
  }

  ### data file ####
  if ($input{'data_file'} ne "") {
    ### internal data file ###
    $parameters .= "-i $input{'data_file'} ";
  } else {
      unless ($input{'data'} =~ /\S/) {
	print &PrintHeader;
	&cgiError("Error: the data box should not be empty.");
      }
      $data_file = "$tmp_file_name.data";
      open DATA, ">$TMP/$data_file";
      print DATA $input{'data'};
      close DATA;
      $parameters .= "-i $TMP/$data_file ";
  }
  
  ### graph file ###
  $graph_file = "$tmp_file_name.gif";
  $parameters .= "-o $TMP/$graph_file ";


  if ($input{'htmap'} eq "yes") {
    $htmap = 1;
    $htmap_file = "$tmp_file_name.html";
    $parameters .= "-htmap ";
    $parameters .= "-htmap > $TMP/$htmap_file ";
  }
  
  ### execute the command ###
  @data_report = `$XYgraph_command $parameters`;

  ### print the result ###
  if ($htmap) {
    print "Location: $WWW_TMP/$htmap_file", "\n\n";
  } else {
    ### display the result ###
	print &PrintHeader;
    print <<End_Header;
<HEADER>
<TITLE>RSA-tools - XY graph result</TITLE>
</HEADER><BODY BGCOLOR="#FFFFFF">
<H3 ALIGN=CENTER><A HREF="$WWW_RSA/RSA_home.cgi">
RSA-tools</A> - XY graph result</H3>
End_Header
    print "<CENTER><IMG SRC=\"$WWW_TMP/$graph_file\"></CENTER><P>\n";
    print "<H4 ALIGN=CENTER>Data report</H4>";
    print "<PRE>";
    print @data_report;
    print "</PRE>";
    print "<HR SIZE = 3>";
    print &HtmlBot;
  }    
  
  DelayedRemoval("$TMP/$graph_file");
  DelayedRemoval("$TMP/$data_file");

  exit(0);
}



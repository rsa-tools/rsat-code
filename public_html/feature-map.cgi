#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib.pl";
require "RSA.cgi.lib.pl";

### intialization
$feature_map_command = "$SCRIPTS/feature-map";
$tmp_file_name = sprintf "feature-map.%s", &AlphaDate;

$features_from_swissprot_cmd = "$SCRIPTS/features-from-swissprot";
$features_from_msf_cmd = "$SCRIPTS/features-from-msf";
$features_from_gibbs_cmd = "$SCRIPTS/features-from-gibbs";
$features_from_fugue_cmd = "$SCRIPTS/features-from-fugue";
$features_from_dssp_cmd = "$SCRIPTS/features-from-dssp";
$features_from_matins_cmd = "$SCRIPTS/features-from-matins";
$features_from_sigscan_cmd = "$SCRIPTS/features-from-sigscan";
$features_from_dnapat_cmd = "$SCRIPTS/features-from-dnapat";
$features_from_tffact_cmd = "$SCRIPTS/features-from-tffact";
$features_from_tfgene_cmd = "$SCRIPTS/features-from-tfgene";
$features_from_patser_cmd = "$SCRIPTS/features-from-patser";


### Read the CGI query
$query = new CGI;

#### this cgi script can perform either of two actions :
### fill the form or execute the action requested by that form
if (($query->param()) && !($query->param('fill_form'))){
  &ExecFeatureMap;
} else {
  &FillFeatureMapForm;
}


exit(0);

######### SUBROUTINE DEFINITION ###############


####### Execute feature-map perl script #########
sub ExecFeatureMap {
  #### execute the script
  $title = "feature-map result";
  
  
  #### update log file ####
  &UpdateLogFile;
  
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
  
  ### horizontal format ###
  if ($query->param('orientation') =~ /vertic/i) {
    $parameters .= " -vertical ";
  } else {
    $parameters .= " -horizontal ";
  }
  
  
  ### scale bar step ###
  if ($query->param('scalebarstep') =~ /\d+/) {
    $parameters .= " -scalestep ".$query->param('scalebarstep');
  }
  
  ### title ###
  if ($query->param('title') ne "") {
    $title = $query->param('title');
    $title =~ s/\'//g;
    $query->param('title',$title);
    $parameters .= " -title \'".$query->param('title'). "\'";
  }
  
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
  
  ### handle ###
  if (lc($query->param('palette')) =~ /mono/i) {
    $parameters .= " -mono ";
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

  ### data file ####
  if ($query->param('feature_file') =~ /\S/) {
    ### file on the server
    $feature_file = $query->param('feature_file');
  } elsif (($query->param('uploaded_file')) ||
	   ($query->param('data') =~ /\S/)) {
    $feature_file = "$TMP/$tmp_file_name.ft";
    ### convert data towards feature-map format
    if ($query->param('format') =~ /swiss/i) {
      open DATA, "| $features_from_swissprot_cmd -o $feature_file";
    } elsif ($query->param('format') =~ /transfac factor/i) {
      open DATA, "| $features_from_tffact_cmd -o $feature_file";
    } elsif ($query->param('format') =~ /transfac gene/i) {
      open DATA, "| $features_from_tfgene_cmd -o $feature_file";
    } elsif ($query->param('format') =~ /msf/i) {
      open DATA, "| $features_from_msf_cmd -o $feature_file";
      $parameters .= " -aacolors ";
      $parameters .= " -horiz ";
    } elsif ($query->param('format') =~ /dssp/i) {
      open DATA, "| $features_from_dssp_cmd -o $feature_file";
    } elsif ($query->param('format') =~ /matins/i) {
      open DATA, "| $features_from_matins_cmd -o $feature_file";
    } elsif ($query->param('format') =~ /dna\-pattern/i) {
      open DATA, "| $features_from_dnapat_cmd -o $feature_file";
    } elsif ($query->param('format') =~ /patser/i) {
      open DATA, "| $features_from_patser_cmd -o $feature_file";
    } elsif ($query->param('format') =~ /signal scan/i) {
      open DATA, "| $features_from_sigscan_cmd -o $feature_file";
    } elsif ($query->param('format') =~ /gibbs/i) {
      open DATA, "| $features_from_gibbs_cmd -o $feature_file";
    } elsif ($query->param('format') =~ /fugue/i) {
      open DATA, "| $features_from_fugue_cmd -o $feature_file";
    } else {
      open DATA, ">$feature_file";
    }
    
    if ($query->param('uploaded_file')) {
      ### upload file from the client
      $fh = $query->param('uploaded_file');
      while (<$fh>) {
	print DATA;
      }
    } else {
      ### data from the textarea
      print DATA $query->param('data');
    }
    close DATA;
  } else {
    print $query->header();
    &cgiError("Error: the feature list should not be empty.<P>Read on-line manual for more information.");
  }
  
  
  $parameters .= " -i $feature_file ";
  &DelayedRemoval($feature_file);
  
  ### map file ###
  $map_file = "$tmp_file_name.gif";
  $html_file = "$tmp_file_name.html";
  $parameters .= " -o $TMP/$map_file > $TMP/$html_file";
  DelayedRemoval("$TMP/$map_file");
  DelayedRemoval("$TMP/$html_file");
  
  ### executre the command
  system "$feature_map_command $parameters ";
  
  ### display the result ###
  if (lc($query->param('htmap')) eq "on") {
    $location = "$WWW_RSA/tmp/$html_file";
  } else {
    $location = "$WWW_RSA/tmp/$map_file";
  }
  print "Location: $location", "\n\n";

  ### debug only
  #print $query->header();
  #print $query->start_html;
  #print "<PRE>command = $feature_map_command $parameters \n</PRE>";
  #&ListParameters;
  #print $query->end_html;


  exit(0);
}



######## generate and fill the feature-map form #########

sub FillFeatureMapForm {

  ### read values to fill the form ###
  $title = $query->param('title');
  $title =~ s/\"//g;


  if (-e $query->param('feature_file')) {
    $file = $query->param('feature_file');
    $data = `cat $file`;
  } else {
    $data = $query->param('data');
  }
  $data =~ s/\"//g;


  if ($query->param('format')) {
    $format = $query->param('format');
  } else {
    $format = 'feature map';
  }

  if ($query->param('from')) {
    $from = $query->param('from');
  } else {
    $from = 'auto';
  }

  if ($query->param('to')) {
    $to = $query->param('to');
  } else {
    $to = 'auto';
  }

  if ($query->param('origin')) {
    $origin = $query->param('origin');
  } else {
    $origin = '0';
  }


  ### print the form ###
  print $query->header;
  print $query->start_html(-title=>'RSA-tools : feature-map',
                            -author=>'jvanheld@ucmb.ulb.ac.be',
                            -base=>'true',
                            -meta=>{'keywords'=>['regulatory', 'sequence', 'analysis', 'feature-map'], 
                                    'copyright'=>'copyright 1999 Jacques van Helden'},
                            -BGCOLOR=>'#FFEEDD');
  print "<CENTER>";
  print "<FONT SIZE=+1><B>RSA-Tools : feature-map</B></FONT><BR>\n";
  print "Generates a physical map of genetic features for one or several sequences<P>\n";
  print "</CENTER>";

  print "<FONT FACE='Helvetica'>";

  print $query->start_multipart_form;
  print "<B>Feature list</B>&nbsp;&nbsp;&nbsp;&nbsp;(";
  print "<A HREF='help.feature-map.html#formats'>";
  print "Format</A>&nbsp;";

  print $query->popup_menu(-name=>'format',
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
			   -default=>"$format");
  print ")<BR>\n";


  print $query->textarea(-name=>'data',
			 -default=>"$data",
			 -rows=>6,
			 -columns=>60);
  print "<BR>\n";

  print "<B><A HREF='help.feature-map.html#file'>";
  print "File ";
  print "</A></B>&nbsp;";
  print $query->filefield(-name=>'uploaded_file',
			  -size=>45);
  

  print "<B><A HREF='help.feature-map.html#title'>";

  print "<P>";
  print "Title";
  print "</A></B>&nbsp;";
  print $query->textfield(-name=>'title',
                            -default=>"$title",
                            -size=>50,
                            -maxlength=>80);

  print "<BR>\n";
  print "<A HREF='help.feature-map.html#legend'>";
  print $query->checkbox(-name=>'legend',
			 -checked=>'checked',
			 -label=>' Legend ');
  print "</A>";
  print "&nbsp;&nbsp;&nbsp;";

  print "<A HREF='help.feature-map.html#scalebar'>";
  print $query->checkbox(-name=>'scalebar',
			 -checked=>'checked',
			 -label=>' Scalebar ');
  print "</A>";

  print "(step ";
  print $query->textfield(-name=>'scalebarstep',
			  -default=>'auto',
			  -size=>5);
  print ")\n";
  print "&nbsp;&nbsp;&nbsp;";

  print "<B><A HREF='help.feature-map.html#orientation'>Orientation</B></A>&nbsp;&nbsp;&nbsp;&nbsp";
  print $query->popup_menu(-name=>'orientation',
			   -Values=>['horizontal','vertical'],
			   -default=>'horizontal');

  print "<BR>\n";
  print "<B><A HREF='help.feature-map.html#limits'>Display limits</A></B>&nbsp;";
  print "&nbsp;From&nbsp;";
  print $query->textfield(-name=>'from',
			  -default=>"$from",
			  -size=>5);
  
  print "&nbsp;To&nbsp;";
  print $query->textfield(-name=>'to',
			  -default=>"$to",
			  -size=>5);
  print "&nbsp;origin&nbsp;";
  print $query->textfield(-name=>'origin',
			  -default=>"$origin",
			  -size=>5);
  



  print "<BR>\n";
  print "<A HREF='help.feature-map.html#dimensions'><B>Map dimensions</B></A>\n";
  print "&nbsp;Length&nbsp;";
  print $query->textfield(-name=>'mlen',
			  -default=>'500',
			  -size=>5);
  
  print "&nbsp;thickness&nbsp;";
  print $query->textfield(-name=>'mapthick',
			  -default=>'auto',
			  -size=>5);
  
  print "&nbsp;spacing&nbsp;";
  print $query->textfield(-name=>'mspacing',
			  -default=>'6',
			  -size=>5);
  
  print "<BR>\n";
  print "<B><A HREF='help.feature-map.html#handle'>Feature handle</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
  print $query->popup_menu(-name=>'handle',
			   -Values=>['color dot','symbol','none'],
			   -default=>'color dot');

  print "<B><A HREF='help.feature-map.html#palette'>Color palette</A></B>&nbsp;&nbsp;&nbsp;&nbsp";
  print $query->popup_menu(-name=>'palette',
			   -Values=>['color','monochrome'],
			   -default=>'color');


  print "<BR>\n";
  print "<B>Feature thickness</B> \n";
  print "max ";
  print $query->textfield(-name=>'maxfthick',
			  -default=>'auto',
			  -size=>5);
  print "&nbsp;&nbsp;&nbsp;\n";

  print "min ";
  print $query->textfield(-name=>'minfthick',
			  -default=>'auto',
			  -size=>5);
  print "&nbsp;&nbsp;&nbsp;\n";

  print "<A HREF='help.feature-map.html#scorethick'>";
  print $query->checkbox(-name=>'scorethick',
			 -checked=>'checked',
			 -label=>' Proportional to score ');
  print "</A>";
  print "</B>";

  

  print "<BR>\n";
  print "<B>";
  print "<A HREF='help.feature-map.html#htmap'>";
  print $query->checkbox(-name=>'htmap',
			 -checked=>'checked',
			 -label=>' Dynamic ');
  print "</A>";
  print "&nbsp;&nbsp;&nbsp;";

  print "<BR>\n";

  print "<B><A HREF='help.feature-map.html#dynamic'>";
  print "Label keys";
  print "</A></B>\n";
  print $query->checkbox(-name=>'label_strand',
			 -label=>' strand ');
  print "&nbsp;&nbsp;&nbsp;";

  print $query->checkbox(-name=>'label_pos',
			 -label=>' position ');
  print "&nbsp;&nbsp;&nbsp;";

  print $query->checkbox(-name=>'label_id',
			 -label=>' identifier ');
  print "&nbsp;&nbsp;&nbsp;";

  print $query->checkbox(-name=>'label_descr',
			 -label=>' description ');
  print "&nbsp;&nbsp;&nbsp;";

  print $query->checkbox(-name=>'label_score',
			 -label=>' score ');
  print "&nbsp;&nbsp;&nbsp;";


  print "<P>";
  
  
  print "<UL><UL><TABLE><TR>";
  print "<TD>", $query->submit(-label=>"GO"), "</TD>\n";
  print"<TD>", $query->reset, "</TD>\n";
  print "<TD><B><A HREF='help.feature-map.html'>MANUAL</A></B></TD>\n";
  print "<TD><B><A HREF='demo.feature-map.html'>DEMO</A></B></TD>\n";
  print "<TD><B><A HREF='mailto:jvanheld\@cifn.unam.mx'>MAIL</A></B></TD>\n";
  print "</TR></TABLE></UL></UL>";

  print "</FONT>";

  print $query->end_form;
  print "<HR>\n";
  print $query->end_html;
}







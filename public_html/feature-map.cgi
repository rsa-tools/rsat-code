#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
#    $ERR_LOG = "$TMP/RSA_ERROR_LOG.txt";
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
#$tmp_file_name = sprintf "feature-map.%s", &AlphaDate;
#$tmp_color_file = sprintf "color-file.%s", &AlphaDate();

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
push @result_files, "Map file ($image_format)", $graph_file;
$htmap_file = $tmp_file_path.".html";
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

### display the result ###
my $graph_URL = $ENV{rsat_www}."/tmp/"; $graph_URL .= &RSAT::util::RelativePath($TMP, $graph_file);
my $html_URL = $ENV{rsat_www}."/tmp/"; $html_URL .= &RSAT::util::RelativePath($TMP, $htmap_file);
my $short_graph_file = &ShortFileName($graph_file);
my $short_feature_file = &ShortFileName($feature_file);

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
  print $htmap_content;
} else {
  print "<center><a href='".$graph_URL."'><img src='".$graph_URL."'></a></center><P>\n";
}

&PrintURLTable(@result_files);


print "<hr>";

# die join "\n",
#   "graph_file = ".$graph_file,
#   "graph_URL = ".$graph_URL,
#   "short_graph_file = ".$short_graph_file;

exit(0);




#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
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


## Draw a heat map for the transition table
my $draw_heatmap = 0;

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("convert-background-model result", 'results');


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = "$SCRIPTS/convert-background-model";
$prefix = "convert-bg";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
@result_files = ();

#### read parameters ####
my $parameters;


################################################################
## Background model method
my $bg_method = $query->param('bg_choose');
if ($bg_method eq "rsat") {
  ## Select pre-computed background file in RSAT genome directory
  #my $bg_taxo = $query->param('bg_taxo');
  #if ($bg_taxo eq "organism"){
  my $organism_name = $query->param("organism");
  $parameters .= " -org ".$organism_name;
  #}

  my $background_model = $query->param("background");
  $parameters .= " -bg ".$background_model;

  my  $markov_order = $query->param('markov_order');
  &RSAT::error::FatalError("Markov model should be a Natural number") unless &IsNatural($markov_order);
  $parameters .= " -markov ".$markov_order;

  ## draw the heatmap
  if ($markov_order <= 3) {
    $draw_heatmap = 1;
  }


  my $strands = $query->param('strands');
  if ($strands=~/single/) {
    $parameters .= " -1str ";
  }

  if ($strands=~/both/) {
    $parameters .= " -2str ";
  }
    
  if ($query->param('noov')) {
    $parameters .= " -noov ";
  }
    
  #    $bg_file = &ExpectedFreqFile($organism_name,
  #				 $oligo_length, $background_model,
  #				 noov=>$noov, str=>$str);
  #    $parameters .= " -bgfile ".$bg_file;
  $parameters .= " -from oligo-analysis";
    
} elsif ($bg_method =~ /upload/i) {
  ## Upload user-specified background file
  my $bgfile = $tmp_file_path."_bgfile.txt";
  push (@result_files, "Input background file", $bgfile);
  my $upload_bgfile = $query->param('upload_bgfile');
  if ($upload_bgfile) {
    if ($upload_bgfile =~ /\.gz$/) {
      $bgfile .= ".gz";
    }
    my $type = $query->uploadInfo($upload_bgfile)->{'Content-Type'};
    open BGFILE, ">$bgfile" ||
      &cgiError("Cannot store background file in temp dir.");
    while (<$upload_bgfile>) {
      print BGFILE;
    }
    close BGFILE;
    $parameters .= " -i $bgfile";
    $parameters .= " -from ".$query->param('bg_format');
  } else {
    &FatalError ("If you want to upload a background model file, you should specify the location of this file on your hard drive with the Browse button");
  }
    
} else {
  &RSAT::error::FatalError($bg_method," is not a valid method for background specification");
}



################################################################
## bg_pseudo
if (&IsReal($query->param('bg_pseudo'))) {
    $parameters .= " -bg_pseudo ".$query->param('bg_pseudo');
}


################################################################
## decimals
if (&IsInteger($query->param('decimals'))) {
    $parameters .= " -decimals ".$query->param('decimals');
} else {
    &FatalError("Decimals should be an integer number");
}


################################################################
## Matrix output format
my $output_format = lc($query->param('output_format'));
$parameters .= " -to ".$output_format;
$draw_heatmap = 0 unless ($output_format eq 'transitions');

## Output file
$result_file = $tmp_file_path.".".$output_format;
push @result_files, "Background model file ($output_format)", $result_file;

&ReportWebCommand($command." ".$parameters);

### execute the command ###
if ($query->param('output') eq "display") {
#    &PipingWarning();

 ## prepare figures
    ### prepare data for piping
    open RESULT, "$command $parameters |";

  ### prepare data for piping
  open RESULT, "$command $parameters |";

  ### open the sequence file on the server
  if (open MIRROR, ">$result_file") {
    $mirror = 1;
    &DelayedRemoval($result_file);
  }

    print "<PRE>";
    while (<RESULT>) {
      print "$_" unless ($query->param('output') =~ /server/i);
      print MIRROR $_ if ($mirror);
      if ($_ =~ /Error<\/h4><blockquote/ ) {
	$error_found = 1;
      }	
    }
    print "</PRE>";
    close RESULT;
    close MIRROR if ($mirror);

    if ($draw_heatmap) {
      my $heatmap_file = $tmp_file_path."_heatmap.png";
      my $command2 = "draw-heatmap -i $result_file";
      $command2 = "cut -f 1-5 ".$result_file; 
      $command2 .= "| ${SCRIPTS}/draw-heatmap -min 0 -max 1 ";
      $command2 .= " -out_format png";
      $command2 .= " -col_width 50 ";
      $command2 .= " -rownames ";
      $command2 .= " -o ".$heatmap_file;
#      $command2 .= " -html ".${TMP}."/".$heatmap_file.".html";
      push @result_files, ("Heatmap", $heatmap_file);

      print "<pre>command: $command2<p>\n</pre>" if ($ENV{rsat_echo} >= 1);
      `$command2`;
      $heatmap_URL = $ENV{rsat_www}."/tmp/".&RSAT::util::RelativePath(&RSAT::util::get_pub_temp(), $heatmap_file);
      print "<center><a href = \"".$heatmap_URL."\"><img src=\"".$heatmap_URL."\"></a>";
      &DelayedRemoval(&RSAT::util::get_pub_temp()."/${heatmap_file}.png");
      &DelayedRemoval(&RSAT::util::get_pub_temp()."/${heatmap_file}.html");
    }

    &PrintURLTable(@result_files);

    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $result_file);
}
print $query->end_html;

exit(0);


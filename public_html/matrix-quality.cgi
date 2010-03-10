#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
#BEGIN {
#    $ERR_LOG = "/dev/null";
##    $ERR_LOG = "$TMP/RSA_ERROR_LOG.txt";
#    use CGI::Carp qw(carpout);
#    open (LOG, ">> $ERR_LOG")
#	|| die "Unable to redirect log\n";
#    carpout(*LOG);
#}
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$command = "$SCRIPTS/matrix-quality";

$ENV{rsat_echo} = 1;

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("matrix-quality result", 'results');

#### update log file ####
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

#### read parameters ####
local $parameters = " -v 1";
################################################################
## File prefix
$tmp_file_name = join( "_", "matrix-quality", &AlphaDate());
$result_subdir = $tmp_file_name;
$result_dir = $TMP."/".$result_subdir;
$result_dir =~ s|\/\/|\/|g;
`mkdir -p $result_dir`;

$file_prefix = $result_dir."/";


################################################################
#### Matrix specification

$matrix_file = $file_prefix."matrix_input";
local $input_format = lc($query->param('matrix_format'));

if ($query->param('matrix')) {
    open MAT, "> $matrix_file";
    print MAT $query->param('matrix');
    close MAT;
    &DelayedRemoval($matrix_file);
      ($input_format) = split (/\s+/, $input_format);
    if (  ( $input_format eq "consensus" ) ||( $input_format eq "meme" ) ||( $input_format eq "infogibbs" ) ){
	$parameters .= " -ms $matrix_file";
    }
    else{
	$parameters .= " -m $matrix_file";
    }
} else {
    &RSAT::error::FatalError('You did not enter any data in the matrix box');
}


$parameters .=  " -matrix_format " . $input_format;
################################################################
## Pseudo-counts
if (&IsReal($query->param('pseudo_counts'))) {
    $parameters .= " -pseudo ".$query->param('pseudo_counts');
} else {
    &FatalError("Pseudo-count should be a real number");
}
if ($query->param('pseudo_distribution') eq "equi_pseudo") {
    $parameters .= " -equi_pseudo ";
}

################################################################
## sequence file
($sequence_file1,$sequence_format1) = &MultiGetSequenceFile(1,$file_prefix."sequence1.input", 1);

if ($query->param('tag1') ){
    $tag1 =$query->param('tag1') ;
}
$parameters .= " -seq ". $tag1 ." ".$sequence_file1 ;
$parameters .= " -seq_format ". $sequence_format1 ;

($sequence_file2) = &MultiGetSequenceFile(2,$file_prefix."sequence2.input", 0);
if ($query->param('tag2') ){
    $tag2 =$query->param('tag2') ;
}
$parameters .= " -seq ". $tag2 ." ".$sequence_file2 ;




################################################################
## permutations
if (&IsInteger($query->param('permutation1'))) {
    $parameters .= " -perm ".$tag1." ".$query->param('permutation1');
}

if (&IsInteger($query->param('permutation2'))) {
    $parameters .= " -perm ".$tag2." ".$query->param('permutation2');
}

################################################################
## Background model method
local $bg_method = $query->param('bg_method');
if ($bg_method eq "from_matrix") {

} elsif ($bg_method eq "bgfile") {
  ## Select pre-computed background file in RSAT genome directory
  local $organism_name = $query->param("organism");
  local $noov = "ovlp";
  local $background_model = $query->param("background");
  local $oligo_length = 1;
  $bg_file = &ExpectedFreqFile($organism_name,
			       $oligo_length, $background_model,
			       noov=>$noov, str=>"-1str");
  $parameters .= " -bgfile ".$bg_file;

} elsif ($bg_method =~ /upload/i) {
  ## Upload user-specified background file
  local $bgfile = "${TMP}/${tmp_file_name}_bgfile.txt";
  local $upload_bgfile = $query->param('upload_bgfile');
  if ($upload_bgfile) {
    if ($upload_bgfile =~ /\.gz$/) {
      $bgfile .= ".gz";
    }
    local $type = $query->uploadInfo($upload_bgfile)->{'Content-Type'};
    open BGFILE, ">$bgfile" ||
      &cgiError("Cannot store background file in temp dir.");
    while (<$upload_bgfile>) {
      print BGFILE;
    }
    close BGFILE;
    $parameters .= " -bgfile $bgfile";
    $parameters .= " -bg_format ".$query->param('bg_format');
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
###############
#output folder

$parameters .= " -o ".$file_prefix ."matrix_quality ";

###########################
#Command

print "<PRE>command: $command $parameters<P>\n</PRE>" if ($ENV{rsat_echo} >= 1);

$index_file = $result_subdir."/matrix_quality_index.html";
my $mail_title = join (" ", "[RSAT]", "matrix-quality",  &AlphaDate());
&EmailTheResult("$command $parameters", $query->param('user_email'), $index_file, title=>$mail_title);

print $query->end_html();


exit(0);

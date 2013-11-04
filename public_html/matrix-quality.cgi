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
require RSAT::util;
require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

#$ENV{rsat_echo} = 1;

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("matrix-quality result", 'results');

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = $SCRIPTS."/matrix-quality";
my $result_dir = &RSAT::util::make_temp_file("","matrix-quality", 1,1);
system("mkdir -p $result_dir ; chmod 755 $result_dir ");

my $file_prefix = "matrix-quality_".&AlphaDate();
my $tmp_file_name = $result_dir."/".$file_prefix;

## We remove the file created by mktemp and create a directory instead
#`rm -f $result_dir; mkdir -p $result_dir; chmod 777 $result_dir`;
&RSAT::message::Info("<br>Temporary file: ", $tmp_file_name,
		     "<br>Result dir: ", $result_dir,
#		     "<br>Result subdir: ", $result_subdir, 
		     "<br>File prefix: ", $file_prefix) if ($ENV{rsat_echo} >= 1);

#### read parameters ####
local $parameters = " -v 0";
################################################################
## File prefix


#####################
#Title specification
if ($query->param('html_title')) {
    $html_title=$query->param('html_title');
    $parameters .= " -html_title ' ".$html_title." ' ";
}

################################################################
#### Matrix specification
$matrix_file = $result_dir."/input_matrix";
local $input_format = lc($query->param('matrix_format'));

if ($query->param('matrix')) {
    open MAT, "> ".$matrix_file;
    print MAT $query->param('matrix');
    close MAT;
    &DelayedRemoval($matrix_file);
    ($input_format) = split (/\s+/, $input_format);
    if (  ( $input_format eq "consensus" ) ||( $input_format eq "meme" ) ||( $input_format eq "infogibbs" ) ||( $input_format eq "transfac" ) ){
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
## k parameter for the k-fold validation
if (&IsInteger($query->param('kfold'))) {
    $parameters .= " -kfold ".$query->param('kfold');
}

################################################################
## First sequence file
my $sequence_format="fasta";
($sequence_file1, $sequence_format1) = &MultiGetSequenceFile(1,$result_dir."/sequence1.fasta", 1);
if ($query->param('tag1') ){
  $tag1 = $query->param('tag1') ;
  $tag1 =~ s|\s|_|g;
  $tag1 =~ s|/|_|g;
  $tag1 =~ s|:|_|g;
}
$parameters .= " -seq ". $tag1 ." ".$sequence_file1 ;
$parameters .= " -seq_format ". "fasta" ;
if (lc($query->param('nwd')) eq "on") {
    
  $parameters .= " -plot ". $tag1." nwd " ;
}


################################################################
## Secod sequence file
($sequence_file2) = &MultiGetSequenceFile(2,$result_dir."/sequence2.fasta", 0);
if ($sequence_file2) {
  if ($query->param('tag2') ){
    $tag2 =$query->param('tag2') ;
    $tag2 =~ s|\s|_|g;
    $tag2 =~ s|/|_|g;
    $tag2 =~ s|:|_|g;
  } else {
    $tag2 = "seq_file2";
  }
  $parameters .= " -seq ". $tag2 ." ".$sequence_file2 ;
  if (lc($query->param('nwd')) eq "on") {
      
      $parameters .= " -plot ". $tag2." nwd " ;
  }

}

################################################################
## Permutations
if (&IsInteger($query->param('permutation1'))) {
    $parameters .= " -perm ".$tag1." ".$query->param('permutation1');
}

if (&IsInteger($query->param('permutation2'))) {
    $parameters .= " -perm ".$tag2." ".$query->param('permutation2')  if  $sequence_file2 ;
}

################################################################
## scan options
if ($query->param('scanopt1') ){
    $scanopt2 =$query->param('scanopt1') ;
}
$parameters .= " -scanopt ".$tag1." ".$scanopt1 if (  $scanopt1 ); 

if ($query->param('scanopt2') ){
    $scanopt1 =$query->param('scanopt2') ;
}
$parameters .= " -scanopt ".$tag2." ".$scanopt2   if ( $sequence_file2 && $scanopt2 ) ;

################################################################
## Markov order
my $markov_order = $query->param('markov_order');
&RSAT::error::FatalError("Markov model should be a Natural number") unless &IsNatural($markov_order);


################################################################
## Background model method
local $bg_method = $query->param('bg_method');
if ($bg_method eq "from_matrix") {
  
} elsif ($bg_method eq "bgfile") {
  ## Select pre-computed background file in RSAT genome directory
  local $organism_name = $query->param("organism");
  local $noov = "ovlp";
  local $background_model = $query->param("background");
  #local $oligo_length = 1;
  local $oligo_length = $markov_order + 1;
  $bg_file = &ExpectedFreqFile($organism_name,
			       $oligo_length, $background_model,
			       noov=>$noov, str=>"-1str");
  $parameters .= " -bgfile ".$bg_file.".gz";
  $parameters .= " -bg_format ".'oligo-analysis';

} elsif ($bg_method =~ /upload/i) {
  ## Upload user-specified background file
  local $bgfile = $tmp_file_name."_bgfile.txt";
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

$parameters .= " -archive -o ".$result_dir."/".$file_prefix ." ";

###########################
#Command

&ReportWebCommand($command." ".$parameters);

## Convert the absolute path of the directory into a path relative to the tmp directory for the Web link
#$result_subdir = $tmp_file_name;
#$result_subdir =~ s/${TMP}//;
#$index_file = $result_subdir."/".$file_prefix."_synthesis.html";
$index_file = $tmp_file_name."_synthesis.html";
my $mail_title = join (" ", "[RSAT]", "matrix-quality",  &AlphaDate());
&EmailTheResult($command." ". $parameters, $query->param('user_email'), $index_file, title=>$mail_title);

print $query->end_html();


exit(0);

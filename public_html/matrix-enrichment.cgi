#!/usr/bin/env perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
#require "cgi-lib.pl";
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
#### redirect error log to a file
#BEGIN {
#    $ERR_LOG = "/dev/null";
##    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
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
&RSA_header("matrix-enrichment result", 'results');

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = $SCRIPTS."/matrix-enrichment -top_matrices 50 ";
my $result_dir = &RSAT::util::make_temp_file("","matrix-enrichment", 1,1);
system("mkdir -p $result_dir ; chmod 755 $result_dir ");

my $file_prefix = "matrix-enrichment_".&AlphaDate();
my $tmp_file_name = $result_dir."/".$file_prefix;

## We remove the file created by mktemp and create a directory instead
#`rm -f $result_dir; mkdir -p $result_dir; chmod 777 $result_dir`;
&RSAT::message::Info("<br>Temporary file: ", $tmp_file_name,
		     "<br>Result dir: ", $result_dir,
#		     "<br>Result subdir: ", $result_subdir, 
		     "<br>File prefix: ", $file_prefix) if ($ENV{rsat_echo} >= 1);


################################################################
## Set parameters
local $parameters = " -v 0";


################################################################
## Title specification
if ($query->param('html_title')) {
    $html_title=$query->param('html_title');
    $parameters .= " -title ' ".$html_title." ' ";
}

################################################################
#### Matrix specification
local $matrix_file = &GetMatrixFile($result_dir."/input_matrix");


local $input_format = lc($query->param('matrix_format'));

($input_format) = split (/\s+/, $input_format);


$parameters .= " -matrix $matrix_file";
$parameters .=  " -matrix_format " . $input_format;

################################################################
## Pseudo-counts
if (&IsReal($query->param('pseudo_counts'))) {
    $parameters .= " -pseudo ".$query->param('pseudo_counts');
} else {
    &FatalError("Pseudo-count should be a real number");
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
} 




################################################################
## Secod sequence file
($sequence_file3) = &MultiGetSequenceFile(3,$result_dir."/sequence3.fasta", 0);
if ($sequence_file3) {
    if ($query->param('tag3') ){
	$tag3 =$query->param('tag3') ;
	$tag3 =~ s|\s|_|g;
	$tag3 =~ s|/|_|g;
	$tag3 =~ s|:|_|g;
    } else {
	$tag3 = "seq_file3";
    }
    $parameters .= " -seq ". $tag3 ." ".$sequence_file3 ;
}


################################################################
## Secod sequence file
($sequence_file4) = &MultiGetSequenceFile(4,$result_dir."/sequence4.fasta", 0);
if ($sequence_file4) {
  if ($query->param('tag4') ){
    $tag4 =$query->param('tag4') ;
    $tag4 =~ s|\s|_|g;
    $tag4 =~ s|/|_|g;
    $tag4 =~ s|:|_|g;
  } else {
    $tag4 = "seq_file4";
  }
  $parameters .= " -seq ". $tag4 ." ".$sequence_file4 ;

}



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
    local $organism_name = $query->param("organism_bg");
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
	
    } 
} elsif ($bg_method =~ /url/i) {
    ## Retrieve user-specified URL for background file
    my $url = $query->param('bgmodel_url');
    
    &RSAT::message::Info("Fetching background from URL ".$url) if ($ENV{rsat_echo} >= 1);
      my $bgmodel = "";
    local $bgfile = $tmp_file_name."_bgfile.txt";
    
    if (open BGM, ">$bgfile") {
	$bg = get($url);
	if ($bg =~ /\S/) {
	    print BGM $bg;
	    close BGM;
	} else {
	    &RSAT::error::FatalError("No background model could be downloaded from the URL ".$url);
	}
	
    }
    
    close BGFILE;
	$parameters .= " -bgfile $bgfile";
	$parameters .= " -bg_format ".$query->param('bg_format');
    
} else {
    &RSAT::error::FatalError($bg_method," is not a valid method for background specification");
}

################################################################
## bg_pseudo This is not an option in the main program we have to code it and reactivate it
#if (&IsReal($query->param('bg_pseudo'))) {
#  $parameters .= " -bg_pseudo ".$query->param('bg_pseudo');
#}
###############
#output folder

## archive is not an option in the program yet we have to code it and reactivate this code
#$parameters .= " -archive -o ".$result_dir."/".$file_prefix ." ";
$parameters .= " -o ".$result_dir."/".$file_prefix ." ";

###########################
#Command

&ReportWebCommand($command." ".$parameters);

$index_file = $tmp_file_name."_report.html";
my $mail_title = join (" ", "[RSAT]", "matrix-enrichment",  &AlphaDate());
#&EmailTheResult($command." ". $parameters, $query->param('user_email'), $index_file, title=>$mail_title);
if ($query->param('output') =~ /display/i) {
  &EmailTheResult("$command $parameters", "nobody@nowhere", $index_file, title=>"$mail_title",no_email=>1);
} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), $index_file, title=>"$mail_title");
}

print $query->end_html();


exit(0);

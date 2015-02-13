#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push @INC, "$`lib/";
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA.disco.lib";
require "RSA2.cgi.lib";

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

#$ENV{rsat_echo}=2;

### Read the CGI query
$query = new CGI;

### Print the header
&RSA_header("variation-scan result", "results");

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

$command = "$SCRIPTS/variation-scan";

#### read parameters ####
$parameters = " -v 1 ";

$prefix = "variation-scan";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1,0); $tmp_file_name = &ShortFileName($tmp_file_path);
@result_files = ();

################################################################
#### Matrix specification
$matrix_file_aux = $tmp_file_path."_input_matrix";
local $input_format = lc($query->param('matrix_format'));

if ($query->param('matrix')) {
    open MAT, "> $matrix_file_aux";
    print MAT $query->param('matrix');
    close MAT;
    &DelayedRemoval($matrix_file_aux);
    ($input_format) = split (/\s+/, $input_format);
} else {
    &RSAT::error::FatalError('You did not enter any data in the matrix box');
}

## variation-scan only takes transfac format matrices as input
## if the matrix in other formating transforme it to trasnfac and store it
my $convert_mtx_cmd;
if (!($input_format eq 'transfac')){
    $matrix_file= $tmp_file_path."_input_matrix_transfac_format";
    $convert_mtx_cmd=" convert-matrix -v 3 -from $input_format ";
    $convert_mtx_cmd.=" -to transfac -i $matrix_file_aux ";
    $convert_mtx_cmd.=" -o $matrix_file ";
    $matrix_covert=1;

    # open MTX_RESULT, "$convert_mtx_cmd |";

    # if (open MTX_FILE, ">$matrix_file") {
    # 	&DelayedRemoval($matrix_file);
    # }
    # while (<MTX_RESULT>) {
    # 	print  MTX_FILE  $_ ;
    # }
    # close MTX_RESULT;
    # close MTX_FILE;   
    
}else{
    $matrix_file=$matrix_file_aux;
}

$parameters .= " -m $matrix_file";

## Get input

unless ($input_seq_file = $query->param('variants_seq_file')){
    
    $input_seq_file= $tmp_file_path."variation-scan_sequence_input";
    
    ### sequence file already on the server side
    ### create a new temporary sequence file
    
    if ($query->param('uploaded_file')) {
	$upload_file = $query->param('uploaded_file');
	if ($upload_file =~ /\.gz$/) {
	    $input_seq_file .= ".gz";
	}
	$type = $query->uploadInfo($upload_file)->{'Content-Type'};
	open INPUT_FILE, ">". $input_seq_file ||
	    &cgiError("Cannot store input file in temporary directory");
	while (<$upload_file>) {
	    print INPUT_FILE;
	}
	close INPUT_FILE;
    } else {
	my $input_var = $query->param('variants_seqs');
	$input_varc =~ s/\r/\n/g;
	my @input_var = split (/[\n\r]/, $input_var);
	if ($input_var =~ /\S/) {
	    open QUERY, ">".$input_seq_file;
	    foreach my $row (@input_var) {
		next unless $row =~ /\S/; ## Skip empty rows
		chomp($row); ## Suppress newline character
		$row =~ s/ +/\t/; ## replace white spaces by a tab for the multiple genomes option
		print QUERY $row, "\n";
	    }
	    close QUERY;
	} else {
	    &cgiError("You should enter at least one variant sequence in rsat format");
	}
    }
    &DelayedRemoval($input_seq_file);
    
}
$parameters .= " -i ".$input_seq_file;
push @result_files ,("input", $input_seq_file) ;

################################################################
## Background model

## Markov order
my $markov_order = $query->param('markov_order');
&RSAT::error::FatalError("Markov model should be a Natural number") unless &IsNatural($markov_order);

## Method for specifyung the background model
my $bg_method = $query->param('bg_method');

if ($bg_method eq "bginput") {
    $parameters .= " -bginput";
    $parameters .= " -markov ".$markov_order;

} elsif  ($bg_method eq "bgfile") {
  ## Select pre-computed background file in RSAT genome directory
  my $organism_name = $query->param("organism");
  @org_name_split=split(" ",$organism_name);
  $species=join("_", $org_name_split[0], $org_name_split[1]);
  $species=~s/_$//;
  $assembly=$org_name_split[2];
  $organism_name_aux=$species;
  my $noov = "ovlp";
  my $background_model = $query->param("background");
  my $oligo_length = $markov_order + 1;
  $bg_file = &ExpectedFreqFile($organism_name_aux,
			       $oligo_length, $background_model,
			       noov=>$noov, str=>"-1str");

  $parameters .= " -bg ".$bg_file.".gz";

} elsif ($bg_method =~ /upload/i) {
    ## Upload user-specified background file
  my $bgfile = $tmp_file_path."_bgfile.txt";
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
      
      $bg_format=$query->param('bg_format');
      ##NEED TO CONVERT BG MODELS IN OTHER FORMAT NOT SUPPORTED
      if (!($bg_format eq 'oligo-analysis')){
	  $bg_file_oligo= $tmp_file_path."input_bgfile_oligoformat";
	  
	  $convert_bg_cmd=" convert-matrix -from $bg_format ";
	  $convert_bg_cmd.=" -to oligo-analysis -i $bg_file ";
	  $convert_bg_cmd.="-o $bg_file_oligo ";
	  $bgfile=$bg_file_oligo ;
	  $bg_convert=1;
      }
      $parameters .= " -bg $bgfile";
      
  } else {
      &FatalError ("If you want to upload a background model file, you should specify the location of this file on your hard drive with the Browse button");
  }
  
} else {
    &RSAT::error::FatalError($bg_method," is not a valid method for background specification");
}


################
##  scanning thresholds
my @l_threshold_fields = qw(score w_diff pval_ratio);
my @u_threshold_fields = qw(pval);
foreach my $field (@l_threshold_fields) {
    my $lth = $query->param("lth_".$field);
    &RSAT::error::FatalError($lth." is not a valid value for the lower $field threshold. Should be a number. ") unless (&IsReal($lth));
    $parameters .= " -lth $field $lth ";
}
foreach my $field (@u_threshold_fields) {    
    my $uth = $query->param("uth_".$field);
    &RSAT::error::FatalError($uth." is not a valid value for the upper $field threshold. Should be a number. ") unless (&IsReal($uth));
    $parameters .= " -uth $field $uth ";
}


if ($matrix_covert){
    $command=$convert_mtx_cmd." ; ".$command ;
}
if ( $bg_conver){
    $command=$convert_bg_cmd. " ; ". $command ;

}

## Report the command
&ReportWebCommand($command." ".$parameters);
$var_scan_file = "$tmp_file_path.variants-seq_result";

#### execute the command #####
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {
    &PipingWarning();

    open RESULT, "$command $parameters |";
    &PrintHtmlTable(RESULT, $var_scan_file, true, 1000);
    ### print the result
    close RESULT;
    

    push @result_files ,("variation-scan results",$var_scan_file ) ;
    &PrintURLTable(@result_files);

    ### prepare data for piping
   # &PipingForm();

    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $var_scan_file);
}



print $query->end_html();

exit(0);


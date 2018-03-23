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

$command = "$C_SCRIPTS/variation-scan1.0";

#### read parameters ####
$parameters = " -v 1 ";

$prefix = "variation-scan";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1,0); $tmp_file_name = &ShortFileName($tmp_file_path);
@result_files = ();

################################################################
#### Matrix specification

#local $matrix_file= &GetMatrixFile($tmp_file_path."_input_matrix");


#local $input_format = "transfac";
local $input_format=lc($query->param('matrix_format'));
($input_format) = split (/\s+/, $input_format);


## Get motif set or input file to scan variants

## MatrixDB selected 
my ($mat_db_params, @selected_db) = &GetMatrixDBchoice_select2("mode"=>"radio", "more"=>1);
if (scalar(@selected_db) > 0) {
  $mat_db_params=~s/-file2 //;
  $mat_db_params=~s/-format2.+//;  
  $parameters .= " -m ".$mat_db_params;
  $parameters .= " -m_format  transfac ";
  
} else {
  ## This should work but for some reason it does not AMR
  #elsif($query->param('custom_motif_db')) {
  my $persomotif_file = $tmp_file_path."variation-scan_sequence_custom_motif_manualinput.tf";
  local $matrix_file =&GetMatrixFile($persomotif_file);
  $parameters .= " -m $matrix_file";
  $parameters .= " -m_format  $input_format ";
}


## Code to get the motif from URL, no longer used
# elsif ($query->param('custom_motif_db_url_txt')=~ /\S/) {
#     my $url = $query->param('custom_motif_db_url_txt');
#     &RSAT::message::Info("Fetching matrices from URL ".$url) if ($ENV{rsat_echo} >= 1);
#     my $motif = "";
#     my $persomotif_file = $tmp_file_path."variation-scan_sequence_custom_motif_db_fromURL.tf";
#     if (open SEQ, ">$persomotif_file") {
# 	$motif = get($url);
# 	if ($motif =~ /\S/) {
# 	    print SEQ $motif;
# 	    close SEQ;
# 	} else {
# 	    &RSAT::error::FatalError("No sequence could be downloaded from the URL ".$url);
# 	}
#     }
    
#     ## Check sequence file
#     my $file_type = `file $persomotif_file`;
#     if ($file_type =~ "gzip") {
# 	&RSAT::message::TimeWarn("Uncompressing motifs file", $persomotif_file);
# 	my $cmd = "mv ".$persomotif_file." ".$persomotif_file.".gz";
# 	$cmd .= " ; gunzip ".$persomotif_file.".gz";
# 	&doit($cmd);
#     }  

#     $parameters .= " -m ".$persomotif_file;    
# }


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
  my $organism_name = $query->param("organism_bg");
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
	&RSAT::error::FatalError ("If you want to upload a background model file, you should specify the location of this file on your hard drive with the Browse button");
    }
    
}elsif ($bg_method =~/url/i) {
    ## Retrieve user-specified URL for background file
    my $url = $query->param('bgmodel_url');
    
    &RSAT::message::Info("Fetching background from URL ".$url) if ($ENV{rsat_echo} >= 1);
      my $bgmodel = "";
    $bgfile = $tmp_file_path."_bgfile.txt";
    
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
	$parameters .= " -bg $bgfile";
	#$parameters .= " -bg_format ".$query->param('bg_format');
    
}else {
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


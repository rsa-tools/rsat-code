#!/usr/bin/env perl
# matrix-clustering.cgi was used as template

############################################
## redirect error log to a file
if ($0 =~ /([^(\/)]+)$/) {
  push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
### redirect error log to a file
BEGIN {
  $ERR_LOG = "/dev/null";
#    $ERR_LOG = "/tmp/RSA_ERROR_LOG.txt";
  use CGI::Carp qw(carpout);
  open (LOG, ">> $ERR_LOG")
      || die "Unable to redirect log\n";
  carpout(*LOG);
}
require "RSA.lib";
require "RSA2.cgi.lib";
require RSAT::util;
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

################################################################
## Result page header

## Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("network-interactions result", "results");
&ListParameters() if ($ENV{rsat_echo} >= 2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();


################################################################
## Output paths
$command = $ENV{RSAT}."/perl-scripts/network-interactions";
#$return_fields = "-return json"; i dont need it

$output_prefix = "";
$output_path = &RSAT::util::make_temp_file("",$output_prefix, 1);

local $dir_name=$query->param('html_title');
$dir_name =~ s/\s/_/g; # correct dir name
$dir_name =~ s/\'//g; # correct dir name
$tmp_dir = &RSAT::util::get_pub_temp();

$output_dir = $output_path."/".$dir_name;


# initialize list of files to display
@result_files = ();
################################################################
## Command line parameters
local $parameters .= " -v 20";

###############
## Add title
local $title = lc($query->param('html_title'));
if($title){
    $title =~ s/\s+/_/g;
    $parameters .= " -title '".$title."'";
}

###############
## Add tf_list
my $tf_list_file = $output_path."_TFsList.tfs";
push @result_files, ("Input TF list",$tf_list_file);

# file option or...
if ($query->param('uploaded_tf_file')) {
  $upload_file = $query->param('uploaded_tf_file');
  if ($upload_file =~ /\.gz$/) {
    $tf_list_file .= ".gz";
  }
  $type = $query->uploadInfo($upload_file)->{'Content-Type'};
  open TF_FILE, ">".$tf_list_file ||
    &cgiError("Cannot store tf list file in temporary directory");
  while (<$upload_file>) {
    print TF_FILE;
  }
  close TF_FILE;

  # ... from text box
} else {
  my $tf_selection = $query->param('tf_selection');
  $tf_selection =~ s/\r/\n/g;
  my @tf_selection = split (/[\n\r]/, $tf_selection);
  if ($tf_selection =~ /\S/) {
    open QUERY, ">".$tf_list_file;
    foreach my $row (@tf_selection) {
        next unless $row =~ /\S/; ## Skip empty rows
        chomp($row); ## Suppress newline character
        $row =~ s/ +/\t/; ## replace white spaces by a tab for the multiple genomes option
        print QUERY $row, "\n";
    }
    close QUERY;
  } else {
    &cgiError("You should enter a transcription factor list.");
  }
}
$parameters .= " -tfs ".$tf_list_file;
&DelayedRemoval($tf_list_file);

###############
## Add cre_bed
my $cre_file = $output_path."_RegSeq.bed";
push @result_files, ("CRE BED file",$cre_file);

if ($query->param('uploaded_cre_file')) {
  $upload_file = $query->param('uploaded_cre_file');
  if ($upload_file =~ /\.gz$/) {
    $cre_file .= ".gz";
  }
  $type = $query->uploadInfo($upload_file)->{'Content-Type'};
  open CRE_FILE, ">".$cre_file ||
    &cgiError("Cannot store BED file in temporary directory");
  while (<$upload_file>) {
    print CRE_FILE;
  }
  close CRE_FILE;
} else {
  my $cre_selection = $query->param('cre_selection');
  $cre_selection =~ s/\r/\n/g;
  my @cre_selection = split (/[\n\r]/, $cre_selection);
  if ($cre_selection =~ /\S/) {
    open QUERY, ">".$cre_file;
    foreach my $row (@cre_selection) {
      next unless $row =~ /\S/; ## Skip empty rows
      chomp($row); ## Suppress newline character
      $row =~ s/ +/\t/; ## replace white spaces by a tab for the multiple genomes option
      print QUERY $row, "\n";
    }
    close QUERY;
  } else {
    &cgiError("You should enter a BED File either by upload or by copy-paste.");
  }
}
$parameters .= " -cre ".$cre_file;
&DelayedRemoval($cre_file);

###############
## Add genome
if ($query->param('genome_v')) {
    ($genome_v) = split " ", $query->param('genome_v'); ### take the first word
    $parameters .= " -genome ".$genome_v;
}

###############
## Add user network
if ($query->param('net_selection') || $query->param('uploaded_net_file')) {

  my $net_file = $output_path."_InputNet.tsv";
  push @result_files, ("Input network",$net_file);

  if ($query->param('uploaded_net_file')) {
    $upload_file = $query->param('uploaded_net_file');
    if ($upload_file =~ /\.gz$/) {
      $net_file .= ".gz";
    }
    $type = $query->uploadInfo($upload_file)->{'Content-Type'};
    open NET_FILE, ">".$net_file ||
      &cgiError("Cannot store network file in temporary directory");
    while (<$upload_file>) {
      print NET_FILE;
    }
    close NET_FILE;

    # ... from text box
  } else {
    my $net_selection = $query->param('net_selection');
    $net_selection =~ s/\r/\n/g;
    my @net_selection = split (/[\n\r]/, $net_selection);
    if ($net_selection =~ /\S/) {
      open QUERY, ">".$net_file;
      foreach my $row (@net_selection) {
          next unless $row =~ /\S/; ## Skip empty rows
          chomp($row); ## Suppress newline character
          $row =~ s/ +/\t/; ## replace white spaces by a tab for the multiple genomes option
          print QUERY $row, "\n";
      }
      close QUERY;
    } else {
      &cgiError("You should at least one interaction.");
    }
  }

  $parameters .= " -report_net -net ".$net_file;
  &DelayedRemoval($net_file);
}

## use user-specified matrixes
## check if user matrix has been provided
if ($query->param('uploaded_file_matrix1') || ($query->param('matrix1') =~ /\S/)) {

  ## Matrix input format
  local $query_matrix_1_format = lc($query->param('matrix_format1'));

  #### Query matrix file
  local $matrix_file_1 = &GetMatrixFile($output_dir."/IntermediateFiles/user_matrices1.".$query_matrix_1_format,1);

  $parameters .= " -m_format ".$query_matrix_1_format." -m ".$matrix_file_1;
}

################################################################
## Advanced options

## matrix-scan

## score lth
local $matscan_score = lc($query->param('matrixscan_score'));
if($matscan_score){
    $matscan_score =~ s/\s+/_/g;
    $parameters .= " -score '".$matscan_score."'";
}

## pval uth
local $matscan_pval = lc($query->param('matrixscan_pval'));
if($matscan_pval){
    $matscan_pval =~ s/\s+/_/g;
    $parameters .= " -pval '".$matscan_pval."'";
}

## Add motif database


################################################################
## Output dir
$parameters .= " -outdir ".$output_dir;

$parameters .= " -html "; # output html SUMMRARY file

$date = &AlphaDate();
$parameters .= " -date ".$date; #date to use to name files

## Add an error-log file for matrix-clustering
#$err_file = $output_path."/".$output_prefix."_err.txt";
#$parameters .= " 2> ".$err_file; i dont know what '2>' does

## Report the full command before executing
&ReportWebCommand($command." ".$parameters, $err_file);


#################################################################
### Display or send result by email
$index_file = $output_dir."/ResultsFiles/";
$index_file .= $output_prefix."_" if ($output_prefix);
$index_file .= "SUMMARY_".$date.".html";

my $mail_title = join (" ", "[RSAT]", "network-interactions", $date);
if ($query->param('output') =~ /display/i) {
  &EmailTheResult("$command $parameters", "nobody@nowhere", "", title=>$mail_title, index=>$index_file, no_email=>1);
} else {
  &EmailTheResult("$command $parameters", $query->param('user_email'), "", title=>$mail_title,index=>$index_file);
}

################################################################
## Result page footer
print $query->end_html;

exit(0);

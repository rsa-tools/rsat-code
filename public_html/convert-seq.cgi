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

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("convert-seq result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

$command = "$SCRIPTS/convert-seq";
$prefix="convert-seq";
$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);
#$tmp_file_name = sprintf "convert-seq.%s", &AlphaDate();
@result_files = ();

#### read parameters ####
$parameters = "";

### sequence file
#($sequence_file,$out_format) = &GetSequenceFile();
#$parameters .= " -i $sequence_file -format $out_format";

##### input sequence file #####
#unless ($query->param('sequence') =~ /\S/) {
#    &cgiError("The sequence box should not be empty.");
#}
#open INSEQ, ">$tmp_file_path";
#print INSEQ $query->param('sequence');
#close INSEQ;
#$parameters .= " -i $tmp_file_path ";
#&DelayedRemoval("$tmp_file_path");

################################################################
## sequence file
($in_sequence_file,$sequence_format) = &GetSequenceFile();
push @result_files, ("Input sequence ($sequence_format)",$in_sequence_file);

#&cgiWarning("in sequence = ".$in_sequence_file);

#### matrix-scan parameters
$parameters .= " -i $in_sequence_file -from $sequence_format";
&DelayedRemoval("$in_sequence_file");

## Short sequence treatment
my $short_action = lc($query->param('short_action'));
if (($short_action eq "mask") || ($short_action eq "skip")) {
  $parameters .= " -".$short_action."_short ";
  if (&IsNatural($query->param('short_size'))) {
    $parameters .= $query->param('short_size');
  } else {
    &FatalError("Minimal sequence size must be a strictly positive Integer number");
  }
}

## add reverse-complement
if ($query->param('addrc')) {
    $parameters .= " -addrc ";
}

## output format
$out_format = lc($query->param('output_format'));
$parameters .= " -to $out_format ";

## line width
if ($query->param('line_width') =~ /\d+/) {
    $parameters .= " -lw ".$query->param('line_width');
}

### open the sequence file on the server
$sequence_file = $tmp_file_path.".".$out_format;
push @result_files, ("Converted sequence ($out_format)",$sequence_file);

&ReportWebCommand($command." ".$parameters);

#### execute the command #####
if (($query->param('output') =~ /display/i) ||
    ($query->param('output') =~ /server/i)) {

#    $ENV{RSA_OUTPUT_CONTEXT} = "text";
    open RESULT, "$command $parameters |";

    ### print the result ### 
    &PipingWarning();
    if ($query->param('output') =~ /server/i) {
	&Info("The result will appear below ...");
    }

    print '<H4>Result</H4>';

    if (open MIRROR, ">$sequence_file") {
	$mirror = 1;
	&DelayedRemoval($sequence_file);
    }

    print "<PRE>";
    while (<RESULT>) {
	print "$_" unless ($query->param('output') =~ /server/i);
	print MIRROR $_ if ($mirror);
    }
    print "</PRE>";
    close RESULT;
    close MIRROR if ($mirror);

    ### prepare data for piping
    &PrintURLTable(@result_files);
    &PipingFormForSequence();
    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'), $sequence_file);
}

exit(0);



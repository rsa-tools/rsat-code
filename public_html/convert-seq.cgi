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
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$command = "$SCRIPTS/convert-seq";
$tmp_file_name = sprintf "convert-seq.%s", &AlphaDate;

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("convert-seq result");
#&ListParameters;

#### update log file ####
&UpdateLogFile;

#### read parameters ####
$parameters = "";

##### input format #####
$input_format = lc($query->param('input_format'));
if ($accepted_input_seq{$input_format}) {
    $parameters .= " -from $input_format ";
} else {
    &cgiError("Invalid format for input sequence: $input_format.");
}

##### output format #####
$output_format = lc($query->param('output_format'));
if ($accepted_output_seq{$output_format}) {
    $parameters .= " -to $output_format ";
} else {
    &cgiError("Invalid format for output sequence: $output_format.");
}

### sequence file
#($sequence_file,$sequence_format) = &GetSequenceFile();
#$parameters .= " -i $sequence_file -format $sequence_format";

##### input sequence file #####
unless ($query->param('sequence') =~ /\S/) {
    &cgiError("The sequence box should not be empty.");
}
open INSEQ, ">$TMP/$tmp_file_name";
print INSEQ $query->param('sequence');
close INSEQ;
$parameters .= " -i $TMP/$tmp_file_name ";
&DelayedRemoval("$TMP/$tmp_file_name");

##### add reverse-complement #####
if ($query->param('addrc') eq "yes") {
    $parameters .= " -addrc ";
}

##### line width #####
if ($query->param('line_width') =~ /\d+/) {
    $parameters .= " -lw ".$query->param('line_width');
}

### execute the command ###
if ($query->param('output') eq "display") {
    ### Print the result on Web page
    open RESULT, "$command $parameters  & |";
    
    print "<PRE>$command $parameters</PRE> " if ($ECHO >= 1);
    
    print "<PRE>";
    while (<RESULT>) {
	print "$_";
    }
    print "</PRE>";
    close RESULT;
    
    print "<HR SIZE = 3>";
} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult( "$command $parameters", , $query->param('user_email'));
}
print $query->end_html;

exit(0);



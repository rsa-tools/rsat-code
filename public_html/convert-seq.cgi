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
$command = "$SCRIPTS/convert-seq";
$tmp_file_name = sprintf "convert-seq.%s", &AlphaDate();

### Read the CGI query
$query = new CGI;

### print the result page
&RSA_header("convert-seq result", "results");
&ListParameters() if ($ENV{rsat_echo} >=2);

#### update log file ####
&UpdateLogFile();

#### read parameters ####
$parameters = "";

### sequence file
#($sequence_file,$out_format) = &GetSequenceFile();
#$parameters .= " -i $sequence_file -format $out_format";

##### input sequence file #####
#unless ($query->param('sequence') =~ /\S/) {
#    &cgiError("The sequence box should not be empty.");
#}
#open INSEQ, ">$TMP/$tmp_file_name";
#print INSEQ $query->param('sequence');
#close INSEQ;
#$parameters .= " -i $TMP/$tmp_file_name ";
#&DelayedRemoval("$TMP/$tmp_file_name");

################################################################
## sequence file
($in_sequence_file,$sequence_format) = &GetSequenceFile();

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

print "<PRE>command: $command $parameters<P>\n</PRE>" if ($ENV{rsat_echo} >= 1);

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

    ### open the sequence file on the server
    $sequence_file = "$TMP/$tmp_file_name.res";
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

    if ($query->param('output') =~ /server/i) {
	$result_URL = "$ENV{rsat_www}/tmp/${tmp_file_name}.res";
	print ("The result is available at the following URL: ", "\n<br>",
	       "<a href=${result_URL}>${result_URL}</a>",
	       "<p>\n");
    }

    ### prepare data for piping
    &PipingFormForSequence();
    print "<HR SIZE = 3>";

#} elsif 
#    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}

# ### execute the command ###
# if ($query->param('output') eq "display") {
#     ### Print the result on Web page
#     open RESULT, "$command $parameters  & |";
    
#     print "<PRE>$command $parameters</PRE> " if ($ENV{rsat_echo} >= 1);
    
#     &PipingWarning();
#     print "<PRE>";
#     while (<RESULT>) {
# 	print "$_";
#     }
#     print "</PRE>";
#     close RESULT;
#     print "<HR SIZE = 3>";
# } elsif ($query->param('output') =~ /server/i) {
#     &ServerOutput("$command $parameters", $query->param('user_email'));
# } else {
#     &EmailTheResult( "$command $parameters", , $query->param('user_email'));
# }
# print $query->end_html;

exit(0);



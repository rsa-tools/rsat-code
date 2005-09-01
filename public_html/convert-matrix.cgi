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
#    $ERR_LOG = "$TMP/RSA_ERROR_LOG.txt";
    use CGI::Carp qw(carpout);
    open (LOG, ">> $ERR_LOG")
	|| die "Unable to redirect log\n";
    carpout(*LOG);
}
require "RSA.lib";
require "RSA.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$command = "$SCRIPTS/convert-matrix";
$tmp_file_name = sprintf "convert-matrix.%s", &AlphaDate();
$result_file = "$TMP/$tmp_file_name.res";

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("convert-matrix result");

#### update log file ####
&UpdateLogFile();

&ListParameters() if ($ECHO >= 2);

#### read parameters ####
$parameters = " -v 1 ";


################################################################
#### Matrix specification
$matrix_file = "$TMP/$tmp_file_name.input";
open MAT, "> $matrix_file";
print MAT $query->param('matrix');
close MAT;
&DelayedRemoval($matrix_file);

$parameters .= " -i $matrix_file";

################################################################
#### Matrix format
my $matrix_format = lc($query->param('matrix_format'));
$matrix_format =~ s/(\S+)/$1/; ## Only retain the first word
$parameters .= " -in_format $matrix_format";

## Return fields
my @return_fields = ();
foreach my $stat qw (counts frequencies weights information consensus parameters profile margins) {
    if ($query->param($stat)) {
	push @return_fields, $stat;
    }
}
$parameters .= " -return ";
$parameters .= join ",", @return_fields;

print "<PRE>command: $command $parameters<P>\n</PRE>" if ($ECHO >= 1);

### execute the command ###
if ($query->param('output') eq "display") {
#    &PipingWarning();

    ### prepare data for piping
    open RESULT, "$command $parameters |";

    print '<H2>Result</H2>';
    print '<PRE>';
    while (<RESULT>) {
	s|${TMP}/||g;
	s|${BIN}/||g;
	print $_;
	$genes .= $_;
    }
    print '</PRE>';
    close(RESULT);

#    &PipingForm();

    print "<HR SIZE = 3>";

} elsif ($query->param('output') =~ /server/i) {
    &ServerOutput("$command $parameters", $query->param('user_email'));
} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}
print $query->end_html;

exit(0);


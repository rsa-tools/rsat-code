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
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";
$command = "$SCRIPTS/convert-features -v 1";
$tmp_file_name = sprintf "convert-features.%s", &AlphaDate();
$result_file = "$TMP/$tmp_file_name.res";
$ENV{rsat_echo} = 1;

### Read the CGI query
$query = new CGI;

### print the header
&RSA_header("convert-features result", 'results');

#### update log file ####
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

#### read parameters ####
my $parameters;

################################################################
#### Feature from input box
my $feature_file = "$TMP/$tmp_file_name.input";
if ($query->param('feature')){
	open FEAT, "> $feature_file";
	print FEAT $query->param('feature');
	close FEAT;
} else  {
    ## Upload user-specified  file
    my $upload_file = $query->param('uploaded_file');
    if ($upload_file) {
	if ($upload_file =~ /\.gz$/) {
	    $feature_file .= ".gz";
	}
	my $type = $query->uploadInfo($upload_file)->{'Content-Type'};
	open FEAT, ">$feature_file" ||
	    &cgiError("Cannot store feature file in temp dir.");
	while (<$upload_bgfile>) {
	    print FEAT;
	}
	close FEAT;
    } else {
	&FatalError ("If you want to upload a file, you should specify the location of this file on your hard drive with the Browse button");
    }

}

&DelayedRemoval($feature_file);
$parameters .= " -i $feature_file";


################################################################
## feature input format

my $input_format = lc($query->param('feature_format'));
$parameters .= " -from ".$input_format;


################################################################
## feature output format
my $output_format = lc($query->param('output_format'));
$parameters .= " -to ".$output_format;


print "<PRE>command: $command $parameters<P>\n</PRE>" if ($ENV{rsat_echo} >= 1);

### execute the command ###
if ($query->param('output') eq "display") {
#    &PipingWarning();

 ## prepare figures
    ### prepare data for piping
    open RESULT, "$command $parameters |";
    
    print '<H4>Result</H4>';
    print '<PRE>';
    while (<RESULT>) {
		print $_;
    }
    print '</PRE>';
    close(RESULT);

#    &PipingForm();

    print "<HR SIZE = 3>";

} else {
    &EmailTheResult("$command $parameters", $query->param('user_email'));
}
print $query->end_html;

exit(0);


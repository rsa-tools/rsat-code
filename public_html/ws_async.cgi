#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;

$| = 1;

require "RSA.lib";
require "RSA2.cgi.lib";

use SOAP::WSDL; ## Requires version 2.0 or later of SOAP::WSDL
use lib '../ws_clients/perl_clients/RSATWS';
use MyInterfaces::RSATWebServices::RSATWSPortType;

use File::Basename;

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$query = new CGI;

my $status = &JobStatus($query->param('ticket'));

$sleepTime;
if ($query->param('sleep_time')) {
	$sleepTime = $query->param('sleep_time');
	if ($sleepTime > 100){
		$sleepTime = 100;
	}
} else {
	$sleepTime = 10;
}

if ($status eq 'Running') {

	### Print the waiting page
	my $css_body_class = "form";
	my $title = $query->param('title');
	print $query->header();
	print $query->start_html(-title=>"RSA-tools : $title",
			   -class => "$css_body_class",
			   -author=>'jvanheld@bigre.ulb.ac.be',
			   -style => { 	-src => "$ENV{rsat_www}/main.css",
                             	       	-type => 'text/css',
                             		-media => 'screen' });
	
	&PrintPage();

	sleep($sleepTime);

	print '<script type="text/javascript" language="JavaScript">
	//submit form
	document.form1.submit();
	</script>';
	print $query->end_html . "\n";       

} else {

	### Print the result page
	### print the header
	my $header = $query->param('title')." result";
	&RSA_header($header, 'results');

	&ListParameters() if ($ENV{rsat_echo} >= 2);

	&GetResult($query->param('ticket'));
	print $query->end_html;

	exit(0);
}

sub JobStatus {
    my $ticket = shift;
    my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();

    my %args = (
	'ticket' => $ticket
	);

    my $som = $soap->monitor({'request' => \%args});

    unless ($som) {
	printf "A fault (%s) occured: %s\n", $som->get_faultcode(), $som->get_faultstring();
    } else {
	my $response = $som->get_response();

	## Report the status
	my $status = $response -> get_status();
	return $status;
    }
}

sub GetResult {
    my $soap=MyInterfaces::RSATWebServices::RSATWSPortType->new();
    my $ticket = shift;

    my %args = (
	'ticket' => $ticket
	);

    my $som = $soap->get_result({'request' => \%args});

    ### print the result ### 
    &PipingWarning();

    ## Get the result
    unless ($som) {
	print '<H4>Error</H4>';
	my $error_message = $som->get_faultstring();
	$error_message =~ s/command/\<BR\>Command/;
	print $error_message;
    } else {
	print '<H4>Result</H4>';
	my $results = $som->get_response();

	## Report the remote command
	my $command = $query->param('command'), "\n";
	$command =~ s|$ENV{RSAT}/perl-scripts/||;
	print "Command used on the server: ".$command, "\n";

	## Report the result
	$result = $results -> get_client();
	if ($ENV{rsat_ws_tmp} =~ /ulb\.ac\.be/) {
#	    $server_file = $results -> get_server();
	    $sequence_file = $results -> get_server();     ## variable name has to be '$sequence_file' for PipingFormForSequence subroutine
	} else {
	    $sequence_file = `mktemp $ENV{rsat_tmp}/retrieve-ensembl-seq.XXXXXXXXXX`;
	    open TMP_IN, ">".$sequence_file or die "cannot open temp file ".$sequence_file."\n";
	    print TMP_IN $result;
	    close TMP_IN;
	}
    }

    if ($query->param('output') =~ /server/i) {
	my $server_file_name = basename($sequence_file);
	$result_URL = "$ENV{rsat_ws_tmp}/$server_file_name";
	print ("<p>The result is available at the following URL: ", "\n",
	       "<a href=${result_URL}>${result_URL}</a>",
	       "<p>\n");
    } elsif ($query->param('output') =~ /display/i) {
	print "<PRE>";
	print $result;
	print "</PRE>";
    }

    ### prepare data for piping
    $out_format = 'fasta';    
    &PipingFormForSequence();

    print "<HR SIZE = 3>";

}

sub PrintPage {
    my @params = $query->param();

    print $query->h2('Processing your ' . $query->param('title') .' request'). "\n";
    print "<hr/>";
    print '<TABLE>' . "\n";
    print "<tr><th>Job ID</th><td>" . $query->param('ticket') . "</td></tr>\n";
    print "<tr><th>Status</th><td>Running</td><td><img src='images/loader.gif'/></tr>\n";
    print "<tr><th>Submitted at </th><td>".$query->param('submit_time')."</td></tr>\n";
    print "<tr><th>Current time </th><td>".&current_time()."</td></tr>\n";
    print "</TABLE>\n"; 
    print "<p>This page will be automatically updated in ".$sleepTime." seconds</p>";
    print "<hr/>";
    print "<h2>The RSAT servers are working for you, take a break with the latest strips from 
     <a href='http://www.phdcomics.com/' target=_blank> PHD Comics </a> !</h2>";
    print '<iframe width="800" height="400" src="http://www.rss-info.com/rss2.php?integration=if&windowopen=1&rss=http%3A%2F%2Fwww.phdcomics.com%2Fgradfeed_justcomics.php&number=10&width=800&ifbgcol=FFFFFF&bordercol=D0D0D0&textbgcol=F0F0F0&rssbgcol=F0F0F0&showrsstitle=1&showtext=1" frameborder=0></iframe>';

    print '<FORM name="form1" id="form1" method="POST" action="ws_async.cgi">';
    foreach my $parameter (sort @params) {
        print ' <INPUT type="hidden" name="'.$parameter.'" value="'.$query->param($parameter).'">';  
    }
    $sleepTime += 5;
    print ' <INPUT type="hidden" name="sleep_time" value="'.$sleepTime.'">';
    print '</FORM>';
}

sub current_time {
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	if ($second <10) {
		$second = "0".$second;
	}
	if ($minute <10) {
		$minute = "0".$minute;
	}
	return("$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year");
}

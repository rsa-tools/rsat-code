#!/usr/bin/perl
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;

#### redirect error log to a file
BEGIN {
    $ERR_LOG = "/dev/null";
    #    $ERR_LOG = &RSAT::util::get_pub_temp()."/RSA_ERROR_LOG.txt";
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

### print the header
&RSA_header("Supported orthologs", "results");


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);
my $group = "";
if ($ENV{group_specificity}) {
    $group = $ENV{group_specificity};
}

my @selected_organisms = ();
push @selected_organisms, &GetOrganismsForTaxon("Bacteria")
if (($group_specificity eq "Bacteria") ||
($group_specificity eq "Prokaryotes"));
push @selected_organisms, &GetOrganismsForTaxon("Archaea")
if (($group_specificity eq "Archaea") ||
($group_specificity eq "Prokaryotes"));
@selected_organisms = sort(@selected_organisms);


################################################################
## Print general information about this RSAT instance
print "<h2>RSAT instance: ", $ENV{rsat_site}, "</h2>\n";

print "<p><b>Organisms supported: </b>", scalar @selected_organisms, "</p>\n";

if ($group) {
    print "<p><b>Group specificity: </b>", $group, "</p>\n";
}

foreach $_ (@selected_organisms){
    print $_ . "<br/>";
}

print '<hr size=3>';
print "</div>";
print "</div>";

print $query->end_html;

exit(0);


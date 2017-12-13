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
&RSA_header("Supported taxonomy", "results");


## Check security issues
&CheckWebInput($query);

## update log file
&UpdateLogFile();

&ListParameters() if ($ENV{rsat_echo} >= 2);

@taxon = &get_taxons_web("all");
##$tmp_file_name = sprintf "supported-organisms.%s", &AlphaDate();
#$prefix = "supported-organisms";
#$tmp_file_path = &RSAT::util::make_temp_file("",$prefix, 1); ($tmp_file_dir, $tmp_file_name) = &SplitFileName($tmp_file_path);

$font{variable} = 1;

################################################################
## treat taxon specificity of the server if required
my $group = "";
if ($ENV{group_specificity}) {
    $group = $ENV{group_specificity};
}
################################################################
## Print general information about this RSAT instance
print "<h2>RSAT instance: ", $ENV{rsat_site}, "</h2>\n";

print "<p><b>Taxonomy supported: </b>", scalar @taxon, "</p>\n";

if ($group) {
    print "<p><b>Group specificity: </b>", $group, "</p>\n";
}


foreach $t (@taxon){
    print $t . "<br/>";
}
print '<hr size=3>';
print "</div>";
print "</div>";

print $query->end_html;

exit(0);

#!/usr/bin/env perl
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

my @supported_orthologs = ();
my $fp_org_file = $ENV{RSAT}."/public_html/data/supported_blast_tables.tab";
&RSAT::message::Debug("supported BLAST tables", $fp_org_file) if ($main::verbose >= 5);
if (-e $fp_org_file) {
    my ($fp_org) = &OpenInputFile($fp_org_file);
    while (<$fp_org>) {
        next unless (/\S/) ; # skip empty rows
        next if (/^;/); # skip comment lines
        next if (/^\#/); # Skip header line
        
        my @split_line = split(/\t/, $_);
        my $org = $split_line[0];
        my $taxa = $split_line[1];
        
        if ($org) {
            
            #	      if (defined($supported_organism{$org})) {## This control seems
            #	      not to work, probably because the library is read before the
            #	      list of supported organisms.
            $supported_orthologs{$org} = 1;
            $supported_genome_blast{$org}{$taxa} = 1;
            #	      }
        }
    }
    @supported_orthologs = sort keys %supported_orthologs;
} else {
    &RSAT::error::FatalError("Missing file", $fp_org_file);
}
################################################################
## Print general information about this RSAT instance
print "<h2>RSAT instance: ", $ENV{rsat_site}, "</h2>\n";

print "<p><b>Organisms supported: </b>", scalar @supported_orthologs, "</p>\n";

if ($group) {
    print "<p><b>Group specificity: </b>", $group, "</p>\n";
}

foreach $_ (@supported_orthologs){
    print $_ . "<br/>";
}

print '<hr size=3>';
print "</div>";
print "</div>";

print $query->end_html;

exit(0);


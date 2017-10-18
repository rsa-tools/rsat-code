#!/usr/bin/perl
############################################################
#
# $Id: getMatrix.cgi: get the matrix from the motif name
#
############################################################
#### this cgi script fills the HTML form for the program dna-pattern
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$query = new CGI;

print "Content-type:txt/html\n\n";
print "<html>";
print "<body>";

$db_choice = $query->param("db_choice");

my %matrix_db = &RSAT::server::supported_motif_databases();
foreach my $db (sort keys %matrix_db){
    my $db_name = $matrix_db{$db}->{name};
    my $format = $matrix_db{$db}->{format};
    my $file = $matrix_db{$db}->{file};

    if($db_choice eq $db_name){
        $file_name = $ENV{RSAT} . "/public_html/motif_databases/" . $file;
        print "$format</format>";
        open($fh, "<", $file_name) or die "Cannot open file $!";
        while(my $row = <$fh>){
           print $row;
        }
    }
}

print "</body></html>";
exit(0);

#!/usr/bin/perl
############################################################
#
# $Id: getMatrixIds.cgi: get the identifiants from the motif file
#
############################################################
#### this cgi script fills the HTML form for the program dna-pattern
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
use CGI;
use JSON;
use CGI::Carp qw/fatalsToBrowser/;
require "RSA.lib";
require "RSA2.cgi.lib";

$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$query = new CGI;

print "Content-type: application/json; charset=iso-8859-1\n\n";

$db_choice = $query->param("db_choice");

my %matrix_db = &RSAT::server::supported_motif_databases();
my @ids;

my %db = %{$matrix_db{$db_choice}};
my $format = $db{format};
my $file = $db{file};

$file_name = $ENV{RSAT} . "/public_html/motif_databases/" . $file;
open($fh, "<", $file_name) or die "Cannot open file $!";
my $id = "";
while(my $row = <$fh>){
    if($row =~ /^AC\s+/){
        my @f = split(/\s+/, $row);
        $id = $f[1];
    }
    if($row =~ /^ID\s+/){
        my @f = split(/\s+/, $row);
        my $idac = "";
        if($f[1] ne $id){
            $idac = $f[1] . " - " . $id;
        }else{
            $idac = $id;
        }
        push @ids, { "id" => $id, "idac" => $idac };
    }
    if($row =~ /^OS\s+/){
        push @ids, { "id" => $id, "idac" => $id };
    }
}

print JSON::encode_json( {entries => \@ids} );
exit(0);

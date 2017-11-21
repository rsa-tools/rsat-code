#!/usr/bin/perl
############################################################
#
# $Id: getOrganisms.cgi: get the supported organisms
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

my $term = $query->param("term");
my @selected_organisms = &RSAT::OrganismManager::get_supported_organisms_web();

my @organisms_popup = ();
foreach my $org (@selected_organisms){
    my $name = $org;
    $name =~ s/\_/ /g;
    if($name =~ /$term/i){
        push @organisms_popup, { "value" => $org, "label" => $name };
    }
}
if(scalar @organisms_popup == 0){
    push @organisms_popup, {"value" => "null", "label" => "null"};
}
print JSON::to_json( \@organisms_popup );


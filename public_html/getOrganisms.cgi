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

my $supported_orgs = $query->param("supported");
my $term = $query->param("term");
my @terms = split(" ", $term);

#### Get taxons
if($query->param("taxon") eq "yes"){
    my $taxon_nodetype = $query->param("nodetype");
    @taxa = &get_taxons_web($taxon_nodetype);
    
    my @taxonomy_popup = ();
   
    foreach my $taxon (@taxa) {
        my $name = $taxon;
        $name =~ s/_/ /g;
        my $reg = "";
        foreach $t (@terms){
            $reg .= $t . "(.*)";
        }
        if($name =~ /$reg/i){
            push @taxonomy_popup, { "value" => $taxon, "label" => $name };
        }
    }
    if(scalar @taxonomy_popup == 0){
        push @taxonomy_popup, {"value" => "null", "label" => "null"};
    }
    print JSON::to_json( \@taxonomy_popup );
}
#### Get organisms
else{
    my @selected_organisms = ();
    #### PrintOrthoSelectionSection
    if($supported_orgs eq "orthologs"){
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
            @selected_organisms = sort keys %supported_orthologs;
        } else {
            &RSAT::error::FatalError("Missing file", $fp_org_file);
        }
    }
    #### infer-operons
    elsif($supported_orgs eq "infer"){
        push @selected_organisms, &GetOrganismsForTaxon("Bacteria")
        if (($group_specificity eq "Bacteria") ||
        ($group_specificity eq "Prokaryotes"));
        push @selected_organisms, &GetOrganismsForTaxon("Archaea")
        if (($group_specificity eq "Archaea") ||
        ($group_specificity eq "Prokaryotes"));
        @selected_organisms = sort(@selected_organisms);
    }else{
        @selected_organisms = &RSAT::OrganismManager::get_supported_organisms_web();
    }
    my @organisms_popup = ();
    foreach my $org (@selected_organisms){
        my $name = $org;
        $name =~ s/\_/ /g;
        my $reg = "";
        foreach $t (@terms){
            $reg .= $t . "(.*)";
        }
        if($name =~ /$reg/i){
            push @organisms_popup, { "value" => $org, "label" => $name };
        }
    }
    if(scalar @organisms_popup == 0){
        push @organisms_popup, {"value" => "null", "label" => "null"};
    }
    print JSON::to_json( \@organisms_popup );

}



#!/usr/bin/env perl
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

my $group = "";
if($query->param("get") eq 'json'){
    print "Content-type: application/json; charset=iso-8859-1\n\n";
}elsif($query->param("get") eq 'display'){
    ### print the header
    &RSA_header("Supported " . $query->param("supported"), "results");
    
    ## Check security issues
    &CheckWebInput($query);
    ## update log file
    &UpdateLogFile();
    &ListParameters() if ($ENV{rsat_echo} >= 2);
    
    if ($ENV{group_specificity}) {
        $group = $ENV{group_specificity};
    }
}

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
    
    if($query->param("get") eq "json"){
        print JSON::to_json( \@taxonomy_popup );
    }else{
        ## Print general information about this RSAT instance
        print "<h2>RSAT instance: ", $ENV{rsat_site}, "</h2>\n";
        
        print "<p><b>Taxon supported: </b>", scalar @taxa, "</p>\n";
        
        if ($group) {
            print "<p><b>Group specificity: </b>", $group, "</p>\n";
        }
        
        foreach $_ (@taxa){
            print $_ . "<br/>";
        }
        
        print '<hr size=3>';
        print "</div>";
        print "</div>";
        
        print $query->end_html;
        
        exit(0);
    }
}
#### Get organisms
else{
    my @selected_organisms = ();
    #### PrintOrthoSelectionSection
    #### infer-operons
    if($supported_orgs eq "infer"){
        push @selected_organisms, &GetOrganismsForTaxon("Bacteria")
        if (($group_specificity eq "Bacteria") ||
        ($group_specificity eq "Prokaryotes"));
        
        push @selected_organisms, &GetOrganismsForTaxon("Archaea")
        if (($group_specificity eq "Archaea") ||
        ($group_specificity eq "Prokaryotes"));
        
        @selected_organisms = sort(@selected_organisms);
    }
    ###### variations
    else{
        my @selected_organisms_ = &RSAT::OrganismManager::get_supported_organisms_web();
        @selected_organisms = ();
        if($supported_orgs eq "orthologs"){
            foreach my $org (@selected_organisms_){
                if($main::supported_organism{$org}->{"blast_available"} eq "1"){
                    push @selected_organisms, $org;
                }
            }
        }elsif($supported_orgs eq "variations"){
            foreach my $org (@selected_organisms_){
                if($main::supported_organism{$org}->{"variant_available"} eq "1"){
                    push @selected_organisms, $org;
                }
            }
        }else{
            @selected_organisms = @selected_organisms_;
        }
    }
    
    if($query->param("get") eq "json"){
    
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
    }else{
        ## Print general information about this RSAT instance
        print "<h2>RSAT instance: ", $ENV{rsat_site}, "</h2>\n";
        $orglist = "";
        if($supported_orgs eq "variations"){
            $orglist .= " with variations";
        }
        print "<p><b>Organisms $orglist supported: </b>", scalar @selected_organisms, "</p>\n";
        
        if ($group) {
            print "<p><b>Group specificity: </b>", $group, "</p>\n";
        }
        
        foreach $org (@selected_organisms){
            $org =~ s/\_/ /g;
            print $org . "<br/>";
        }
        print '<hr size=3>';
        print "</div>";
        print "</div>";
        
        print $query->end_html;
        
        exit(0);
    }

}



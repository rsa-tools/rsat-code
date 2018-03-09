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
    }
    ###### variations
    elsif($supported_orgs eq "variations"){
        my $data_rsat=join("/",$ENV{RSAT},"data") ;
        my $supported_variation_organims_file=join ("/",$data_rsat,"supported_organisms_variation.tab");
        my ($var_org) = &OpenInputFile($supported_variation_organims_file);
        while(<$var_org>){
            chomp;
            next unless (/\S/) ; # skip empty rows
            next if (/^;/); # skip comment lines
            next if (/^\#/); # Skip header line
            my $org=$_ ;
            push (@selected_organisms, $org) ;
        }
    }
    else{
        @selected_organisms = &RSAT::OrganismManager::get_supported_organisms_web();
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



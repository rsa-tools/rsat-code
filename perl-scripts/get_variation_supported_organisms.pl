#!/usr/bin/env perl

################################################################
## This script retrieves the list of organisms supporting
## variation tools
require "RSA.lib";
require RSAT::Tree;
require RSAT::TreeNode;
require RSAT::OrganismManager;
require RSAT::organism;

if  (scalar(@ARGV) > 0) {
    if ( $ARGV[0] == "-h"){
	&PrintHelp
    }
    else {
	print "Only option avialable is -h";
	die;
    }
}

package main;
{
    my $data_rsat=join("/",$ENV{RSAT},"data") ;
    print $data_rsat;
    my $supported_variation_organims_file=join ("/",$data_rsat,"supported_organisms_variation.tab");
    open VARORG, ">", $supported_variation_organims_file  or die $!;

    ## Select organims to retrieve variants sequences from
    ## Get supported organims
    my @installed_organisms = &RSAT::OrganismManager::get_supported_organisms();
    ##my @installed_organisms = &RSAT::OrganismManager::get_supported_organisms_with_variations();
    ## Intialize array to store organisms with variation files
    my @org_variations=(); 
    
    foreach my $org_aux  ( @installed_organisms){
	## Check by organims if there is variation file installed
	my $org_var=&RSAT::organism::has_variations($org_aux);
	if ($org_var){
	    #print $org_var."++";
	    push (@org_variations, $org_aux);
	}
    }
    
    print VARORG join("\n", @org_variations);
    exit(0);
}


sub PrintHelp {
  print <<End_of_help;

This script allows to create and update the file containing
the organims supporting variation tools

Author: Jacques.van-Helden\@univ-amu.fr

Usage: 
   perl perl-scripts/get_variation_supported_organisms.pl

The option "-h" givs this help description.
End_of_help
  exit(0);
}

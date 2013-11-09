#!/usr/bin/perl
############################################################
#
# $Id: parse_genes.pl,v 1.25 2011/02/17 05:07:46 rsat Exp $
#
# Time-stamp: <2003-07-10 11:53:04 jvanheld>
#
############################################################

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "config.pl";
require "lib/load_classes.pl";
require "lib/util.pl";
require "lib/loading_util.pl"; ### for converting polypeptide IDs into ACs
require "lib/parsing_util.pl";


package KEGG::Gene;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "gene_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     names=>"ARRAY",
			     organism=>"SCALAR",
			     description=>"SCALAR",
			     chrom_position=>"SCALAR",
			     chromosome=>"SCALAR",
			     strand=>"SCALAR",
			     start=>"SCALAR",
			     end=>"SCALAR",
			     source=>"SCALAR",
			     xrefs=>"EXPANDED"
			     );
}

package main ;
{

    ################################################################
    #### initialization

    #### oracle schema
    $schema = "kegg";
    $user = "kegg";
    $password = "kegg";

    
    ### files to parse
    @selected_organisms= ();
    
#    $dir{KEGG} = "${Databases}/kegg.genome.ad.jp/pub"; ### older version (Dec 2002)
    #$dir{genes} = $dir{KEGG}."/genomes/previous_genes/";
#    $dir{genes} = $dir{KEGG}."/genomes/genes_2000/";


    $dir{KEGG} = "${Databases}/ftp.genome.ad.jp/pub/kegg";
    $dir{genes} = $dir{KEGG}."/genomes/genes/";

    foreach $d (keys %dir) {
	unless (-d ($dir{$d}) ) {
	    die "Error: $d directory $dir{$d} does not exist\n";
	}
    }
    
    $kegg_file{yeast} = "S.cerevisiae.ent";
    $kegg_file{human} = "H.sapiens.ent";
    $kegg_file{ecoli} = "E.coli.ent";
    
    @all_files = ();
    push @all_files, glob($dir{genes}."/*\.ent");
    push @all_files, glob($dir{genes}."/*\.ent\.gz");
    foreach $file (@all_files) {
	if ($file =~ /\/([^\/]*)$/) {
	    $file = $1;
	}
	my $organism = $file;
	$organism =~ s/\.ent$//;
	$organism =~ s/\.ent\.gz$//;
	$organism_name{$organism} = $organism;
	$kegg_file{lc($organism)} = $file;
	push @all_organisms, $organism;
    }

#    $organism_name{"S.cerevisiae"} = "Saccharomyces cerevisiae";
#    $organism_name{"E.coli"} = "Escherichia coli";
#    $organism_name{"H.sapiens"} = "Homo sapiens";

#    $organism_name{"yeast"} = "Saccharomyces cerevisiae";
#    $organism_name{"ecoli"} = "Escherichia coli";
#    $organism_name{"human"} = "Homo sapiens";

    $verbose = 0;
    $out_format = "obj";

    #### classes and classholders
    @classes = qw( KEGG::Gene );
    $genes = classes::ClassFactory->new_class(object_type=>"KEGG::Gene",
					   prefix=>"gene_");

    #### default output directory
    $export_subdir = "kegg_genes";
    $dir{output} = "$parsed_data/$export_subdir/$delivery_date";


    #### read command-line arguments
    &ReadArguments;
    
    #### check output directory
    &CheckOutputDir();

    $outfile{error} = "$dir{output}/gene.errors.txt";
    $outfile{stats} = "$dir{output}/gene.stats.txt";
    $outfile{genes} = "$dir{output}/gene.obj" if ($export{obj});
    ### open error report file
    open ERR, ">$outfile{error}" || die "Error: cannot write error file $outfile{error}\n";


    ### default output fields for each class
    if ($rsa) {
        #### specific export format for RSA-tools
	$single_name = 1;
	$genes->set_out_fields(qw( id name chromosome start_pos end_pos strand description chrom_position names xrefs));
    } else {

	$genes->set_out_fields(qw( id
				   source 
				   organism
				   chrom_position
				   chromosome
				   strand
				   start_pos
				   end_pos 
				   description
				   names
				   exons
				   introns
				   xrefs
				   dblinks
				   ECs));

	#@{$out_fields{'KEGG::Gene'}} = qw( id source organism raw_position chromosome strand start_base end_base description names exons );
    }



    ### select all organisms if none was selected (-org)
    unless ($#selected_organisms >= 0) {
	@selected_organisms = @all_organisms;
    }

    foreach $org (@selected_organisms) {
	die "; Fatal Error: no data file for organism $org\n"
	    unless (defined($kegg_file{lc($org)}));
	my $data_file = $dir{genes}.$kegg_file{lc($org)};
	if (-e $data_file) {
	    if ($data_file =~ /\.gz$/) {
		$in_file{$org} = "gunzip -c ${data_file} | ";
	    } else {
		$in_file{$org} = "cat ${data_file} | ";
	    }
	} elsif (-e "${data_file}.gz") {
	    $in_file{$org} = "gunzip -c ${data_file}.gz | ";
	} else {
	    die ("Error: cannot find data file for organism ", $org, "\n",
		 "\t", $data_file{$org}, "\n");
	}
    }

    ### test conditions
    if ($test) {
	warn ";TEST\n" if ($verbose >= 1);
	### fast partial parsing for debugging
	foreach $key (keys %in_file) {
	    $in_file{$key} .= " head -1000 |";
	}
    }

    &DefaultVerbose if ($verbose >= 1);
    warn "; Selected organisms\n;\t", join("\n;\t", @selected_organisms), "\n"
	if ($verbose >= 1);

 
    ### parse data from original files
    foreach $org (@selected_organisms) {
	warn join ("\t", "; Organism name", $org, $organism_name{$org}), "\n" 
	    if ($verbose >= 1);
	&ParseKeggFile($in_file{$org}, 
		       $genes, 
		       organism_name=>$organism_name{$org},
		       source=>"KEGG:".$kegg_file{lc($org)});

    }

    &ParsePositions($genes);

    foreach $gene ($genes->get_objects()) {
	### define a single name  (take the first value in the name list)
	if ($single_name) {
	    if ($name = $gene->get_name()) {
		$gene->set_attribute("name",$name);
	    } else {
		$gene->set_attribute("name",$gene->get_id());
	    }
	}

	#### check for genes without description
	if (($gene->get_attribute("description") eq "") ||
	    ($gene->get_attribute("description") eq $null)) {
	    &ErrorMessage(join ("\t", 
				"Warning",
				"no description for gene", 
				$gene->get_name(),
				), "\n");
	    $gene->set_attribute("description",$gene->get_attribute("organism")." ".$gene->get_name());
	} else {
	    my $description = $gene->get_attribute("description");

	    #### extract the swissprot link from the description
	    if ($description =~ /\[SP\:(\S+)\]/) {
		my $swissprot_id = $1;
		$gene->push_expanded_attribute("xrefs", "swissprot", $swissprot_id);
	    } else {
		&ErrorMessage(join ("\t", "Warning", "no swissprot reference", $gene->get_attribute("id"),  "name", $gene->get_name(), $gene->get_attribute("organism")), "\n");
	    }
	    
	    #### extract EC numbers from the description
	    if ($description =~ /\[EC\:([\d\.\s]+)\]/) {
		my $ec_string = $1;
		my @ecs = split /\s+/, $ec_string;
		foreach $ec (@ecs) {
		    $gene->push_attribute("ECs", $ec);
		}
	    }
	}

	### check organism attribute
	unless ($gene->get_attribute("organism")) {
	    &ErrorMessage("Warning: gene ", $gene->get_attribute("id"), " has no organism attribute\n");
	    $gene->set_attribute("organism",$null);
	    next;
	}
	    
	#### choose the full organism name instead of the abbreviated name
	if ($gene->get_attribute("organism_name")) {
	    $gene->force_attribute("organism",$gene->get_attribute("organism_name"));
	}

	#### transform dblinks (single string) to xrefs 
	#### (expanded attribute with database name + external ID)
	my @dblinks = $gene->get_attribute("dblinks");
	foreach my $xref (@dblinks) {
	    $xref =~ s/: /\t/;
	    $gene->push_expanded_attribute("xrefs", $xref);
	}
	
    }


    ### print result
    &PrintStats($outfile{stats}, @classes);
    $genes->dump_tables();
    $genes->generate_sql(schema=>$schema, 
			 user=>$user,
			 password=>$password,
			 dir=>"$dir{output}/sql_scripts",
			 prefix=>"",
			 dbms=>$dbms
			 );
    &ExportClasses($outfile{genes}, $out_format, @classes) if $export{obj};


    close ERR;

    CompressParsedData();

    ### report execution time
    if ($verbose >= 1) {
	$done_time = &AlphaDate;
	warn ";\n";
	warn "; job started $start_time";
	warn "; job done    $done_time\n";
    }

    exit(0);
}


################################################################
############## subroutines for the main package ################
################################################################

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
	parse_kegg.pl

DESCRIPTION
	Parse genes from a KEGG file (http://www.genome.ad.jp/kegg/). 
	
AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

VERSION
	0.01
	Created		1999/12/16
	Last modified	2000/01/08
	
OPTIONS	
	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-outdir output directory
	-v #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-obj    Export results in "obj" format (human readable)
	-enz	export enzymes only
	-all	export all polypeptides for the selected organisms
		(default)
	-org	select an organism for exportation
		can be used reiteratively in the command line 
		to select several organisms
		by default, all organisms found in the input directory
		are selected
	-name
		exports a name as single value attribute in
		the main table (this is redundant but can be useful)
	-clean	remove all files from the output directory before
		parsing
EndHelp
  close HELP;
}

  



### read arguments from the command line
sub ReadArguments {
    for my $a (0..$#ARGV) {


	&ReadGenericOptions($a);
	
	    
	    ### export single name in main table
	if ($ARGV[$a] eq "-name") {
	    $main::single_name = 1;
	    
	    ### specific export format  for RSA-tools
	} elsif ($ARGV[$a] eq "-rsa") {
	    $main::rsa = 1;
	    
	    ### select enzymes for exportation
	} elsif ($ARGV[$a] =~ /^-org/) {
	    push @selected_organisms, $ARGV[$a+1];
	    
	    ### select enzymes for exportation
	} elsif ($ARGV[$a] =~ /^-enz/) {
	    $main::export{enzymes} = 1;

	} elsif ($ARGV[$a] =~ /^-all/) {
	    $main::export{all} = 1;

	}
	
    }
}



#!/usr/bin/perl
############################################################
#
# $Id: parse_genbank.pl,v 1.4 2002/12/09 00:22:31 jvanheld Exp $
#
# Time-stamp: <2002-10-25 12:05:54 jvanheld>
#
############################################################

### parse_genbank.pl
### type parse_genbank.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "PFBP_classes.pl";
require "PFBP_config.pl";
require "PFBP_util.pl";
require "PFBP_loading_util.pl"; ### for converting polypeptide IDs into ACs
require "PFBP_parsing_util.pl";


package Genbank::Contig; ### for parsing genbank files
{
  @ISA = qw ( PFBP::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "ctg_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     names=>"ARRAY",
			     organism=>"SCALAR",
			     type=>"SCALAR",
			     xrefs=>"EXPANDED"
			     );
}

package Genbank::Feature;
{
  @ISA = qw ( PFBP::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "ft_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     names=>"ARRAY",
			     organism=>"SCALAR",
			     type=>"SCALAR",
			     description=>"SCALAR",
			     chrom_position=>"SCALAR",
			     chrom_position=>"SCALAR",
			     chromosome=>"SCALAR",
			     strand=>"SCALAR",
			     start_pos=>"SCALAR",
			     end_pos=>"SCALAR",
			     source=>"SCALAR",
			     note=>"ARRAY",
			     xrefs=>"EXPANDED"
			     );

}




package main;
{
    
    #### initialization
    $test;
    $suffix = "";
    $verbose = 0;
    $out_format = "obj";
    
    
    #### selected organisms
#    @selected_organisms = ();
#    push @selected_organisms, qw ( Bacteria/Mycoplasma_genitalium ); 
#    push @selected_organisms, qw ( Bacteria/Escherichia_coli_K12 );
#    push @selected_organisms, qw ( Saccharomyces_cerevisiae ); #### problem with multiple chromosomes

    #### input directory
    $dir{genbank} = $Databases."/ftp.ncbi.nih.gov/genomes/";
    
    #### default export directory
    $export_subdir = "genbank";
    $dir{output} = "$parsed_data/${export_subdir}/$delivery_date";

    #### temporary
    if ($hostname =~ /^brol/i) {
	$dir{genbank} = "/lin/genomics/genbank/ftp.ncbi.nih.gov/genomes/";
    }

    #### classes and classholders
    $features = PFBP::ClassFactory->new_class(object_type=>"Genbank::Feature",
					   prefix=>"feature_");

    $contigs = PFBP::ClassFactory->new_class(object_type=>"Genbank::Contig",
					  prefix=>"contig_");
    @classes = qw( Genbank::Feature Genbank::Contig );

    #### read command arguments
    &ReadArguments();

    ################################################################
    #### check arguments


    ### select all organisms if none was selected (-org)
    if (($export{all}) || 
	($#selected_organisms < 0)) {
	@selected_organisms = &SelectAllOrganisms();
    }
    #### check existence of input directories
    foreach my $org (@selected_organisms) {
	my $org_dir = "${org}/";
	my $data_dir = $dir{genbank}."/${org_dir}";
	die ("Error: directory ", $data_dir, " does not exist.\n") 
	    unless (-d $data_dir);
    }


    ### default output fields
    @out_fields = qw( id type name chromosome start_pos end_pos strand description chrom_position names  exons introns db_xref note);
    if ($rsa) {
	#### specific export format for RSA-tools
	warn "; Exporting in special format for rsa-tools.\n" if ($verbose >=1);
	$features->set_out_fields(@out_fields);
	$single_name = 1;
    } else {
	$features->set_out_fields(@out_fields, qw(source organism));
    }
    $contigs->set_out_fields(qw(id	
				organism
				type
				length
				form	
				genbank_file
				
				));
    
    &CheckOutputDir();
    chdir $dir{output};
    $out_file{features} = "$dir{output}/feature.obj" if ($export{obj});
    $out_file{error} = "$dir{output}/feature.errors.txt";
    $out_file{stats} = "$dir{output}/feature.stats.txt";
    
    ### open error report file
    open ERR, ">$out_file{error}" || die "Error: cannot write error file $out_file{error}\n";
    

    #### define organism-specific data directories
    my %data_dir = ();
    foreach my $org (@selected_organisms) {
	my $org_dir = "${org}/";
	$data_dir{$org} = $dir{genbank}."/${org_dir}";
#	die ("Error: directory ", $data_dir, " does not exist.\n") 
#	    unless (-d $data_dir);

#  	@genbank_files = glob($data_dir."/*.gbk");
#  	push @genbank_files, glob($data_dir."/*.gbk.gz"); #### compressed files are supported
	
#  #die join "\n", $#genbank_files, @genbank_files;
#  	if ($#genbank_files <= -1) {
#  	    die ("Error: there is not genbank file in the directory ", 
#  		 $data_dir,
#  		 "\n");
#  	} elsif ($#genbank_files > 0) {
#  	    warn ("Error: there are several genbank files in the directory ", 
#  		 $data_dir, "\n\t",
#  		 join ("\n\t", @genbank_files),
#  		 "\n");
#  	    next;
#  	} else {
#  	    $in_file{$org} = $genbank_files[0];
#  	    push @parsed_organisms, $org;
#  	}

#  	$short_file{$org} = `basename $in_file{$org} .gbk`;
#  	chomp $short_file{$org};
	
#  	if (-e $in_file{$org}) {
#  	    if ($in_file{$org} =~ /\.gz$/) {
#  		$in_stream{$org} = "gunzip -c $in_file{$org} | ";
#  	    } else {
#  		$in_stream{$org} = "cat $in_file{$org} | ";
#  	    }
#  	} elsif (-e "$in_file{$org}.gz") {
#  	    $in_stream{$org} = "gunzip -c $in_file{$org}.gz | ";
#  	} else {
#  	    die ("Error: cannot find data file for organism ", $org, "\n",
#  		 "\t", $in_file{$org}, "\n");
#  	}
       
#  	### test conditions
#  	if ($test) {
#  	    $in_stream{$org} .= " head -1000 |";
#  	}
    }
    
    if ($verbose >=1) {
	warn ";TEST\n" if ($test); 
	
	&DefaultVerbose;

	warn "; Selected organisms\n;\t", join("\n;\t", @selected_organisms), "\n";
	warn "; Parsed organisms\n";
	foreach my $org (@selected_organisms) {
	    warn ";\t$org\tdir\t",$data_dir{$org},"\n";
	}
    }
    
    ### parse data from original files
    foreach $org (@selected_organisms) {
	warn "; Parsing genome data for organism $org\n" if ($verbose >=1);
	chdir ($data_dir{$org});
	my @genbank_files = ();
	push @genbank_files, glob("*.gbk");
	push @genbank_files, glob("*.gbk.gz");
	if ($#genbank_files < 0) {
	    &FatalError("There is no genbank file in the input directory $data_dir{$org}\n");
	} else {
	    warn "; \tGenbank files\n;\t", join("\n;\t", @genbank_files), "\n" if ($verbose >= 1);
	}
	foreach my $file (@genbank_files) {
	    $file{input} = "$data_dir{$org}/$file";
	    
	    #### check whether the file has to be uncompressed on the fly
	    if ($file{input} =~ /.gz/) {
		$file{input} = "gunzip -c $file{input} |";
	    } else {
		$file{input} = "cat $file{input} |";
	    }
	    
	    #### for quick testing; only parse the first 10000 lines
	    if ($test) {
		$file{input} .= " head -10000 | ";
	    }
	    &ParseGenbankFile($file{input}, 
			      $features, 
			      $contigs, 
			      source=>"genbank:".$short_file{$org},
			      no_seq=>1);
	}
	
#  	&ParseGenbankFile($in_stream{$org}, 
#  			  $features, 
#  			  $contigs, 
#  			  source=>"genbank:".$short_file{$org},
#  			  no_seq=>1);
	
    }
    
    &ParsePositions($features);

    &GuessSynonyms($features);
    
    
    foreach $feature ($features->get_objects()) {
	foreach my $name ($feature->get_attribute("gene")) {
	    $feature->push_attribute("names",$name);
	};
	#$feature->push_attribute("names", $feature->get_attribute("gene"));
	
	### define a single name  (take the first value in the name list)
	if ($single_name) {
	    if ($name = $feature->get_name()) {
		$feature->set_attribute("name",$name);
	    } else {
		$feature->set_attribute("name",$feature->get_id());
	    }
	}
	
	#### check for features without description
	if (($feature->get_attribute("description") eq $null) 
	    || ($feature->get_attribute("description") eq "")) {
	    $feature->set_attribute("description",$feature->get_attribute("product"));
	}

	################################################################
	#### cross-references
	my @xrefs = $feature->get_attribute("db_xref");
	my $gi = "";
	foreach my $xref (@xrefs) {
#	    my @fields = split ":", $xref;
	    #### use GI as feature identifier
	    if ($xref =~ /GI:/) {
		$gi = $';
		last;
	    } 
	}

	if ($gi) {
	    $feature->force_attribute("id",$gi);
	} else {
	    &ErrorMessage("; Error\tfeature ".$feature->get_attribute("id")." has no GI.\n"); 
	}

	#### use genbank name as chromosome name
	my $source = $feature->get_attribute("source");
	if ($source =~ /genbank:/) {
	    my $chromosome = $';
	    $chromosome =~ s/\.gz$//;
	    $chromosome =~ s/\.gbk$//;
	    $feature->force_attribute("chromosome",$chromosome);
	}



    }
    
    
    ################################################################
    ### export result in various formats
    chdir $dir{output};
    &PrintStats($out_file{stats}, @classes);
#    $features->dump_tables();
#    $contigs->dump_tables();
    @class_factories = qw (
			   features
			   contigs
			   );
    foreach $class_factory (@class_factories) {
	
	$$class_factory->dump_tables();
	$$class_factory->generate_sql(schema=>$schema,
				      dir=>"$dir{output}/sql_scripts",
				      prefix=>"g_$class_factory_",
				      dbms=>$dbms
				      );
    }
    &ExportClasses($out_file{features}, $out_format, @classes) if $export{obj};


    
    
    ### report execution time
    if ($verbose >= 1) {
	$done_time = &AlphaDate;
	warn ";\n";
	warn "; job started $start_time";
	warn "; job done    $done_time\n";
    }
    
    close ERR;
    
    &CompressParsedData();
    
    exit(0);
    
}

### subroutines for the main package

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
	parse_genbank.pl

DESCRIPTION
	Parse features from a Genbank file 
	
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
	-all	export all polypeptides for the selected organisms
		(default)
	-org	select an organism for exportation
		can be used reiteratively in the command line 
		to select several organisms
		by default, all organisms found in the input directory
		are selected
	-suffix suffix
		add a suffix to output dir name
	-clean	remove all files from the output directory before
		parsing
	-name
		exports a name as single value attribute in
		the main table (this is redundant but can be useful)

EndHelp
  close HELP;
}

  



### read arguments from the command line
sub ReadArguments {
    for my $a (0..$#ARGV) {
	
	### warn level
	if (($ARGV[$a] eq "-v" ) && 
	    ($ARGV[$a+1] =~ /^\d+$/)){
	    $main::verbose = $ARGV[$a+1];
	    
	    ### test run
	} elsif ($ARGV[$a] eq "-test") {
	    $main::test = 1;
	    
	    ### export single name in main table
	} elsif ($ARGV[$a] eq "-name") {
	    $main::single_name = 1;
	    
	    ### specific export format  for RSA-tools
	} elsif ($ARGV[$a] eq "-rsa") {
	    $main::rsa = 1;
	    
	    ### clean
	} elsif ($ARGV[$a] eq "-clean") {
	    $main::clean = 1;
	    
	    ### suffix
	} elsif ($ARGV[$a] eq "-suffix") {
	    $a++;
	    $main::suffix = $ARGV[$a];
	    
	    ### output dir
	} elsif ($ARGV[$a] eq "-outdir") {
	    $a++;
	    $main::dir{output} = $ARGV[$a];
	    
	    ### help
	} elsif (($ARGV[$a] eq "-h") ||
		 ($ARGV[$a] eq "-help")) {
	    &PrintHelp;
	    exit(0);

	    ### select enzymes for exportation
	} elsif ($ARGV[$a] =~ /^-org/) {
	    push @selected_organisms, $ARGV[$a+1];

	    #### export all organisms
	} elsif ($ARGV[$a] =~ /^-all/) {
	    $main::export{all} = 1;

	    #### export object file
	} elsif ($ARGV[$a] =~ /^-obj/) {
	    $main::export{obj} = 1;
	}
	
    }
}


sub SelectAllOrganisms {
    my @organisms = ();
    my $dir = $dir{genbank};
    die "Error: directory $dir does not exist.\n"
	unless (-d $dir);
    chdir $dir;
    @organisms = glob "Bacteria/*_*";
    die "Error: there are no organisms in the directory $dir.\n"
	unless $#organisms >=0;
    return @organisms;
}


#!/usr/bin/perl
############################################################
#
# $Id: parse_genbank.pl,v 1.14 2011/02/17 05:07:46 rsat Exp $
#
# Time-stamp: <2003-10-01 15:41:54 jvanheld>
#
############################################################

### parse_genbank.pl
### type parse_genbank.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "lib/load_classes.pl";
require "config.pl";
require "lib/util.pl";
#require "lib/loading_util.pl"; ### for converting polypeptide IDs into ACs
require "lib/parsing_util.pl";
require "lib/parse_genbank_lib.pl";

require "classes/Genbank_classes.pl";


################################################################
#### main package
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


    #### classes and classholders
    
    $features = classes::ClassFactory->new_class(object_type=>"Genbank::Feature", prefix=>"ft_");

    $genes = classes::ClassFactory->new_class(object_type=>"Genbank::Gene", prefix=>"gn_");

    $mRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::mRNA", prefix=>"mRNA_");

    $scRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::scRNA", prefix=>"scRNA_");

    $tRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::tRNA", prefix=>"tRNA_");
    
    $rRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::rRNA", prefix=>"rRNA");

    $misc_RNAs = classes::ClassFactory->new_class(object_type=>"Genbank::misc_RNA", prefix=>"misc_RNA");

    $misc_features = classes::ClassFactory->new_class(object_type=>"Genbank::misc_feature", prefix=>"misc_feature");

    $CDSs = classes::ClassFactory->new_class(object_type=>"Genbank::CDS", prefix=>"cds_");

    $sources = classes::ClassFactory->new_class(object_type=>"Genbank::Source", prefix=>"src_");

    $contigs = classes::ClassFactory->new_class(object_type=>"Genbank::Contig", prefix=>"ctg_");

    $organisms = classes::ClassFactory->new_class(object_type=>"Genbank::Organism", prefix=>"org_");

#      $features = classes::ClassFactory->new_class(object_type=>"Genbank::Feature",
#  					   prefix=>"ft_");

#      $genes = classes::ClassFactory->new_class(object_type=>"Genbank::Gene",
#  					   prefix=>"gn_");

#      $mRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::mRNA",
#  					   prefix=>"mRNA_");

#      $CDSs = classes::ClassFactory->new_class(object_type=>"Genbank::CDS",
#  					   prefix=>"CDS_");

#      $contigs = classes::ClassFactory->new_class(object_type=>"Genbank::Contig",
#  					  prefix=>"contig_");

#     $organisms = classes::ClassFactory->new_class(object_type=>"Genbank::Organism",
# 					  prefix=>"org_");

    @classes = qw( Genbank::Feature Genbank::Contig Genbank::Organism Genbank::Gene Genbank::CDS Genbank::mRNA Genbank::scRNA);

    #### read command arguments
    &ReadArguments();

    ################################################################
    #### check arguments


    ### select all organisms if none was selected (-org)
    if (($export{all}) || 
	   ($#selected_organisms < 0)) {
	   %data_dir = &SelectAllOrganisms();
	   foreach my $org (keys %data_dir) {
	   		push @selected_organisms, $org;
	   }
    } else {
	    #### check existence of input directories
	    foreach my $org (@selected_organisms) {
		my $org_dir = "${org}/";
		my $data_dir = $dir{genbank}."${org_dir}";
		die ("Error: directory ", $data_dir, " does not exist.\n") 
		    unless (-d $data_dir);
		$data_dir{$org} = $data_dir;
    }


    ################################################################
    ### default output fields

    ### output fields for features
    @feature_out_fields = qw( id type name chromosome start_pos end_pos strand description chrom_position names  exons introns db_xref note);
    if ($rsa) {
	#### specific export format for RSA-tools
	warn "; Exporting in special format for rsa-tools.\n" if ($verbose >=1);
	$features->set_out_fields(@feature_out_fields);
	$single_name = 1;
    } else {
	$features->set_out_fields(@feature_out_fields, qw(source organism));
    }

    ### output fields for contigs
    $contigs->set_out_fields(qw(id
				organism
				type
				length
				form
				genbank_file
			       ));

    #### output directory
    &CheckOutputDir();
    chdir $dir{output};
    $outfile{features} = "$dir{output}/genbank.obj" if ($export{obj});
    $outfile{error} = "$dir{output}/genbank.errors.txt";
    $outfile{stats} = "$dir{output}/genbank.stats.txt";
    
    ### open error report file
    open ERR, ">$outfile{error}" || die "Error: cannot write error file $outfile{error}\n";
    

    #### define organism-specific data directories
#    foreach my $org (@selected_organisms) {
#	my $org_dir = "${org}/";
#	$data_dir{$org} = $dir{genbank}."${org_dir}";
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
	
	&DefaultVerbose();

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
			  $genes,
			  $mRNAs,
			  $scRNAs,
			  $tRNAs,
			  $rRNAs,
			  $misc_RNAs,
			  $misc_features,
			  $CDSs,
			  $contigs, 
			  $organisms, 
			  $sources,
#			  source=>"Genbank", 
#			  seq_dir=>$dir{output}
			  );
#  	    &ParseGenbankFile($file{input}, 
#  			      $features,
#  			      $genes,
#  			      $mRNAs,
#  			      $scRNAs,
#  			      $CDSs,
#  			      $contigs, 
#  			      $organisms, 
#  			      source=>$contig, 
#  			      seq_dir=>$dir{output});
#  	    &ParseGenbankFile($file{input}, 
#  			      $features, 
#  			      $contigs, 
#  			      $organisms, 
#  			      source=>"genbank:".$short_file{$org},
#  			      no_seq=>1);
	}
	
	
    }
    
    

    #### index names
    $features->index_names();
    $genes->index_names();
    $mRNAs->index_names();
    $scRNAs->index_names();
    $CDSs->index_names();


    #### parse chromosomal poitions
    &ParsePositions($features);
    &ParsePositions($genes);
    &ParsePositions($mRNAs);
    &ParsePositions($scRNAs);
    &ParsePositions($CDSs);
    
    #### Create features from CDSs and mRNAs
    &CreateGenbankFeatures($features, $genes, $mRNAs, $scRNAs, $CDSs);


    ################################################################
    #### export result 
    chdir $dir{output};

    #### print parsing statistics
    &PrintStats($outfile{stats}, @classes);

    @class_factories = qw (
			   organisms 
			   contigs
			   features
			   genes
			   mRNAs
			   scRNAs
			   CDSs
			   );

    foreach $class_factory (@class_factories) {	
	$$class_factory->dump_tables();
	$$class_factory->generate_sql(schema=>$schema,
				      dir=>"$dir{output}/sql_scripts",
				      prefix=>"g_$class_factory_",
				      dbms=>$dbms
				      );
    }
    &ExportClasses($outfile{features}, $out_format, @classes) if $export{obj};


    
    
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
	Parse one or sveral Genbank files for extracting genome
	information.
	
	Genbank genomes can be retrieved by anonymous ftp :
		ftp://ftp.ncbi.nlm.nih.gov/genomes

	Genome sequence and annotations are parsed from the genbank
	flat files (extension .gbk).  

	Each directory contains the genome of one organism.  Note that

	  - a single directoy can contain several files, if the
	    organism has several chromosomes.  

          - a single .gbk file can contain several contigs, if the
	    genome is not fully assembled for example.

	Parsed data is exported in tab-delimited format, according to
	the normalization rules for relational databases : one main
	table regroups all the single-value attributes(one attribute
	per column), and each multi-value attribute comes in a
	separate table with the foregin key in the first column.

	The program also exports some parsing statistics, an error
	log, and a set of SQL scripts for creating a relational
	database (supported standards: mysql, postgresql, oracle).


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
	-indir  input directory (with genbank .gbk files)
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
	    
	    ### input dir (with genbank files)
	} elsif ($ARGV[$a] eq "-indir") {
	    $a++;
	    $main::dir{genbank} = $ARGV[$a];
	    
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

## find all gbk files in the genbank directory, and deduce the directories containing an organism
sub SelectAllOrganisms {
    my @organisms = ();
    my $dir = $dir{genbank};
    die "Error: directory $dir does not exist.\n" unless (-d $dir);
#    chdir $dir;
#    @organisms = glob "Bacteria/*_*";   ####problem: recognizes all_bacteria files (not directories!)
	$gbk_dirs = `find $dir -name '*.gbk*' -maxdepth 2 -exec dirname {} \\; | sort -u`;
	$gbk_dirs .= `find $dir/Bacteria -name '*.gbk*' -maxdepth 2 -exec dirname {} \\; | sort -u`;
	
	@gbk_dirs = split "\n", $gbk_dirs;
	my %gbk_dirs = ();
	foreach my $d (@gbk_dirs) {	
		my $organism_name = `basename $d`;
		chomp($organism_name);
		$gbk_dirs{$organism_name} = $d;
		warn join( "\t", $organism_name, $d), "\n" if ($verbose >= 2);
	}
	
	@organisms = keys %gbk_dirs;
	unless ($#organisms >=0) {
	    die "Error: could not identify the organisms in the directory $dir.\n" ;
	}
    return %gbk_dirs;
}


#!/usr/bin/perl 
############################################################
#
# $Id: parse-genbank.pl,v 1.14 2003/12/17 18:09:18 jvanheld Exp $
#
# Time-stamp: <2003-10-01 16:17:10 jvanheld>
#
############################################################
#use strict;;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "RSA.lib";
push @INC, "$RSA/perl-scripts/parsers/";
require "lib/load_classes.pl";
require "lib/util.pl";
require "lib/parsing_util.pl";
require "classes/Genbank_classes.pl";


################################################################
#### initialization
$no_suffix=1;
$host= $default{'host'};
$schema="genbank";
$user="genbank";
$password="genbank";

$test_lines = 10000;

################################################################
#### main package
package main;
{
    #### initialise parameters ####
    my $start_time = &AlphaDate;
    
    local %infile = ();
    local %outfile = ();
    
    local $verbose = 0;
    local $in = STDIN;
    local $out = STDOUT;

    local $single_name = 1;
    
    $features = classes::ClassFactory->new_class(object_type=>"Genbank::Feature", prefix=>"ft_");

    $genes = classes::ClassFactory->new_class(object_type=>"Genbank::Gene", prefix=>"gn_");

    $mRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::mRNA", prefix=>"mRNA_");

    $tRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::tRNA", prefix=>"tRNA_");
    
    $rRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::rRNA", prefix=>"rRNA");

    $misc_RNAs = classes::ClassFactory->new_class(object_type=>"Genbank::misc_RNA", prefix=>"misc_RNA");

    $misc_features = classes::ClassFactory->new_class(object_type=>"Genbank::misc_feature", prefix=>"misc_feature");

    $CDSs = classes::ClassFactory->new_class(object_type=>"Genbank::CDS", prefix=>"cds_");

    $sources = classes::ClassFactory->new_class(object_type=>"Genbank::Source", prefix=>"src_");

    $contigs = classes::ClassFactory->new_class(object_type=>"Genbank::Contig", prefix=>"ctg_");

    $organisms = classes::ClassFactory->new_class(object_type=>"Genbank::Organism", prefix=>"org_");

    $features->set_out_fields(qw(id 
				 type
				 name
				 contig
				 start_pos
				 end_pos
				 strand
				 description
				 chrom_position
				 organism

				 names
				 db_xref
				 introns
				 exons
				 EC_number
				 ));
    
    @classes = qw( Genbank::Feature
		   Genbank::Contig
		   Genbank::Organism
		   Genbank::Gene
		   Genbank::CDS
		   Genbank::mRNA 
		   Genbank::tRNA 
		   Genbank::rRNA 
		   Genbank::misc_RNA 
		   Genbank::misc_feature 
		   Genbank::Source);
    
    &ReadArguments();
    
    ################################################################
    #### check argument values ####

    #### input directory
    unless (defined($dir{input})) {
	&FatalError("You must specify the input directory.\n");
    }
    unless (-d $dir{input}) {
	&FatalError("Input directory '$dir{input}' does not exist.\n");
    }

    #### organism name
    unless (defined($org)) {
	$org = `basename $dir{input}`;
	chomp($org);
	warn "; Auto selection of organism name\t$org\n" if ($verbose >= 1);
    }

    ################################################################
    #### find genbank files in the input directory
    $dir{main} = `pwd`; #### remember working directory
    chomp($dir{main});
    chdir ($dir{input});
    push @genbank_files, glob("*.gbk");
    push @genbank_files, glob("*.gbk.gz");
    push @genbank_files, glob("*.genomic.gbff");
    push @genbank_files, glob("*.genomic.gbff.gz");
    if ($#genbank_files < 0) {
	&FatalError("There is no genbank file in the input directory $dir{input}\n");
    } else {
	warn "; Genbank files\n;\t", join("\n;\t", @genbank_files), "\n" if ($verbose >= 1);
    }
    chdir($dir{main});     #### come back to the starting directory
    


    #### output directory
    unless (defined($dir{output})) {
	if ($dir{input} =~ /refseq/) {
	    $dir{output} = $parsed_data."/refseq/gbff/".$delivery_date;
	    warn "; Auto selection of output dir\t$dir{output}\n" if ($verbose >= 1);
	} else {
	    $dir{output} = "$RSA/data/genomes/".$org."/genome";
	    warn "; Auto selection of output dir\t$dir{output}\n" if ($verbose >= 1);
	}
    }
    &CheckOutputDir();
    chdir $dir{main};
    $out_file{features} = "$dir{output}/genbank.obj" if ($export{obj});
    $out_file{error} = "$dir{output}/genbank.errors.txt";
    $out_file{stats} = "$dir{output}/genbank.stats.txt";

    ### open error report file
    open ERR, ">$out_file{error}" 
	|| die "Error: cannot write error file $out_file{error}\n";

    #### verbose ####
    &Verbose() if ($verbose);

    ################################################################
    #### parse the genbank files
    chdir $dir{main};
    foreach my $file (@genbank_files) {
	$file{input} = "$dir{input}/$file";

	#### check whether the file has to be uncompressed on the fly
	if ($file{input} =~ /.gz/) {
	    $file{input} = "gunzip -c $file{input} |";
	} else {
	    $file{input} = "cat $file{input} |";
	}
	
	#### for quick testing; only parse the first  lines
	if ($test) {
	    $file{input} .= " head -$test_lines | ";
	}

	&ParseGenbankFile($file{input}, 
			  $features,
			  $genes,
			  $mRNAs,
			  $tRNAs,
			  $rRNAs,
			  $misc_RNAs,
			  $misc_features,
			  $CDSs,
			  $contigs, 
			  $organisms, 
			  $sources,
#			  source=>"Genbank", 
			  seq_dir=>$dir{output}
			  );
    }

    #### write the chromosome file
    chdir $dir{main};
    $chrom = &OpenOutputFile("$dir{output}/Contigs.txt"); # file with chromosome IDs
    foreach my $contig ($contigs->get_objects()) {
	print $chrom join ("\t", 
			   $contig->get_attribute("file"),
			   $contig->get_attribute("id"),
			   $contig->get_attribute("form")), "\n";
    }
    close $chrom;
    
    #### index names
    $features->index_names();
    $genes->index_names();
    $mRNAs->index_names();
    $tRNAs->index_names();
    $rRNAs->index_names();
    $misc_RNAs->index_names();
    $misc_features->index_names();
    $CDSs->index_names();


    #### Create features from CDSs and RNAs
    &CreateGenbankFeatures($features, $genes, $mRNAs, $tRNAs, $rRNAs, $misc_RNAs, $misc_features, $CDSs, $sources);

    #### parse chromosomal positions
    &ParsePositions($features);
    &ParsePositions($genes);
    &ParsePositions($mRNAs);
    &ParsePositions($tRNAs);
    &ParsePositions($rRNAs);
    &ParsePositions($misc_RNAs);
    &ParsePositions($misc_features);
    &ParsePositions($CDSs);

    ################################################################
    ### export result in various formats
#    chdir $dir{output};

    #### print parsing statistics
    &PrintStats($out_file{stats}, @classes);

    @class_factories = qw (
			   organisms 
			   contigs
			   features
			   genes
			   mRNAs
			   tRNAs
			   rRNAs
			   misc_RNAs
			   misc_features
			   CDSs
			   sources
			   );

    foreach my $factory_name (@class_factories) {
	my $class_factory = $$factory_name;
	warn "; Dumping class $factory_name $class_factory\n" if ($verbose >= 1);
	$suffix = "_$org" unless ($no_suffix);
	$class_factory->dump_tables($suffix);
	$class_factory->generate_sql(dir=>"$dir{output}/sql_scripts",
				     prefix=>"$class_factory_",
				     schema=>$schema,
				     host=>$host,
				     dbms=>$dbms,
				     user=>$user,
				     password=>$password
				     );
    }
    &ExportMakefile(@classes);
    &ExportClasses($out_file{features}, $out_format, @classes) if $export{obj};
    
    ###### verbose ######
    if ($verbose) {
	my $done_time = &AlphaDate();
	print $out "; Job started $start_time\n";
	print $out "; Job done    $done_time\n";
    }
    
    
    ###### close output file ######
    close $out if ($outfile{output});
    
    
    exit(0);
}


################################################################
#################### subroutine definition #####################
################################################################

################################################################
#### display full help message #####
sub PrintHelp {
  open HELP, "| more";
  print HELP <<End_of_help;
NAME
	parse-genbank

        2001 by Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)
	
USAGE
        parse-genbank [-dir input_dir][-i inputfile] [-o outputfile] [-v]

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

CATEGORY
	parser

OPTIONS
	-h	(must be first argument) display full help message
	-help	(must be first argument) display options
	-v	verbose
	-i	input directory
		input directory. This directory must contain one or
		several genbank files (extension .gbk). 
	-o	output directory
		The parsing result will be saved in this directory. If
		the directory does not exist, it will be created.
	-test #	quick test (for debugging): only parse the # first
		lines of each Genabnk file (default $test_lines).

   Options for the automaticaly generated SQL scripts
	-schema database schema (default: $schema)
	-host	database host (efault: $host)
	-user	database user (efault: $user)
	-password	
		database password (default: $password)
End_of_help
  close HELP;
  exit;
}


################################################################
#### display short help message #####
sub PrintOptions {
  open HELP, "| more";
  print HELP <<End_short_help;
parse-genbank options
----------------
-h	(must be first argument) display full help message
-help	(must be first argument) display options
-i	input dir
-o	output dir
-v	verbose
-test #	quick test (for debugging)
-schema database schema (default: $schema)
-host	database host (default: $host)
-user	database user (default: $user)
-password	database password (default: $password)
End_short_help
  close HELP;
  exit(0);
}


################################################################
#### read arguments 
sub ReadArguments {
    foreach my $a (0..$#ARGV) {
	### verbose ###
	if ($ARGV[$a] eq "-v") {
	    if (&IsNatural($ARGV[$a+1])) {
		$verbose = $ARGV[$a+1];
	    } else {
		$verbose = 1;
	    }
	    
	    ### detailed help
	} elsif ($ARGV[$a] eq "-h") {
	    &PrintHelp();
	    
	    ### list of options
	} elsif ($ARGV[$a] eq "-help") {
	    &PrintOptions();
	    
	    ### input file ###
	} elsif ($ARGV[$a] eq "-i") {
	    $dir{input} = $ARGV[$a+1];

	    ### output file ###
	} elsif ($ARGV[$a] eq "-o") {
	    $dir{output} = $ARGV[$a+1];

	    ### quick test
	} elsif ($ARGV[$a] eq "-test") {
	    $test = 1;
	    if (&IsNatural($ARGV[$a+1])) {
		$test_lines = $ARGV[$a+1];
	    }

	    ### database schema
	} elsif ($ARGV[$a] eq "-schema") {
	    $schema = $ARGV[$a+1];
	    
	    ### database host
	} elsif ($ARGV[$a] eq "-host") {
	    $host = $ARGV[$a+1];
	    
	    ### database user
	} elsif ($ARGV[$a] eq "-user") {
	    $user = $ARGV[$a+1];
	    
	    ### password 
	} elsif ($ARGV[$a] eq "-password") {
	    $password = $ARGV[$a+1];
	    
	}
    }
}


################################################################
#### report parameters
sub Verbose {
    print $out "; parse-genbank ";
    &PrintArguments($out);
    if (defined(%dir)) {
	print $out "; Directories\n";
	while (($key,$value) = each %dir) {
	    print $out ";\t$key\t$value\n";
	}
    }
    if (defined(%infile)) {
	print $out "; Input files\n";
	while (($key,$value) = each %infile) {
	    print $out ";\t$key\t$value\n";
	}
    }
    if (defined(%outfile)) {
	print $out "; Output files\n";
	while (($key,$value) = each %outfile) {
	    print $out ";\t$key\t$value\n";
	}
    }
    printf $out "; %-29s\t%s\n", "organism", $org;
    printf $out "; %-29s\t%s\n", "genbank files", join (" ", @genbank_files);
}


#!/usr/bin/perl 
############################################################
#
# $Id: parse-genbank.pl,v 1.6 2003/06/05 22:08:18 jvanheld Exp $
#
# Time-stamp: <2003-06-05 01:14:53 jvanheld>
#
############################################################
#use strict;;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "RSA.lib";
push @INC, $ENV{PARSER};
require "PFBP_classes.pl";
require "Genbank_classes.pl";
require "PFBP_parsing_util.pl";


################################################################
#### initialization
$dbms = "mysql";


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
    
    $features = PFBP::ClassFactory->new_class(object_type=>"Genbank::Feature",
						 prefix=>"ft_");

    $genes = PFBP::ClassFactory->new_class(object_type=>"Genbank::Gene",
					   prefix=>"gn_");

    $mRNAs = PFBP::ClassFactory->new_class(object_type=>"Genbank::mRNA",
					   prefix=>"mRNA_");

    $CDSs = PFBP::ClassFactory->new_class(object_type=>"Genbank::CDS",
					   prefix=>"CDS_");

    $contigs = PFBP::ClassFactory->new_class(object_type=>"Genbank::Contig",
					     prefix=>"ctg_");
    $features->set_out_fields(qw(id 
				 type
				 name
				 contig
				 start_pos
				 end_pos
				 strand
				 description
				 chrom_position

				 names
				 db_xref
				 introns
				 exons
				 EC_number
				 ));
    
    $organisms = PFBP::ClassFactory->new_class(object_type=>"Genbank::Organism",
					  prefix=>"org_");
    @classes = qw( Genbank::Feature Genbank::Contig Genbank::Organism );
    
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
    #### output directory
    unless (defined($dir{output})) {
	$dir{output} = "$RSA/data/genomes/".$org."/genome";
	warn "; Auto selection of output dir\t$dir{output}\n" if ($verbose >= 1);
    }
    &CheckOutputDir();
    
    
    ################################################################
    #### find genbank files in the input directory
    $wd = `pwd`; #### remember working directory
    chomp($wd);
    chdir ($dir{input});
    push @genbank_files, glob("*.gbk");
    push @genbank_files, glob("*.gbk.gz");
    if ($#genbank_files < 0) {
	&FatalError("There is no genbank file in the input directory $dir{input}\n");
    } else {
	warn "; Genbank files\n;\t", join("\n;\t", @genbank_files), "\n" if ($verbose >= 1);
    }
    chdir($wd);     #### come back to the starting directory
    

    #### verbose ####
    &Verbose if ($verbose);
    
    ################################################################
    #### parse the genbank files
    foreach my $file (@genbank_files) {
	$file{input} = "$dir{input}/$file";

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
			  $CDSs,
			  $contigs, 
			  $organisms, 
			  source=>"Genbank", 
#			  source=>$contig, 
			  seq_dir=>$dir{output});
    }

    #### write the chromosome file
    $chrom = &OpenOutputFile("$dir{output}/Contigs_${org}.txt"); # file with chromosome IDs
    foreach my $contig ($contigs->get_objects()) {
	print $chrom join ("\t", 
			   $contig->get_attribute("file"),
			   $contig->get_attribute("id"),
			   $contig->get_attribute("form")), "\n";
    }
    close $chrom;
    
    #### parse feature positions
#    die "HELLO\t$ENV{PARSER}\t$0\t@INC\n";   

    #### parse chromosomal poitions
    &ParsePositions($features);
#    &ParsePositions($genes);
#    &ParsePositions($mRNAs);
#    &ParsePositions($CDSs);

    #### check feature names, description, cross-references
#    &CheckGenbankFeatures($features);
    &CheckGenbankFeatures($features, $genes, $mRNAs, $CDSs);
#    &GuessSynonyms($features);


    ################################################################
    ### export result in various formats
    chdir $dir{output};
    &PrintStats($out_file{stats}, @classes);

    @class_factories = qw (organisms contigs features genes mRNAs CDSs);

    foreach my $factory_name (@class_factories) {
	my $class_factory = $$factory_name;
	warn "; Dumping class $factory_name $class_factory\n" if ($verbose >= 1);
	if ($factory_name eq "features") {
	    $suffix = "_$org";
	} else {
	    $suffix = "";
	}
	$class_factory->dump_tables($suffix);
	$class_factory->generate_sql(schema=>$schema,
				     dir=>"$dir{output}/sql_scripts",
				     prefix=>"$class_factory_",
				     dbms=>$dbms
				     );
    }
    &ExportClasses($out_file{features}, $out_format, @classes) if $export{obj};



#      #### print result
#      chdir $dir{output};
#      #&PrintStats($out_file{stats}, @classes);
#      $features->dump_tables("_$org");
#      $genes->dump_tables("_$org");
#      $mRNAs->dump_tables("_$org");
#      $CDSs->dump_tables("_$org");
#      $contigs->dump_tables("_$org");
#      $organisms->dump_tables("_$org");

#      $features->generate_sql();
#      $genes->generate_sql("_$org");
#      $mRNAs->generate_sql("_$org");
#      $CDSs->generate_sql("_$org");
#      $contigs->generate_sql("_$org");
#      $organisms->generate_sql("_$org");

#      #&ExportClasses($out_file{features}, $out_format, Genbank::Feature) if $export{obj};
    
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
End_short_help
  close HELP;
  exit;
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
	    &PrintHelp;
	    
	    ### list of options
	} elsif ($ARGV[0] eq "-help") {
	    &PrintOptions;
	    
	    ### input file ###
	} elsif ($ARGV[$a] eq "-i") {
	    $dir{input} = $ARGV[$a+1];

	    ### output file ###
	} elsif ($ARGV[$a] eq "-o") {
	    $dir{output} = $ARGV[$a+1];

	    ### quick test
	} elsif ($ARGV[$a] eq "-test") {
	    $test = 1;
	    
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


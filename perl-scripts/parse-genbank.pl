#!/usr/bin/perl 
############################################################
#
# $Id: parse-genbank.pl,v 1.3 2002/05/10 21:03:28 jvanheld Exp $
#
# Time-stamp: <2002-05-10 23:03:14 jvanheld>
#
############################################################
#use strict;;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}
require "RSA.lib";
push @INC, $ENV{PARSER};
require "PFBP_classes.pl";
require "PFBP_parsing_util.pl";


package PFBP::Feature;
{
  @ISA = qw ( PFBP::BiochemicalEntity );
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
			     type=>"SCALAR",
			     description=>"SCALAR",
			     position=>"SCALAR",
			     chromosome=>"SCALAR",
			     strand=>"SCALAR",
			     start=>"SCALAR",
			     end=>"SCALAR",
			     source=>"SCALAR",
			     xrefs=>"EXPANDED"
			     );
}

package main;
{
    #### initialise parameters ####
    my $start_time = &AlphaDate;
    
    local %infile = ();
    local %outfile = ();
    
    local $verbose = 0;
    local $in = STDIN;
    local $out = STDOUT;
    
    my $features = PFBP::ClassFactory->new_class(object_type=>"PFBP::Feature",
						 prefix=>"ft_");
    my $contigs = PFBP::ClassFactory->new_class(object_type=>"PFBP::Contig",
					     prefix=>"ctg_");
    $features->set_out_fields(qw( id type name contig start end strand description position names db_xref introns exons EC_number));
    
    #### working directory
    $wd = `pwd`;
    chomp($wd);
    
    &ReadArguments();
    
    
    #### check argument values ####
    
    #### input directory
    unless (defined($dir{input})) {
	&FatalError("You must specify the input directory.\n");
    }
    unless (-d $dir{input}) {
	&FatalError("Input directory '$dir{input}' does not exist.\n");
    }
    unless (defined($org)) {
	$org = `basename $dir{input}`;
	chomp($org);
	warn "; Auto selection of organism name\t$org\n" if ($verbose >= 1);
    }
    
    #### find genbank files in the input directory
    chdir ($dir{input});
    push @genbank_files, glob("*.gbk");
    push @genbank_files, glob("*.gbk.gz");
    if ($#genbank_files < 0) {
	&FatalError("There is no genbank file in the input directory $dir{input}\n");
    } else {
	warn "; Genbank files\n;\t", join("\n;\t", @genbank_files), "\n" if ($verbose >= 1);
    }
    
    #### come back to the starting directory
    chdir($wd);
    
    #### output directory
    unless (defined($dir{output})) {
	$dir{output} = "$RSA/data/$org/genome";
	warn "; Auto selection of output dir\t$dir{output}\n" if ($verbose >= 1);
    }
    unless (-d $dir{output}) {
	#### create output directory if necessary
	warn "; Creating output dir $dir{output}\n" if ($verbose >= 1);
	system "mkdir -p $dir{output}";
	unless (-d $dir{output}) {
	    &FatalError("Could not create output directory $dir{output}\n");
	}
    }
    
    #### verbose ####
    &Verbose if ($verbose);
    
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
			  $contigs, 
			  source=>$contig, 
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
    &ParsePositions($features);
    
    #### check some feature attributes (name, description, ...)
    foreach $feature ($features->get_objects()) {
	foreach my $name ($feature->get_attribute("gene")) {
	    $feature->push_attribute("names",$name);
	};
	#$feature->push_attribute("names", $feature->get_attribute("gene"));
	
	### define a single name  (take the first value in the name list)
	if ($name = $feature->get_name()) {
	    $feature->set_attribute("name",$name);
	} else {
	    $feature->set_attribute("name",$feature->get_id());
	}
	
	#### check for features without description
	if (($feature->get_attribute("description") eq $null) 
	    || ($feature->get_attribute("description") eq "")) {
	    $feature->set_attribute("description",
				    join("; ", $feature->get_attribute("product"),
					 $feature->get_attribute("note")));
	}
	
	#### use GI as feature identifier
	my @xrefs = $feature->get_attribute("db_xref");
	my $gi = "";
	foreach my $xref (@xrefs) {
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
    
    #### print result
    chdir $dir{output};
    #&PrintStats($out_file{stats}, @classes);
    $features->dump_tables("_$org");
    $contigs->dump_tables("_$org");

#    $features->generate_sql();
#    $contigs->generate_sql();
    #&ExportClasses($out_file{features}, $out_format, PFBP::Feature) if $export{obj};
    
    ###### verbose ######
    if ($verbose) {
	my $done_time = &AlphaDate;
	print $out "; Job started $start_time\n";
	print $out "; Job done    $done_time\n";
    }
    
    
    ###### close output file ######
    close $out if ($outfile{output});
    
    
    exit(0);
}

########################## subroutine definition ############################

sub PrintHelp {
#### display full help message #####
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
		ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes

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

sub PrintOptions {
#### display short help message #####
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


sub ReadArguments {
#### read arguments ####
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

	    ### output file ###
	} elsif ($ARGV[$a] eq "-test") {
	    $test = 1;
	    
	}
    }
}

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


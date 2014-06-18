#!/usr/bin/perl

#############################################################
# $Id: parse-genbank.pl,v 1.69 2012/09/23 05:34:51 jvanheld Exp $
#
# Time-stamp: <2003-10-01 16:17:10 jvanheld>
#
############################################################
#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}


require "RSA.lib";
push @INC, "$ENV{RSAT}/perl-scripts/parsers/";
require "lib/load_classes.pl";
require "lib/util.pl";
require "lib/parsing_util.pl";
require "lib/parse_genbank_lib.pl";
require "classes/Genbank_classes.pl";


################################################################
#### main package
package main;
{
    #### initialise parameters ####
    local $start_time = &RSAT::util::StartScript();
    local %infile = ();
    local %outfile = ();
    local $verbose = 0;
    local $in = STDIN;
    local $out = STDOUT;
    local $single_name = 1;

    ################################################################
    #### initialization
    $null = "<NULL_VALUE>";
    $data_source = "NCBI";
    $ext = "gbk";
    $data_type = "gbk";
    $no_suffix=1;
    $host= $default{'host'};
    $schema="rsat";
    $user="rsat";
    $password="rsat";
    $full_path = 0;
    $test = 0;
    $test_files = 2; ## Maximal number of genbank files to parse for a given organism (there is generally one contig per chromosome)
    $test_lines = 10000; ## macimal number of lines to parse per file

    %supported_prefid_key = (cds=>1,
			     mrna=>1,
			     gene=>1,
			     trna=>1,
			     srna=>1,
			     rrna=>1,
			     misc_rna=>1,
			     scrna=>1,
	);
    $supported_prefid_keys = join(",", keys(%supported_prefid_key));

    %supported_prefid_value = (
	GeneID=>1,
	transcript_id=>1,
	protein_id=>1,
	locus_tag=>1,
	gene=>1, ## Obsolete
	GI=>1, ## Obsolete
	);
    $supported_prefid_values = join(",", keys(%supported_prefid_value));

    ## Default IDs for different feature types
    %preferred_id = (
		     cds=>protein_id,
		     mrna=>transcript_id,
		     gene=>GeneID,
		     trna=>locus_tag,
		     srna=>locus_tag,
		     rrna=>locus_tag,
		     misc_rna=>locus_tag,
		     scrna=>locus_tag,
		    );

    %org_specific_preferred_id = ('Arabidopsis_thaliana'=>'logus_tag');

    ## Initialize class factories for each feature type
    $features = classes::ClassFactory->new_class(object_type=>"Genbank::Feature", prefix=>"ft_");
    $genes = classes::ClassFactory->new_class(object_type=>"Genbank::Gene", prefix=>"gn_");
    $mRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::mRNA", prefix=>"mRNA_");
    $scRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::scRNA", prefix=>"scRNA_");
    $tRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::tRNA", prefix=>"tRNA_");
    $rRNAs = classes::ClassFactory->new_class(object_type=>"Genbank::rRNA", prefix=>"rRNA");
    $repeat_regions = classes::ClassFactory->new_class(object_type=>"Genbank::repeat_region", prefix=>"rep_");
    $misc_RNAs = classes::ClassFactory->new_class(object_type=>"Genbank::misc_RNA", prefix=>"misc_RNA");
    $misc_features = classes::ClassFactory->new_class(object_type=>"Genbank::misc_feature", prefix=>"misc_feature");
    $CDSs = classes::ClassFactory->new_class(object_type=>"Genbank::CDS", prefix=>"cds_");
    $sources = classes::ClassFactory->new_class(object_type=>"Genbank::Source", prefix=>"src_");
    $contigs = classes::ClassFactory->new_class(object_type=>"Genbank::Contig", prefix=>"ctg_");
    $organisms = classes::ClassFactory->new_class(object_type=>"Genbank::Organism", prefix=>"org_");

    @classes = qw( Genbank::Feature
		   Genbank::Contig
		   Genbank::Organism
		   Genbank::Gene
		   Genbank::CDS
		   Genbank::mRNA
		   Genbank::scRNA
		   Genbank::tRNA
		   Genbank::rRNA
		   Genbank::repeat_region
		   Genbank::misc_RNA
		   Genbank::misc_feature
		   Genbank::Source);

    ## Default output fields
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
				 gene_id
				 GeneID

				 names
				 db_xref
				 introns
				 exons
				 EC_number
				 ));
    $organisms->set_out_fields(qw(id taxonomy source));


    ################################################################
    ## Read arguments of the command line
    &ReadArguments();

    ################################################################
    #### Check argument values ####

    #### input directory
    unless (($inputfiles) || (defined($dir{input}))) {
	&RSAT::error::FatalError("You must specify an input directory or input file(s).\n");
    }
    if ($dir{input}) {
	unless (-d $dir{input}) {
	    &RSAT::error::FatalError("Input directory '$dir{input}' does not exist.\n");
	}
    }

    #### organism name
    unless ($org) {
	$org = `basename $dir{input}`;
	chomp($org);
	warn "; Auto selection of organism name\t$org\n" if ($verbose >= 2);
    }

    ## treat organism-specific preferred IDs
    if (defined($org_specific_preferred_id{$org})) {
      %preferred_id = (
		       cds=>$org_specific_preferred_id{$org},
		       mrna=>$org_specific_preferred_id{$org},
		       gene=>GeneID,
		       trna=>$org_specific_preferred_id{$org},
		       srna=>$org_specific_preferred_id{$org},
		       rrna=>$org_specific_preferred_id{$org},
		       misc_rna=>$org_specific_preferred_id{$org},
		       scrna=>$org_specific_preferred_id{$org},
		      );
    }


    ### initial directory
    $dir{main} = `pwd`; #### remember working directory
    chomp($dir{main});

    ################################################################
    #### find genbank files in the input directory
    unless ($inputfiles) {
#	chdir ($dir{input});
	push @genbank_files, glob($dir{input}."/*.${ext}");
	push @genbank_files, glob($dir{input}."/*.${ext}.gz");
	if (scalar(@genbank_files) > 1) {
	    &RSAT::message::Info(scalar(@genbank_files)." Genbank files found in the directory");
	}

	## If no genbank files were found in the main directory (Bacteria),
	## search in CHR_* directories (eukaryotes)
	if (scalar(@genbank_files) ==0 ) {
	    push @genbank_files, glob($dir{input}."/CHR*/*.${ext}");
	    push @genbank_files, glob($dir{input}."/CHR*/*.${ext}.gz");
	    if (scalar(@genbank_files) > 1) {
		&RSAT::message::Info(scalar(@genbank_files)." Genbank files found in the CHR_* subdirectories");
	    }
	}

	## Remove alternative haplotypes (e.g. in Apis_mellifera).
	my @noalt = ();
	foreach my $file (@genbank_files) {
	    if ($file =~ /_alt_/) {
		&RSAT::message::Warning("Ignoring alternative haplotype file", $file);
	    } else {
		push @noalt, $file;
	    }
	}
	@genbank_files = @noalt;
    }

    if ($#genbank_files < 0) {
	system "ls -l";
	&RSAT::error::FatalError("There is no genbank file in the input directory $dir{input}\n");
    } else {
	&RSAT::message::Info(scalar(@genbank_files)." Genbank files") if ($verbose >= 1);
	&RSAT::message::Info("Genbank file names\n;\t", join("\n;\t", @genbank_files)) if ($verbose >= 2);
    }

    my @filtered_genbank_files = grep (!/\/CHR_Un\//, @genbank_files);
    if (scalar(@filtered_genbank_files) < scalar(@genbank_files)) {
      @genbank_files = @filtered_genbank_files;
      &RSAT::message::Warning("Filtering out CHR_Un folder", scalar(@genbank_files). " remaining files");
    }

    #### come back to the starting directory
    chdir($dir{main});

    #### output directory
    unless (defined($dir{output})) {
	if (($dir{input} =~ /refseq/) ||
	    ($data_type eq "refseq")) {
	    $dir{output} = $parsed_data."/refseq/".$ext."/".$delivery_date;
	    warn "; Auto selection of output dir\t$dir{output}\n" if ($verbose >= 2);
	} else {
	    $dir{output} = "$ENV{RSAT}/public_html/data/genomes/".$org."/genome";
	    warn "; Auto selection of output dir\t$dir{output}\n" if ($verbose >= 2);
	}
    }
    chdir $dir{main};
    &CheckOutputDir($dir{output});
    $out_file{features} = "$dir{output}/genbank.obj" if ($export{obj});
    $out_file{error} = "$dir{output}/genbank.errors.txt";
    $out_file{stats} = "$dir{output}/genbank.stats.txt";

    ### Sequence directory
    unless ($noseq) {
	$dir{sequences} = $dir{output};
    }

    ### open error report file
    open ERR, ">$out_file{error}"
	|| die "Error: cannot write error file $out_file{error}\n";

    #### verbose ####
    &Verbose() if ($verbose);

    ## Restrict the number of genbank file for quick tests
    if ($test) {
      @genbank_files = @genbank_files[0..(${test_files})];
    }

    ## Parse the genbank files
#    chdir $dir{input};

    &ParseAllGenbankFiles(@genbank_files);

    ## Export masked sequences
    my @repeats = $repeat_regions->get_objects();
    if (scalar(@repeats) > 1) {
	&ExportMaskedSequences() unless ($noseq);
    }

    #### write the contig file
    chdir $dir{main};
    $chrom = &OpenOutputFile("$dir{output}/contigs.txt"); # file with contig IDs
    foreach my $contig ($contigs->get_objects()) {
	print $chrom join ("\t",
			   $contig->get_attribute("file"),
			   $contig->get_attribute("id"),
			   $contig->get_attribute("form")), "\n";
    }
    close $chrom;

#    $features->index_names();
#    $genes->index_names();
#    $mRNAs->index_names();
#    $scRNAs->index_names();
#    $tRNAs->index_names();
#    $rRNAs->index_names();
#    $misc_RNAs->index_names();
#    $misc_features->index_names();
#    $CDSs->index_names();

    ################################################################
    ### export result in various formats
#    chdir $dir{output};

    #### print parsing statistics
    &PrintStats($out_file{stats}, @classes);

    @class_factories = qw (
			   organisms
			   features
			   contigs
			   genes
			   mRNAs
			   scRNAs
			   tRNAs
			   rRNAs
			   repeat_regions
			   misc_RNAs
			   misc_features
			   CDSs
			   sources
			   );

    foreach my $factory_name (@class_factories) {
	my $class_factory = $$factory_name;
	&RSAT::message::TimeWarn("Dumping class $factory_name")
	  if ($verbose >= 1);
	$suffix = "_$org" unless ($no_suffix);
	$class_factory->dump_tables($suffix);
	$class_factory->generate_sql(dir=>"$dir{output}/sql_scripts",
				     prefix=>"$class_factory_",
				     schema=>$schema,
				     host=>$host,
				     dbms=>$dbms,
				     user=>$user,
				     password=>$password,
				     full_path=>$full_path
				     );
    }
    &ExportMakefile(@classes);
    &ExportClasses($out_file{features}, $out_format, @classes) if $export{obj};

    ## Export protein sequences
    &ExportProteinSequences($CDSs,$org);

    ### Report the output directory
    &RSAT::message::Info(join("\t", "Output directory", $dir{output}));

    ###### close output file ######
    my $exec_time = &RSAT::util::ReportExecutionTime($start_time);
    print $main::out $exec_time if ($main::verbose >= 1);
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
	parse-genbank.pl

        2001 by Jacques van Helden (jvanheld\@bigre.ulb.ac.be)

USAGE
        parse-genbank.pl [-dir input_dir][-i inputfile] [-o outputfile] [-v]

DESCRIPTION
	Parse one or sveral Genbank files for extracting genome
	information.

	Genbank genomes can be retrieved by anonymous ftp :
		ftp://ftp.ncbi.nlm.nih.gov/genomes

	Genome sequence and annotations are parsed from the genbank
	flat files (extension .gbk or .gb).

	Each directory contains the genome of one organism.  Note that

	  - a single directoy can contain several files, if the
	    organism has several contigs.

          - a single .gbk file can contain several contigs, if the
	    genome is not fully assembled for example.

	Parsed data is exported in tab-delimited format, according to
	the normalization rules for relational databases : one main
	table regroups all the single-value attributes(one attribute
	per column), and each multi-value attribute comes in a
	separate table with the foreign key in the first column.

	The program also exports some parsing statistics, an error
	log, and a set of SQL scripts for creating a relational
	database (supported standards: mysql, postgresql, oracle).

CATEGORY
	parser

OPTIONS
	-h	(must be first argument) display full help message

	-help	(must be first argument) display options

	-v	verbose

	-f      input file (mutually exclusive with -i)
		Specify one file to be parsed.
		Can be used iteratively to specify several files.
		 -f file1 -f file2 -f ...

	-i	input directory
		input directory. This directory must contain one or
		several genbank files (extension .gbk by default).

	-source	data source (default: $data_source).

	-ext    extension to be found in the directory specified
	        with the option -i. Default: $ext

	-org    organism name (you should replace spaces by
		underscores to avoid problems)

		This name is used for creating the organism path in
		the export directory.

		If not provided, the basename of the genbank directory
		is used as organism name. It can be convenient to
		specify it manually, when the Genbank directory is not
		the complete organism name (example: H_sapiens instead
		of Homo_sapiens).

	-o	output directory
		The parsing result will be saved in this directory. If
		the directory does not exist, it will be created.

	-test #	quick test (for debugging): only parse the # first
		lines of each Genabnk file (default $test_lines).

	-refseq	input files are refseq entries

	-noseq  do not export sequences in .raw files

	-prefid feattype idname

		Specify the preferred ID for a given feature type
		(feattype). In the original NCBI files, each feature
		(gene, CDS, mRNA,tRNA, ...) is specified by a unique
		identifier (GI), which is in principle reserved for
		internal use at NCBI. However, other databases use
		distinct IDs.

		For example,
		     gene 	   GeneID
		     gene	   locus_tag
		     CDS	   protein_id
		     CDS	   locus_tag
		     mRNA	   transcript_id
		     mRNA	   locus_tag

		In order to cross-link RSAT to other databases, it can
		be useful to impose one or another ID for a given
		feature type.

		A typical use of this option is for genome communities
		having a strong usage of the locus_tag, for historical
		reasons (yeast, Arabidopsis). For the yeast
		Saccharomyces cerevisiae, we use the following option:

		  parse-genbank -i [...] -prefid cds locus_tag


		Note that GeneID and locus_tag are attributes of CDS
		and mRNA, but they are merely cross-references to
		their gene. It is not a good idea to use gene IDS to
		specify CDSs or mRNA, because this raises a confusion
		for genes associated with multiple transcripts and CDS
		(e.g. alternative splicing).

		Examples of utilization
			 -prefid cds protein_id
			 -prefid cds transcript_id
			 -prefid cds locus_tag

   Options for the automaticaly generated SQL scripts
	-schema database schema (default: $schema)
	-host	database host (efault: $host)
	-user	database user (efault: $user)
	-password
		database password (default: $password)
	-fullpath
		Specify the full table path
		(e.g. ${schema}.cds_names) in the SQL scripts, instead
		of the simple table name (e.g. cds_names).

End_of_help
  close HELP;
  exit;
}


################################################################
#### display short help message #####
sub PrintOptions {
  open HELP, "| more";
  print HELP <<End_short_help;
parse-genbank.pl options
----------------
-h		(must be first argument) display full help message
-help		(must be first argument) display options
-f		input file
-i		input dir
-source		data source (default: $data_source).
-ext    	extension of the input files (default: $ext).
-org		organism name (you should replace spaces by underscores)
-refseq		input files are refseq entries
-noseq  	do not export sequences in .raw files
-o		output dir
-v		verbose
-test #		quick test (for debugging)
-prefid		Specify the preferred ID for a given feature type
-schema 	database schema (default: $schema)
-host		database host (default: $host)
-user		database user (default: $user)
-password	database password (default: $password)
-fullpath	Specify the full table path in the SQL scripts.
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

	    ### input file
	} elsif ($ARGV[$a] eq "-f") {
	    $inputfiles++;
	    push @genbank_files, $ARGV[$a+1];

	    ### input directory
	} elsif ($ARGV[$a] eq "-i") {
	    &RSAT::error::FatalError("Option -i is incompatible with option -f") if ($inputfiles);
	    $dir{input} = $ARGV[$a+1];

	    ### extension to be searched in the input directory
	} elsif ($ARGV[$a] eq "-ext") {
	    $ext = $ARGV[$a+1];

	    ### extension to be searched in the input directory
	} elsif ($ARGV[$a] eq "-source") {
	    $data_source = $ARGV[$a+1];

	    ### organism name
	} elsif ($ARGV[$a] eq "-org") {
	    $org = $ARGV[$a+1];

	    ### refseq
	} elsif ($ARGV[$a] eq "-refseq") {
	    $data_type = "refseq";

	    ### do not export sequences
	} elsif ($ARGV[$a] eq "-noseq") {
	    $noseq = 1;

	    ### output file ###
	} elsif ($ARGV[$a] eq "-o") {
	    $dir{output} = $ARGV[$a+1];

	    ### Preferred ID
	} elsif ($ARGV[$a] eq "-prefid") {
	    my $key = lc($ARGV[$a+1]);
	    &RSAT::error::FatalError($key, " is not a valid key for preferred ID. Supported: $supported_prefid_keys")
		unless ($supported_prefid_key{$key});
	    my $value = $ARGV[$a+2];
	    &RSAT::error::FatalError($value, " is not a valid value for preferred ID. Supported: $supported_prefid_values")
		unless ($supported_prefid_value{$value});
	    $preferred_id{$key} = $value;

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

	    ### full path
	} elsif ($ARGV[$a] eq "-fullpath") {
	    $full_path  = 1;

	}
    }
}


################################################################
#### report parameters
sub Verbose {
    print $out "; parse-genbank.pl ";
    &PrintArguments($out);
    if (%dir) {
	print $out "; Directories\n";
	while (($key,$value) = each %dir) {
	    print $out ";\t$key\t$value\n";
	}
    }
    if (%main::infile) {
	print $out "; Input files\n";
	while (($key,$value) = each %infile) {
	    print $out ";\t$key\t$value\n";
	}
    }
    if (%main::outfile) {
	print $out "; Output files\n";
	while (($key,$value) = each %outfile) {
	    print $out ";\t$key\t$value\n";
	}
    }
    printf $out "; %-29s\t%s\n", "organism", $org;
    printf $out "; %-29s\t%d\n", "genbank files", scalar (@genbank_files);
    printf $out "; %-29s\t%s\n", "genbank files", join (" ", @genbank_files) if ($main::verbose >= 2);
}

################################################################
## Specific postprocessing for Refseq data
sub RefseqPostProcessing {
    warn "; Postprocessing proteins parsed from Refseq\n" if ($main::verbose >= 2);
    foreach my $protein ($contigs->get_objects()) {
	my $id = $protein->get_attribute("id");
	my $version = $protein->get_attribute("version");
	my $comment = $protein->get_attribute("comment");
	my ($version_id,$gi) = split /\s+/, $version;
	my ($comment_1, $comment_2) = split /The reference sequence was derived from /, $comment;
	my ($embl_xref, $comment_3) = split /\./, $comment_2;
	if ($gi =~ /^GI:/) {
	    $gi = "$'";
	} else {
	    &ErrorMessage(join ("\t", "Protein", $id, "version", $version, "Invalid GI"));
	}

	$protein->set_attribute("GI", $gi);
	$protein->set_attribute("version_id", $version_id);
	$protein->set_attribute("EMBl_XRef", $embl_xref);
	warn join ("\t", $protein->get_attribute('id'),$version,
		  $version_id,
		  $gi, $embl_xref), "\n" if ($main::verbose >= 3);
    }

}


################################################################
## Export protein sequences in fasta format
sub ExportProteinSequences {
    my ($CDSs, $org) = @_;
    my $separator="; ";
    $out_file{pp} = $dir{output}."/".$org."_aa.fasta";

    &RSAT::message::TimeWarn("Exporting translated sequences to file", $out_file{pp})
	if ($main::verbose >= 2);

    open PP, ">$out_file{pp}";
    foreach my $cds ($CDSs->get_objects()) {
	next unless ($cds);
	my ($translation) = $cds->get_attribute("translation");
	next unless ($translation =~ /\S+/);
	my $id = $cds->get_attribute("id");
	my $gene = $cds->get_attribute("gene");
	if (!($gene) || ($gene eq $null)) {
	    $gene = $id;
	}
	my $pp_id = $id;

        ## Get CDS description
        my $description = join($separator, $org,$id,$gene);
	$description .= $separator;
	if ($cds->get_attribute("description")) {
          $description .= $cds->get_attribute("description");
	} else {
          $description .= join ("; ", $org, $gene, $cds->get_attribute("note"));
       }

        my $pp_description;
        if ($description) {
          $pp_description .= $description;
        }
        $pp_description .= "; ".join ("|", $cds->get_attribute("names"));

        print PP $header, "\n";
#        &PrintNextSequence(PP,"fasta",60,$translation,$pp_id, $pp_description);
        &PrintNextSequence(PP,"fasta",60,$translation,$pp_id);
    }
    close PP;
}


################################################################
## Export masked sequences
sub ExportMaskedSequences {
    chdir $dir{sequences};
    foreach my $contig ($contigs->get_objects()) {
	my $file =    $contig->get_attribute("file");
	my $contig_id = $contig->get_attribute("id");
	my $sequence = `cat $file`;
	&RSAT::message::TimeWarn(join("\t", "Masking sequence",
				      "contig",
				      $contig_id,
				      "length", $contig_len)) if ($main::verbose >= 2);

	foreach my $region ($repeat_regions->get_objects) {
	    my $contig_len = length($sequence);
	    my $region_contig = $region->get_attribute("contig");
	    next unless ($region_contig eq $contig_id);
	    my $start = $region->get_attribute("start_pos");
	    $start =~ s/\<//g;
	    $start =~ s/\>//g;
	    my $offset = $start -1;
	    my $end = $region->get_attribute("end_pos");
	    $end =~ s/\<//g;
	    $end =~ s/\>//g;
	    next unless ((&IsReal($start)) && (&IsReal($end)));
	    my $len = $end - $start + 1;

#	    my $before = substr($sequence, $offset - 3, $len + 6);
#	    my $after = substr($sequence, $offset - 3, $len + 6);
# 	    &RSAT::message::Debug("Masking region",
# 				  $region_contig,
# 				  "contig len=".$contig_len,
# 				  "start=".$start,
# 				  "end=".$end,
# 				  "offset=".$offset,
# 				  "len=".$len,
# ##				  $before,
# ##				  $after,
# 				  ) if ($main::verbose >= 20);

	    substr($sequence, $offset, $len) = "n"x$len;
	}
	my $masked_sequence_file = $contig_id."_repeat_masked.raw";
	my ($masked, $dir) = &OpenOutputFile($masked_sequence_file);
	print $masked $sequence;
	close $masked;
	&RSAT::message::TimeWarn(join("\t", "Exported masked sequence",
				      "contig",
				      $contig_id,
				      "dir",
				      $dir{sequence},
				      "file",
				      $masked_sequence_file)
				 ) if ($main::verbose >= 2);

    }
}

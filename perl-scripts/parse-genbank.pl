#!/usr/bin/perl 
############################################################
#
# $Id: parse-genbank.pl,v 1.2 2002/01/07 00:30:04 jvanheld Exp $
#
# Time-stamp: <2002-01-07 01:26:55 jvanheld>
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

#### initialise parameters ####
my $start_time = &AlphaDate;

local %infile = ();
local %outfile = ();

local $verbose = 0;
local $in = STDIN;
local $out = STDOUT;

my $genes = PFBP::ClassFactory->new_class(object_type=>"PFBP::Gene",
					  prefix=>"gene_");
$genes->set_out_fields(qw( id type name chromosome start end strand description position names db_xref introns exons ));

#### working directory
$wd = `pwd`;
chomp($wd);

&ReadArguments;


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
@genbank_files = glob("*.gbk");
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
$chrom = &OpenOutputFile("$dir{output}/Chromosomes_${org}.txt"); # file with chromosome IDs
foreach my $file (@genbank_files) {
    my $contig = `basename $file .gbk`;
    chomp($contig);
    my $sequence = &ParseGenbankFile("$dir{input}/$file", $genes, source=>$contig);
    open RAW, ">$dir{output}/${contig}.raw";
    &PrintNextSequence(RAW, "raw", 0, $sequence, $contig);
    close RAW;
    print $chrom "$contig.raw", "\t", $contig, "\n";
}
close $chrom;

#### parse gene positions
&ParsePositions($genes);

#### check some gene attributes (name, description, ...)
foreach $gene ($genes->get_objects()) {
    foreach my $name ($gene->get_attribute("gene")) {
	$gene->push_attribute("names",$name);
    };
    #$gene->push_attribute("names", $gene->get_attribute("gene"));
    
    ### define a single name  (take the first value in the name list)
    if ($name = $gene->get_name()) {
	$gene->set_attribute("name",$name);
    } else {
	$gene->set_attribute("name",$gene->get_id());
    }
    
    #### check for genes without description
    if (($gene->get_attribute("description") eq $null) 
	|| ($gene->get_attribute("description") eq "")) {
	$gene->set_attribute("description",
			     join("; ", $gene->get_attribute("product"),
				  $gene->get_attribute("note")));
    }
    
    #### use GI as gene identifier
    my @xrefs = $gene->get_attribute("db_xref");
    my $gi = "";
    foreach my $xref (@xrefs) {
	if ($xref =~ /GI:/) {
	    $gi = $';
	    last;
	} 
    }
    
      if ($gi) {
  	$gene->force_attribute("id",$gi);
      } else {
  	&ErrorMessage("; Error\tgene ".$gene->get_attribute("id")." has no GI.\n"); 
      }
    
    #### use genbank name as chromosome name
    my $source = $gene->get_attribute("source");
    if ($source =~ /genbank:/) {
	my $chromosome = $';
	$chromosome =~ s/\.gz$//;
	$chromosome =~ s/\.gbk$//;
	$gene->force_attribute("chromosome",$chromosome);
    }
}

### print result
chdir $dir{output};
#&PrintStats($out_file{stats}, @classes);
$genes->dump_tables("_$org");
#&ExportClasses($out_file{genes}, $out_format, PFBP::Gene) if $export{obj};

###### verbose ######
if ($verbose) {
    my $done_time = &AlphaDate;
    print $out "; Job started $start_time\n";
    print $out "; Job done    $done_time\n";
}


###### close output file ######
close $out if ($outfile{output});


exit(0);


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


#  sub ParseGenbankFile {
#      my ($input_file, $class_holder, %args) = @_;
#      warn ";\n; Parsing file $input_file.\n" if ($verbose >= 1);
    
#      $gbk = &OpenInputFile($input_file);
#  #    open GBK, $input_file 
#  #	|| die "Error: cannot open input file $input_file.\n";
#      my $l = 0;
#      my $in_features = 0;
#      my $in_feature = 0;
#      my $in_gene = 0;
#      my $in_cds = 0;
#      my $current_gene = null;
#      my $organism_name = "";
#      while (my $line = <$gbk>) {
#  	$l++;
#  	print STDERR $line if ($verbose >= 10);
#  	chomp $line;
#  	next unless ($line =~ /\S/);
#  	unless ($in_features) {
#  	    if ($line =~ /^\s+ORGANISM\s+/) {
#  		$organism_name = $';
#  		warn "; Organism name\t\t$organism_name\n" if ($verbose >= 1);
#  	    }
#  	    if ($line =~ /^FEATURES/) {
#  		$in_features = 1 ;
#  		warn "; Reading features\n" if ($verbose >= 1);
#  	    }
#  	    next;
#  	}
#  	if ($line =~ /^     CDS\s+(.*)/) {
#  	    $in_gene = 0;
#  	    $in_cds = 1;
#  	    $position = &trim($1);
#  	    $current_gene = $class_holder->new_object(%args);
#  	    $current_gene->set_attribute("type","CDS");
#  	    $current_gene->set_attribute("organism",$organism_name);
#  	    $current_gene->set_attribute("position",$position);
#  	    warn ";\tline $l\tnew CDS\n" if ($verbose >= 2);
#  	} elsif ($line =~ /     gene\s+(.*)/) {
#  	    $in_cds = 0;
#  	    $in_gene  = 1;
#  	    $position = $1;
#  	    &trim($position);
#  	} elsif ($in_cds) {
#  	    unless ($current_gene) {
#  		die "Error: $file $input_file\tline $l\tno gene defined\n.";
#  	    }
#  	    if ($line =~ / +\/(\S+)\=(\d+)/) {
#  		#### numerical value
#  		$key = $1;
#  		$value = $2;
#  		$current_gene->new_attribute_value($key,$value);
		
#  	    } elsif ($line =~ / +\/(\S+)\=\"(.+)\"/) {
#  		#### short string
#  		$key = $1;
#  		$value = $2;
#  		$current_gene->new_attribute_value($key,$value);

#  	    } elsif ($line =~ / +\/(\S+)\=\"(.+)/) {
#  		#### long string
#  		$key = $1;
#  		$value = $2;
#  		$in_feature = 1;
		
#  	    } elsif ($in_feature) {
#  		$to_add = &trim($line);
#  		if ($to_add =~ /\"$/) {
#  		    $value .= " " unless ($key eq "translation");
#  		    $value .= $`;
#  #		    die "HELLO\n$to_add\n$`\n$value\n";
#  		    $current_gene->new_attribute_value($key,$value);
#  		    $key = "";
#  		    $value  = "";
#  		    $in_feature = 0;
#  		} else {
#  		    $value .= " " unless ($key eq "translation");
#  		    $value .= $to_add;
#  		}
#  	    } elsif ($in_gene) {
#  		#### genes are not parsed for the time being (which means that tRNA are not parsed, since there is no corresponding CDS)
#  	    }
#  	} else {
#  	    warn ("file ".$input_file."\tline $l\tnot parsed\t$line\n") if ($verbose >= 2);
#  	}
#      }
#      close $gbk;
	    
#  }


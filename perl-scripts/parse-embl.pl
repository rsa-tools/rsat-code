#!/usr/bin/perl 
############################################################
#
# $Id: parse-embl.pl,v 1.2 2002/06/06 11:23:37 jvanheld Exp $
#
# Time-stamp: <2002-06-06 13:23:32 jvanheld>
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

#### find embl files in the input directory
chdir ($dir{input});
@embl_files = glob("*.contig");
if ($#embl_files < 0) {
    &FatalError("There is no embl file in the input directory $dir{input}\n");
} else {
    warn "; EMBL files\n;\t", join("\n;\t", @embl_files), "\n" if ($verbose >= 1);
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

#### parse the embl files
$chrom = &OpenOutputFile("$dir{output}/Contigs_${org}.txt"); # file with chromosome IDs
foreach my $file (@embl_files) {
    my $contig = `basename $file .contig`;
    chomp($contig);
    my $sequence = &ParseEMBLFile("$dir{input}/$file", $genes, source=>$contig);
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
    
    #### use embl name as chromosome name
    my $source = $gene->get_attribute("source");
    if ($source =~ /embl:/) {
	my $chromosome = $';
	$chromosome =~ s/\.gz$//;
	$chromosome =~ s/\.contig$//;
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
	parse-embl

        2001 by Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)
	
USAGE
        parse-embl [-dir input_dir][-i inputfile] [-o outputfile] [-v]

DESCRIPTION
	Parse one or sveral EMBL files for extracting genome
	information.

CATEGORY
	parser

OPTIONS
	-h	(must be first argument) display full help message
	-help	(must be first argument) display options
	-v	verbose
	-i	input directory
		input directory. This directory must contain one or
		several embl files (extension .contig). 
	-o	output directory
		The parsing result will be saved in this directory. If
		the directory does not exist, it will be created.
	-org	organism name
End_of_help
  close HELP;
  exit;
}

sub PrintOptions {
#### display short help message #####
  open HELP, "| more";
  print HELP <<End_short_help;
parse-embl options
----------------
-h	(must be first argument) display full help message
-help	(must be first argument) display options
-i	input dir
-o	output dir
-v	verbose
-org	organism name
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
	    
	    ### organism ###
	} elsif ($ARGV[$a] eq "-org") {
	    $org = $ARGV[$a+1];
	    $org =~ s/\s+/_/g;
	}
    }
}

sub Verbose {
    print $out "; parse-embl ";
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
    printf $out "; %-29s\t%s\n", "embl files", join (" ", @embl_files);
}


sub ParseEMBLFile {
    my ($input_file, $class_holder, %args) = @_;
    warn ";\n; Parsing file $input_file.\n" if ($verbose >= 1);
    
    $embl = &OpenInputFile($input_file);
#    open EMBL, $input_file 
#	|| die "Error: cannot open input file $input_file.\n";
    my $l = 0;
    my $in_features = 0;
    my $in_feature = 0;
    my $in_gene = 0;
    my $in_cds = 0;
    my $current_gene = null;
    my $organism_name = "";
    my $sequence = "";
    while (my $line = <$embl>) {
	$l++;
	print STDERR $line if ($verbose >= 10);
	chomp $line;
	next unless ($line =~ /\S/);

	#### read the full sequence
	if  ($line =~ /^SQ/) {
	    $in_featuress = 0;
	    $in_sequence = 1;
	    while (my $line = <$embl>) {
		if ($line =~ /\d+\s*$/) {
		    $sequence .= $`;
		} elsif ($line =~ /^\/\/$/) {
		    $in_sequence = 0;
		}
	    }
	    
        }


	unless ($in_features) {
#	    if ($line =~ /^\s+ORGANISM\s+/) {
#		$organism_name = $';
#		warn "; Organism name\t\t$organism_name\n" if ($verbose >= 1);
#	    }
	    if ($line =~ /^FT\s+/) {
		$in_features = 1 ;
		warn "; Reading features\n" if ($verbose >= 1);
	    }
	    next;
	}
	if ($line =~ /^FT   CDS\s+(.*)/) {
	    $in_gene = 0;
	    $in_cds = 1;
	    $position = &trim($1);
	    if ($position =~ /join\(/){
		### check that the position is complete
		my $start_line = $l;
		unless ($' =~ /\)/) {
		    do {
			die "Error: position starting at line $l is not terminated properly.\n"
			    unless $position_suite = <$embl>;
			$position_suite =~ s/^FT//;
			$position .= &trim($position_suite);
		    } until ($position =~ /\)/);
		}
	    }
	    $current_gene = $class_holder->new_object(%args);
	    $current_gene->set_attribute("type","CDS");
	    $current_gene->set_attribute("organism",$organism_name);
	    $current_gene->set_attribute("position",$position);
	    warn ";\tline $l\tnew CDS\n" if ($verbose >= 2);
	} elsif ($line =~ /^FT   \S+/) {
	    $in_cds = 0;
#	} elsif ($line =~ /     gene\s+(.*)/) {
#	    $in_cds = 0;
#	    $in_gene  = 1;
#	    $position = $1;
#	    &trim($position);
	} elsif ($in_cds) {
	    unless ($current_gene) {
		die "Error: $file $input_file\tline $l\tno gene defined\n.";
	    }
	    if ($line =~ /^FT +\/(\S+)\=(\d+)/) {
		#### numerical value
		$key = $1;
		$value = $2;
		$current_gene->new_attribute_value($key,$value);
		
	    } elsif ($line =~ /^FT +\/(\S+)\=\"(.+)\"/) {
		#### short string
		$key = $1;
		$value = $2;
		$current_gene->new_attribute_value($key,$value);

	    } elsif ($line =~ /^FT +\/(\S+)\=\"(.+)/) {
		#### long string
		$key = $1;
		$value = $2;
		$in_feature = 1;
		
	    } elsif ($in_feature) {
		if ($line =~ /^FT +/) {
		    $to_add = &trim($');
		    if ($to_add =~ /\"$/) {
			$value .= " " unless ($key eq "translation");
			$value .= $`;
#		    die "HELLO\n$to_add\n$`\n$value\n";
			$current_gene->new_attribute_value($key,$value);
			$key = "";
			$value  = "";
			$in_feature = 0;
		    } else {
			$value .= " " unless ($key eq "translation");
			$value .= $to_add;
		    }
		} else {
		    warn ("file ".$input_file."\tline $l\tnot parsed\t$line\n") if ($verbose >= 2);
		}
	    } elsif ($in_gene) {
		#### genes are not parsed for the time being (which means that tRNA are not parsed, since there is no corresponding CDS)
	    }
	} else {
	    warn ("file ".$input_file."\tline $l\tnot parsed\t$line\n") if ($verbose >= 2);
	}
    }
    close $embl;
    return $sequence;
}


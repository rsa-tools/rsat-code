#!/usr/bin/perl
############################################################
#
# $Id: parse_swissprot.pl,v 1.10 2002/03/19 11:43:38 jvanheld Exp $
#
# Time-stamp: <2002-03-19 12:41:45 jvanheld>
#
############################################################

### add the program's directory to the lib path
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	$script_dir = $` || ".";
	push @INC, "$script_dir";
	push @INC, "$script_dir/perllib/";
    }
}
require "PFBP_config.pl";
require "PFBP_classes.pl";
require "PFBP_util.pl";
require "PFBP_parsing_util.pl";
require "PFBP_loading_util.pl";
use SWISS::Entry;

package main;
{
    
    
    ### organism selection
    @selected_organisms = ();
    $full_name{yeast} = "Saccharomyces cerevisiae (Baker's yeast)";
    $full_name{ecoli} = "Escherichia coli";
    $full_name{human} = "Homo sapiens (Human)";

    $dir{delivery} = "/rubens/dsk3/genomics/delivery/internal/swissprot_parsed";

    ### input directories and files
#    $dir{input} = "/win/databases/downloads/ftp.ebi.ac.uk/pub/databases/";
    $dir{input} = "/win/databases/downloads/ftp.expasy.org/databases/sp_tr_nrdb";
#    $dir{swissprot} = $dir{input};
#    $dir{swissprot} = $dir{input}."/swiss-prot/release_compressed";
#    $dir{trembl} = $dir{input};

    $source{swissprot} = "sprot";
    $source{trembl} = "trembl";
    $source{trembl_new} = "trembl_new";
#    $source{swissprot} = "sprot40";
#    $source{yeast} = "sptrembl/fun";
#    $source{human} = "sptrembl/hum";
#    $source{ecoli} = "sptrembl/pro";
#  &OpenIndex("ECSet");
    
    ### export classes
    push @classes, "PFBP::Polypeptide";
#    @{$out_fields{'PFBP::Polypeptide'}} = qw( id name description gene organisms swissprot_acs swissprot_ids names );
    
    &ReadArguments;
    
    #### output directory
    $out_format = "obj";
    unless (defined($dir{output})) {
	$dir = $parsed_data."/swissprot_parsed";
	unless (-d $dir) {
	    warn "Creating output dir $dir", "\n";
	    mkdir $dir, 0775 || die "Error: cannot create directory $dir\n";
	}
	$dir{output} = $dir."/".$delivery_date;
    }
    unless (-d $dir{output}) {
	warn "Creating output dir $dir{output}\n";
	mkdir $dir{output}, 0775 || die "Error: cannot create directory $dir\n";
    }
    die unless chdir $dir{output};
    if ($clean) {
	system "\\rm -f $dir{output}/*";
    }

    #### output file names depend on the organisms to parse and source 
    #### -> defined after reading arguments
#  $suffix = "";
#  foreach $organism (@selected_organisms) {
#      $suffix .= "_$organism";
#  }
#  if ($export{enzymes}) {
#      $suffix .= "_enz";
#  }
#  $suffix .= "_test" if ($test);
    
    $out_file{polypeptides} = $dir{output}."/Polypeptide".$suffix.".obj";
##  $out_file{mldbm_polypeptides} = $dir{output}."/Polypeptide".$suffix.".mldbm";
#  $out_file{mldbm_polypeptide_index} = $dir{output}."/Polypeptide".$suffix."__name_index.mldbm";
#  $out_file{catalyses} = $dir{output}."/Catalysis".$suffix.".obj";
#  $out_file{mldbm_catalyses} = $dir{output}."/Catalysis".$suffix.".mldbm";
    $out_file{errors} = $dir{output}."/swissprot".$suffix.".errors.txt";
    $out_file{stats} = $dir{output}."/swissprot".$suffix.".stats.txt";
    
    ### open error report file
    open ERR, ">$out_file{errors}" || 
	die "Error: cannot write error file $out_file{errors}\n";
    
    ### select three predefined organisms if none was selected (-org)
    unless (($#selected_organisms >= 0) || 
	    ($export{allorg})){
	push @selected_organisms, "ecoli";
	push @selected_organisms, "human";
	push @selected_organisms, "yeast";
    }

    ### define regular expressions to match the selected organisms
    foreach $organism (@selected_organisms) {
	$regexp = $full_name{$organism};
	$regexp = $` if ($regexp =~ /\(/); 
	$regexp = $` if ($regexp =~ /\)/);
	$regexp = $` if ($regexp =~ /\'/);
	$selected_organism{uc($full_name{$organism})} = $regexp;
	push @regexps, qr /$regexp/; ### list of regular expressions
    }
    

    #### read a list of selected ACCESSION NUMBERS
    if ($in_file{acs}) {
	warn "; Reading Accession Number list from file $in_file{acs}\n" 
	    if ($warn_level >=1);
	unless (-e $in_file{acs}) {
	    die "Accession number file $in_file{acs} does not exist.\n";
	}
	unless (-r $in_file{acs}) {
	    die "Cannot read aSccession number file $in_file{acs}\n";
	}
	open ACS, $in_file{acs} || die 
	    "Error: cannot open file $in_file{acs}\n";
	while (<ACS>) {
	    chomp;
	    my @fields = split /\s+/;
	    my $ac = $fields[0];
	    $selected_acs{$ac}++;
	}
	close ACS;
	warn "; Selected ACs\n;\t",  join ("\n;\t", keys %selected_acs), "\n" if ($warn_level >= 0);
    }

    ### by default, parse swissprot only
    unless (defined(%data_sources)) {
	$data_sources{swissprot} = 1;
    }
#    if ($data_sources{swissprot}) {
#	$in_file{$source{swissprot}} = "gunzip -c ".$dir{swissprot}."/".$source{swissprot}.".dat.gz | ";
#	$files_to_parse{$source{swissprot}} = 1;
#    }
    foreach $db (keys %data_sources) {
	if ($data_sources{$db}) {
	    $in_file{$source{$db}} = "gunzip -c ".$dir{input}."/".$source{$db}.".dat.gz | ";
	    $files_to_parse{$source{$db}} = 1;
	}
	warn $db, "\t", $source{$db}, "\t", $in_file{$source{$db}}, "\n";
    }

    ### create class holders
    $polypeptides = PFBP::ClassFactory->new_class(object_type=>"PFBP::Polypeptide",
						  prefix=>"spp_");
    $polypeptides->set_out_fields( qw( id 
				       source
#				       name 
				       description
				       gene
				       names
				       organisms
				       features
				       comments
				       swissprot_acs
				       swissprot_ids
				       ECs
				       ));
    $polypeptides->set_attribute_header("features", join ("\t", "Feature_key", "from", "to", "description") );
    $polypeptides->set_attribute_header("comments", join ("\t", "topic", "comment") );

    ### testing mode
    if ($test) {
	warn ";TEST\n" if ($warn_level >= 1);
	### fast partial parsing for debugging
	foreach $key (keys %in_file) {
	    next if ($key eq 'acs');
	    $in_file{$key} .= " head -20000 |";
	}
    }


    ### default verbose message
    if ($warn_level >= 1) {
	warn "; Selected organisms\n;\t", join("\n;\t", @selected_organisms), "\n";
	warn "; Polypeptide classes\n;\t", join("\n;\t", keys %export), "\n";
	warn "; Data sources\n;\t", join("\n;\t",  keys %data_sources), "\n";
	&DefaultVerbose;
    }


    ### parse data from original files
    foreach $source (keys %files_to_parse) {
	&ParseSwissprot($in_file{$source}, "SWISSPROT:".$source);
    }
#  &CreateCatalyses;

    #### use first name as primary name
#  foreach $polypeptide ($polypeptides->get_objects()) {
#      if ($name = $polypeptide->get_name()) {
#	  $polypeptide->set_attribute("primary_name",$name);
#      }
#  }


#### print the result
    &PrintStats($out_file{stats}, @classes);
    $polypeptides->dump_tables($suffix);
    #  $catalyses->dump_tables($suffix);
    &ExportClasses($out_file{polypeptides}, $out_format,PFBP::Polypeptide) if ($export{obj});
#  &ExportClasses($out_file{catalyses}, $out_format, PFBP::Catalysis);
#  $catalyses->export('MLDBM',$out_file{mldbm_catalyses});
#  $polypeptides->export("MLDBM",$out_file{mldbm_polypeptides});
#  $polypeptides->export_name_index("MLDBM",$out_file{mldbm_polypeptide_index});


#if ($out_file{stats}) {
#  open STDOUT, ">>$out_file{stats}"
#    || die "Error: cannot write stat file $out_file{stats}\n";
#  &SwissStats;
#  close $out if ($out_file{stats});
#}


### report execution time
    if ($warn_level >= 1) {
	$done_time = &AlphaDate;
	warn (";\n",
	      "; job started $start_time",
	      "; job done    $done_time\n")
	}

#  &CloseIndex("ECSet");
    close ERR;

    system "gzip -f $dir{output}/*.tab $dir{output}/*.txt";
    system "gzip -f $dir{output}/*.obj" if ($export{obj});

    exit(0);
}

### subroutines for the main package

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
	parse_swissprot.pl

DESCRIPTION
	Parse polypeptides from a Swissprot flat file

AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

VERSION
	0.01
	Created		2000/01/11
	Last modified	2000/12/11
	
OPTIONS
	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-indir	input directory. 
		This directory should contain a download of the
		non-redundant swiss-prot/trembl ftp site :
		       ftp://ftp.expasy.org/databases/sp_tr_nrdb
	-outdir	output directory.
		The parsed data will be stored in this directory.
	-w #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-enz	export enzymes only
	-org	select an organism for exportation
		can be used reiteratively in the command line 
		to select several organisms
		   Supported organisms :
			ecoli
			human
			yeast
		by default, these three organisms are selected
	-allorg export all organisms
	-data	database to be parsed.
		   Valid data sources:
			swissprot (default)
			trembl
			trembl_new
		This argument can be used iteratively to parse several
		databases. E.g.
			   -data swissprot -data trembl
	-obj	export data in object format (.obj file)
		which are human-readable (with some patience 
		and a good cup of coffee)
	-acs	accession file
		The parsing will be restricted to the swissprot
		accession numbers (AC field) specified in the
		file. Each accession number must come as the first
		word of a new line (all subsequent words are ignored).
	-clean	remove all files from the output directory before
		parsing
EXAMPLE
	parse_polypeptides.pl -w 2 -org ecoli -data swissprot -enz
EndHelp
  close HELP;
}

### read arguments from the command line
sub ReadArguments {
    for my $a (0..$#ARGV) {
	
	### warn level
	if (($ARGV[$a] eq "-w" ) && 
	    ($ARGV[$a+1] =~ /^\d+$/)){
	    $main::warn_level = $ARGV[$a+1];
	    
	    ### test run
	} elsif ($ARGV[$a] eq "-test") {
	    $main::test = 1;
	    
	    ### clean
	} elsif ($ARGV[$a] eq "-clean") {
	    $main::clean = 1;
	    
	    ### output file
	} elsif ($ARGV[$a] eq "-obj") {
	    $a++;
	    $main::export{obj} = 1;
	    
	    ### help
	} elsif (($ARGV[$a] eq "-h") ||
		 ($ARGV[$a] eq "-help")) {
	    &PrintHelp;
	    exit(0);

	    ### input directory
	} elsif ($ARGV[$a] =~ /^-indir/) {
	    $dir{input} = $ARGV[$a+1];

	    ### input file with a list of Swissprot AC
	} elsif ($ARGV[$a] =~ /^-acs/) {
	    $in_file{acs} = $ARGV[$a+1];

	    ### output directory
	} elsif ($ARGV[$a] =~ /^-outdir/) {
	    $dir{output} = $ARGV[$a+1];

	    ### select organisms for exportation
	} elsif ($ARGV[$a] =~ /^-org/) {
	    push @selected_organisms, lc($ARGV[$a+1]);
	    #### export all organisms
	} elsif ($ARGV[$a] =~ /^-allorg/) {
	    $main::export{allorg} = 1;
	    
	    ### select enzymes for exportation
	} elsif ($ARGV[$a] =~ /^-data/) {
	    $data_sources{lc($ARGV[$a+1])} = 1;
	    

	    ### select enzymes for exportation
	} elsif ($ARGV[$a] =~ /^-enz/) {
	    $main::export{enzymes} = 1;
	}


    }

}

#### this is the routine for actually parsing the swissprot or trembl file
sub ParseSwissprot {
    my ($input_file, $source) = @_;
    $source = $input_file unless ($source);
#    my %attrib_keys = (
#		       id=>"ID",
#		       ac=>"AC",
#		       gn=>"Gene",
#		       ft=>"Features",
#		       cc=>"Comments"
#		       );
#    my $class = PFBP::Polypeptide;
    
    warn (";\n; ", &AlphaDate,  " parsing polypeptides from $input_file\n")
	if ($warn_level >= 1);
    
    open DATA, $input_file || 
	die "Error: cannot open data file $input_file\n";
    
    # Read an entire record at a time
    $/ = "\/\/\n";
    
    my $entries = 0;
    while ($text_entry = <DATA>){
	$entries++;
	
	my $export = 0;
	
	### check that the entry matches organism name 
	### before converting it to an object
	my $parse = 1;
	unless ($export{allorg}) {
	    $parse = 0;
	    foreach $regexp (@regexps) {
		if ($text_entry =~ /$regexp/i) {
		    $parse = 1;
		    last;
		}
	    }
	}
	next unless $parse;
	
	# Read the entry
	my $object_entry = SWISS::Entry->fromText($text_entry);
	
	### check whether the polypeptide organism(s)
	### match the organism selection
	my @organisms = $object_entry->OSs->elements;
	unless ($export{allorg}) {
	    foreach $organism (@organisms) {
		if ($selected_organism{uc($organism)}) {
		    $export = 1;
		    last;
		}
	    }
	}
	
	# get the polypeptide accession number
	my $swissprot_ac = $object_entry->AC;

	#### check whether the accession number was specified for export 
	$export = 1 if ($selected_acs{$swissprot_ac});
	
	next unless $export;
	my @swissprot_ids = $object_entry->IDs->elements;
	my @swissprot_acs = $object_entry->ACs->elements;
	my $descr = $object_entry->DEs->text;
	my $geneNames = $object_entry->GNs->text;
	my @geneNames = split " OR ", $geneNames;
	### extract name from description
	my @names = ();
	if ($descr =~ / \(/) {
	    push @names, lc($`);
	    #if ($descr =~ /^([^( \()]+)/) {
	    #push @names, lc($1);
	} else {
	    push @names, lc($descr);
	}
	

	### extract ECs from description
	my @ECs = ();
	my $tmp = $descr;
	while ($tmp =~ /\(EC ([^\(]+)\)/) {
	    $ec_number = $1;
	    $ec_number =~ s/^ //g;
	    $ec_number =~ s/ $//g;
	    $tmp = $';
	    if ($ec_number =~ /^\S+\.\S+\.\S+\.\S+$/ ) {
		push @ECs, $ec_number;
	    } else {
		&ErrorMessage("Error: $swissprot_ac invalid EC number $ec_number\n");
		if ($ec_number =~ /^(\S+\.\S+\.\S+\.\S+)/ ) { #### try to recuperate th EC number if it is credible
		    push @ECs, $1;
		}
	    }
	}
	if (($export{enzymes})  &&
	    ($#ECs < 0)){
	    next;
	}
	
	### create a new polypeptide
	warn "$source\tentry $entries\t$swissprot_ids[0]\t$names[0]\n"
	    if ($warn_level >= 2); 
	my $polypeptide = $polypeptides->new_object(id=>$swissprot_ac,
						    source=>$source);
	if (defined($geneNames[0])) {
	    $polypeptide->set_attribute("gene",$geneNames[0]);
	} else {
	    $polypeptide->set_attribute("gene",$null);
	}

	my %already_assigned = ();
#	foreach my $name (@names, @geneNames, @swissprot_ids, @swissprot_acs) {
	foreach my $name (@names, @geneNames) {
	    $polypeptide->push_attribute("names",$name) 
		unless $already_assigned{uc($name)};
	    $already_assigned{uc($name)}++;  #### prevent assigning twice the same name
	}

	
	my @features = $object_entry->FTs->elements;
	foreach my $ft (@features) {
	    $polypeptide->push_expanded_attribute("features", @{$ft});
	}
	my @comments = $object_entry->CCs->elements;
	foreach my $cc (@comments) {
	    $polypeptide->push_expanded_attribute("comments", @{$cc});
	}

	foreach my $id (@swissprot_ids) {
#print STDERR "id\t$id\n";
	    $polypeptide->push_attribute("swissprot_ids",$id);
	}
	foreach my $ac (@swissprot_acs) {
	    $polypeptide->push_attribute("swissprot_acs",$ac);
	}
	$polypeptide->set_attribute("description", $descr);

	my $pp_id = $polypeptide->get_id();
	foreach my $ec (@ECs) {
#	    my $catalysis = $catalyses->new_object(source=>      $source,
#						   catalyst=>    $pp_id,
#						   catalyzed=>   $ec);
	    $polypeptide->push_attribute("ECs",$ec);
	}
	foreach my $organism (@organisms) {
	    $polypeptide->push_attribute("organisms",$organism);
	}
	
	
    }
    close DATA;
    
#    $polypeptides->index_names();
    return 1;
}



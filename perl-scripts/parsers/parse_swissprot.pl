#!/usr/bin/perl
############################################################
#
# $Id: parse_swissprot.pl,v 1.37 2011/02/17 05:07:46 rsat Exp $
#
# Time-stamp: <2003-07-10 11:52:54 jvanheld>
#
############################################################

### add the program's directory to the lib path
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	$script_dir = $` || ".";
	push @INC, "$script_dir";
	push @INC, "$script_dir/perllib";
    }
}
require "config.pl";
require "lib/load_classes.pl";
require "lib/util.pl";
require "lib/parsing_util.pl";
require "lib/loading_util.pl";
use SWISS::Entry;

package main;
{
    
    ################################################################
    #### initialize parameters
    $distrib = "uniprot" ; # can be uniprot or swissprot


    #### output fields
    @out_fields = qw( id 
		      source
		      description
		      gene
		      names
		      organisms
		      features
		      db_refs
		      comments
		      swissprot_acs
		      swissprot_ids
		      ECs
		      taxids
		      );

    #### organism selection
    @selected_organisms = ();
    $full_name{yeast} = "Saccharomyces cerevisiae (Baker's yeast)";
    $full_name{ecoli} = "Escherichia coli";
    $full_name{human} = "Homo sapiens (Human)";
    $full_name{mouse} = "Mus musculus";

    #### input directories and files for the old swissprot database
    if ($distrib eq "swissprot") {
	$dir{input} = "$Databases/ftp.expasy.org/databases/sp_tr_nrdb";
	$source{swissprot} = "sprot";
	$source{trembl} = "trembl";
	$source{trembl_new} = "trembl_new";
	$export_subdir = "swissprot";
    } else {
	#### input directories and files for the uniprot database
	$dir{input} = "$Databases/ftp.uniprot.org/pub/databases/uniprot/knowledgebase";
	$source{swissprot} = "uniprot_sprot";
	$source{trembl} = "uniprot_trembl";
	$source{trembl_new} = "uniprot_trembl_new";
	$export_subdir = "uniprot";
    }

    #### default output directory
    $dir{output} = "$parsed_data/$export_subdir/$delivery_date";
    
    #### export classes
    push @classes, "classes::Polypeptide";

    #### oracle schema
    $schema = "swissprot";
    $user = "swissprot";
    $password = "swissprot";

    #### pseudo-serialized format 
    $out_format = "obj";
    
    &ReadArguments();
    
    #### add optional fields
    push @out_fields, keys %out_fields;

    #### check option compatibility
    if ($in_file{acs} && ($#selected_organisms >= 0)) {
	die "Error: options -acs and -org are incompatible\n";
    }

    #### read a list of selected ACCESSION NUMBERS
    if ($in_file{acs}) {
	warn "; Reading Accession Number list from file $in_file{acs}\n" 
	    if ($verbose >=1);
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
	    my $ac = &trim($fields[0]);
	    next unless ($ac);
	    $selected_acs{uc($ac)}++; #### case-insensitive
	}
	close ACS;
	warn "; Selected ACs\n;\t",  join ("\n;\t", sort (keys %selected_acs)), "\n" if ($verbose >= 0);
    }

    #### output directory
    &CheckOutputDir();

    #### output file names
    $outfile{polypeptides} = $dir{output}."/Polypeptide".$suffix.".obj";
    $outfile{errors} = $dir{output}."/swissprot".$suffix.".errors.txt";
    $outfile{stats} = $dir{output}."/swissprot".$suffix.".stats.txt";
    
    #### clean output directory
    if ($clean) {
	system "\\rm -f $dir{output}/swissprot*.txt";
	system "\\rm -f $dir{output}/swissprot*.txt.gz";
	system "\\rm -f $dir{output}/Polypeptide*.tab";
	system "\\rm -f $dir{output}/Polypeptide*.tab.gz";
	system "\\rm -f $dir{output}/Polypeptide*.gz";
	system "\\rm -f $dir{output}/Polypeptide*.obj.gz";
	system "\\rm -f $outfile{stats}";
	system "\\rm -f $outfile{errors}";
    }

    #### open error report file
    open ERR, ">$outfile{errors}" || 
	die "Error: cannot write error file $outfile{errors}\n";
#    ERR = STDERR;
   
    #### select three predefined organisms if none was selected (-org)
    unless (($in_file{acs}) ||
	    ($#selected_organisms >= 0) || 
	    ($export{allorg})){
	push @selected_organisms, "ecoli";
	push @selected_organisms, "human";
	push @selected_organisms, "mouse";
	push @selected_organisms, "yeast";
    }

    #### define regular expressions to match the selected organisms
    foreach $organism (@selected_organisms) {
	$regexp = $full_name{$organism};
	$regexp = $` if ($regexp =~ /\(/); 
	$regexp = $` if ($regexp =~ /\)/);
	$regexp = $` if ($regexp =~ /\'/);
	$selected_organism{uc($full_name{$organism})} = $regexp;
	warn join ("\t", "; Regular expression for organism", $organism, $regexp), "\n" if ($verbose >= 2);
	push @regexps, qr /$regexp/; #### list of regular expressions
    }
    

    #### by default, parse swissprot only
    unless (defined(%data_sources)) {
	$data_sources{swissprot} = 1;
    }

    foreach $db (keys %data_sources) {
	if ($data_sources{$db}) {
	    #### the file can be compressed or not
	    if (-e "$dir{input}/$source{$db}.dat") {
		$in_file{$source{$db}} = "cat $dir{input}/$source{$db}.dat | ";
	    } elsif (-e "$dir{input}/$source{$db}.dat.gz") {
		$in_file{$source{$db}} = "gunzip -c $dir{input}/$source{$db}.dat.gz | ";
	    } elsif (-e "$dir{input}/$source{$db}.dat.Z") {
		$in_file{$source{$db}} = "uncompress -c $dir{input}/$source{$db}.dat.Z | ";
	    } else {
		die "Error: file $dir{input}/$source{$db}.dat does not exist\n";
	    }
	    $files_to_parse{$source{$db}} = 1;
	}
	warn ";\tdata source\t", $db, "\t", $in_file{$source{$db}}, "\n" if ($verbose >= 2);
    }

    #### create class holders
    $polypeptides = classes::ClassFactory->new_class(object_type=>"classes::Polypeptide",
						     prefix=>"spp_");
    $polypeptides->set_out_fields(@out_fields);
    $polypeptides->set_attribute_header("features", join ("\t", 
							  "Feature_key", 
							  "start_pos", 
							  "end_pos", 
							  "description", 
							  "method", 
							  "FTId", 
							  ) );
#    die join( "\t",  "BOUM",  $polypeptides->get_object_type(),
#	      $polypeptides->get_object_type()->get_attribute_header("features")), "\n";

    $polypeptides->set_attribute_header("db_refs", join ("\t", 
							 "DB_identifier", 
							 "primary_key", 
							 "secondary_key", 
							 "description"));
#    $polypeptides->set_attribute_header("comments", join ("\t", "topic", "comments") );

    #### testing mode
    if ($test) {
	warn ";TEST\n" if ($verbose >= 1);
	### fast partial parsing for debugging
	foreach $key (keys %in_file) {
	    next if ($key eq 'acs');
	    $in_file{$key} .= " head -20000 |";
	}
    }


    #### default verbose message
    if ($verbose >= 1) {
	warn "; Selected organisms\n;\t", join("\n;\t", @selected_organisms), "\n";
	warn "; Polypeptide classes\n;\t", join("\n;\t", keys %export), "\n";
	warn "; Data sources\n;\t", join("\n;\t",  keys %data_sources), "\n";
	warn "; Output fields\n;\t", join("\n;\t",  @out_fields), "\n";
	&DefaultVerbose();
    }

    #### parse data from original files
    foreach $source (keys %files_to_parse) {
	&ParseSwissprot($in_file{$source}, $source);
    }

    #### check that all the selected acs have been found
    &CheckSelectedACs() if ($in_file{acs});
    
    #### print the result
    &PrintStats($outfile{stats}, @classes);
    $polypeptides->dump_tables($suffix, 0, $dir{output});

    #### special field formats for SQL
    $special_field_size{features} = 2047;
#    $special_field_size{db_refs} = 2047;
    $special_field_size{description} = 1023;
    $special_field_size{bioseq} = 10000; 
#    $special_field_size{comments} = 10000;

    ### generate SQL scripts for loading the data
    $polypeptides->generate_sql(schema=>$main::schema, 
				user=>$main::user,
				password=>$main::password,
				dir=>"$dir{output}/sql_scripts",
				prefix=>"",
				host=>$main::host,
				);
    &ExportMakefile(@classes);
    &ExportClasses($outfile{polypeptides}, $out_format,classes::Polypeptide) if ($export{obj});
#    &CompressParsedData();
    close ERR;

    ### report execution time
    if ($verbose >= 1) {
	$done_time = &AlphaDate();
	warn (";\n",
	      "; job started $start_time",
	      "; job done    $done_time\n");
    }
    
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
$generic_option_message
	-indir	input directory. 
		This directory should contain a download of the
		non-redundant swiss-prot/trembl ftp site :
		       ftp://ftp.expasy.org/databases/sp_tr_nrdb
	-outdir	output directory.
		The parsed data will be stored in this directory.
	-enz	export enzymes only
	-org	select an organism for exportation
		can be used reiteratively in the command line 
		to select several organisms
		   Supported organisms :
			ecoli
			human
			mouse
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
	-acs	accession file
		The parsing will be restricted to the swissprot
		accession numbers (AC field) or identifers (ID)
		specified in the file. Each accession number or
		identifier must come as the first word of a new line
		(all subsequent words are ignored).  
		This option is incompatible with the option -org.
	-fields
		Fields to be exported. Several fields can be provided,
		separated by commas. Example:
			  -fields id,names,ene,ECs,swissprot_ids
	-seq    return polypeptide bioseqs
EXAMPLE
	parse_swissprot.pl -v 2 -allorg -data swissprot 
EndHelp
  close HELP;
}

### read arguments from the command line
sub ReadArguments {

    for my $a (0..$#ARGV) {

	&ReadGenericOptions($a);
	    
	    #### input directory
	if ($ARGV[$a] =~ /^-indir/) {
	    $dir{input} = $ARGV[$a+1];

	    #### input file with a list of Swissprot AC
	} elsif ($ARGV[$a] =~ /^-acs/) {
	    $in_file{acs} = $ARGV[$a+1];

	    #### select organisms for exportation
	} elsif ($ARGV[$a] =~ /^-org/) {
	    push @selected_organisms, lc($ARGV[$a+1]);

	    #### export all organisms
	} elsif ($ARGV[$a] =~ /^-allorg/) {
	    $main::export{allorg} = 1;
	    
	    #### select enzymes for exportation
	} elsif ($ARGV[$a] =~ /^-data/) {
	    $data_sources{lc($ARGV[$a+1])} = 1;
	    

	    #### select enzymes for exportation
	} elsif ($ARGV[$a] =~ /^-enz/) {
	    $main::export{enzymes} = 1;

	    #### output fields
	} elsif ($ARGV[$a] =~ /^-field/) {
	    @main::out_fields = split ",", $ARGV[$a+1];

	    #### output fields
	} elsif ($ARGV[$a] =~ /^-seq/) {
	    $main::out_fields{bioseq} = 1;

	}


    }

}

################################################################
#### Parse one  swissprot or trembl file
sub ParseSwissprot {
    my ($input_file, $source) = @_;
    $source = $input_file unless ($source);

    ### regular expression for matching organism names
    my $regexp = "";
    if ($#selected_organisms >= 0) {
	$regexp = "(?:".join ("|", @regexps).")";
	warn "; Selecting organisms with regular expression\n;\t", $regexp, "\n" if ($verbose >= 1);
    }

    warn (";\n; ", &AlphaDate,  " parsing polypeptides from $input_file\n")
	if ($verbose >= 1);
    
    open DATA, $input_file || 
	die "Error: cannot open data file $input_file\n";
    
    #### Read an entire record at a time
    $/ = "\/\/\n";
    
    my $entries = 0;
    while ($text_entry = <DATA>){
	$entries++;
	warn "; Entry\t$entries\t",&AlphaDate(), "\n" if (($verbose >= 2) && ($entries%100 == 0));
	warn "$text_entry\n" if ($verbose >= 4);
	
	my $parse = 1;

	################################################################ 
	# Pre-parsing checks, to save memory and time
        #


	#### check that the entry contains a selected AC or ID before
#	if ($in_file{acs}) {
#	    $regexp = "(?:".join ("|", keys %selected_acs).")";
#	}
	if ($regexp) {
	    unless ($text_entry =~ /$regexp/i) {
		$parse = 0;
		warn("; Skipping entry ", $entries, " (no match with selected organisms regexp)\n") if ($verbose >= 3); 
		next;
	    }
	}

	#### convert the entry into an object
	my $object_entry = SWISS::Entry->fromText($text_entry);

	#### get the polypeptide accession number
	my $swissprot_ac = $object_entry->AC;
	my $swissprot_id = $object_entry->ID;

	warn join ("\t", "; Parsing entry", $entries, $swissprot_ac, $swissprot_id), "\n" if ($verbose >= 3);

	#### initialize the export flag
	my $export = 0;

	#### Organisms
	my @organism_objects = $object_entry->OSs->elements();
	my @organisms = ();
	foreach my $organism_object (@organism_objects) {
	    push @organisms, $organism_object->toText();
	}
	warn "; \tOrganisms\t", join ("; ", @organisms), "\n" if ($verbose >= 3);

	#### check whether the accession number was specified for export 
	if ($in_file{acs}) {

	    #### check if the current polypeptide was part of the selection
	    if ($selected_acs{uc($swissprot_ac)}) {
		$found_acs{$swissprot_ac}++;
		delete $selected_acs{uc($swissprot_ac)};
		$export = 1;
	    } elsif ($selected_acs{uc($swissprot_id)}) {
		$found_acs{$swissprot_id}++;
		delete $selected_acs{uc($swissprot_id)};
		$export = 1;
	    }

	} elsif ($export{allorg}) {
	    $export = 1;
	} else {
	    #### check whether the polypeptide organism(s)
	    #### match the organism selection
	    foreach $organism (@organisms) {
		if ($selected_organism{uc($organism)}) {
		    $export = 1;
		    last;
		}
	    }
	}
	

	unless ($export) {
	    warn(";\tSkipping entry ", $entries, " (field organisms does not match selected organisms )\n") if ($verbose >= 3); 
	    next;
	}

#	warn ";\tExporting polypeptide $swissprot_ac\t$swissprot_id\n" if ($verbose >= 3);

	#### Swissprot identifiers
	my @swissprot_ids = $object_entry->IDs->elements;

	#### Swissprot accession numbers
	my @swissprot_acs = $object_entry->ACs->elements;

	#### Description
	my $descr = $object_entry->DEs->text;

	### Extract polypeptide name from description
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
	    $tmp = "$'";
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
	warn join ("\t", ";", "Creating polypeptide",$source, $entries, $swissprot_acs[0], $swissprot_ids[0], $names[0]), "\n"
	    if ($verbose >= 3); 
	my $polypeptide = $polypeptides->new_object(id=>$swissprot_ac,
						    source=>$source);

	#### Gene names
	my $geneNames = $object_entry->GNs->text;
	if ($distrib eq "swissprot") {
	    @geneNames = split " OR ", $geneNames; ## SUPPORT FOR THE OLD FORMAT
	} 

	## parse qualified gene names from uniprot
	my @geneNames = ();
	my %qualifier = (); # names with qualifies
	my @fields = split "\;", $geneNames;
	foreach my $g (@fields) {
	    $g = &trim($g);
	    if ($g =~ /\=/) {
		$name = $'; ####';
		$qualifier = $`;
	    } else {
		$name = $g;
		$qualifier = 'alternate';
	    }
	    push @geneNames, $name;
	    $qualifier{$name} = $qualifier; # index gene names
	    warn join ("\t", ";", "gene name", $name, "qualifier", $qualifier), "\n" if ($main::verbose >= 4);

	}
	## Set single-value "gene" attribute
	if (defined($geneNames[0])) {
	    $polypeptide->set_attribute("gene",$geneNames[0]);
	} else {
	    $polypeptide->set_attribute("gene",$null);
	}

	my %already_assigned = ();
#	foreach my $name (@names, @geneNames, @swissprot_ids, @swissprot_acs) {
	if ($distrib eq "swissprot") {
	    foreach my $name (@names, @geneNames) {
		$polypeptide->push_attribute("names",$name) 
		    unless $already_assigned{uc($name)};
		$already_assigned{uc($name)}++;  #### prevent assigning twice the same name
	    }
	} else {
	    foreach my $name (@names, @geneNames) {
		## Define name qualifier
		my $qualifier = "alternate";
		if (defined($geneNames{$name})) {
		    $qualifier = geneNames{$name};
		}
		unless ($already_assigned{uc($name)}) { #### prevent assigning twice the same name
		    $polypeptide->push_attribute("names",$name) ;
		    $polypeptide->push_expanded_attribute("names_qualified",$name, $qualifier) ;
		    $already_assigned{uc($name)}++;  
		}
	    }
	}
	
	#### Features
	my @features = $object_entry->FTs->elements;
	foreach my $ft (@features) {
#	    warn join "'|'", "HELLO ft", $#{$ft}, @{$ft} , "'\n" if ($verbose >= 10);
	    $polypeptide->push_expanded_attribute("features", @{$ft}[0..5]);
	}
	my @db_refs = $object_entry->DRs->elements;
	foreach my $dr (@db_refs) {
	    $polypeptide->push_expanded_attribute("db_refs", @{$dr});
	}
#	my @comments = $object_entry->CCs->elements;
#	foreach my $cc (@comments) {
#	    $polypeptide->push_expanded_attribute("comments", @{$cc});
#	}

	#### Swissprot IDs
	foreach my $id (@swissprot_ids) {
#print STDERR "id\t$id\n";
	    $polypeptide->push_attribute("swissprot_ids",$id);
	}

	#### Swissprot accession numbers
	foreach my $ac (@swissprot_acs) {
	    $polypeptide->push_attribute("swissprot_acs",$ac);
	}

	#### Description
	$polypeptide->set_attribute("description", $descr);

	#### Sequence
        if ($out_fields{bioseq}) {
	    my $bioseq = $object_entry->SQs->seq;
	    $polypeptide->set_attribute("bioseq", $bioseq);
#	    die join "\t",  $polypeptide->get_attribute("id"), $polypeptide->get_attribute("bioseq"), "\n";
	}

	my $pp_id = $polypeptide->get_id();

	#### Taxids
	my @taxid_objects = $object_entry->OXs->{NCBI_TaxID}->elements();
	my @taxids = ();
	foreach my $taxid_object (@taxid_objects) {
	    my $taxid = $taxid_object->toText();
	    push @taxids, $taxid;
	    $polypeptide->push_attribute("taxids", $taxid);
	}
	warn "; \tTaxids\t", join ("; ", @taxids), "\n" if ($verbose >= 3);


	#### EC numbers
	foreach my $ec (@ECs) {
#	    my $catalysis = $catalyses->new_object(source=>      $source,
#						   catalyst=>    $pp_id,
#						   catalyzed=>   $ec);
	    $polypeptide->push_attribute("ECs",$ec);
	}

	#### Organisms
	foreach my $organism (@organisms) {
	    $polypeptide->push_attribute("organisms",$organism);
	}


#### in construction 
#  	@multi_attributes = qw (
#  				 IDs 
#  				 ACs 
#  				);
#  	@single_attributes = qw (
#  				 DTs 
#  				 OSs 
#  				 DEs 
#  				 GNs 
#  				 OCs 
#  				 Refs 
#  				 CCs 
#  				 KWs 
#  				 DRs 
#  				 FTs 
#  				 Stars 
#  				 SQs 
#  				 );
#  	foreach my $attribute (@multi_attributes) {
#  	    print "\n", $attribute, "\t";
#  	    if ($attr = $object_entry->{$attribute}) {
#  		foreach my $value ($attr->elements) {
#  		    $polypeptide->push_attribute($attribute, $value);
#  		}
#  		print join ";", $object_entry->{$attribute}->elements;
#  	    };
#  	}
#  	foreach my $attribute (@single_attributes) {
#  	    print "\n", $attribute, "\t";
#  	    if ($attr = $object_entry->{$attribute}) {
#  		$polypeptide->set_attribute($attribute, $value);
#  		print join ";", $object_entry->{$attribute}->toText;
#  	    };
#  	}

	#### check how many polypeptides remain to be found
	if ($in_file{acs}) {
	    my $remaining = scalar(keys %selected_acs);
	    warn ";\tfound\t$swissprot_ac\t$swissprot_id\tremaning ACs\t$remaining\n" if ($verbose >= 2);
	    last if ($remaining == 0);
	}
	
    }
    close DATA;

    
#    $polypeptides->index_names();
    return 1;
}


################################################################
# Check that all the selected ACs/IDs have been found
#
sub CheckSelectedACs {
    #### report the selected ACs/IDs which were not found
    foreach my $ac (keys %selected_acs) {
	&ErrorMessage ("Did not find any polypeptide with accession or identifier\t", $ac, "\n");
    }
}

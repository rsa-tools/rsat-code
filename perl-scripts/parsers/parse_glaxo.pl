#!/usr/bin/perl
############################################################
#
# $Id: parse_glaxo.pl,v 1.3 2001/05/30 08:54:55 jvanheld Exp $
#
# Time-stamp: <2001-05-30 10:51:30 jvanheld>
#
############################################################
### parse_ligand.plt
### type parse_ligand.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "PFBP_classes.pl";
require "PFBP_config.pl";
require "PFBP_util.pl";
require "PFBP_parsing_util.pl";

### specific classes for Glaxo parsing
package PFBP::Smiles;
{
    @ISA = qw ( PFBP::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "gsm_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       names=>"ARRAY",
			       SMILES=>"SCALAR",
			       source=>"SCALAR");
}

package PFBP::ECreference;
{
    @ISA = qw ( PFBP::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "grf_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       source=>"SCALAR");
}

package PFBP::Biosequence;
{
    @ISA = qw ( PFBP::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "grf_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       source=>"SCALAR");
}


package main;

### initialization
$start_time = `date +%Y-%m-%d.%H%M%S`;
$delivery_date = `date +%Y%m%d`;
chomp $delivery_date;
$warn_level = 0;

### files to parse
$file_name = "didier.fdt";
$source = "GLAXO".$file_name;
$in_file{glaxo} = "gunzip -c  $data_glaxo/".$file_name.".gz | ";

$dir{output} = "$parsed_data/glaxo_parsed/$delivery_date";
unless (-d $dir{output}) {
    warn "Creating output dir $dir{output}";
    mkdir $dir{output}, 0775 || die "Error: cannot create directory $dir\n";
}
chdir $dir{output};
$out_file{glaxo} = "$dir{output}/glaxo.obj";
$out_file{stats} = "$dir{output}/glaxo.stats.txt";
$out_file{errors} = "$dir{output}/glaxo.errors.txt";

open ERR, ">$out_file{errors}" || die "Error: cannot write error report fle $out_file{errors}\n";

$out_format = "obj";

push @classes, "PFBP::Smiles";
push @classes, "PFBP::Compound";
push @classes, "PFBP::Biosequence";
push @classes, "PFBP::ECreference";

&ReadArguments;

### testing mode
if ($test) {
    warn ";TEST\n" if ($warn_level >= 1);
    ### fast partial parsing for debugging
    $in_file{glaxo} .= " head -10000 |";
}

### default verbose message
&DefaultVerbose if ($warn_level >= 1);

### parse data from original files
$smiles = PFBP::ClassFactory->new_class(object_type=>"PFBP::Smiles",
					prefix=>"gsm_");
$brenda_names = PFBP::ClassFactory->new_class(object_type=>"PFBP::Compound",
					      prefix=>"gbn_");
$biosequences = PFBP::ClassFactory->new_class(object_type=>"PFBP::Biosequence",
					      prefix=>"gbs_");
$ec_references = PFBP::ClassFactory->new_class(object_type=>"PFBP::ECreference",
					       prefix=>"grf_");
&ParseGlaxoFile();
#$smiles->index_names();
#$brenda_names->index_names();

### define a primary name (take the first value in the name list)
foreach $smile ($smiles->get_objects()) {
    @BRENDA_names = $smile->get_attribute("BRENDA_name");
    if ($#BRENDA_names >= 0) {
	$smile->set_attribute("primary_name",$BRENDA_names[0]);
    } elsif ($name = $smile->get_name()) {
	$smile->set_attribute("primary_name",$name);
    }
}

foreach $brenda_name ($brenda_names->get_objects()) {
    if ($name = $brenda_name->get_name()) {
	$brenda_name->set_attribute("primary_name",$name);
    }
}

### print result
$smiles->dump_tables();
$brenda_names->dump_tables();
$biosequences->dump_tables();
$ec_references->dump_tables();
&ExportClasses($out_file{glaxo}, $out_format, @classes);

### print some stats after parsing
&PrintStats($out_file{stats}, @classes);

### report execution time
if ($warn_level >= 1) {
    $done_time = `date +%Y-%m-%d.%H%M%S`;
    warn ";\n";
    warn "; job started $start_time";
    warn "; job done    $done_time";
}

close ERR;

warn "compressing parsed files\n" 
    if ($warn_level >= 1);
system "gzip -f $dir{output}/*.tab $dir{output}/*.obj $dir{output}/*.txt";

if ($warn_level >= 1) {
    my @commands = ("ls -l '$dir{output}'",
		    "gunzip -c '${out_file{stats}}.gz'");
    foreach my $command (@commands) {
	warn "\t$command\n";
#	system $command;
    }
}

exit(0);

### subroutines for the main package

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
    open HELP, "| more";
    print <<EndHelp;
    NAME
        parse_ligand.pl

    DESCRIPTION

	Parse smiles from Glaxo.  The daylight file is parsed in two
	separate files Smiles one object per SMILES entry in the Glaxo
	file Compounds one object per "BRENDA name" entry in the Glaxo
	file

					
	I think that BRENDA names are nothing else than synonyms for
	the compound associated to each SMILES. If this happens to be
	true, all attributes of a BRENDA name should become attributes
	of the parent SMILES. This needs to be checked with
	Ceara. inbetween, I prefer to create separate objects. It will
	be easy to merge separate objects than to split a composite
	object.

    AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

    VERSION
	0.01
	Created		2000/10/25
	Last modified	2000/10/25

    SYNOPSIS	
	parse_ligand.pl [-v] [-vv] [-i infile] [-format output_format]
	    [-o outfile]

    OPTIONS	
	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-w #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-i	input file
		If ommited, the defult file is used (specified in the
		initialization of the code).
	-o	output file
EndHelp
    close HELP;
}


### read arguments from the command line
sub ReadArguments {
    my $a = "";
    for $a (0..$#ARGV) {

	### warn level
	if (($ARGV[$a] eq "-w" ) && 
	    ($ARGV[$a+1] =~ /^\d+$/)){
	    $main::warn_level = $ARGV[$a+1];
	    
	    ### fast test
	} elsif ($ARGV[$a] eq "-test") {
	    $test = 1;

	    ### input file
	} elsif ($ARGV[$a] eq "-i") {
	    $a++;
	    $in_file{glaxo} = $ARGV[$a];
	    
	    ### output file
	} elsif ($ARGV[$a] eq "-o") {
	    $a++;
	    $out_file{glaxo} = $ARGV[$a];

	    ### help
	} elsif (($ARGV[$a] eq "-h") ||
		 ($ARGV[$a] eq "-help")) {
	    &PrintHelp;
	    exit(0);
	}
    }
}


sub ParseGlaxoFile{
    ### open the input stream
    if ($in_file{glaxo}) {
	warn (";\n; ", &AlphaDate,
	      " parsing class $class from $in_file{glaxo}\n")
	    if ($warn_level >= 1);
	open COMP, $in_file{glaxo} 
	|| die "Error: cannot $in_file{glaxo}\n";
	$in = COMP;
    } else {
	warn (";\n; ", &AlphaDate,
	      " parsing class $class from STDIN\n")
	    if ($warn_level >= 1);
	$in = STDIN;
    } 
    

    my $current_smiles;
    my %smiles_multi_value_attributes;
    $smiles_multi_value_attributes{"Remark"} = 1;
    $smiles_multi_value_attributes{"BRENDA name"} = 1;
    $smiles_multi_value_attributes{"Biosequence"} = 1;

    my @smiles_single_value_attributes;
    $smiles_single_value_attributes{"Parent Avg molecular weight"} = 1;
    $smiles_single_value_attributes{"Parent Formula"} = 1;
#  $smiles_single_value_attributes{"Biosequence"} = 1;
    $smiles_single_value_attributes{"LIGAND AccNo"} = 1;
    
    my @smiles_single_value_attributes;
    foreach $attr (keys %smiles_multi_value_attributes, keys %smiles_single_value_attributes) {
	$field_name{$attr} = $attr;
	$field_name{$attr} =~ s/\s/_/g;
    }

    my $current_child;
    my %child_multi_value_attributes;
    $child_multi_value_attributes{"Remark"} = 1;
    $child_single_value_attributes{"Isomer"} = 1;
    $child_multi_value_attributes{"Local Name"} = 1;
    $child_multi_value_attributes{"Avg molecular weight"} = 1;
    $child_single_value_attributes{"Mol Formula"} = 1;
    $child_multi_value_attributes{"EC Number"} = 1;
    $child_single_value_attributes{"Parent Avg molecular weight"} = 1;
    $child_single_value_attributes{"Parent Formula"} = 1;
    
    my @child_single_value_attributes;
    foreach $attr (keys %child_multi_value_attributes, keys %child_single_value_attributes) {
	$field_name{$attr} = $attr;
	$field_name{$attr} =~ s/\s/_/g;
    }

    #### some specific names
    $field_name{"Local Name"} = "names";

    my $l = 0;
    ### read data
    while (<$in>) {
	$l++;
	chomp;
	if (/^SMILES\s*/) { ### new record
	    my $smiles_string = $';
	    undef($current_child);

	    ### instantiate new object
	    warn ";creating object\t$current_id in class\t$class\n" 
		if ($warn_level >= 2);
	    warn ("\t", $smiles, "\n")
		if ($warn_level >= 3);
	    $current_smiles = $smiles->new_object(source=>$source);
	    $current_smiles->set_attribute("SMILES",$smiles_string);
	    
	    ### reeinitialize the arguments
	    my %attributes = ();
	    
	    ### read next line
	    next;
	}

	### SMILES attributes (single \t at the beginning of the line)
	my $found = 0;
	foreach $attr (keys %smiles_multi_value_attributes, keys %smiles_single_value_attributes) {
	    if (/^\t$attr\s+/) {
		$value = $';
		undef($current_chlid);
		if (defined($smiles_multi_value_attributes{$attr})) {
		    $current_smiles->new_attribute_value($field_name{$attr},$value);
		    if ($attr eq "BRENDA name") {
			$current_child = $brenda_names->new_object(source=>$source);
			$current_child->set_attribute("BRENDA_name",$value);
			$current_child->set_attribute("parent_smiles",$current_smiles->get_id());
		    } elsif ($attr eq "Biosequence") {
			$current_child = $biosequences->new_object(source=>$source);
			$current_child->set_attribute("biosequence",$value);
			$current_child->set_attribute("parent_smiles",$current_smiles->get_id());
		    }
		    $found = 1;
		    last;
		} elsif (defined($smiles_single_value_attributes{$attr})) {
		    $current_smiles->set_attribute($field_name{$attr},$value);
		    $found = 1;
		    last;
		}
	    }
	}
	next if ($found);
	if (/^\t\S/) {
	    print ERR ("Error line $l\tunknown attribute for SMILES, not parsed\t", $_, "\n");
	    next;
	}

	### treatment of the attributes starting with \t\t these are
	### actually attributes of the "BRENDA name" However my
	### understanding is that all the BRENDA names included in a
	### SMILE are just synonyms for the same compound. If this is
	### true (should be checked with Ceara), all the attributes
	### belong in fact to the unique compound (the SMILES entry in
	### Glaxo file). Without confirmation, I prefer to hold the
	### information in separate objects (it will be easier to merge
	### then afterwards)
	
	my $found = 0;
	foreach $attr (keys %child_multi_value_attributes, keys %child_single_value_attributes) {
	    if (/^\t\t$attr\s+/) {
		$value = $';
		if (defined($current_child)) {
		    if ($attr eq "EC Number") { #### specific reatment for EC references, which are actually whole objects
			$current_ECref = $ec_references->new_object(source=>$source);
			@fields = split "\t", $value;
			$current_ECref->set_attribute("compound", $current_child->get_id());
			$current_ECref->set_attribute("EC", $fields[0]);
			$current_ECref->set_attribute("action", $fields[2]);
			$current_ECref->set_attribute("reference", $fields[4]);
			$found = 1;
			last;
		    } elsif (defined($child_multi_value_attributes{$attr})) {
			$current_child->new_attribute_value($field_name{$attr},$value);
			$found = 1;
			last;
		    } elsif (defined($child_single_value_attributes{$attr})) {
			$current_child->set_attribute($field_name{$attr},$value);
			$found = 1;
			last;
		    }
		} else {
		    print ERR "Error line $l\tchild undefined\t$_\n";
		}
	    }
	}
	next if ($found);
	
	### if we reached here, there is a parsing problem. report the line
	print ERR ("error line $l\t", $_, "\n");
    }
    
    close $in if ($in_file{glaxo});
}

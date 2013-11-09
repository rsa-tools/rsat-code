#!/usr/bin/perl
############################################################
#
# $Id: parse_glaxo.pl,v 1.7 2011/02/17 05:07:46 rsat Exp $
#
# Time-stamp: <2003-07-10 11:53:01 jvanheld>
#
############################################################
### parse_ligand.plt
### type parse_ligand.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "lib/load_classes.pl";
require "config.pl";
require "lib/util.pl";
require "lib/parsing_util.pl";

### specific classes for Glaxo parsing
package classes::Smiles;
{
    @ISA = qw ( classes::DatabaseObject );
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
			       Biosequence=>"SCALAR",
			       Parent_Formula=>"ARRAY",
			       source=>"SCALAR");
}

package classes::ECreference;
{
    @ISA = qw ( classes::DatabaseObject );
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

package classes::BrendaName;
{
    @ISA = qw ( classes::Compound );
    ### class attributes
    $_count = 0;
    $_prefix = "gbn_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       names=>"ARRAY",
			       source=>"SCALAR",
			       parent_smiles=>"SCALAR",
			       BRENDA_name=>"SCALAR",
			       Mol_Formula=>"ARRAY",
			       Isomer=>"ARRAY");
}

package classes::Biosequence;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "gbs_";
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
$verbose = 0;
$clean = 0;

### files to parse
$file_name = "didier.fdt";
$source = "GLAXO:".$file_name;
$data_glaxo = "/win/amaze/Databases/Glaxo/glaxo_compounds/2000_09_glaxo_compounds";
$dir{input} = $data_glaxo;
$in_file{glaxo} = "gunzip -c ".$dir{input}."/".$file_name.".gz | ";


open ERR, ">$outfile{errors}" || die "Error: cannot write error report fle $outfile{errors}\n";

$out_format = "obj";

#### default export directory
$export_subdir = "glaxo";
$dir{output} = "$parsed_data/${export_subdir}/$delivery_date";

push @classes, "classes::Smiles";
push @classes, "classes::BrendaName";
push @classes, "classes::Biosequence";
push @classes, "classes::ECreference";

&ReadArguments;


#### output directory
&CheckOutputDir();

$outfile{glaxo} = "$dir{output}/glaxo.obj";
$outfile{stats} = "$dir{output}/glaxo.stats.txt";
$outfile{errors} = "$dir{output}/glaxo.errors.txt";

#### clean the output directory
if ($clean) {
    warn "Cleaning output directory $dir{output}\n" if ($verbose >=1);
    system "\\rm -f $dir{output}/*.tab.gz $dir{output}/*.txt.gz $dir{output}/*.obj.gz" ;
#    system "\\rm -rf $dir{output}/mirror";
}


### test mode
if ($test) {
    warn ";TEST\n" if ($verbose >= 1);
    ### fast partial parsing for debugging
    $in_file{glaxo} .= " head -5000 |";
}



### default verbose message
&DefaultVerbose if ($verbose >= 1);

### parse data from original files
$smiles = classes::ClassFactory->new_class(object_type=>"classes::Smiles",
					prefix=>"gsm_");
@smiles_out_fields = qw (id
			 source
			 Parent_Avg_molecular_weight	
			 LIGAND_AccNo	
			 primary_name	
			 SMILES	
			 Parent_Formula
			 BRENDA_name
			 Biosequence
			 Remark
#			 names
			 );
if ($export{analysis}) {
    push @smiles_out_fields, qw( gbn_count gbs_count );
}
$smiles->set_out_fields(@smiles_out_fields);

$brenda_names = classes::ClassFactory->new_class(object_type=>"classes::BrendaName",
					      prefix=>"gbn_");

$biosequences = classes::ClassFactory->new_class(object_type=>"classes::Biosequence",
					      prefix=>"gbs_");
$ec_references = classes::ClassFactory->new_class(object_type=>"classes::ECreference",
					       prefix=>"grf_");
&ParseGlaxoFile();
#$smiles->index_names();
#$brenda_names->index_names();

### define a primary name (take the first value in the name list)
foreach $smile ($smiles->get_objects()) {
    @BRENDA_names = $smile->get_attribute("BRENDA_name");
#    if ($#BRENDA_names >= 0) {
#	$smile->set_attribute("primary_name",$BRENDA_names[0]);
#    } elsif ($name = $smile->get_name()) {
#	$smile->set_attribute("primary_name",$name);
#    }
}

foreach $brenda_name ($brenda_names->get_objects()) {
    if ($name = $brenda_name->get_name()) {
	$brenda_name->set_attribute("primary_name",$name);
    }
}

#### restore the biosequence attributes to their parent SMILES
foreach my $gbs_obj ($biosequences->get_objects()) {
    my $gsm_id = $gbs_obj->get_attribute("parent_smiles");
    my $gsm_obj = $smiles->get_object($gsm_id);

    foreach my $paf ($gbs_obj->get_attribute("Parent Formula")) {
	$gsm_obj->push_attribute("Parent Formula", $paf);
    }

#    print STDERR join ("\t", "HELLO", 
#		       $gsm_id, 
#		       $gbs_obj->get_attribute("biosequence"), 
#		       $gbs_obj->get_attribute("Parent_Formula")), "\n";
    

    foreach my $pamw ($gbs_obj->get_attribute("Parent_Avg_molecular_weight")) {
	$gsm_obj->push_attribute("Parent_Avg_molecular_weight", $pamw);
#	print STDERR "PAMW\t", $gsm_id, "\t", $pamw, "\n";
    }
    foreach my $remark ($gbs_obj->get_attribute("Remark")) {
	$gsm_obj->push_attribute("Remark", $remark);
	#print STDERR "REMARK\t", $gsm_id, "\t", $remark, "\n";
    }
}


#### analyze the result
if ($export{analysis}) {
    warn "; Analyzing parsing result\n" if ($verbose >=1);
    foreach my $gbn_obj ($brenda_names->get_objects()) {
	my $gsm_id = $gbn_obj->get_attribute("parent_smiles");
	$gbn_count{$gsm_id}++;
    }
    foreach my $gbs_obj ($biosequences->get_objects()) {
	my $gsm_id = $gbs_obj->get_attribute("parent_smiles");
	$gbs_count{$gsm_id}++;
    }
    foreach my $gsm_obj ($smiles->get_objects()) {
	my $gsm_id = $gsm_obj->get_attribute("id");
	$gsm_obj->set_attribute("gbn_count", $gbn_count{$gsm_id} || 0);
	$gsm_obj->set_attribute("gbs_count", $gbs_count{$gsm_id} || 0);
    }
}

### export the result
warn "Dumping objects to tab files\n" if ($verbose >=1);
$smiles->dump_tables();
$brenda_names->dump_tables();
$biosequences->dump_tables();
$ec_references->dump_tables();

if ($export{obj}) {
    warn "Exporting objects in .obj format\n" if ($verbose >=1);
    &ExportClasses($outfile{glaxo}, $out_format, @classes) 
}


### print some stats after parsing
warn "Printing parsing statitsics\n" if ($verbose >=1);
&PrintStats($outfile{stats}, @classes);

### report execution time
if ($verbose >= 1) {
    $done_time = `date +%Y-%m-%d.%H%M%S`;
    warn ";\n";
    warn "; job started $start_time";
    warn "; job done    $done_time";
}

close ERR;

&CompressParsedData();

if ($verbose >= 1) {
    my @commands = ("ls -l '$dir{output}'",
		    "gunzip -c '${outfile{stats}}.gz'");
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
	separate files :
	- Smiles:      one object per SMILES entry in the Glaxo file 
	- BrendaName:  one object per "BRENDA name" entry in the 
		       Glaxo file
					
	I think that BRENDA names are nothing else than synonyms for
	the compound associated to each SMILES. If this happens to be
	true, all attributes of a BRENDA name should become attributes
	of the parent SMILES. This needs to be checked with
	Ceara. inbetween, I prefer to create separate objects. It will
	be easy to merge separate objects than to split a composite
	object.

	Note added 2001/09/10 
	=====================
	Apparently some SMILES correspond to generic compounds, which
	are parent for a list of well distinct compounds in KEGG. See
	for example OC1COC(O)C(O)C1O.

    AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

    VERSION
	0.01
	Created		2000/10/25
	Last modified	2000/10/25

    OPTIONS	
	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-v #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-obj	export data in object format (.obj file)
		which are human-readable
	-analysis
	        analyze the data and report the result. 
	-clean remove all files from the outpu directory before
	       parsing.
EndHelp
    close HELP;
}


### read arguments from the command line
sub ReadArguments {
    my $a = "";
    for $a (0..$#ARGV) {

	### warn level
	if (($ARGV[$a] eq "-v" ) && 
	    ($ARGV[$a+1] =~ /^\d+$/)){
	    $main::verbose = $ARGV[$a+1];
	    
	    ### fast test
	} elsif ($ARGV[$a] eq "-test") {
	    $test = 1;
	    
	    ### export in object format
 	} elsif ($ARGV[$a] eq "-obj") {
	    $a++;
	    $main::export{obj} = 1;
	    
	    ### analyze the data
 	} elsif ($ARGV[$a] eq "-analysis") {
	    $a++;
	    $main::export{analysis} = 1;
	    
	    ### clean
	} elsif ($ARGV[$a] eq "-clean") {
	    $main::clean = 1;

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
	    if ($verbose >= 1);
	open COMP, $in_file{glaxo} 
	|| die "Error: cannot $in_file{glaxo}\n";
	$in = COMP;
    } else {
	warn (";\n; ", &AlphaDate,
	      " parsing class $class from STDIN\n")
	    if ($verbose >= 1);
	$in = STDIN;
    } 
    

    my $current_smiles;
    my %smiles_multi_value_attributes;
    $smiles_multi_value_attributes{"Remark"} = 1;
    $smiles_multi_value_attributes{"BRENDA name"} = 1;
    $smiles_single_value_attributes{"Biosequence"} = 1;
    $smiles_multi_value_attributes{"LIGAND AccNo"} = 1;
    $smiles_multi_value_attributes{"Parent Avg molecular weight"} = 1;

    my @smiles_single_value_attributes;
    $smiles_multi_value_attributes{"Parent Formula"} = 1;
#  $smiles_single_value_attributes{"Biosequence"} = 1;
    
    my @smiles_single_value_attributes;
    foreach $attr (keys %smiles_multi_value_attributes, keys %smiles_single_value_attributes) {
	$field_name{$attr} = $attr;
	$field_name{$attr} =~ s/\s/_/g;
    }

    my $current_child;
    my %child_multi_value_attributes;
    $child_multi_value_attributes{"Remark"} = 1;
    $child_multi_value_attributes{"Local Name"} = 1;
    $child_multi_value_attributes{"Avg molecular weight"} = 1;
    $child_multi_value_attributes{"EC Number"} = 1;
    $child_multi_value_attributes{"Isomer"} = 1; #### THERE IS A SINGLE CASE WITH 2 VALUES IN THE WHOLE FILE ! IS THIS AN ERROR ? 

    $child_multi_value_attributes{"Mol Formula"} = 1; #### there is at least one case at the 14921th BRENDA name
    $child_multi_value_attributes{"Parent Avg molecular weight"} = 1; ### this is for biosequence
    $child_multi_value_attributes{"Parent Formula"} = 1; ### this is for biosequence
    
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
	    warn ("\tsmiles string", $smiles_string, "\n")
		if ($verbose >= 3);
	    $current_smiles = $smiles->new_object(source=>$source);
	    $current_smiles->set_attribute("SMILES",$smiles_string);
	    $smiles_count++;
	    if (($verbose >=1) and ($smiles_count % 100 == 0)) {
		warn "; Treated\t", $smiles_count, " Smiles and ", "\t", $brenda_names_count, " BrendaNames", "\n";
	    }
	    warn ";created object\t", $current_smiles->get_attribute("id"), "\n" 
		if ($verbose >= 2);
	    
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
			$brenda_names_count++;
			$current_child->set_attribute("BRENDA_name",$value);
			$current_child->set_attribute("parent_smiles",$current_smiles->get_id());
		    }
		    $found = 1;
		    last;
		} elsif (defined($smiles_single_value_attributes{$attr})) {
		    $current_smiles->set_attribute($field_name{$attr},$value);
		    if ($attr eq "Biosequence") {
			$current_child = $biosequences->new_object(source=>$source);
			$current_child->set_attribute("biosequence",$value);
			$current_child->set_attribute("parent_smiles",$current_smiles->get_id());
		    }
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

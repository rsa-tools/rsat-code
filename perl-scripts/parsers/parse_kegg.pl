#!/usr/bin/perl
############################################################
#
# $Id: parse_kegg.pl,v 1.2 2000/12/14 20:55:49 jvanheld Exp $
#
# Time-stamp: <2000-12-14 12:44:45 jvanheld>
#
############################################################

### type parse_kegg.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "PFBP_classes.pl";
require "PFBP_config.pl";
require "PFBP_util.pl";
require "PFBP_parsing_util.pl";

package PFBP::GenericPathway;
### A class to treat EC numenclature
{
  @ISA = qw ( PFBP::ObjectSet );
  ### class attributes
  $_count = 0;
  $_prefix = "pthw_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      names=>"ARRAY",
		      reactions=>"ARRAY",
		      ECs=>"ARRAY",
		     );
}

package PFBP::Pathway;
### A class to treat EC numenclature
{
  @ISA = qw ( PFBP::ObjectSet );
  ### class attributes
  $_count = 0;
  $_prefix = "pthw_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      organism=>"SCALAR",
		      parent=>"SCALAR",
		      source=>"SCALAR",
		      names=>"ARRAY",
		      genes=>"ARRAY"
		     );
}


package main;

################################################################
#
# files to parse
#

#### ligand
$in_file{compounds} = "gunzip -c ".$dir{KEGG}."/molecules/ligand/compound.gz | ";
$in_file{reactions} = "gunzip -c ".$dir{KEGG}."/molecules/ligand/reaction.lst.gz | ";
$in_file{ec} = "gunzip -c ".$dir{KEGG}."/molecules/ligand/ECtable.gz | ";
$in_file{reaction2ec} = "uncompress -c ".$dir{KEGG}."/molecules/ligand/reaction.tar.Z | tar -xpOf - | cut -d ':' -f 1,2 | sort -u |";

#### pathways
$dir{pathway_reactions} = $dir{KEGG}."/molecules/ligand/reaction/";
$dir{eco} = $dir{KEGG}."/pathways/eco/";
$dir{sce} = $dir{KEGG}."/pathways/sce/";
$dir{hsa} = $dir{KEGG}."/pathways/hsa/";
$in_file{pathway_names} = $dir{KEGG}."/pathways/map_title.tab";

################################################################
#
# output files
#
$dir{output} = "$parsed_data/kegg_parsed/$delivery_date";
unless (-d $dir{output}) {
    warn "Creating output dir $dir{output}";
    mkdir $dir{output}, 0775 || die "Error: cannot create directory $dir\n";
}
chdir $dir{output};
$out_file{kegg} = "$dir{output}/kegg.obj";
$out_file{connectivity} = "$dir{output}/kegg.connectivity.txt";
$out_file{stats} = "$dir{output}/kegg.stats.txt";
$out_file{errors} = "$dir{output}/kegg.errors.txt";

open ERR, ">$out_file{errors}" || die "Error: cannot write error report fle $$out_file{errors}\n";

$out_format = "obj";

push @classes, ("PFBP::Compound");
push @classes, ("PFBP::Reaction");
push @classes, ("PFBP::Reactant");
push @classes, ("PFBP::ECSet");
push @classes, ("PFBP::KeggPathway");

### default output fields for each class
@{$out_fields{'PFBP::Compound'}} = qw( id names formula source );
@{$out_fields{'PFBP::Reaction'}} = qw( id  source);
@{$out_fields{'PFBP::Reactant'}} = qw( id reactant_type reaction_id compound_id stoichio valid_interm );
@{$out_fields{'PFBP::ECSet'}} = qw( id names parent reactions );
@{$out_fields{'PFBP::Pathway'}} = qw( id parent organism source names reactions ECs genes );
@{$out_fields{'PFBP::GenericPathway'}} = qw( id source names reactions ECs );

&ReadArguments;

### testing mode
if ($test) {
    warn ";TEST\n" if ($warn_level >= 1);
    ### fast partial parsing for debugging
    $in_file{compounds} .= " head -5000 |";
    $in_file{reactions} .= " head -50 |";
    $in_file{ec} .= " head -20 |";
    $in_file{reaction2ec} .= " head -200 |";
}

### default verbose message
&DefaultVerbose if ($warn_level >= 1);

### instantiate class factories
$compounds = PFBP::ClassFactory->new_class(object_type=>"PFBP::Compound",
					   prefix=>"comp_");
$reactions = PFBP::ClassFactory->new_class(object_type=>"PFBP::Reaction",
					   prefix=>"rctn_");
$reactants = PFBP::ClassFactory->new_class(object_type=>"PFBP::Reactant",
					   prefix=>"rctt_");
$ecs = PFBP::ClassFactory->new_class(object_type=>"PFBP::ECSet",
				     prefix=>"ec_");
$pathways = PFBP::ClassFactory->new_class(object_type=>"PFBP::Pathway",
					  prefix=>"pthw_");
$genericPathways = PFBP::ClassFactory->new_class(object_type=>"PFBP::GenericPathway",
						 prefix=>"gptw_");

#### parse compounds
&ParseKeggFile($in_file{compounds}, $compounds, source=>'KEGG:compound');
$compounds->index_names();


#### parse reactions
&ParseReactions($in_file{reactions}, $reactions);

### parse EC numbers
&ParseEC($in_file{ec}, $ecs);
&ParseReactionEC($in_file{reaction2ec});
&ParseKeggPathways();

foreach $compound ($compounds->get_objects()) {
    ### define a primary name (take the first value in the name list)
#    if ($name = $compound->get_name()) {
#	$compound->set_attribute("primary_name",$name);
#    }
    #### check for compounds without formula
    if ($compound->get_attribute("formula") eq "<UNDEF>") {
	$compound->set_attribute("formula","<NULL>");
    }
}



### print result
$compounds->dump_tables();
$reactions->dump_tables();
$reactants->dump_tables();
$ecs->dump_tables();
$pathways->dump_tables();
$genericPathways->dump_tables();
&ExportClasses($out_file{kegg}, $out_format, @classes) if ($export{obj});


### print some stats after parsing
&PrintStats($out_file{stats}, @classes);
&Connectivity($out_file{connectivity});
#&AppendECstats($out_file{stats});

### report execution time
if ($warn_level >= 1) {
    $done_time = `date +%Y-%m-%d.%H%M%S`;
    warn ";\n";
    warn "; job started $start_time";
    warn "; job done    $done_time";
}

close ERR;

system "gzip -f $dir{output}/*.tab $dir{output}/*.txt";
system "gzip -f $dir{output}/*.obj" if ($export{obj});

exit(0);

### subroutines for the main package

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
        parse_kegg.pl

DESCRIPTION
	Parse	compounds from KEGG (http://www.genome.ad.jp/kegg/). 
	
AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

VERSION
	0.01
	Created		1999/12/16
	Last modified	2000/10/25
	
OPTIONS	
	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-w #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-obj	export data in object format (.obj file)
		which is human-readable (with some patience and
		a good cup of coffee)
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

	    ### export .obj file
	} elsif ($ARGV[$a] eq "-obj") {
	    $export{obj} = 1;

	    ### help
	} elsif (($ARGV[$a] eq "-h") ||
		 ($ARGV[$a] eq "-help")) {
	    &PrintHelp;
	    exit(0);
	}
    }
}



sub TrivialCompounds {
    $trivial{"C00001"} = 1; ## C00001	660	284	944	H2O
    $trivial{"C00002"} = 1; ## C00002	326	3	329	ATP
    $trivial{"C00007"} = 1; ## C00007	303	2	305	Oxygen
    $trivial{"C00003"} = 1; ## C00003	208	49	257	NAD+
    $trivial{"C00004"} = 1; ## C00004	47	200	247	NADH
    $trivial{"C00008"} = 1; ## C00008	6	235	241	ADP
    $trivial{"C00005"} = 1; ## C00005	82	155	237	NADPH
    $trivial{"C00006"} = 1; ## C00006	154	83	237	NADP+
    $trivial{"C00009"} = 1; ## C00009	36	194	230	Orthophosphate
    $trivial{"C00010"} = 1; ## C00010	67	139	206	CoA
    $trivial{"C00011"} = 1; ## C00011	12	191	203	CO2
    $trivial{"C00013"} = 1; ## C00013	20	161	181	Pyrophosphate
    $trivial{"C00015"} = 1; ## C00015	0	135	135	UDP
    $trivial{"C00028"} = 1; ## C00028	76	32	108	Acceptor
    $trivial{"C00030"} = 1; ## C00030	32	73	105	Reduced acceptor
    $trivial{"C00020"} = 1; ## C00020	11	92	103	AMP
    $trivial{"C00027"} = 1; ## C00027	12	72	84	H2O2
}


sub ParseReactions {
    ### read the reaction file from KEGG
    my ($input_file, $class_holder) = @_;

    &TrivialCompounds;
    
    if ($warn_level >= 1) {
	warn (";\n; ", 
	      &AlphaDate, 
	      "\tparsing class ", 
	      $class_holder->get_object_type(),
	      " from ");
	if ($input_file) {
	    warn "$input_file\n";
	} else {
	    warn "STDIN\n";
	}
    }
    ### open the input stream
    if ($input_file) {
	open RXN, $input_file 
	    || die "Error: cannot read $input_file\n";
	$rxn_handle = RXN;
    } else {
	$rxn_handle = STDIN;
    } 
    
    ### read data
    while (<$rxn_handle>) {
	if (/(R\d+):\s*(.*)\s*$/) {
	    $reaction_id = $1;
	    $reaction = $reactions->new_object(id=>$reaction_id,
					       source=>'KEGG:reaction.lst');
	    $equation = $2;
	    $reaction->set_attribute("equation",$equation);
	    
	    if ($equation =~ /^\s*(.*)\s*\<\=\>\s*(.*)\s*$/) {
		$left = $1;
		$right = $2;

		### temporary
#		$reaction->set_attribute("left",$left);
#		$reaction->set_attribute("right",$right);

		### split the equation into reactants
		%substrate = &SplitReactants($left);
		%product = &SplitReactants($right);
		
		foreach $reactant_type ("substrate","product") {
		    while (($compound, $stoichio) = each %{$reactant_type}) {
			### instantiate a new reactant
			$reactant = $reactants->new_object();
			$reactant_id = $reactant->get_attribute("id");
			$reactant->set_attribute("compound_id",$compound);
			$reactant->set_attribute("reaction_id",$reaction_id);
			$reactant->set_attribute("reactant_type",$reactant_type);
			$reactant->new_attribute_value("stoichio",$stoichio);
			
			### validity of the reactant as intermediate
			### between two successive reactions in a pathway
			if ($trivial{$compound}) {
			    $validity = "0";
			} else {
			    $validity = "1";
			}
			$reactant->set_attribute("valid_interm", $validity);
#			printf STDERR ";\t%s\t%d\n", $compound, 1 - $trivial{$compound};
			### assign the reactant to the reaction
#			$reaction->new_attribute_value($reactant_type,$reactant_id);
		    }
		}
	    }
	}
    }
    close $rxn_handle if ($input_file);
}


sub ParseEC {
    my ($file, $class_holder, %args) = @_;
    warn (";\n; ", &AlphaDate, "; parsing class ",
	  $class_holder->get_object_type(),
	  " from $file\n")
	if ($warn_level >= 1);
    
    ### create the class for 0.0.0.0
    $ec_object = $ecs->new_object(%args,id=>"0.0.0.0");
    $ec_object->push_attribute("names", "non-enzymatic or not clearly enzymatic");
    $ec_object->set_attribute("parent","<NULL>");

    ### create the class for -.-.-.-
    $ec_object = $ecs->new_object(%args,id=>"-.-.-.-");
    $ec_object->push_attribute("names", "non-assigned EC number");
    $ec_object->set_attribute("parent","<NULL>");

    ### open the input stream
    if ($file) {
	open EC, $file 
	    || die "Error: cannot $kegg_file\n";
	$in = EC;
    } else { ### use standard input
	$in = STDIN;
    } 
    
    ### read data
    while (<$in>) {
	chomp;
	if (/^\s*((\d+\.*)+)\s+/) {
	    $ec_number = $1;
	    $names = $';
	    
	    $ec_number =~ s/\.$//; ### suppress the trailing dot from the ec number
	    @names = split /\; /, $names;
	    
	    $ec_object = $ecs->new_object(%args,id=>$ec_number);
	    foreach $name (@names) {
		$ec_object->push_attribute("names",$name);
	    }
	    
	    ### partially specified ecs
	    unless ($ec_number =~ /^\d+\.\d+\.\d+\.\d+$/) {
		if ($ec_number =~ /^\d+\.\d+\.\d+$/) { 
		    $uncomplete_id = $ec_number.".-";
		} elsif ($ec_number =~ /^\d+\.\d+$/) { 
		    $uncomplete_id = $ec_number.".-.-";
		} elsif ($ec_number =~ /^\d+$/) { 
		    $uncomplete_id = $ec_number.".-.-.-";
		} else {
		    $uncomplete_id = undef;
		    &Errormessage("Error: non-standard EC number\n");
		}
		if ($uncomplete_id) {
		    $uncomplete_object = $ecs->new_object(%args,id=>$uncomplete_id);
		    foreach $name (@names) {
			$uncomplete_object->push_attribute("names",$name);
		    }
#		    $ec_object->push_attribute("subsets",$uncomplete_id);
		    $uncomplete_object->set_attribute("parent",$ec_number);
		}
	    }
	} 
    }
    close $in if ($file);
    
    
    ### create indexes
    $ecs->index_ids();
    $ecs->index_names();
    
    ### create taxonomic relationships between EC numbersc
    #my %objects = $class->get_id_index();
    foreach my $ec_object ($ecs->get_objects()) {
	if (my $ec_number = $ec_object->get_attribute("id")) {
	    if ($ec_number =~ /\.\d+$/) {
		$parent_ec = $`;
		if ($parent_object = $ecs->get_object($parent_ec)) {
#		    $superset->push_attribute("subsets",$object_ec);
		    $ec_object->set_attribute("parent",$parent_ec);
		} else {
		    &ErrorMessage("Error: parent $parent_ec not defined for child $ec_number\n");
		}
	    } elsif ($ec_number =~ /^\d+$/) {
		    $ec_object->set_attribute("parent","<NULL>");
	    }
	}
    }
}


sub ParseReactionEC {
    my ($file) = @_;
    
    warn (";\n; ", &AlphaDate, " parsing reaction ECs from $file\n")
	if ($warn_level >= 1);
    
    open EC_REACT, $file || die "Error: cannot read ec2reaction file $file\n";
    
    while (<EC_REACT>) {
	chomp;
	if (/(\S+)\:(.+)\s*/) {
	    $reaction_id = $1;
	    $ec_number = $2;
#      $ec_number =~ s/\.\-//g; ### remove the .- used for partly specified EC
	    $ec_number =~ s/\.\s*$//g; ### remove trailing dot
	    $ec_number =~ s/^ +//g; ### remove leading spaces
	    $ec_number =~ s/ +$//g; ### remove trailing spaces
#print STDERR "trying to get EC number '", $ec_number, "' for reaction $reaction_id\n";
	    unless ($reaction_object = $reactions->get_object($reaction_id)) {
		&ErrorMessage(";WARNING: reaction '$reaction_id' has not been defined, cannot assign EC number $ec_number to it\n");
		next;
	    }
	    if ($ec_object = $ecs->get_object($ec_number)) {
		$ec_object->push_attribute("elements",$reaction_id);
	    } else {
		&ErrorMessage(";WARNING: could not identify the set '$ec_number' for reaction '$reaction_id'\n");
	    }
	} else {
	    &ErrorMessage(";WARNING: invalid reaction line\t$_\n");
	}
    }
    

    close EC_REACT;
}

#### calculate the number of reactions in which each compound is involved
#### and the number of substrates/products for each reactions
sub Connectivity {
    my ($file) = @_;
    warn (";\n; ", &AlphaDate, "; printing connectivity to file $file\n")
	if ($warn_level >= 1);

    open CONNECTIVITY, ">$file" 
	|| die "Error: cannot write file $file\n";
    
    my %reaction_connectivity;
    my %compound_connecivity;
    foreach my $reactant ($reactants->get_objects()) {
	
	my $compound_id = $reactant->get_attribute("compound_id");
	my $reaction_id = $reactant->get_attribute("reaction_id");
	my $valid = $reactant->get_attribute("valid_interm");
	my $reactant_type = $reactant->get_attribute("reactant_type");
	$reaction_connectivity{$reaction_id}->{$reactant_type}++;
	$compound_connectivity{$compound_id}->{$reactant_type}++;
	if ($valid > 0) {
	    $compound_connectivity{$compound_id}->{"valid_".$reactant_type}++;
	    $reaction_connectivity{$reaction_id}->{"valid_".$reactant_type}++;
	}
    }

    printf CONNECTIVITY (";%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			 "id",
			 "node.type",
			 "in",
			 "out",
			 "in.ok",
			 "out.ok",
			 "description",
			 );
    foreach my $compound_id (sort keys %compound_connectivity) {
	my $compound = $compounds->get_object($compound_id);
	my $compound_name = "";
	if ($compound) {
	    $compound_name = $compound->get_name();
	} 

#	warn $compound_id, "\t", $compound, "\t", $compound_name, "\n";
	
	printf CONNECTIVITY ("%s\t%s\t%d\t%d\t%d\t%d\t%s\n",
			     $compound_id,
			     "compound",
			     $compound_connectivity{$compound_id}->{"product"},
			     $compound_connectivity{$compound_id}->{"substrate"},
			     $compound_connectivity{$compound_id}->{"valid_product"},
			     $compound_connectivity{$compound_id}->{"valid_substrate"},
			     $compound_name,
			     );
    }

    foreach my $reaction_id (sort keys %reaction_connectivity) {
	my $reaction = $reactions->get_object($reaction_id);
	my $reaction_equation = "";
	if ($reaction) {
	    $reaction_equation = $reaction->get_attribute("equation");
	} 
	printf CONNECTIVITY ("%s\t%s\t%d\t%d\t%d\t%d\t%s\n",
			     $reaction_id,
			     "reaction",
			     $reaction_connectivity{$reaction_id}->{"substrate"},
			     $reaction_connectivity{$reaction_id}->{"product"},
			     $reaction_connectivity{$reaction_id}->{"valid_substrate"},
			     $reaction_connectivity{$reaction_id}->{"valid_product"},
			     $reaction_equation
			     );
    }

    close CONNECTIVITY;

}

#  sub AppendECstats {
#      my ($file) = @_;
#      my @roots = qw( 1 2 3 4 5 6 0.0.0.0 ); ### for tree expansion

#      if ($file) {
#  	open STDOUT, ">>$out_file{stats}" 
#  	    || die "Error: cannot write $file\n";
#      }
    
#      ### WARNING : there is a problem: I have more objects in 
#      ### total than expected. TO CHECK
#      if ($warn_level >= 1) {
#  	warn (";\n; ", &AlphaDate, " writing statistics\n");
#  	$total = 0;
#  	foreach $root (@roots) {
#  	    #warn "; root : $root\n";
#  	    if ($object = PFBP::ECSet->get_object($root)) {
#  		if ($count = $object->count_elements()) {
#  		    printf STDERR "; ECset %10s\t%6s reactions\n", $root, $count;
#  		    $total += $count;
#  		}
#  	    }
#  	}
#  	warn "; total\t$total\n";
#      }
    
#      foreach $root (@roots) {
#  	if ($object = PFBP::ECSet->get_object($root)) {
#  	    $object->expand(count_elements=>1);
#  	}
	
#  	close STDOUT if ($file);
#      }
#  }
    
    

sub ParseKeggPathways {
    ### full names of organisms
    $full_name{eco} = "Escherichia coli";
    $full_name{sce} = "Saccharomyces cerevisiae";
    $full_name{hsa} = "Homo sapiens";

    ### read pathway names
    %pathway_names = ();
    if (-r $in_file{pathway_names}) {
	open NAMES, $in_file{pathway_names};
	while (<NAMES>) {
	    chomp;
	    my ($map_id, $name) = split "\t";
	    $pathway_names{"map".$map_id} = $name;
	}
	close NAMES;
    } else {
	&ErrorMessage("Cannot read the file with pathway names '",
		      $in_file{pathway_names},
		      "'\n");
    }
    
    
    ### input files
    @in_files = glob($dir{pathway_reactions}."/map*.rea*");
    if ($test) {
	warn ";TEST\n" if ($warn_level >= 1);
	### fast partial parsing for debugging
	@in_files = shift @in_files;
    }
    
    ################################################################
    # parse generic pathways (lists of reactions and ECs)
    #
    foreach $filename (@in_files) {
	$short_file_name = $filename;
	$short_file_name =~ s|.*/||g;
	if ($short_file_name =~ /(map\d+)\.rea/) {
	    $map_id = $1;
	}

	$pathway = $genericPathways->new_object(id=>$map_id,
						source=>"KEGG:$short_file_name");
	if (defined($pathway_names{$map_id})) {
	    $pathway->new_attribute_value("names",$pathway_names{$map_id});
	} else {
	    $pathway->new_attribute_value("names",<NULL>);
	    &ErrorMessage("WARNING: no name for pathway $map_id\n");
	}
	
	warn ("; reading reactions",
	      "\tfile: ", $short_file_name,
	      "\tid: ", $pathway->get_attribute("id"),
	      "\tname: ", $pathway->get_name(),
	      "\n")
	    if ($warn_level >= 1);
	&ParsePathwayReactions($pathway,$filename);
    }

    ################################################################
    # parse organism-specific pathways (lists of genes)
    #
    foreach $org (keys %full_name) {
	$organism = $full_name{$org} || $org;
	@in_files = glob($dir{$org}."/*.gene*");
	@in_files = shift @in_files if ($test);
	
	foreach $filename (@in_files) {
	    $short_file_name = $filename;
	    $short_file_name =~ s|.*/||g;
	    if ($short_file_name =~ /(\d+)\.gene/) {
		$map_id = "map".$1;
	    }
	    
	    $pathway = $pathways->new_object(id=>$short_file_name,
					     source=>"KEGG:$short_file_name");
	    
	    ### include the generic pathway as a sub-pathway
	    if ($sub_pathway = $genericPathways->get_object($map_id)) {
		$parent_id = $sub_pathway->get_attribute("id");
		$pathway->set_attribute("parent", $parent_id);
	    } else {
		&ErrorMessage("Cannot identify generic pathway for id $map_id\n");
		$pathway->set_attribute("parent","<NULL>");
	    }

	    ### pathway name
	    if (defined($pathway_names{$map_id})) {
		$pathway->new_attribute_value("names",$pathway_names{$map_id});
	    } else {
		$pathway->new_attribute_value("names","<NULL>");
		&ErrorMessage("WARNING: no name for pathway $map_id\n");
	    }
	    $pathway->new_attribute_value("organism",$organism);

	    warn (";reading genes",
		  "\tfile: ", $short_file_name,
		  "\tid: ", $pathway->get_attribute("id"),
		  "\tname: ", $pathway->get_name(),
		  "\n")
		if ($warn_level >= 1);
	    &ParseGenes($pathway,$filename);
	}
	
    }
}

sub ParsePathwayReactions {
  my ($pathway, $input_file) = @_;
  
  if ($input_file =~ /\.gz$/) {
    $input_stream = "gunzip -c $input_file |";
  }

  open PATH_REACT, $input_file ||
      die "Error: cannot read file $input_file\n";
  my %pathway_reaction = ();
  my %pathway_ec = ();
  while ($line = <PATH_REACT>) {
    chomp $line;
    my ($reaction_id, $ec_id, $equation)  = split ":", $line;
    $ec_id = &trim($ec_id);
    $ec_id =~ s/\.$/\.-/;
    
    $pathway_reaction{$reaction_id}++;
    $pathway_ec{$ec_id}++;
  }

  #### assign EC numbers as elements
  foreach $ec_id (sort keys %pathway_ec) {
      if ($ecs->get_object($ec_id)) {
	  $pathway->new_attribute_value("ECs", $ec_id);
      } else {
	  &ErrorMessage("Error in file $input_file: unknown EC number '$ec_id'\n$line\n");
      }
  }
  
  #### assign reactions as elements
  foreach $reaction_id (sort keys %pathway_reaction) {
      if ($reactions->get_object($reaction_id)) {
	  $pathway->new_attribute_value("reactions", $reaction_id);
      } else {
	  &ErrorMessage("Error in file $input_file: unknown reaction '$reaction_id'\n");
	  next;
      }
  }

  close PATH_REACT;
  return;
}


sub ParseGenes {
  my ($pathway, $input_file) = @_;
  
  if ($input_file =~ /\.gz$/) {
    $input_file = "gunzip -c $input_file |";
  }
  open GENES, $input_file ||
      die "Error: cannot read file $input_file\n";
  while ($line = <GENES>) {
    chomp $line;
#    my $expr_pointer = undef;
    my ($Gene, $description)  = split "\t", $line;
    
    if ($org eq "hsa") {
      $Gene = "KEGG_gn_".$Gene."_CDS_H.sapiens";
    }

#    my $gene_id = undef;
#    unless ($gene_id = $database->get_id("PFBP::Gene::".$Gene)) {
#      &ErrorMessage("line $l file $in_file{pathway_index}: unknown gene $Gene\n", $line, "\n");
#      next;
#    }
#    unless (($expr_pointer) = $database->get_objects_by_input("PFBP::Expression","PFBP::Gene::".$gene_id)) {
#      &ErrorMessage("line $l file $in_file{pathway_index}: no expression for gene $Gene\n");
#      next;
#    }
    
#    my $expression = $database->get_object($expr_pointer);
#    $pathway->new_attribute_value(elements=>"PFBP::Expression::".$expression->get_attribute("id"));
    $pathway->new_attribute_value(genes=>$Gene);
    
#    $pathway->new_attribute_value("elements",$gene_id);
  }

  close PATH_REACT;
  return;
}

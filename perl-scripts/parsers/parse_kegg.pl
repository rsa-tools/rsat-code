#!/usr/bin/perl
############################################################
#
# $Id: parse_kegg.pl,v 1.19 2011/02/17 05:07:46 rsat Exp $
#
# Time-stamp: <2003-07-10 11:53:00 jvanheld>
#
############################################################

################################################################
## THIS SCRIPT IS NOT FUNCTIONAL ANYMORE? IT HAS TO BE CHECKED

### type parse_kegg.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "lib/load_classes.pl";
require "config.pl";
require "lib/util.pl";
require "lib/parsing_util.pl";

package KEGG::GenericPathway;
### A class to treat EC numenclature
{
  @ISA = qw ( classes::ObjectSet );
  ### class attributes
  $_count = 0;
  $_prefix = "pthw_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
		      names=>"ARRAY",
		      reactions=>"EXPANDED",
		      ECs=>"ARRAY",
		     );
}

package KEGG::Pathway;
### A class to treat EC numenclature
{
  @ISA = qw ( classes::ObjectSet );
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

package KEGG::Reaction;
{
  @ISA = qw ( classes::HasReactants );
  ### class attributes
  $_count = 0;
  $_prefix = "rctn_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     equation=>"SCALAR",
			     description=>"SCALAR",
			     enzyme=>"SCALAR",
			     pathway=>"ARRAY",
#		      substrates=>"ARRAY",
#		      products=>"ARRAY",
			     source=>"SCALAR");
}


### reactant is just a class to store a hash attribute for Reaction
### it describes 
### - the compound ID, 
### - the stoichiometry 
### - the validity of the reactant as intermediate betweeen 2 successive reactions in a pathway
package KEGG::Reactant;
{
  @ISA = qw ( classes::DatabaseObject );
  ### class attributes
  $_count = 0;
  $_prefix = "rctt_";
  @_objects = ();
  %_name_index = ();
  %_id_index = ();
  %_attribute_count = ();
  %_attribute_cardinality = (id=>"SCALAR",
			     compound_id=>"SCALAR",
			     stoichio=>"SCALAR",
			     valid_interm=>"SCALAR"
			     );
}

package main;

################################################################
# initialisation

#### oracle schema
$schema = "kegg";
$user = "kegg";
$password = "kegg";


################################################################
#
# files to parse
#

#### ligand
$dir{KEGG} = "${Databases}/ftp.genome.ad.jp/pub/kegg";
$data_file{compound} = $dir{KEGG}."/ligand/compound";
$data_file{reaction} = $dir{KEGG}."/ligand/reaction";
#$data_file{reaction} = $dir{KEGG}."/ligand/reaction.lst";
#$data_file{reaction_name} = $dir{KEGG}."/ligand/reaction_name.lst";
$data_file{ec} = $dir{KEGG}."/ligand/ECtable";

#### pathways
#$dir{pathway_reactions} = $dir{KEGG}."/ligand/reaction.main/";
$dir{eco} = $dir{KEGG}."/pathways/eco/";
$dir{sce} = $dir{KEGG}."/pathways/sce/";
$dir{hsa} = $dir{KEGG}."/pathways/hsa/";
$in_file{pathway_names} = $dir{KEGG}."/pathways/map_title.tab";

################################################################
#
# output files
#

#### default export directory
$export_subdir = "kegg_ligand";
$dir{output} = "$parsed_data/${export_subdir}/$delivery_date";

$out_format = "obj";

push @classes, ("classes::Compound");
push @classes, ("KEGG::Reaction");
push @classes, ("KEGG::Reactant");
push @classes, ("classes::ECSet");
push @classes, ("KEGG::GenericPathway");

&ReadArguments();


#### input directories
#unless (-d $dir{pathway_reactions}) {
#    warn "Uncompressing reactions from KEGG archive\n";
#    system "cd $dir{KEGG}/ligand/; uncompress -c reaction.main.tar.Z | tar -xvf -";
#}

#### output directory
&CheckOutputDir();

#### output files
$outfile{kegg} = "$dir{output}/kegg.obj";
$outfile{connectivity} = "$dir{output}/kegg.connectivity.txt";
$outfile{stats} = "$dir{output}/kegg.stats.txt";
$outfile{errors} = "$dir{output}/kegg.errors.txt";
open ERR, ">$outfile{errors}" || die "Error: cannot write error report fle $$outfile{errors}\n";

#### check input files
foreach $key (keys %data_file) {
    if (-e $data_file{$key}) {
	$in_file{$key} = $data_file{$key};
	if ($test) {
	    $in_file{$key} = "head -1000 $in_file{$key} | ";
	}
    } elsif (-e "$data_file{$key}.gz") {
	$in_file{$key} = "gunzip -c $data_file{$key}.gz | ";
	if ($test) {
	    $in_file{$key} = "| head -1000 ";
	}
    } elsif (-e "$data_file{$key}.Z") {
	$in_file{$key} = "uncompress -c $data_file{$key}.Z | ";
	if ($test) {
	    $in_file{$key} = "| head -1000 ";
	}
    } else {
	warn "Warning: $key data file $data_file{$key} does not exist\n";
    }
}

################################################################
#### KEGG reactions come in a compressed archive
$in_file{reaction2ec} = "uncompress -c ".$dir{KEGG}."/ligand/reaction.tar.Z | tar -xpOf - | cut -d ':' -f 1,2 | sort -u |";

### default verbose message
&DefaultVerbose() if ($verbose >= 1);

### instantiate class factories
$compounds = classes::ClassFactory->new_class(object_type=>"classes::Compound",
					      prefix=>"comp_");
$reactions = classes::ClassFactory->new_class(object_type=>"KEGG::Reaction",
					      prefix=>"rctn_");
$reactants = classes::ClassFactory->new_class(object_type=>"KEGG::Reactant",
					      prefix=>"rctt_");
$ecs = classes::ClassFactory->new_class(object_type=>"classes::ECSet",
					prefix=>"ec_");
$pathways = classes::ClassFactory->new_class(object_type=>"KEGG::Pathway",
					     prefix=>"pthw_");
$genericPathways = classes::ClassFactory->new_class(object_type=>"KEGG::GenericPathway",
						    prefix=>"gptw_");

#  $compounds->generate_sql(schema=>$schema,
#  			 user=>$user,
#  			 password=>$password,
#  #				  dir=>"$dir{output}/sql_scripts",
#  			 prefix=>"$class_factory_",
#  			 dbms=>$dbms
#  			 );

### default output fields for each class
$compounds->set_out_fields(qw( id names formula description source ));
$reactions->set_out_fields(qw( id  source equation description ecs pathways
			       name pathway enzyme  left right)); ## for parsing only
$reactants->set_out_fields(qw( id reactant_type reaction_id compound_id stoichio valid_interm ));
$ecs->set_out_fields(qw( id parent description names  ));
$pathways->set_out_fields(qw( id parent organism source names reactions ECs genes ));
$genericPathways->set_out_fields(qw( id source names reactions ECs ));

$genericPathways->set_attribute_header("reactions", join ("\t", "reaction", "direction") );


#### parse compounds
&ParseKeggFile($in_file{compound}, $compounds, source=>'KEGG:compound');
if ($in_file{additional_compounds}) {
    &ParseKeggFile($in_file{additional_compounds}, $compounds, source=>"amaze");
}


$compounds->index_names();
foreach $compound ($compounds->get_objects()) {
    #### check for compounds without formula
    if ($compound->get_attribute("formula") eq "") {
	$compound->set_attribute("formula",$null);
    }

    #### use the formula as dscription
    $compound->set_attribute("description", $compound->get_attribute("formula"));
}


##### parse reactions
&ParseReactions($in_file{reaction}, $reactions);

$reactions->set_out_fields(qw( id  source equation description ecs pathways)); ### remove temporary fields, which were necessary for parsing but not for export

#### parse EC numbers
&ParseEC($in_file{ec}, $ecs);

#### define a description for EC numbers
foreach my $ec ($ecs->get_objects()) {
    my $name = "";
    my @names = $ec->get_attribute("names");
    if ($#names > 0) {
	$name = $names[1];
    } else {
	$name = $names[0];
    }
    $ec->set_attribute('description', $name)
}

#### parse pathways
&ParseKeggPathways();



################################################################
### export parsing result
@class_factories = qw (
		       compounds
		       reactions
		       reactants
		       ecs
		       pathways
		       genericPathways
		       );
foreach $class_factory (@class_factories) {
    $$class_factory->dump_tables();
    $$class_factory->generate_sql(schema=>$schema,
				  user=>$user,
				  password=>$password,
				  prefix=>"$class_factory_",
				  dbms=>$dbms
				  );
}

&ExportClasses($outfile{kegg}, $out_format, @classes) if ($export{obj});


### print some stats after parsing
&PrintStats($outfile{stats}, @classes);
&Connectivity($outfile{connectivity});
#&AppendECstats($outfile{stats});

### report execution time
if ($verbose >= 1) {
    $done_time = `date +%Y-%m-%d.%H%M%S`;
    warn ";\n";
    warn "; job started $start_time";
    warn "; job done    $done_time";
}

close ERR;

&CompressParsedData();

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
	Parse compounds from KEGG (http://www.genome.ad.jp/kegg/). 
	
AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

OPTIONS	
$generic_option_message
	-add_react additional_reaction_file
	       Specify a file containing additional reactions. This
	       file must be exactly in the same format as the flat
	       file containing KEGG reactions.

	-add_comp additional_compound_file
	       Specify a file containing additional compounds. This
	       file must be exactly in the same format as the flat
	       file containing KEGG compounds.

EXAMPLE
	perl parse_kegg.pl -v 1 -clean -dbms postgresql -schema atest

EndHelp
  close HELP;
}

  


################################################################
### read arguments from the command line
sub ReadArguments {
    my $a = "";
    for $a (0..$#ARGV) {
	&ReadGenericOptions($a);

	### Manual specification of the reaction file reactions
	if ($ARGV[$a] eq "-react") {
	    $main::data_file{reactions} = $ARGV[$a+1];
	}

	### Manual specification of the compounds
	if ($ARGV[$a] eq "-comp") {
	    $main::data_file{compounds} = $ARGV[$a+1];
	}

	### additional reactions
	if ($ARGV[$a] eq "-add_react") {
	    $main::in_file{additional_reactions} = $ARGV[$a+1];
	}

	### additional compounds
	if ($ARGV[$a] eq "-add_comp") {
	    $main::in_file{additional_compounds} = $ARGV[$a+1];
	}
    }
}


################################################################
#### define a list of compounds which cannot be used as intermediate between two reactions
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

################################################################
#### parse reactions
sub ParseReactions {
    ### read the reaction file from KEGG
    my ($reaction_file, $class_holder, $source) = @_;
    $source = $source || "KEGG:reaction";

    &TrivialCompounds();
    
#     warn (";\n; ", 
# 	  &AlphaDate, 
# 	  "\tparsing class ", 
# 	  $class_holder->get_object_type(),
# 	  " from ",
# 	  $reaction_file,
# 	  "\n")     if ($verbose >= 1);

    &ParseKeggFile($reaction_file,$class_holder, source=>$source);
    if ($in_file{additional_reactions}) {
	&ParseKeggFile($in_file{additional_reactions}, $reactions, source=>$source);
    }

    ################################################################
    #### post-processing of KEGG information to better fit with our objects
    foreach my $reaction ($reactions->get_objects()) {

	my $reaction_id = $reaction->get_attribute("id");

	################################################################
	#### parse reaction equations to create reactants
	my $equation = $reaction->get_attribute("equation");
	if ($equation =~ /^\s*(.*)\s*\<\=\>\s*(.*)\s*$/) {
	    $left = $1;
	    $right = $2;
	    
	    ### temporary
	    $reaction->set_attribute("left",$left);
	    $reaction->set_attribute("right",$right);
	    
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
		}
	    }
	} else {
	    &ErrorMessage("invalid equation for reaction $reaction_id\t$equation");
	}

	################################################################
	#### parse pathways in which each reaction is involved
	foreach my $pathway ($reaction->get_attribute("pathway")) {
	    if ($pathway =~ /PATH: (MAP\d+)/) {
		my $pathway_id = $1;
		$reaction->push_attribute("pathways", $pathway_id);
	    } else {
		&ErrorMessage("reaction $reaction_id\tinvalid pathway\t$pathway");
	    }
	}


	################################################################
	#### parse ECs associated to each reaction
	foreach my $enzyme_string ($reaction->get_attribute("enzyme")) {

	    if ($enzyme_string eq $null) {
		warn "; Reaction $reaction_id has no enzyme\n" if ($verbose >= 2);
		next;
	    }
	    my @ecs = split /\s+/, $enzyme_string;
	    foreach my $ec (@ecs) {
		$reaction->push_attribute("ecs", $ec);
	    }
	}

    }
}


################################################################
#### Parse EC numbers
sub ParseEC {
    my ($file, $class_holder, %args) = @_;
    warn (";\n; ", &AlphaDate, "; parsing class ",
	  $class_holder->get_object_type(),
	  " from $file\n")
	if ($verbose >= 1);
    
    ### create the class for 0.0.0.0
    $ec_object = $ecs->new_object(%args,id=>"0.0.0.0");
#    $ec_object = $ecs->new_object(%args);
    $ec_object->push_attribute("names", "0.0.0.0");
    $ec_object->push_attribute("names", "non-enzymatic or not clearly enzymatic");
#    $ec_object->set_attribute("ec", "0.0.0.0");
    $ec_object->set_attribute("parent",$null);

    ### create the class for -.-.-.-
    $ec_object = $ecs->new_object(%args,id=>"-.-.-.-");
#    $ec_object = $ecs->new_object(%args);
    $ec_object->push_attribute("names", "-.-.-.-");
    $ec_object->push_attribute("names", "non-assigned EC number");
#    $ec_object->set_attribute("ec", "-.-.-.-");
    $ec_object->set_attribute("parent",$null);

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
	    $names = "$'";
	    
	    $ec_number =~ s/\.$//; ### suppress the trailing dot from the ec number
	    @names = split /\; /, $names;
	    
	    $ec_object = $ecs->new_object(%args,id=>$ec_number);
#	    $ec_object = $ecs->new_object(%args);
	    foreach $name (@names) {
		$ec_object->push_attribute("names",$ec_number);
		$ec_object->push_attribute("names",$name);
#		$ec_object->force_attribute("ec",$ec_number);
	    }
	    
	    ### partially specified ecs
	    unless ($ec_number =~ /^\d+\.\d+\.\d+\.\d+$/) {
		if ($ec_number =~ /^\d+\.\d+\.\d+$/) { 
		    $partial_ec = $ec_number.".-";
		} elsif ($ec_number =~ /^\d+\.\d+$/) { 
		    $partial_ec = $ec_number.".-.-";
		} elsif ($ec_number =~ /^\d+$/) { 
		    $partial_ec = $ec_number.".-.-.-";
		} else {
		    $partial_ec = undef;
		    &Errormessage("Error: non-standard EC number\n");
		}
		if ($partial_ec) {
		    $uncomplete_object = $ecs->new_object(%args,id=>$partial_ec);
#		    $uncomplete_object = $ecs->new_object(%args);
		    foreach $name (@names) {
			$uncomplete_object->push_attribute("names",$partial_ec);
			$uncomplete_object->push_attribute("names",$name);
		    }
#		    $ec_object->push_attribute("subsets",$partial_ec);
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
		    $ec_object->set_attribute("parent",$null);
	    }
	}
    }
}


################################################################
#### calculate the number of reactions in which each compound is
#### involved and the number of substrates/products for each 
#### reaction
sub Connectivity {
    my ($file) = @_;
    warn (";\n; ", &AlphaDate, "; printing connectivity to file $file\n")
	if ($verbose >= 1);

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

################################################################
#### Parse pathways
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

    #### generic pathways
#    &ParseGenericPathways();

    #### organism-specific pathways
    &ParseSpecificPathways(); 
}

################################################################
#### parse links between reactions and pathways
sub ParsePathwayReactions {
    my ($pathway, $input_file) = @_;
    
    if ($input_file =~ /\.gz$/) {
	$input_stream = "gunzip -c $input_file |";
    } else {
	$input_stream = $input_file;
    }

    open PATH_REACT, $input_stream ||
	die "Error: cannot read file $input_file\n";
    my %pathway_reaction = ();
    my %pathway_ec = ();
    while ($line = <PATH_REACT>) {
	my $reaction_direction = $null;
	chomp $line;
	my ($reaction_id, $ec_id, $equation)  = split ":", $line;
	$ec_id = &trim($ec_id);
	$ec_id =~ s/\.$/\.-/;
	
	if ($equation =~ / <=> /) {
	    $reaction_direction = "both";
	} elsif ($equation =~ / => /) {
	    $reaction_direction = "forward";
	} elsif ($equation =~ / <= /) {
	    $reaction_direction = "reverse";
	}
	$pathway_reaction{$reaction_id} = $reaction_direction;
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
	my $reaction_direction = $pathway_reaction{$reaction_id};
	if ($reactions->get_object($reaction_id)) {
	    $pathway->push_expanded_attribute("reactions", $reaction_id, $reaction_direction);
	} else {
	    &ErrorMessage("Error in file $input_file: unknown reaction '$reaction_id'\n");
	    next;
	}
    }

    close PATH_REACT;
    return;
}

################################################################
#### parse links between genes
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
#    unless ($gene_id = $database->get_id("classes::Gene::".$Gene)) {
#      &ErrorMessage("line $l file $in_file{pathway_index}: unknown gene $Gene\n", $line, "\n");
#      next;
#    }
#    unless (($expr_pointer) = $database->get_objects_by_input("classes::Expression","classes::Gene::".$gene_id)) {
#      &ErrorMessage("line $l file $in_file{pathway_index}: no expression for gene $Gene\n");
#      next;
#    }
	
#    my $expression = $database->get_object($expr_pointer);
#    $pathway->new_attribute_value(elements=>"classes::Expression::".$expression->get_attribute("id"));
	$pathway->new_attribute_value(genes=>$Gene);
	
#    $pathway->new_attribute_value("elements",$gene_id);
    }

    close PATH_REACT;
    return;
}


################################################################
# parse organism-specific pathways (lists of genes)
#
sub ParseSpecificPathways {
  foreach my $org (keys %full_name) {
    warn "; Parsing specific pathways for organism $org\n" if ($verbose >=1); 
    
    unless (-d $dir{$org}) {
      die "Error: directory $dir{$org} does not exist\n";
    }
    my $organism = $full_name{$org} || $org;
    my @in_files = glob($dir{$org}."/*.gene*");
    @in_files = shift @in_files if ($test);
    
    warn join ("\n;\t", "; List of files to parse", @in_files), "\n" if ($verbose >=2); 
    
    foreach my $filename (@in_files) {
      my $short_file_name = $filename;
      $short_file_name =~ s|.*/||g;
      if ($short_file_name =~ /(\d+)\.gene/) {
	$map_id = "map".$1;
      }
      
      #### create a new pathway
      $pathway = $pathways->new_object(id=>$short_file_name,
				       source=>"KEGG:$short_file_name");
      
      ### include the generic pathway as a sub-pathway
      if ($sub_pathway = $genericPathways->get_object($map_id)) {
	$parent_id = $sub_pathway->get_attribute("id");
	$pathway->set_attribute("parent", $parent_id);
      } else {
	&ErrorMessage("Cannot identify generic pathway for id $map_id\n");
	$pathway->set_attribute("parent",$null);
      }
      
      ### pathway name
      if (defined($pathway_names{$map_id})) {
	$pathway->new_attribute_value("names",$pathway_names{$map_id});
      } else {
	$pathway->new_attribute_value("names",$null);
	&ErrorMessage("WARNING: no name for pathway $map_id\n");
      }
      $pathway->new_attribute_value("organism",$organism);
      
      warn ("; Reading genes",
	    "\tfrom file: ", $short_file_name,
	    "\tid: ", $pathway->get_attribute("id"),
	    "\tname: ", $pathway->get_name(),
	    "\n")
	if ($verbose >= 1);
      &ParseGenes($pathway,$filename);
    }
    
  }
}


################################################################
# parse generic pathways (lists of reactions and ECs)
#
sub ParseGenericPathways {
    warn "; Parsing generic pathways\n" if ($verbose >= 1);
    ### input files
    warn "; Getting generic pathways from $dir{pathway_reactions}\n" if ($verbose >= 1);
    @in_files = glob($dir{pathway_reactions}."/map*.rea*");
    if ($#in_files < 0) {
	die "Error: cannot find generic pathways in directory $dir{pathway_reactions}\n";
    }
    if ($test) {
	warn ";TEST\n" if ($verbose >= 1);
	### fast partial parsing for debugging
	@in_files = ($in_files[0], $in_files[1], $in_files[2]);
    }

    warn join ("\n;\t", "; List of files to parse", @in_files), "\n" if ($verbose >=2);
    foreach $filename (@in_files) {
	$short_file_name = $filename;
	$short_file_name =~ s|\.*/||g;
	if ($short_file_name =~ /(map\d+)\.rea/) {
	    $map_id = $1;
	}
	
	$pathway = $genericPathways->new_object(id=>$map_id,
						source=>"KEGG:$short_file_name");
	if (defined($pathway_names{$map_id})) {
	    $pathway->new_attribute_value("names",$pathway_names{$map_id});
	} else {
	    $pathway->new_attribute_value("names",$null);
	    &ErrorMessage("WARNING: no name for pathway $map_id\n");
	}
	
	warn ("; reading reactions",
	      "\tfile: ", $filename,
	      "\tfile: ", $short_file_name,
	      "\tid: ", $pathway->get_attribute("id"),
	      "\tname: ", $pathway->get_name(),
	      "\n")
	    if ($verbose >= 1);
	&ParsePathwayReactions($pathway,$filename);
    }
}

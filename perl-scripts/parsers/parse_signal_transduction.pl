#!/usr/bin/perl

### add the program's directory to the lib path
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`"); ### add the program's directory to the lib path
    }
}
require "PFBP_config.pl";
require "PFBP_classes.pl";
require "PFBP_parsing_util.pl";

#use strict;
#no strict "refs";
#no strict "vars";

package PFBP::MoleculeInState;
{
    @ISA = qw ( PFBP::BiochemicalEntity );
    %_attribute_cardinality = (state=>"SCALAR",
			       parent=>"SCALAR");
}

package PFBP::Interaction;
{
    @ISA = qw ( PFBP::DatabaseObject );
    %_attribute_cardinality = (
			       id=>"SCALAR",
			       type=>"SCALAR",
			       description=>"SCALAR",
			       inputes=>"EXPANDED",
			       outputs=>"EXPANDED"
			       );

    #### collect IDs of the input objects  by taking the first field of the inputs expanded attribute
    sub get_input_IDs {
	my ($self) = @_;
	my @objects = ();
	my @reactants = $self->get_attribute("inputs");
	foreach my $r (0..$#reactants) {
	    @reactant = @{$reactants[$r]};
	    push @objects, @reactant[0];
	}
	return @objects;
    }

    #### collect IDs of the output objects  by taking the first field of the outputs expanded attribute
    sub get_output_IDs {
	my ($self) = @_;
	my @objects = ();
	my @reactants = $self->get_attribute("outputs");
	foreach my $r (0..$#reactants) {
	    @reactant = @{$reactants[$r]};
	    push @objects, @reactant[0];
	}
	return @objects;
    }
}


package main;
{

    ### initialization
    $null = "<NULL>";
    $start_time = `date +%Y-%m-%d.%H%M%S`;
    $clean = 1;
    
    $out_format = "obj";
    @classes = qw( PFBP::BiochemicalEntity PFBP::MoleculeInState PFBP::Interaction PFBP::Pathway );
    
    #### class factory for entities
    $entities = PFBP::ClassFactory->new_class(object_type=>"PFBP::BiochemicalEntity",
					      prefix=>"ent_");
    $entities->set_out_fields(qw( id type primary_name names states ));

    #### class for entities in state
    $entities_in_state = PFBP::ClassFactory->new_class(object_type=>"PFBP::MoleculeInState",
						       prefix=>"ens_");
    $entities_in_state->set_out_fields(qw( id type primary_name names state location parent ));

    #### class factory for interactions
    $interactions = PFBP::ClassFactory->new_class(object_type=>"PFBP::Interaction",
						  prefix=>"int_");
    $interactions->set_attribute_header("inputs", "object\tsubunit\tstate\tstoichiometry\tlocation" );
    $interactions->set_attribute_header("outputs", "object\tsubunit\tstate\tstoichiometry\tlocation" );
    $interactions->set_out_fields(qw( id type description inputs outputs  ));


    #### class factory  for pathways
    $pathways = PFBP::ClassFactory->new_class(object_type=>"PFBP::Pathway",
					      prefix=>"pth_");
    $pathways->set_out_fields(qw( id names description entities interactions subpathways ));

    #### class factory  for diagrams
    $diagrams = PFBP::ClassFactory->new_class(object_type=>"PFBP::PathwayDiagram",
					     prefix=>"dgm_");

    ### old initialisation
    $path_element_count = 0;
    $Step = 0;

    ### directory containing the input files ###
    $dir{input} = "/win/amaze/amaze_team/Sandra/export_20011208";
    unless (-d $dir{input}) {
	die "Error : input dir $dir{input} does not exist\n";
    }

    ### entities (entities)
    $in_file{entities} = $dir{input}."/description_of_molecules.txt";
    $in_file{molecule_state} = $dir{input}."/molecule_state.txt";
    $in_file{synonyms} = $dir{input}."/synonyms.txt";

    ### interactions (interactions)
    $in_file{interactions} = $dir{input}."/description_of_interactions.txt";
    $in_file{interaction_input} = $dir{input}."/interaction_source.txt";
    $in_file{interaction_output} = $dir{input}."/interaction_target.txt";

    ### pathways
    $in_file{pathway_description} = $dir{input}."/pathway_description.txt";
    $in_file{pathway_subpathway} = $dir{input}."/pathway_subpathway.txt";
    $in_file{pathway_interaction} = $dir{input}."/pathway_interaction.txt";
    $in_file{pathway_entity} = $dir{input}."/pathway_molecule.txt";

    #### check for the existence of all the input files
    foreach my $file (keys %in_file) {
	&checkfile($in_file{$file});
    }

    ### output dir
    $dir{output} = "${parsed_data}/signal_transduction/$delivery_date";
    unless (-d $dir{output}) {
	warn "Creating output dir $dir{output}\n"  if ($warn_level >= 1);
	`mkdir -p $dir{output}`;
	die "Error: cannot create directory $dir\n" 
	    unless (-d $dir{output});
    }
    if ($clean) {
	warn "Cleaning output directory\n" if ($warn_level >= 1);
	system "\\rm -rf $dir{output}/*";
    }

    ### diagram dir
    $dir{diagrams} = $dir{output}."/diagrams";
    unless (-d $dir{diagrams}) {
	warn "Creating diagram dir $dir{diagrams}\n"  if ($warn_level >= 1);
	`mkdir -p $dir{diagrams}`;
	die "Error: cannot create directory $dir\n" 
	    unless (-d $dir{diagrams});
    }

    chdir $dir{output};

    ### object files
    $out_file{signal_transduction} = $dir{output}."/signal_transduction.obj";
    $out_file{graph} = $dir{diagrams}."/sigtrans_graph.tdd";
    $out_file{stats} = $dir{output}."/sigtrans_parsing_stats";

    ### reports
    $out_file{errors} = $dir{output}."/sigtrans_parsing_errors.txt";

    &ReadArguments;

    &DefaultVerbose() if ($warn_level >= 1);

    open ERR, "> $out_file{errors}" || die "Error : cannot write file $out_file{errors}\n";

    &ReadEntities();
    &ReadMoleculeStates();
    &ReadInteractions();
    &CheckInteractions();
    &ReadPathways();

    ### print each pathway in a separate file
    foreach $pathway ($pathways->get_objects()) {
	&PathwayToDiagram($pathway);
    }

    #####################
    ### print objects ###
    #####################
    $entities->dump_tables();
    $entities_in_state->dump_tables();
    $interactions->dump_tables();
    $pathways->dump_tables();

    &ExportClasses($out_file{signal_transduction}, $out_format, @classes)  if ($export{obj});

    ###################
    ### print stats ###
    ###################
    &PrintStats($out_file{stats}, @classes);

    warn "; Done\t", `date` if ($warn_level >= 1);

    close ERR;
    exit(0);
}


############# SUBROUTINE DEFINITION ################

### read entity information
sub ReadEntities {
    warn "; Reading entities\n" if ($warn_level >= 1);
    undef %col;
    $col{swissprot_ac} = 0;
    $col{swissprot_id} = 1;
    $col{name} = 2;
    $col{type} = 3;
    $col{descr} = 4;
    $entity_nb = -1;

    ### read entity description
    open ENT, "$in_file{entities}" || die "Error : cannot read file $in_file{entities}\n";
    $header = <ENT>;
    $line_nb = 0;
    while (<ENT>) {
	$line_nb++;
	my @fields = &MySplit;
	warn "$_\n" if ($warn_level >= 4);

	### read the fields ###
	$name = $fields[$col{name}];
	$description = $fields[$col{descr}];
	
	### check entity type ###
	if ($fields[$col{type}] =~ /\S/) {
	    $type = lc($fields[$col{type}]);
	    $type =~ s/ /_/g;
	    #### set first letter to uppercase, for compatibility with aMAZE styles
	    $type = ucfirst($type);
	} else {
	    $type = "undef";
	    print ERR "$in_file{entities}, line $line_nb: entity $name has no type\n";
	    print ERR "\t$_\n";
	}
	$entityType{$type} = 1;
	
	### instantiate a new BiochemicalEntity
	my $entity = $entities->new_object();
	$ac = $entity->get_attribute("id");
	$entity->push_attribute("names",$name);
	$entity->set_attribute("description",$description);
	$entity->set_attribute("type",$type);
	warn (join ("\t", "new entity",
		    $entity->get_attribute("id"),
		    $entity->get_attribute("type"),
		    $entity->get_attribute("names")
		   ), "\n")
	    if ($warn_level >= 2);
	
	### remind accession number ###
	$lc_name = lc($entity{$ac}->{name});
	$AC{$lc_name} = $ac;
    } 
    close ENT;

    #### set a primary nme
    foreach my $entity ($entities->get_objects()) {
	$entity->set_attribute("primary_name", $entity->get_name());
    }

    ### create indexes
    $entities->index_ids();
    $entities->index_names();

    ### read alternative names
    undef %col;
    $col{name} = 0;
    $col{altName} = 1;
    open NAMES, "$in_file{synonyms}"  || 
	die "Eannot read file $MolAltNameFile\n";
    $header = <NAMES>;
    $line_nb = 0;
    while (<NAMES>) {
	$line_nb++;
	my @fields = &MySplit;
	
	$name = $fields[$col{name}];
	$altName = $fields[$col{altName}];

	### add new name to the object
	if ($object = $entities->get_object($name)) {
	    $object->push_attribute("names", $altName);
	} else {      
	    print  ERR "$in_file{synonyms} line $line_nb: $name had not been previously declared in $in_file{entities}\n";
	}

    }
    close NAMES;
}

#### import molecule states
sub ReadMoleculeStates {
    my $l = 0;
    open STATES, $in_file{molecule_state} ||
	die "Error: cannot read molecule state file $in_file{molecule_state}\n";
    my $header = <STATES>;
    while (<STATES>) {
	$l++;
	my @fields = &MySplit();
	my $molecule = $fields[0];
	my $state = $fields[1];
	warn ";\tMolecule state\t$molecule\t$state\n" if ($warn_level >= 2);
	if (my $entity = $entities->get_object($molecule)) {
#	    $entity->push_attribute("states", $state);
  	    unless (my $entity_in_state = $entities_in_state->get_object("$state $molecule")) {
  		my $entity_in_state = $entities_in_state->new_object();
  		$entity_in_state->push_attribute("names","$state $molecule");
  		$entity_in_state->push_attribute("names","$molecule");
  		$entity_in_state->set_attribute("state", $state);
  		$entity_in_state->set_attribute("parent", $entity->get_attribute("id"));
  	    }
	} else {
	    &ErrorMessage("$in_file{molecule_state}\t$l\tundefined molecule\t$molecule\n");
	}
    }
    close STATES

}


#### check consistency of interactions
sub CheckInteractions {
    foreach my $interaction ($interactions->get_objects()) {
	my @inputs = $interaction->get_attribute("inputs");
	my $id = $interaction->get_attribute("id");
	my $input_nb = $#inputs +1;
	if ($input_nb < 1) {
	    &ErrorMessage(join ("\t", 
				"$in_file{interactions}",
				"interaction", 
				$id,
				"has no input"), "\n");
	}

	my @outputs = $interaction->get_output_IDs();
	my $id = $interaction->get_attribute("id");
	my $output_nb = $#outputs +1;
	if ($output_nb < 1) {
	    &ErrorMessage(join ("\t", 
				"$in_file{interactions}",
				"interaction", 
				$id,
				"has no output"), "\n");
	}

	warn join ("\t", 
		   ";\tInteraction $id",
		   "inputs: $input_nb",
		   "outputs: $output_nb"
		   ), "\n" if ($warn_level >= 2);
    }
}

### read description of interact<ion ###
sub ReadInteractions {
    warn "; Reading interactions\n" if ($warn_level >= 1);
    undef %cnol;
    $col{ac} = 0;
    $col{type} = 1;
    $col{descr} = 2;
    $assoc_nb = -1;

    open INTERACT, "$in_file{interactions}"  || die "Error : cannot read file $in_file{interactions}\n";
    $header = <INTERACT>;
    $line_nb = 0;
    while (<INTERACT>) {
	$line_nb++;
	my $line = $_;
	my @fields = &MySplit;  
	
	$assoc_nb++;
	$ac = $fields[$col{ac}];
	$type = lc($fields[$col{type}]);
	$type =~ s/^\'//;
	$type =~ s/\'$//;
	$type =~ s/ /_/g;
	unless ($type =~ /\S/) {
	    $type = "undef";
	    &ErrorMessage( "$in_file{interactions} line $line_nb: interaction type is not specified\n");
	    &ErrorMessage( "\t$_\n");
	}
	#$interactionType{$type} = 1;
	
	$name = "$type";
	$description = $fields[$col{descr}];
	$input = "";
	$output = "";
	
	my $interaction = $interactions->new_object(id=>$ac);
	$interaction->set_attribute("type",$type);
	$interaction->set_attribute("description",$description);
 	warn (join ("\t", "new interaction",
		    $interaction->get_attribute("id"),
		    $interaction->get_attribute("type")
		   ), "\n")
	    if ($warn_level >= 2);
	

	&CreateAssociation($ac,$type,$name,$description,$input,$output);
    }
    close INTERACT;

    ### create indexes
    $interactions->index_ids();
    $interactions->index_names();


    ### read interaction inputs
    warn "; Reading interaction inputs\n" if ($warn_level >= 1);

    ### initialisation
    undef %col;
    $col{inter} = 0;
    $col{molec} = 1;
    $col{subunit} = 2;
    $col{state_before} = 3;
    $col{state_after} = 4;
    $col{stoeichiometry} = 5;
    $col{location_before} = 6;
    $col{location_after} = 7;

    open INPUT, "$in_file{interaction_input}"   || die "Error : cannot read file $in_file{interaction_input}\n";
    $header = <INPUT>;
    $line_nb = 0;
    while (<INPUT>) {
	$line_nb++;
	my @fields = &MySplit;  
	$line = $_;
	
	my $inter = $fields[$col{inter}];
	my $molec = $fields[$col{molec}];
	my $subunit = $fields[$col{subunit}] || $null;
	my $state_before = $fields[$col{state_before}] || $null;
	my $state_after = $fields[$col{state_after}] || $null;
	my $location_before = $fields[$col{location_before}] || $null;
	my $location_after = $fields[$col{location_after}] || $null;
	
	### identify the interaction
	my $inter_object = $interactions->get_object($inter);
	unless ($inter_object) {
	    #unless (defined($association{$inter})) {
	    &ErrorMessage( "$in_file{interaction_input} line $line_nb: $inter has not been defined in $InteractionFile\n");
	    &ErrorMessage( "\t$_\n");
	    next;
	}

	### identify the entity
	if ($entity_object = $entities->get_object($molec)) {
	    $ent_id = $entity_object->get_attribute("id");
	} else {
	    &ErrorMessage( "$in_file{interaction_input} line $line_nb: $molecAC has not been defined in $in_file{entities}\n");
	    &ErrorMessage( "\t$_\n");
	    next;
	}


	#### TEMPORARY : for the time being, if before or after is missing, replace it by after or before resp.
	if (($state_before ne $null) && ($state_after eq $null)) {
	    $state_after = $state_before;
	} elsif (($state_after ne $null) && ($state_before eq $null)) {
	    $state_before = $state_after;
	}
	if (($location_before ne $null) && ($location_after eq $null)) {
	    $location_after = $location_before;
	} elsif (($location_after ne $null) && ($location_before eq $null)) {
	    $location_before = $location_after;
	}

	#### compare state and location before and after the interaction
	if (($state_before eq $state_after) && ($location_before eq $location_after)) {
	    $inter_object->push_expanded_attribute("inputs",$ent_id, $subunit, $state, $stoichiometry, $location);
	} else {
	    #### before remains input but after becomes output
	    $inter_object->push_expanded_attribute("inputs",$ent_id, $subunit, $state_before, $stoichiometry, $location_before);
	    $inter_object->push_expanded_attribute("outputs",$ent_id, $subunit, $state_after, $stoichiometry, $location_after);

	    #### check the reason why location before and after differ
	    unless ($location_before eq $location_after) {
		if ($inter_object->get_attribute("type") eq "translocation") {
		    #### specific treatment for translocations : there is no "target" attribute, but the output is in the "location_after" of the "source"
		    warn "Translocation treated\n" if ($warn_level >= 2);
		} else {
		    #### report unexpected cases
		    &ErrorMessage ("Unknown location modification for interaction input", 
				   "\t", $inter_object->get_attribute("type"),
				   "\t", $location_before, 
				   "\t", $location_after, 
				   "\n", $line, "\n");
		}
	    }

	    #### check the reason why state before and after differ
	    unless ($state_before eq $state_after) {
		if ($inter_object->get_attribute("type") eq "HELLO") {
		} else {
		    #### report unexpected cases
		    &ErrorMessage ("Unknown state modification for interaction input", 
				   "\t", $inter_object->get_attribute("type"),
				   "\t", $state_before, 
				   "\t", $state_after, 
				   "\n", $line, "\n");
		}
	    }
	}
    }
    close INPUT;

    open OUTPUT, "$in_file{interaction_output}"   || die "Error : cannot read file $in_file{interaction_output}\n";
    $header = <OUTPUT>;
    $line_nb = 0;
    while (<OUTPUT>) {
	$line_nb++;
	my @fields = &MySplit;  
	$line = $_;
	
	my $inter = $fields[$col{inter}];
	my $molec = $fields[$col{molec}];
	my $subunit = $fields[$col{subunit}] || $null;
	my $state_before = $fields[$col{state_before}] || $null;
	my $state_after = $fields[$col{state_after}] || $null;
	my $location_before = $fields[$col{location_before}] || $null;
	my $location_after = $fields[$col{location_after}] || $null;
	
	### identify the interaction
	my $inter_object = $interactions->get_object($inter);
	unless ($inter_object) {
	    #unless (defined($association{$inter})) {
	    &ErrorMessage( "$in_file{interaction_output} line $line_nb: $inter has not been defined in $InteractionFile\n");
	    &ErrorMessage( "\t$_\n");
	    next;
	}

	### identify the entity
	if ($entity_object = $entities->get_object($molec)) {
	    $ent_id = $entity_object->get_attribute("id");
	} else {
	    &ErrorMessage( "$in_file{interaction_output} line $line_nb: $molecAC has not been defined in $in_file{entities}\n");
	    &ErrorMessage( "\t$_\n");
	    next;
	}


	#### TEMPORARY : for the time being, if before or after is missing, replace it by after or before resp.
	if (($state_before ne $null) && ($state_after eq $null)) {
	    $state_after = $state_before;
	} elsif (($state_after ne $null) && ($state_before eq $null)) {
	    $state_before = $state_after;
	}
	if (($location_before ne $null) && ($location_after eq $null)) {
	    $location_after = $location_before;
	} elsif (($location_after ne $null) && ($location_before eq $null)) {
	    $location_before = $location_after;
	}

	#### compare state and location before and after the interaction
	if (($state_before eq $state_after) && ($location_before eq $location_after)) {
	    $inter_object->push_expanded_attribute("outputs",$ent_id, $subunit, $state, $stoichiometry, $location);
	} else {
	    #### before remains input but after becomes output
	    $inter_object->push_expanded_attribute("inputs",$ent_id, $subunit, $state_before, $stoichiometry, $location_before);
	    $inter_object->push_expanded_attribute("outputs",$ent_id, $subunit, $state_after, $stoichiometry, $location_after);

	    #### check the reason why location before and after differ
	    unless ($location_before eq $location_after) {
		if ($inter_object->get_attribute("type") eq "HELLO") {
		} else {
		    #### report unexpected cases
		    &ErrorMessage ("Unknown location modification for interaction output", 
				   "\t", $inter_object->get_attribute("type"),
				   "\t", $location_before, 
				   "\t", $location_after, 
				   "\n", $line, "\n");
		}
	    }

	    #### check the reason why state before and after differ
	    unless ($state_before eq $state_after) {
		if ($inter_object->get_attribute("type") eq "HELLO") {
		} else {
		    #### report unexpected cases
		    &ErrorMessage ("Unknown state modification for interaction output", 
				   "\t", $inter_object->get_attribute("type"),
				   "\t", $state_before, 
				   "\t", $state_after, 
				   "\n", $line, "\n");
		}
	    }
	}
    }
    close OUTPUT;


#      ### locations (temporary)
#      undef %col;
#      $col{molec} = 0;
#      $col{inter} = 1;
#      $col{location_before} = 2;
#      $col{location_after} = 3;
    
#      open LOCATIONS, "$LocationFile";
#      $header = <LOCATIONS>;
#      $line_nb = 0;
#      while (<LOCATIONS>) {
#  	$line_nb++;
#  	my @fields = &MySplit;  
#  	$inter = $fields[$col{inter}];
#  	unless (defined($association{$inter})) {
#  	    &ErrorMessage( "Error in $LocationFile line $line_nb: $inter has not been defined in $InteractionFile\n");
#  	    &ErrorMessage( "\t$_\n");
#  	    next;
#  	}
	
#  	$molec = &PrologString($fields[$col{molec}]);
#  	$lc_molec = lc($molec);
#  	if (defined ($AC{$lc_molec})) {
#  	    $molecAC = $AC{$lc_molec};
#  	} else {
#  	    &ErrorMessage( "Error in $LocationFile line $line_nb: $molec has not been defined in $in_file{entities}\n");
#  	    &ErrorMessage( "\t$_\n");
#  	    next;
#  	}
#  	$OK = 0;
#  	if (defined(@{$input_records{$inter}{$molecAC}})) {
#  	    foreach $line_nb (@{$input_records{$inter}{$molecAC}}) {
#  		${$input_attributes}[$line_nb]->{location_before} = $fields[$col{location_before}];
#      	        ${$input_attributes}[$line_nb]->{location_after} = $fields[$col{location_after}];
#  	        $OK = 1;
#              }
#          }
#          if (defined(@{$output_records{$inter}{$molecAC}})) {
#              foreach $line_nb (@{$output_records{$inter}{$molecAC}}) {
#  	        ${$output_attributes}[$line_nb]->{location_before} = $fields[$col{location_before}];
#                  ${$output_attributes}[$line_nb]->{location_after} = $fields[$col{location_after}];
#                  $OK = 1;
#              }
#          }
#          unless ($OK) {
#              &ErrorMessage( "Error in $LocationFile line $line_nb: ");
#              &ErrorMessage( "$molec has not been defined as input or output for interaction $inter\n");
#          }
#      }
#      close LOCATIONS;
}


### read pathways
sub ReadPathways {
    #### create a pathway to store all interactions and entities
    $complete_pathway = $pathways->new_object();
    $complete_pathway->push_attribute("names","Signal transduction pathways");

    ### read pathway description
    open PATH_DESC, $in_file{pathway_description} || die ";Error: cannot read pathway description file $in_file{pathway_description}\n";
    $header = <PATH_DESC>; ### skip header line
    while (<PATH_DESC>) {
	my @fields = &MySplit;
	my $name = $fields[0];
	my $description = $fields[1];
	
	my $pathway = $pathways->new_object(names=>$name);
	if ($description) {
	    $pathway->set_attribute("description", $description);
	} else {
	    $pathway->set_attribute("description", $name);
	}
	warn (join ("\t", "new pathway",
		    $pathway->get_attribute("id"),
		    $pathway->get_attribute("names")
		   ), "\n")
	    if ($warn_level >= 2);
	
    }
    close PATH_DESC;
    $pathways->index_ids();
    $pathways->index_names();

    ### read pathway entities
    open PATH_ENT, $in_file{pathway_entity} || die ";Error: cannot read pathway entity file $in_file{pathway_entity}\n";
    $header = <PATH_ENT>; ### skip header line
    my $line_count = 1;
    while (<PATH_ENT>) {
	$line_count++;
	my @fields = &MySplit;
	my $path_id = $fields[0];
	my $entity = $fields[1];
	warn ";\t$path_id\t$entity\n" if ($warn_level >= 2);

	#### identify the pathway
	unless ($pathway = $pathways->get_object($path_id)) {
	    &ErrorMessage( ";ERROR: file $in_file{pathway_entity} line $line_count: pathway $path_id has not been found in $in_file{pathway_descriptions}\n");
	    next;
	}

	#### identify the entity
	if ($entity_object = $entities->get_object($entity)) {
	    $ent_id = $entity_object->get_attribute("id");
	    $pathway->push_attribute("entities",$ent_id);
	    $complete_pathway->push_attribute("entities",$ent_id);
	} else {
	    &ErrorMessage( "Error in $in_file{pathwy_entity} line $line_nb: $entity has not been defined in $in_file{entities}\n");
	    &ErrorMessage( "\t$_\n");
	    next;
	}

    }
    close PATH_ENT;

    ### read pathway interactions
    open PATH_INT, $in_file{pathway_interaction} || die ";Error: cannot read pathway interaction file $in_file{pathway_interaction}\n";
    $header = <PATH_INT>; ### skip header line
    my $line_count = 1;
    while (<PATH_INT>) {
	$line_count++;
	my @fields = &MySplit;
	my $path_id = $fields[0];
	my $interaction = $fields[1];
	warn ";\t$path_id\t$interaction\n" if ($warn_level >= 2);

	#### identify the pathway
	unless ($pathway = $pathways->get_object($path_id)) {
	    &ErrorMessage( ";ERROR: file $in_file{pathway_interaction} line $line_count: pathway $path_id has not been found in $in_file{pathway_descriptions}\n");
	    next;
	}

	#### identify the interaction
	if ($interaction_object = $interactions->get_object($interaction)) {
	    $int_id = $interaction_object->get_attribute("id");
	    $pathway->push_attribute("interactions",$int_id);
	    $complete_pathway->push_attribute("interactions",$int_id);
	} else {
	    &ErrorMessage( "Error in $in_file{pathway_interaction} line $line_nb: $interaction has not been defined in $in_file{interactions}\n");
	    &ErrorMessage( "\t$_\n");
	    next;
	}

    }
    close PATH_INT;

    ### read pathway subpathways
    open PATH_INT, $in_file{pathway_subpathway} || die ";Error: cannot read pathway subpathway file $in_file{pathway_subpathway}\n";
    $header = <PATH_INT>; ### skip header line
    my $line_count = 1;
    while (<PATH_INT>) {
	$line_count++;
	my @fields = &MySplit;
	my $path_id = $fields[0];
	my $subpathway = $fields[1];
	warn ";\t$path_id\t$subpathway\n" if ($warn_level >= 2);

	#### identify the pathway
	unless ($pathway = $pathways->get_object($path_id)) {
	    &ErrorMessage( ";ERROR: file $in_file{pathway_subpathway} line $line_count: pathway $path_id has not been found in $in_file{pathway_descriptions}\n");
	    next;
	}

	#### identify the subpathway
	if ($subpathway_object = $pathways->get_object($subpathway)) {
	    $sub_id = $subpathway_object->get_attribute("id");
	    $pathway->push_attribute("subpathways",$sub_id);
	} else {
	    &ErrorMessage( "Error in $in_file{pathway_subpathway} line $line_nb: $subpathway has not been defined in $in_file{subpathways}\n");
	    &ErrorMessage( "\t$_\n");
	    next;
	}

    }
    close PATH_INT;

}

#  sub GetNextPwelAC {
#      $path_element_count++;
#      local($PwelAC) = sprintf "pwel_%5d", $path_element_count;
#      $PwelAC =~ s/ /0/g;
#      return $PwelAC;
#  }

#  sub CreatePathwayEntity {
#      local($PwelAC) = &GetNextPwelAC;
#      local($PathwayAC) = $_[0];
#      local($EntityAC) = $_[1];
#      local($Type) = $_[2];
#      local($Label) = $_[3];
#      local($Step) = 0;  
#      local($Xpos) = "_";
#      local($Ypos) = "_";

#      $pwel{$PwelAC}->{PwelAC} = $PwelAC;
#      $pwel{$PwelAC}->{PathwayAC} = $PathwayAC;
#      $pwel{$PwelAC}->{Step} = $Step;
#      $pwel{$PwelAC}->{Type} = $Type;
#      $pwel{$PwelAC}->{EntityAC} = $EntityAC;
#      $pwel{$PwelAC}->{Label} = $Label;
#      $pwel{$PwelAC}->{Xpos} = $Xpos;
#      $pwel{$PwelAC}->{Ypos} = $Ypos;
    
#      print PATHEL "pathway_element($PwelAC,$PathwayAC,$Step,entity,$Type,$EntityAC,$Label,$Xpos,$Ypos).\n";
#      ### remind the pertainance of this entity to this pathway
#      $member{$EntityAC}{$PathwayAC} = $PwelAC;

#  }

#  sub CreatePathwayAssociation {
#      print PATHEL "pathway_element($PwelAC,$PathwayAC,$Step,association,$Type,$AssocAC,$FromAC,$ToAC,$Label,$Xpos,$Ypos).\n";
#  }

### create a new entity
sub CreateEntity {
    local($lac) = $_[0];
    local($ltype) = $_[1];
    local($lname) = $_[2];
    local($ldescription) = $_[3];

    $entity{$ac}->{ac} = $lac;
    $entity{$ac}->{type} = $ltype;
    $entity{$ac}->{name} = $lname;
    $entity{$ac}->{description} = $ldescription;
    $entity_count{$ltype}++;
}

sub CreateAssociation {
    local($lac) = $_[0];
    local($ltype) = $_[1];
    local($lname) = $_[2];
    local($ldescription) = $_[3];
    local($linput) = $_[4];
    local($loutput) = $_[5];

    $association{$lac}->{ac} = $lac;
    $association{$lac}->{type} = $ltype;
    $association{$lac}->{name} = $lname;
    $association{$lac}->{description} = $ldescription;
    $association{$lac}->{input} = $linput;
    $association{$lac}->{output} = $loutput;
    $assoc_count{$ltype}++;
}

### print entity file
sub PrintEntities {
    open ENTITY, "> $EntityFile" || die "Error : cannot write file $EntityFile";;
    print ENTITY ":- multifile entity/4.\n";
    print ENTITY "% entity(Type,AC,Name,Description).\n";

    foreach $ac (sort keys %entity) {
	$type = $entity{$ac}->{type};
	$name = $entity{$ac}->{name};
	#    $description = "$entity{$ac}->{description}";
	$description = "[";
	$description .= $entity{$ac}->{name};
	foreach $altName (@{$altName{$name}}) {
	    $description .= ",$altName";
	}
	$description .= "]";
	
	print ENTITY "entity($type,$ac,$name,$description).\n";
	print NODEARC "node($ac, $type, $name).\n";
    }
    close ENTITY;
}

### print a summary of the parsing
sub PrintReport {
    @entity_list = sort keys %entity;
    $entity_nb = $#entity_list +1;
    print REPORT "$entity_nb Entities\n";
    print REPORT "\tTypes of Entities\n";
    foreach $type (sort keys %entityType) {
	print REPORT "\t\t",$entity_count{$type},"\t$type\n";
    }

    @assoc_list = sort keys %association;
    $assoc_nb = $#assoc_list +1;
    print REPORT "$assoc_nb Interactions\n";
    print REPORT "\tTypes of Interactions\n";
    foreach $type (sort keys %interactionType) {
	print REPORT "\t\t",$assoc_count{$type},"\t$type\n";
    }
}

sub PrintAssociations {
    open ASSOC, "> $AssocFile" || die "Error : cannot write file $AssocFile";
    print ASSOC ":- multifile association/5.\n";
    print ASSOC "% association(Type,AC,From,To,Label).\n";
    
    for $ac (sort keys %association) {
	### entity-association format
	$type = lc($association{$ac}->{type});
	$label = $association{$ac}->{type};
	
	$inputs = "";
	foreach $input (@{$association{$ac}->{inputs}}) {
	    $inputs .= $input.",";
	}
	$inputs =~ s/,$//;
	
	$outputs = "";
	foreach $output (@{$association{$ac}->{outputs}}) {
	    $outputs .= $output.",";
	}
	$outputs =~ s/,$//;
	
	print ASSOC "association($type,$ac,\[$inputs\],\[$outputs\],$label).\n";

	### node-arc format
	### export only associations that have at least one input and one output
	$to_export = 1;
	if ($#{$association{$ac}->{inputs}} < 0) {
	    &ErrorMessage( "Error: association $ac has no input\n");
	    $to_export = 0;
	}
	if ($#{$association{$ac}->{outputs}} < 0) {
	    &ErrorMessage( "Error: association $ac has no output\n");
	    $to_export = 0;
	}
	if ($to_export) {
	    foreach $input (@{$association{$ac}->{inputs}}) {
		print NODEARC "arc($input,$ac,+,'${label}_input').\n";
	    }
	    foreach $output (@{$association{$ac}->{outputs}}) {
		print NODEARC "arc($ac,$output,+,'${label}_output').\n";
	    }
	    print NODEARC "node($ac, $type, '${label}_$ac').\n";
	}

    }
    close ASSOC;
}

sub PrintInputsOutputs {
### just temporary : reprint input and output attributes
### and merge the info from the location table
    undef @col;
    $col[0] = inter;
    $col[1] = molec;
    $col[2] = subunit;
    $col[3] = state_before;
    $col[4] = state_after;
    $col[5] = stoeichiometry;
    $col[6] = location_before;
    $col[7] = location_after;

    open UPDATED_INPUTS, "> $Updatedin_file{interaction_input}" || die "Error: cannot write file $Updatedin_file{interaction_input}\n";
    print UPDATED_INPUTS $col[0];
    foreach $c (1..$#col) {
	print UPDATED_INPUTS "\t", $col[$c];
    }
    print UPDATED_INPUTS "\n";
    foreach $line_nb (1..$#{$input_attributes}) {
	print UPDATED_INPUTS ${$input_attributes}[$line_nb]->{$col[0]};
    foreach $c (1..$#col) {
	print UPDATED_INPUTS "\t", ${$input_attributes}[$line_nb]->{$col[$c]};
}
print UPDATED_INPUTS "\n";
}
close UPDATED_INPUTS;

open UPDATED_OUTPUTS, "> $Updatedin_file{interaction_output}" || die "Error: cannot write file $Updatedin_file{interaction_output}\n";
print UPDATED_OUTPUTS $col[0];
foreach $c (1..$#col) {
    print UPDATED_OUTPUTS "\t", $col[$c];
}
print UPDATED_OUTPUTS "\n";
foreach $line_nb (1..$#{$output_attributes}) {
    print UPDATED_OUTPUTS ${$output_attributes}[$line_nb]->{$col[0]};
foreach $c (1..$#col) {
    print UPDATED_OUTPUTS "\t", ${$output_attributes}[$line_nb]->{$col[$c]};
}
print UPDATED_OUTPUTS "\n";
}
close UPDATED_OUTPUTS;
}


### read arguments from the command line
sub ReadArguments {
    my $a = "";
    for $a (0..$#ARGV) {

	### warn level
	if (($ARGV[$a] eq "-v" ) && 
	    ($ARGV[$a+1] =~ /^\d+$/)){
	    $main::warn_level = $ARGV[$a+1];
	    $a++;

	    ### clean
	} elsif ($ARGV[$a] eq "-clean") {
	    $main::clean = 1;
	    

	    ### output file
 	} elsif ($ARGV[$a] eq "-obj") {
	    $a++;
	    $main::export{obj} = 1;

	    #### help
	} elsif (($ARGV[$a] eq "-h") ||
		 ($ARGV[$a] eq "-help")) {
	    &PrintHelp;
	    exit(0);
	    
	}
    }
}


### export a pathway diagram
sub PathwayToDiagram {
    ### usage : &PathwayToDiagram($pathway);
    ### collects all inputs and outputs of these interactions
    ### and prints the diagramin text format
    my ($pathway) = @_;
    my $arc_count = 0;
    my %linked_entities = ();

    my $xsize = 600;
    my $ysize = 600;
    srand (time);
    

    my @entity_ids = $pathway->get_attribute("entities");
    my @interaction_ids = $pathway->get_attribute("interactions");
    my @subpathway_ids = $pathway->get_attribute("subpathways");
    my $name = $pathway->get_name();

    my $filename = &my_trim($pathway->get_name());
    $filename =~ s/ /_/g;
    my $diagram_file = $dir{diagrams}."/$filename.tdd";

    warn ("; ",join ("\t", 
		     "entities", $#entity_ids+1,
		     "interactions", $#interaction_ids+1,
		     "subpathways", $#subpathay_ids+1,
		     "pathway", $pathway->get_attribute("names")
		     ), "\n") if ($warn_level >= 1);
    warn "; Creating diagram\t$name\n" if ($warn_level >= 1); 
    $diagram = $diagrams->new_object();
    $diagram->push_attribute("names", $name);
    $diagram->set_attribute("description", "Saccharomyces cerevisiae - $name");
    $diagram->set_attribute("type", "signal transduction");

    ### collect inputs and outputs for each interaction
    foreach my $interaction_id (@interaction_ids) {
	#### identify the interaction
	if ($interaction = $interactions->get_object($interaction_id)) {
	    foreach my $input ($interaction->get_input_IDs()) {
		$linked_entities{$input}++;
	    }
	    foreach my $output ($interaction->get_output_IDs()) {
		$linked_entities{$output}++;
	    }
	}
    }

    ### create nodes corresponding to subpathway
    warn "; Creating nodes for subpathways\n" if ($main::warn_level >= 2);
    foreach my $sub_id (@subpathway_ids) {
	$subpathway = $pathways->get_object($sub_id);
	unless ($diagram->get_node($sub_id)) {
	    #### create a node for the entity if necessary
	    my $node = $diagram->add_node(id=>$sub_id);
	    $node->set_attribute("label", $subpathway->get_name());
	    $node->set_attribute("type", "pathway");
	    $node->set_attribute("xpos", int(rand $xsize));
	    $node->set_attribute("ypos", int(rand $ysize));
	}

    }
    
    ### create nodes corresponding to interaction inputs/output
    warn "; Creating nodes for interaction inputs/outputs\n" if ($main::warn_level >= 2);
    foreach my $ent_id (keys %linked_entities) {
	$entity = $entities->get_object($ent_id);

	if ($entity) {
#	warn join("\t", "getting node", $ent_id, $diagram->get_node($ent_id), "\n");
	    unless ($diagram->get_node($ent_id)) {
		#### create a node for the entity if necessary
		my $node = $diagram->add_node(id=>$ent_id);
		$node->set_attribute("label", $entity->get_name());
		$node->set_attribute("type", $entity->get_attribute("type"));
		$node->set_attribute("xpos", int(rand $xsize));
		$node->set_attribute("ypos", int(rand $ysize));
	    }
	} else {
	    &ErrorMessage("Cannot identify entity\t$ent_id\n");
	}
    }
    
    #### create nodes for the explicitly specified entities
    warn "; Creating nodes for entities\n" if ($main::warn_level >= 2);
    foreach my $ent_id (@entity_ids) {
	$entity = $entities->get_object($ent_id);

#	warn join("\t", "getting node", $ent_id, $diagram->get_node($ent_id), "\n");
	unless ($diagram->get_node($ent_id)) {
	    #### create a node for the entity if necessary
	    my $node = $diagram->add_node(id=>$ent_id);
	    $node->set_attribute("label", $entity->get_name());
	    $node->set_attribute("type", $entity->get_attribute("type"));
	    $node->set_attribute("xpos", int(rand $xsize));
	    $node->set_attribute("ypos", int(rand $ysize));
	}
    }
    
    ### create a node for each interaction
    warn "; Creating nodes for interactions\n" if ($main::warn_level >= 2);
    foreach my $int_id (@interaction_ids) {
	my $interaction = $interactions->get_object($int_id);
	my $type = $interaction->get_attribute("type");

	#### create a node for the entity
	$abbrev{association} = "asm";
	$abbrev{phosphorylation} = "pho";
	$abbrev{dephosphorylation} = "dpho";
	$abbrev{inhibition} = "inh";
	$abbrev{activation} = "atc";
	$abbrev{transactivation} = "trac";
	$abbrev{transrepression} = "trep";

	my $node = $diagram->add_node(id=>$int_id);
	my $label = "";
	if ($abbrev{$type}) {
	    $label = $abbrev{$type};
	} else {
	    $label = $type;
	}
	$label .= $int_id;
	$node->set_attribute("label", $label);
#	$node->set_attribute("label", substr($label,0,3));
	$node->set_attribute("type", $type);
	$node->set_attribute("xpos", int(rand $xsize));
	$node->set_attribute("ypos", int(rand $ysize));
	
	### print arc between the interaction node and each input
	foreach my $input ($interaction->get_input_IDs()) {
	    my $from = $input;
	    my $to = $int_id; #### interaction id

	    #### create a new arc
	    my $arc = $diagram->add_arc(from=>$from, to=>$to);
	    $arc->set_attribute("type", $interaction->get_attribute("type")."_input");
	}

	### print arc between the interaction node and each output
	foreach my $output ($interaction->get_output_IDs()) {
	    my $from = $int_id; #### interaction id
	    my $to = $output; 

	    #### create a new arc
	    my $arc = $diagram->add_arc(from=>$from, to=>$to);
	    $arc->set_attribute("type", $interaction->get_attribute("type")."_output");
	}
    }

    #### export the diagram in text format
    warn "; Exporting diagram\t$diagram_file\n" if ($warn_level >= 1);
    $diagram->print("tdd", $diagram_file);
}


#  ### print entities and interactions in the tabular format for diagram description
#  sub PrintGraph {
#      ### usage : &PrintGraph($graph_file, @interactions);
#      ### collects all inputs and outputs of these interactions
#      ### and prints the graphs in .tdf format
#      my ($graph_file, @interactions) = @_;
#      my $title = "Roche T cell signalling data";
#      my $xsize = 600;
#      my $ysize = 600;
#      my $arc_count = 0;
#      my %linked_entities = ();
#      srand (time);
    
#      warn "; Printing graph\t$graph_file\n" if ($warn_level >= 1);
    
#      open GRAPH, ">$graph_file" 
#  	|| die "Error: cannot write $out_file{graph}\n";
    
#      print GRAPH "ftype\ttfd\n";
#      print GRAPH "title\t$title\n";
#      print GRAPH "xsize\t$xsize\n";
#      print GRAPH "ysize\t$ysize\n";
#      print GRAPH "creator\troche_parser\n";
#      print GRAPH "cdate\t", `date`;

#      ### collect inputs and outputs for each interaction
#      foreach my $interaction (@interactions) {
#  	foreach $input ($interaction->get_input_IDs()) {
#  	    $linked_entities{$input}++;
#  	}
#  	foreach $output ($interaction->get_output_IDs()) {
#  	    $linked_entities{$output}++;
#  	}
#      }

#      ### print nodes corresponding to these entities
#      foreach $id (keys %linked_entities) {
#  	$entity = $entities->get_object($id);
#  	$xpos = int(rand $xsize);
#  	$ypos = int(rand $ysize);
#  	$label = $entity->get_name;
#  	$type = $entity->get_attribute("type");

#  	&PrintNode;
#      }
    
#      ### print a node for each interaction
#      foreach my $interaction (@interactions) {
#  	$id = $interaction->get_attribute("id");
#  	$xpos = int(rand $xsize);
#  	$ypos = int(rand $ysize);
#  	$type = $interaction->get_attribute("type");
#  	$label = $type."_".$id;
	
#  	&PrintNode;

#  	### print arc between the interaction node and each input
#  	foreach $input ($interaction->get_input_IDs()) {
#  	    $arc_id = "arc_".$arc_count++;
#  	    $from = $input;
#  	    $to = $id; #### interaction id
#  	    $label = $type = "input";
#  	    &PrintArc;
#  	}

#  	### print arc between the interaction node and each output
#  	foreach $output ($interaction->get_output_IDs()) {
#  	    $arc_id = "arc_".$arc_count++;
#  	    $from = $id; #### interaction id
#  	    $to = $output; 
#  	    $label = $type = "output";
#  	    &PrintArc;
#  	}
#      }

#      close GRAPH;
#  }

#  sub PrintArc {
#      print GRAPH "arc";
#      print GRAPH "\t", $arc_id;
#      print GRAPH "\t", $from;
#      print GRAPH "\t", $to;
#      print GRAPH "\t", $label;
#      print GRAPH "\t", $type;
#      print GRAPH "\n";
#  }

#  sub PrintNode {
#      print GRAPH "node";
#      print GRAPH "\t", $id;
#      print GRAPH "\t", $xpos;
#      print GRAPH "\t", $ypos;
#      print GRAPH "\t", $label;
#      print GRAPH "\t", $type;
#      print GRAPH "\n";
#  }


sub LinkEntities {
    ### usage : @interactions = &LinkEntities(@entities);
    ### given a set of entities, return all interactions 
    ### that have at one of these as input and one as output
    ### this is a subgraph extraction with a maximum link length of 1 arc
    my @entities = @_;  ### a list of entity IDs
    my @link_interactions = ();
    my %seed = ();

    ### put entity ids in a hash
    foreach my $entity (@entities) {
	$seed{uc($entity)} = 1;
    }

    #### collect interactions
    foreach my $interaction ($interactions->get_objects()) {
	my $in = 0;
	my $out = 0;
	foreach $input ($interaction->get_input_IDs()) {
	    $in++ if ($seed{uc($input)});
	}
	foreach $output ($interaction->get_output_IDs()) {
	    $out++ if ($seed{uc($output)});
	}
	if (($in) and ($out)) {
	    push @link_interactions, $interaction;
	}
    }

    return @link_interactions;
}

### print the help message 
### when the program is called with the -h or -help option
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
        parse_signal_transduction.pl

DESCRIPTION

	Parse signal transduction data, exported in tab-delimited
	format from the Access database.

AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

OPTIONS	
	-h	detailed help
	-help	short list of options
	-v #	warn level
		Warn level 1 corresponds to a restricted verbosity
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-obj	export the data in .obj format
	-clean	remove all files from the output directory before
		parsing
EndHelp
  close HELP;
}

  


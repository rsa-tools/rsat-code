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

package main;
{

    ### initialization
    $start_time = `date +%Y-%m-%d.%H%M%S`;
    $clean = 1;
    
    $out_format = "obj";
    @classes = qw( PFBP::BiochemicalEntity PFBP::BiochemicalActivity PFBP::Pathway );
    
    #### class factory for entities
    $entities = PFBP::ClassFactory->new_class(object_type=>"PFBP::BiochemicalEntity",
					      prefix=>"ent_");
    $entities->set_out_fields(qw( id type primary_name names ));

    #### class factory for interactions
    $interactions = PFBP::ClassFactory->new_class(object_type=>"PFBP::BiochemicalActivity",
						  prefix=>"int_");
    $interactions->set_out_fields(qw( id type description inputs outputs ));

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
    $in_file{synonyms} = $dir{input}."/synonyms.txt";

    ### interactions (interactions)
    $in_file{interactions} = $dir{input}."/description_of_interactions.txt";
    $in_file{interaction_source} = $dir{input}."/interaction_source.txt";
    $in_file{interaction_target} = $dir{input}."/interaction_target.txt";

    ### pathways
    $in_file{pathway_description} = $dir{input}."/pathway_description.txt";
    $in_file{pathway_subpathway} = $dir{input}."/pathway_subpathway.txt";
    $in_file{pathway_interaction} = $dir{input}."/pathway_interaction.txt";
    $in_file{pathway_entity} = $dir{input}."/pathway_molecule.txt";

    #### check for the existence of all the input files
    foreach my $file (keys %infile) {
	&checkfile($file);
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
	warn "Cleaning output directory\n" if ($verbose >= 1);
#	die "HELLO";
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

    &DefaultVerbose() if ($verbose >= 1);

    open ERR, "> $out_file{errors}" || die "Error : cannot write file $out_file{errors}\n";

    &ReadEntities();
    &ReadInteractions();
    &ReadPathways();

    ### print each pathway in a separate file
    foreach $pathway ($pathways->get_objects()) {
	&ExportPathway($pathway);
    }

    #####################
    ### print objects ###
    #####################
    $entities->dump_tables();
    $interactions->dump_tables();
    $pathways->dump_tables();

    &ExportClasses($out_file{signal_transduction}, $out_format, @classes)  if ($export{obj});

    ###################
    ### print stats ###
    ###################
    &PrintStats($out_file{stats}, @classes);

    warn "; Done\t", `date` if ($verbose >= 1);

    close ERR;
    exit(0);
}


############# SUBROUTINE DEFINITION ################

### read entity information
sub ReadEntities {
    warn "; Reading entities\n" if ($verbose >= 1);
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
	&MySplit;
	warn "$_\n" if ($verbose >= 4);

	### read the fields ###
#	$entity_nb++;
#	$ac = sprintf "ent%5d", $entity_nb;
#	$ac =~ s/ /0/g;
	$name = $fields[$col{name}];
	$description = $fields[$col{descr}];
	
	### check entity type ###
	if ($fields[$col{type}] =~ /\S/) {
	    $type = lc($fields[$col{type}]);
	    $type =~ s/ /_/g;
	} else {
	    $type = "undef";
	    print ERR ";Error in $in_file{entities}, line $line_nb: entity $name has no type\n";
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
	    if ($verbose >= 2);
	
	### old routine to export in prolog
	&CreateEntity($ac,&PrologString($type),&PrologString($name),&PrologString($description));
	
	
	### remind accession number ###
	$lc_name = lc($entity{$ac}->{name});
	$AC{$lc_name} = $ac;
#print "$ac\t$lc_name\n";
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
    open NAMES, "$in_file{synonyms}"  || print "Warning : cannot read file $MolAltNameFile\n";
    $header = <NAMES>;
    $line_nb = 0;
    while (<NAMES>) {
	$line_nb++;
	&MySplit;
	
	$name = $fields[$col{name}];
	$altName = $fields[$col{altName}];
	#$lc_name = lc($name);
	#if (defined($AC{$lc_name})) {
	#  $AC{lc($altName)} =  $AC{$lc_name};
	#  push @{$altName{$name}}, $altName; 
	#} else {
	#  print  ERR ";Error in $in_file{synonyms} line $line_nb: $name had not been previously declared in $in_file{entities}\n";
	#}

	### add new name to the object
	if ($object = $entities->get_object($name)) {
	    $object->push_attribute("names", $altName);
	} else {      
	    print  ERR ";Error in $in_file{synonyms} line $line_nb: $name had not been previously declared in $in_file{entities}\n";
	}

    }
    close NAMES;
}

### read description of interact<ion ###
sub ReadInteractions {
    warn "; Reading interactions\n" if ($verbose >= 1);
    undef %col;
    $col{ac} = 0;
    $col{type} = 1;
    $col{descr} = 2;
    $assoc_nb = -1;

    open INTERACT, "$in_file{interactions}"  || die "Error : cannot read file $in_file{interactions}\n";
    $header = <INTERACT>;
    $line_nb = 0;
    while (<INTERACT>) {
	$line_nb++;
	&MySplit;  
	
	$assoc_nb++;
	$ac = $fields[$col{ac}];
	$type = lc($fields[$col{type}]);
	$type =~ s/^\'//;
	$type =~ s/\'$//;
	$type =~ s/ /_/g;
	unless ($type =~ /\S/) {
	    $type = "undef";
	    print ERR ";Error in $in_file{interactions} line $line_nb: interaction type is not specified\n";
	    print ERR "\t$_\n";
	}
	#$interactionType{$type} = 1;
	
	$name = "$type";
	$description = $fields[$col{descr}];
	$source = "";
	$target = "";
	
	my $interaction = $interactions->new_object(id=>$ac);
	$interaction->set_attribute("type",$type);
	$interaction->set_attribute("description",$description);
 	warn (join ("\t", "new interaction",
		    $interaction->get_attribute("id"),
		    $interaction->get_attribute("type")
		   ), "\n")
	    if ($verbose >= 2);
	

	&CreateAssociation($ac,$type,$name,$description,$source,$target);
    }
    close INTERACT;

    ### create indexes
    $interactions->index_ids();
    $interactions->index_names();


    ### read interaction sources
    warn "; Reading interaction sources\n" if ($verbose >= 1);

    ### initialisation
    undef %col;
    $col{inter} = 0;
    $col{molec} = 1;
    $col{subunit} = 2;
    $col{state_before} = 3;
    $col{state_after} = 4;
    $col{stoeichiometry} = 5;

    open SOURCE, "$in_file{interaction_source}"   || die "Error : cannot read file $in_file{interaction_source}\n";
    $header = <SOURCE>;
    $line_nb = 0;
    while (<SOURCE>) {
	$line_nb++;
	&MySplit;  

	### identify the interaction
	$inter = $fields[$col{inter}];
	unless ($activ_object = $interactions->get_object($inter)) {
	    #unless (defined($association{$inter})) {
	    print ERR "Error in $in_file{interaction_source} line $line_nb: $inter has not been defined in $InteractionFile\n";
	    print ERR "\t$_\n";
	    next;
	}

	### identify the entity
	$molec = $fields[$col{molec}];
	if ($entity_object = $entities->get_object($molec)) {
	    $ent_id = $entity_object->get_attribute("id");
	} else {
	    print ERR "Error in $in_file{interaction_source} line $line_nb: $molecAC has not been defined in $in_file{entities}\n";
	    print ERR "\t$_\n";
	    next;
	}
	
	### additional attributes
	$subunit = $fields[$col{subunit}];
	$state_before = $fields[$col{state_before}];
	$state_after = $fields[$col{state_after}];
	
	$activ_object->push_attribute("inputs",$ent_id);
	
	#$lc_molec = lc($molec);
	#if (defined ($AC{$lc_molec})) {
	#  $molecAC = $AC{$lc_molec};
	#} else {
	#  print ERR "Error in $in_file{interaction_source} line $line_nb: $molecAC has not been defined in $in_file{entities}\n";
	#  print ERR "\t$_\n";
#   #   next;
	#}
	#push @{$association{$inter}->{sources}}, $molecAC;

	#foreach $attribute (keys %col) {
	#  ${$source_attributes}[$line_nb]->{$attribute} = $fields[$col{$attribute}];
	#  push @{$source_records{$inter}{$molecAC}}, $line_nb;
	#}
#print "$inter\t", ${$source_attributes}[$line_nb]->{molec}, "\t", ${$source_attributes}[$line_nb]->{stoeichiometry}, "\n";
    }
    close SOURCE;

    ### interaction targets
    warn "; Reading interaction targets\n" if ($verbose >= 1);
    open TARGET, "$in_file{interaction_target}"   || die "Error : cannot read file $in_file{interaction_target}\n";
    $header = <TARGET>;
    $line_nb = 0;
    while (<TARGET>) {
	$line_nb++;
	&MySplit;  
	
	
	### identify the interaction
	$inter = $fields[$col{inter}];
	unless ($activ_object = $interactions->get_object($inter)) {
	    #unless (defined($association{$inter})) {
	    print ERR "Error in $in_file{interaction_source} line $line_nb: $inter has not been defined in $InteractionFile\n";
	    print ERR "\t$_\n";
	    next;
	}
	
	### identify the entity
	$molec = $fields[$col{molec}];
	if ($entity_object = $entities->get_object($molec)) {
	    $ent_id = $entity_object->get_attribute("id");
	} else {
	    print ERR "Error in $in_file{interaction_source} line $line_nb: $molecAC has not been defined in $in_file{entities}\n";
	    print ERR "\t$_\n";
	    next;
	}
	
	### additional attributes
	$subunit = $fields[$col{subunit}];
	$state_before = $fields[$col{state_before}];
	$state_after = $fields[$col{state_after}];
	
	$activ_object->push_attribute("outputs",$ent_id);


	#$inter = $fields[$col{inter}];
	#unless (defined($association{$inter})) {
	#  print ERR "Error in $in_file{interaction_target} line $line_nb: $inter has not been defined in $InteractionFile\n";
	#  print ERR "\t$_\n";
	#  next;
	#}
	#$molec = &PrologString($fields[$col{molec}]);
	#$lc_molec = lc($molec);
	#if (defined ($AC{$lc_molec})) {
	#  $molecAC = $AC{$lc_molec};
	#} else {
	#  print ERR "Error in $in_file{interaction_target} line $line_nb: $molecAC has not been defined in $in_file{entities}\n";
	#  print ERR "\t$_\n";
	#  next;
	#}

	#$subunit = $fields[$col{subunit}];
	#$state_before = $fields[$col{state_before}];
	#$state_after = $fields[$col{state_after}];
#
	#push @{$association{$inter}->{targets}}, $molecAC;
	#foreach $attribute (keys %col) {
	#  ${$target_attributes}[$line_nb]->{$attribute} = $fields[$col{$attribute}];
	#  push @{$target_records{$inter}{$molecAC}}, $line_nb;
	#}
    }
    close TARGET;

    ### locations (temporary)
    undef %col;
    $col{molec} = 0;
    $col{inter} = 1;
    $col{location_before} = 2;
    $col{location_after} = 3;
    
    open LOCATIONS, "$LocationFile";
    $header = <LOCATIONS>;
    $line_nb = 0;
    while (<LOCATIONS>) {
	$line_nb++;
	&MySplit;  
	$inter = $fields[$col{inter}];
	unless (defined($association{$inter})) {
	    print ERR "Error in $LocationFile line $line_nb: $inter has not been defined in $InteractionFile\n";
	    print ERR "\t$_\n";
	    next;
	}
	
	$molec = &PrologString($fields[$col{molec}]);
	$lc_molec = lc($molec);
	if (defined ($AC{$lc_molec})) {
	    $molecAC = $AC{$lc_molec};
	} else {
	    print ERR "Error in $LocationFile line $line_nb: $molec has not been defined in $in_file{entities}\n";
	    print ERR "\t$_\n";
	    next;
	}
	$OK = 0;
	if (defined(@{$source_records{$inter}{$molecAC}})) {
	    foreach $line_nb (@{$source_records{$inter}{$molecAC}}) {
		${$source_attributes}[$line_nb]->{location_before} = $fields[$col{location_before}];
	    ${$source_attributes}[$line_nb]->{location_after} = $fields[$col{location_after}];
	$OK = 1;
    }
}
if (defined(@{$target_records{$inter}{$molecAC}})) {
    foreach $line_nb (@{$target_records{$inter}{$molecAC}}) {
	${$target_attributes}[$line_nb]->{location_before} = $fields[$col{location_before}];
    ${$target_attributes}[$line_nb]->{location_after} = $fields[$col{location_after}];
$OK = 1;
}
}
unless ($OK) {
    print ERR "Error in $LocationFile line $line_nb: ";
    print ERR "$molec has not been defined as source or target for interaction $inter\n";
}
}
close LOCATIONS;
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
	&MySplit;
	my $name = $fields[0];
	my $description = $fields[1];
#  	my $id = $name;
#  	$id =~ s/^\'//g;
#  	$id =~ s/\'$//g;
#  	$id =~ s/\'/prime/g;
#  	$id =~ s/\s/_/g;
#  	$id =~ s/\-/_/g;
#  	$id = "path_".$id;
	
	my $pathway = $pathways->new_object(names=>$name);
#	$pathway->push_attribute("names", $name);
	if ($description) {
	    $pathway->set_attribute("description", $description);
	} else {
	    $pathway->set_attribute("description", $name);
	}
	warn (join ("\t", "new pathway",
		    $pathway->get_attribute("id"),
		    $pathway->get_attribute("names")
		   ), "\n")
	    if ($verbose >= 2);
	
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
	&MySplit;
	my $path_id = $fields[0];
	my $entity = $fields[1];
	warn ";\t$path_id\t$entity\n" if ($verbose >= 2);

	#### identify the pathway
	unless ($pathway = $pathways->get_object($path_id)) {
	    print ERR ";ERROR: file $in_file{pathway_entity} line $line_count: pathway $path_id has not been found in $in_file{pathway_descriptions}\n";
	    next;
	}

	#### identify the entity
	if ($entity_object = $entities->get_object($entity)) {
	    $ent_id = $entity_object->get_attribute("id");
	    $pathway->push_attribute("entities",$ent_id);
	    $complete_pathway->push_attribute("entities",$ent_id);
	} else {
	    print ERR "Error in $in_file{pathwy_entity} line $line_nb: $entity has not been defined in $in_file{entities}\n";
	    print ERR "\t$_\n";
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
	&MySplit;
	my $path_id = $fields[0];
	my $interaction = $fields[1];
	warn ";\t$path_id\t$interaction\n" if ($verbose >= 2);

	#### identify the pathway
	unless ($pathway = $pathways->get_object($path_id)) {
	    print ERR ";ERROR: file $in_file{pathway_interaction} line $line_count: pathway $path_id has not been found in $in_file{pathway_descriptions}\n";
	    next;
	}

	#### identify the interaction
	if ($interaction_object = $interactions->get_object($interaction)) {
	    $int_id = $interaction_object->get_attribute("id");
	    $pathway->push_attribute("interactions",$int_id);
	    $complete_pathway->push_attribute("interactions",$int_id);
	} else {
	    print ERR "Error in $in_file{pathway_interaction} line $line_nb: $interaction has not been defined in $in_file{interactions}\n";
	    print ERR "\t$_\n";
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
	&MySplit;
	my $path_id = $fields[0];
	my $subpathway = $fields[1];
	warn ";\t$path_id\t$subpathway\n" if ($verbose >= 2);

	#### identify the pathway
	unless ($pathway = $pathways->get_object($path_id)) {
	    print ERR ";ERROR: file $in_file{pathway_subpathway} line $line_count: pathway $path_id has not been found in $in_file{pathway_descriptions}\n";
	    next;
	}

	#### identify the subpathway
	if ($subpathway_object = $pathways->get_object($subpathway)) {
	    $sub_id = $subpathway_object->get_attribute("id");
	    $pathway->push_attribute("subpathways",$sub_id);
	} else {
	    print ERR "Error in $in_file{pathway_subpathway} line $line_nb: $subpathway has not been defined in $in_file{subpathways}\n";
	    print ERR "\t$_\n";
	    next;
	}

    }
    close PATH_INT;

}

#  sub oldReadPathways {
#      open PATHEL, "> $PathwayElementFile" || die  "Error: cannot write pathway element file $PathwayElementFile\n";
#      open PATH, "$PathwayFile" || die "Error: cannot read pathway file $PathwayFile\n";
#      while (<PATH>) {
#  	&MySplit;
#  	### pathway name
#  	$PathwayName = $fields[0];
#  	$PathwayName = "undef" if ($PathwayName eq "");
#  	### pathway accession number
#  	$PathwayAC = "$PathwayName";
#  	$PathwayAC =~ s/^\'//g;
#  	$PathwayAC =~ s/\'$//g;
#  	$PathwayAC =~ s/\'/prime/g;
#  	$PathwayAC =~ s/\s/_/g;
#  	$PathwayAC =~ s/\-/_/g;
#  	$PathwayAC = "path_".$PathwayAC;
#  	$pathway{$PathwayAC}->{ac} = $PathwayAC;
#  	$pathway{$PathwayAC}->{name} = $PathwayName;
#  	### entities
#  	$EntityName = $fields[1];
#  	$lc_name = lc($EntityName);
#  	if ($AC{$lc_name} ne "") {
#  	    $EntityAC = $AC{$lc_name};
#  	    $Type = $entity{$EntityAC}->{type};
#  	} else {
#  	    print  ERR ";Error in $PathwayFile: $name had not been previously declared in $in_file{entities}\n";
#  	}
#  	&CreatePathwayEntity($PathwayAC,$EntityAC,$Type,$EntityName);
#      }
#      close PATH;
#      ### convert associations into pathway elements
#      foreach $AssocAC (keys %association) {
#  	$AssocAC = $association{$AssocAC}->{ac};
#  	$AssocType = $association{$AssocAC}->{type};
#  	$Label = $association{$AssocAC}->{name};
#  	$source = $association{$AssocAC}->{source};
#  	$target = $association{$AssocAC}->{target};
#  	foreach $PathwayAC (keys %pathway) {
#  	    next if ($PathwayAC eq "path_undef");
#  	    if (($member{$source}{$PathwayAC}) && ($member{$target}{$PathwayAC})) {
#  		$PwelAC = &GetNextPwelAC;
#  		$FromAC = $member{$source}{$PathwayAC};
#  		$ToAC = $member{$target}{$PathwayAC};
#  		print PATHEL "pathway_element($PwelAC,$PathwayAC,$Step,association,$AssocType,$AssocAC,$FromAC,$ToAC,$Label,_,_).\n";
#  	    }
#  	}
#      }
#      close PATHEL;
#      ### split pathway_element by pathway
#      foreach $PathwayAC (keys %pathway) {
#  	next if ($PathwayAC eq "path_undef");
#  	$command = "grep $PathwayAC $PathwayElementFile > ${PathwayAC}.pl";
#  	system $command;
#      }
#  }

sub GetNextPwelAC {
    $path_element_count++;
    local($PwelAC) = sprintf "pwel_%5d", $path_element_count;
    $PwelAC =~ s/ /0/g;
    return $PwelAC;
}

sub CreatePathwayEntity {
    local($PwelAC) = &GetNextPwelAC;
    local($PathwayAC) = $_[0];
    local($EntityAC) = $_[1];
    local($Type) = $_[2];
    local($Label) = $_[3];
    local($Step) = 0;  
    local($Xpos) = "_";
    local($Ypos) = "_";

    $pwel{$PwelAC}->{PwelAC} = $PwelAC;
    $pwel{$PwelAC}->{PathwayAC} = $PathwayAC;
    $pwel{$PwelAC}->{Step} = $Step;
    $pwel{$PwelAC}->{Type} = $Type;
    $pwel{$PwelAC}->{EntityAC} = $EntityAC;
    $pwel{$PwelAC}->{Label} = $Label;
    $pwel{$PwelAC}->{Xpos} = $Xpos;
    $pwel{$PwelAC}->{Ypos} = $Ypos;
    
    print PATHEL "pathway_element($PwelAC,$PathwayAC,$Step,entity,$Type,$EntityAC,$Label,$Xpos,$Ypos).\n";
    ### remind the pertainance of this entity to this pathway
    $member{$EntityAC}{$PathwayAC} = $PwelAC;

}

sub CreatePathwayAssociation {
    print PATHEL "pathway_element($PwelAC,$PathwayAC,$Step,association,$Type,$AssocAC,$FromAC,$ToAC,$Label,$Xpos,$Ypos).\n";
}

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
    local($lsource) = $_[4];
    local($ltarget) = $_[5];

    $association{$lac}->{ac} = $lac;
    $association{$lac}->{type} = $ltype;
    $association{$lac}->{name} = $lname;
    $association{$lac}->{description} = $ldescription;
    $association{$lac}->{source} = $lsource;
    $association{$lac}->{target} = $ltarget;
    $assoc_count{$ltype}++;
#print "ac\t$lac\t", $association{$lac}->{ac}, "\n";
#print "type\t$ltype\t", $association{$ltype}->{ac}, "\n";
#print "name\t$lname\t", $association{$lname}->{ac}, "\n";
#print "description\t$ldescription\t", $association{$ldescription}->{ac}, "\n";
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
	
	$sources = "";
	foreach $source (@{$association{$ac}->{sources}}) {
	    $sources .= $source.",";
	}
	$sources =~ s/,$//;
	
	$targets = "";
	foreach $target (@{$association{$ac}->{targets}}) {
	    $targets .= $target.",";
	}
	$targets =~ s/,$//;
	
	print ASSOC "association($type,$ac,\[$sources\],\[$targets\],$label).\n";

	### node-arc format
	### export only associations that have at least one source and one target
	$to_export = 1;
	if ($#{$association{$ac}->{sources}} < 0) {
	    print ERR "Error: association $ac has no source\n";
	    $to_export = 0;
	}
	if ($#{$association{$ac}->{targets}} < 0) {
	    print ERR "Error: association $ac has no target\n";
	    $to_export = 0;
	}
	if ($to_export) {
	    foreach $source (@{$association{$ac}->{sources}}) {
		print NODEARC "arc($source,$ac,+,'${label}_source').\n";
	    }
	    foreach $target (@{$association{$ac}->{targets}}) {
		print NODEARC "arc($ac,$target,+,'${label}_target').\n";
	    }
	    print NODEARC "node($ac, $type, '${label}_$ac').\n";
	}

    }
    close ASSOC;
}

sub PrintSourcesTargets {
### just temporary : reprint source and target attributes
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

    open UPDATED_SOURCES, "> $Updatedin_file{interaction_source}" || die "Error: cannot write file $Updatedin_file{interaction_source}\n";
    print UPDATED_SOURCES $col[0];
    foreach $c (1..$#col) {
	print UPDATED_SOURCES "\t", $col[$c];
    }
    print UPDATED_SOURCES "\n";
    foreach $line_nb (1..$#{$source_attributes}) {
	print UPDATED_SOURCES ${$source_attributes}[$line_nb]->{$col[0]};
    foreach $c (1..$#col) {
	print UPDATED_SOURCES "\t", ${$source_attributes}[$line_nb]->{$col[$c]};
}
print UPDATED_SOURCES "\n";
}
close UPDATED_SOURCES;

open UPDATED_TARGETS, "> $Updatedin_file{interaction_target}" || die "Error: cannot write file $Updatedin_file{interaction_target}\n";
print UPDATED_TARGETS $col[0];
foreach $c (1..$#col) {
    print UPDATED_TARGETS "\t", $col[$c];
}
print UPDATED_TARGETS "\n";
foreach $line_nb (1..$#{$target_attributes}) {
    print UPDATED_TARGETS ${$target_attributes}[$line_nb]->{$col[0]};
foreach $c (1..$#col) {
    print UPDATED_TARGETS "\t", ${$target_attributes}[$line_nb]->{$col[$c]};
}
print UPDATED_TARGETS "\n";
}
close UPDATED_TARGETS;
}


### read arguments from the command line
sub ReadArguments {
    my $a = "";
    for $a (0..$#ARGV) {

	### warn level
	if (($ARGV[$a] eq "-v" ) && 
	    ($ARGV[$a+1] =~ /^\d+$/)){
	    $main::verbose = $ARGV[$a+1];
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
sub ExportPathway {
    ### usage : &ExportPathway($pathway);
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
    my $diagram_file = $dir{diagrams}."/".$pathway->get_attribute("id").".tdd";

    warn ("; ",join ("\t", 
		     "entities", $#entity_ids+1,
		     "interactions", $#interaction_ids+1,
		     "subpathways", $#subpathay_ids+1,
		     "pathway", $pathway->get_attribute("names")
		     ), "\n") if ($verbose >= 1);
    warn "; Creating diagram\t$name\n" if ($verbose >= 1); 
    $diagram = $diagrams->new_object();
    $diagram->push_attribute("names", $name);
    $diagram->set_attribute("description", "Saccharomyces cerevisiae - $name");
    $diagram->set_attribute("type", "signal transduction");

    ### collect inputs and outputs for each interaction
    foreach my $interaction_id (@interaction_ids) {
	#### identify the interaction
	if ($interaction = $interactions->get_object($interaction_id)) {
	    foreach my $input ($interaction->get_attribute("inputs")) {
		$linked_entities{$input}++;
	    }
	    foreach my $output ($interaction->get_attribute("outputs")) {
		$linked_entities{$output}++;
	    }
	}
    }

    ### print nodes corresponding to these entities
    foreach my $ent_id (keys %linked_entities) {
	$entity = $entities->get_object($ent_id);

	unless ($diagram->get_node("ent_id")) {
	    #### create a node for the entity if necessary
	    my $node = $diagram->add_node(id=>$ent_id);
	    $node->set_attribute("label", $entity->get_name());
	    $node->set_attribute("type", $entity->get_attribute("type"));
	    $node->set_attribute("xpos", int(rand $xsize));
	    $node->set_attribute("ypos", int(rand $ysize));
	}
    }
    
    ### print a node for each interaction
    foreach my $int_id (@interaction_ids) {
	my $interaction = $interactions->get_object($int_id);
	my $type = $interaction->get_attribute("type");

	#### create a node for the entity
	$node = $diagram->add_node(id=>$int_id);
	$node->set_attribute("label", $type."_".$int_id);
	$node->set_attribute("type", $type);
	$node->set_attribute("xpos", int(rand $xsize));
	$node->set_attribute("ypos", int(rand $ysize));
	
	### print arc between the interaction node and each input
	foreach my $input ($interaction->get_attribute("inputs")) {
	    my $from = $input;
	    my $to = $int_id; #### interaction id

	    #### create a new arc
	    my $arc = $diagram->add_arc(from=>$from, to=>$to);
	    $arc->set_attribute("type", $interaction->get_attribute("type")."_input");
	}

	### print arc between the interaction node and each output
	foreach my $output ($interaction->get_attribute("outputs")) {
	    my $from = $int_id; #### interaction id
	    my $to = $output; 

	    #### create a new arc
	    my $arc = $diagram->add_arc(from=>$from, to=>$to);
	    $arc->set_attribute("type", $interaction->get_attribute("type")."_output");
	}
    }

    #### export the diagram in text format
    warn "; Exporting diagram\t$diagram_file\n" if ($verbose >= 1);
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
    
#      warn "; Printing graph\t$graph_file\n" if ($verbose >= 1);
    
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
#  	foreach $input ($interaction->get_attribute("inputs")) {
#  	    $linked_entities{$input}++;
#  	}
#  	foreach $output ($interaction->get_attribute("outputs")) {
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
#  	foreach $input ($interaction->get_attribute("inputs")) {
#  	    $arc_id = "arc_".$arc_count++;
#  	    $from = $input;
#  	    $to = $id; #### interaction id
#  	    $label = $type = "input";
#  	    &PrintArc;
#  	}

#  	### print arc between the interaction node and each output
#  	foreach $output ($interaction->get_attribute("outputs")) {
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
    ### usage : @activites = &LinkEntities(@entities);
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
	foreach $input ($interaction->get_attribute("inputs")) {
	    $in++ if ($seed{uc($input)});
	}
	foreach $output ($interaction->get_attribute("outputs")) {
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
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-obj	export the data in .obj format
	-clean	remove all files from the output directory before
		parsing
EndHelp
  close HELP;
}

  


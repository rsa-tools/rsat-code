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

### initialization
$start_time = `date +%Y-%m-%d.%H%M%S`;

$out_format = "obj";
@classes = qw( PFBP::BiochemicalEntity PFBP::BiochemicalActivity );

### default output fields for each class
@{$out_fields{'PFBP::BiochemicalEntity'}} = qw( id type names );
@{$out_fields{'PFBP::BiochemicalActivity'}} = qw( id type description inputs outputs );
@{$out_fields{'PFBP::Pathway'}} = qw( id names description entities );

### old initialisation
$path_element_count = 0;
$Step = 0;

### directory containing the input files ###
$dir{input} = "/win/amaze/amaze_team/Sandra/export";
#$dir{input} = "/home/jvanheld/Databases/Roche/roche_991118/tab_export/";
unless (-e $dir{input}) {
    die "Error : cannot open input dir $dir{input}\n";
}

### molecules (entities)
$in_file{molecules} = $dir{input}."/description_of_molecules.txt";
$in_file{synonyms} = $dir{input}."/synonyms.txt";

### interactions (activities)
$in_file{interactions} = $dir{input}."/description_of_interactions.txt";
$in_file{interaction_source} = $dir{input}."/interaction_source.txt";
$in_file{interaction_target} = $dir{input}."/interaction_target.txt";

### pathways
$in_file{pathway_description} = $dir{input}."/pathway_description.txt";
$in_file{pathway_molecule} = $dir{input}."/pathway_molecule.txt";

#$PathwayFile = $dir{input}."/pathway_description.txt";
#$PathwayInteractionFile = $dir{input}."/pathway_interaction.txt";
#$PathwayMoleculeFile = $dir{input}."/pathway_molecule.txt";
#$PathwaySubPathwayFile = $dir{input}."/pathway_subpathway.txt";

#$LocationFile = $dir{input}."/location.txt";
#$PathwayFile = $dir{input}."/pathway_description.txt";

### output dir
$dir{output} = "${parsed_data}/signal_transduction/$delivery_date";
unless (-d $dir{output}) {
    warn "Creating output dir $dir{output}\n"  if ($warn_level >= 1);
    `mkdir -p $dir{output}`;
    die "Error: cannot create directory $dir\n" 
	unless (-d $dir{output});
}
chdir $dir{output};

### object files
$out_file{entities} = $dir{output}."/sigtrans_entities";
$out_file{activities} = $dir{output}."/sigtrans_activities";
$out_file{pathways} = $dir{output}."/sigtrans_pathways";
$out_file{graph} = $dir{output}."/sigtrans_graph.tdf";
$out_file{stats} = $dir{output}."/sigtrans_parsing_stats";

### reports
$out_file{parsing_errors} = $dir{output}."/sigtrans_parsing_errors.txt";
#$out_file{parsing_report} = $dir{output}."/sigtrans_parsing_report.txt";

### prolog files
#$EntityFile = $dir{output}."/db_sigtrans_entity.pl";
#$AssocFile = $dir{output}."/db_sigtrans_assoc.pl";
#$PathwayElementFile = $dir{output}."/path_sigtrans_pathel.pl";
#$NodeArcFile = $dir{output}."/int_sigtrans.pl";
#$Updatedin_file{interaction_target} = $dir{output}."/updated_targets.txt";
#$Updatedin_file{interaction_source} = $dir{output}."/updated_sources.txt";


&ReadArguments;

### actualize parameters
if ($test) {
    warn ";TEST\n" if ($verbose);
    ### fast partial parsing for debugging
    $in_file{molecules} = " head -20 $in_file{molecules} |";
}

&DefaultVerbose();
if ($verbose) {
    print "; input dir\t$dir{input}\n";
    print "; result dir\t$dir{output}\n";
    print "; input files\n";
    while (($key, $file) = each (%in_file)) {
	print ";\t$key\t$file\n";
    }
    print "; output files\n";
    while (($key, $file) = each (%out_file)) {
	print ";\t$key\t$file\n";
    }
}



#open REPORT, "> $out_file{parsing_report}"  || die "Error : cannot write file $out_file{parsing_report}";;
open ERR, "> $out_file{parsing_errors}" || die "Error : cannot write file $out_file{parsing_errors}\n";
#open NODEARC, "> $NodeArcFile" || die "Error : cannot write file $NodeArcFile\n";

&ReadMolecules;
&ReadInteractions;
&ReadPathways;

#&PrintEntities;
#&PrintAssociations;
#&PrintSourcesTargets;
#&PrintReport;

### print the complete graph of interactions
&PrintGraph($out_file{graph},PFBP::BiochemicalActivity->get_objects());

### print each pathway in a separate file
foreach $pathway (PFBP::Pathway->get_objects()) {
    my @entities = $pathway->get_attribute("entities");
    my @activities = &LinkEntities(@entities);
    my $out_file = "pathway_".$pathway->get_attribute("id").".tdf";
    warn ";",$#entities," entities\n";
    warn ";",$#activities," activities\n";
    &PrintGraph($out_file,@activities);
}

close ERR;
#close REPORT;
#close NODEARC;

#####################
### print objects ###
#####################

#&ExportClasses($out_file{entities}, $out_format,"PFBP::BiochemicalEntity");
#&ExportClasses($out_file{activities}, $out_format,"PFBP::BiochemicalActivity");
#&ExportClasses($out_file{pathways}, $out_format,"PFBP::Pathway");

###################
### print stats ###
###################
&PrintStats($out_file{stats});


warn "; done\t", `date` if ($verbose);

exit(0);

############# SUBROUTINE DEFINITION ################

### read molecule information
sub ReadMolecules {
    warn "; reading molecules\n" if ($verbose);
    undef %col;
    $col{swissprot_ac} = 0;
    $col{swissprot_id} = 1;
    $col{name} = 2;
    $col{type} = 3;
    $col{descr} = 4;
    $entity_nb = -1;
    $class = "PFBP::BiochemicalEntity";  

    ### read molecule description
    open MOL, "$in_file{molecules}" || die "Error : cannot read file $in_file{molecules}\n";
    $header = <MOL>;
    $line_nb = 0;
    while (<MOL>) {
	$line_nb++;
	&MySplit;
	
	### read the fields ###
	$entity_nb++;
	$ac = sprintf "ent%5d", $entity_nb;
	$ac =~ s/ /0/g;
	$name = $fields[$col{name}];
	$description = $fields[$col{descr}];
	
	### check molecule type ###
	if ($fields[$col{type}] =~ /\S/) {
	    $type = lc($fields[$col{type}]);
	    $type =~ s/ /_/g;
	} else {
	    $type = "undef";
	    print ERR ";Error in $in_file{molecules}, line $line_nb: molecule $entity{$ac}->{name} has no type\n";
	    print ERR "\t$_\n";
	}
	$moleculeType{$type} = 1;
	
	### instantiate a new BiochemicalEntity
	my $BiochemicalEntity = new $class(id=>$ac);
	warn "new entity\t$ac\t$BiochemicalEntity\n" if ($hyperverbose);
	$BiochemicalEntity->push_attribute("names",$name);
	$BiochemicalEntity->set_attribute("description",$description);
	$BiochemicalEntity->set_attribute("type",$type);
	
	### old routine to export in prolog
	&CreateEntity($ac,&PrologString($type),&PrologString($name),&PrologString($description));
	
	
	### remind accession number ###
	$lc_name = lc($entity{$ac}->{name});
	$AC{$lc_name} = $ac;
#print "$ac\t$lc_name\n";
    } 
    close MOL;
    
    ### create indexes
    $class->index_ids();
    $class->index_names();
    

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
	#  print  ERR ";Error in $in_file{synonyms} line $line_nb: $name had not been previously declared in $in_file{molecules}\n";
	#}

	### add new name to the object
	if ($object = $class->get_object($name)) {
	    $object->push_attribute("names", $altName);
	} else {      
	    print  ERR ";Error in $in_file{synonyms} line $line_nb: $name had not been previously declared in $in_file{molecules}\n";
	}

    }
    close NAMES;
}

### read description of interact<ion ###
sub ReadInteractions {
    warn "; reading interactions\n" if ($verbose);
    undef %col;
    $col{ac} = 0;
    $col{type} = 1;
    $col{descr} = 2;
    $assoc_nb = -1;

    $class = "PFBP::BiochemicalActivity";
    
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
	
	my $BiochemicalActivity = new $class(id=>$ac);
	$BiochemicalActivity->set_attribute("type",$type);
	$BiochemicalActivity->set_attribute("description",$description);

	&CreateAssociation($ac,$type,$name,$description,$source,$target);
    }
    close INTERACT;

    ### create indexes
    $class->index_ids();
    $class->index_names();


    ### read interaction sources
    warn "; reading interaction sources\n" if ($verbose);

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
	unless ($inter_object = $class->get_object($inter)) {
	    #unless (defined($association{$inter})) {
	    print ERR "Error in $in_file{interaction_source} line $line_nb: $inter has not been defined in $InteractionFile\n";
	    print ERR "\t$_\n";
	    next;
	}

	### identify the molecule
	$molec = $fields[$col{molec}];
	if ($molec_object = PFBP::BiochemicalEntity->get_object($molec)) {
	    $molec_id = $molec_object->get_attribute("id");
	} else {
	    print ERR "Error in $in_file{interaction_source} line $line_nb: $molecAC has not been defined in $in_file{molecules}\n";
	    print ERR "\t$_\n";
	    next;
	}
	
	### additional attributes
	$subunit = $fields[$col{subunit}];
	$state_before = $fields[$col{state_before}];
	$state_after = $fields[$col{state_after}];
	
	$inter_object->push_attribute("inputs",$molec_id);
	
	#$lc_molec = lc($molec);
	#if (defined ($AC{$lc_molec})) {
	#  $molecAC = $AC{$lc_molec};
	#} else {
	#  print ERR "Error in $in_file{interaction_source} line $line_nb: $molecAC has not been defined in $in_file{molecules}\n";
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
    warn "; reading interaction targets\n" if ($verbose);
    open TARGET, "$in_file{interaction_target}"   || die "Error : cannot read file $in_file{interaction_target}\n";
    $header = <TARGET>;
    $line_nb = 0;
    while (<TARGET>) {
	$line_nb++;
	&MySplit;  
	
	
	### identify the interaction
	$inter = $fields[$col{inter}];
	unless ($inter_object = $class->get_object($inter)) {
	    #unless (defined($association{$inter})) {
	    print ERR "Error in $in_file{interaction_source} line $line_nb: $inter has not been defined in $InteractionFile\n";
	    print ERR "\t$_\n";
	    next;
	}
	
	### identify the molecule
	$molec = $fields[$col{molec}];
	if ($molec_object = PFBP::BiochemicalEntity->get_object($molec)) {
	    $molec_id = $molec_object->get_attribute("id");
	} else {
	    print ERR "Error in $in_file{interaction_source} line $line_nb: $molecAC has not been defined in $in_file{molecules}\n";
	    print ERR "\t$_\n";
	    next;
	}
	
	### additional attributes
	$subunit = $fields[$col{subunit}];
	$state_before = $fields[$col{state_before}];
	$state_after = $fields[$col{state_after}];
	
	$inter_object->push_attribute("outputs",$molec_id);
	

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
	#  print ERR "Error in $in_file{interaction_target} line $line_nb: $molecAC has not been defined in $in_file{molecules}\n";
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
	    print ERR "Error in $LocationFile line $line_nb: $molec has not been defined in $in_file{molecules}\n";
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
    ### read pathway description
    open PATH_DESC, $in_file{pathway_description} || die ";Error: cannot read pathway description file $in_file{pathway_description}\n";
    $header = <PATH_DESC>; ### skip header line
    while (<PATH_DESC>) {
	&MySplit;
	my $name = $fields[0];
	my $description = $fields[1];
	my $id = $name;
	$id =~ s/^\'//g;
	$id =~ s/\'$//g;
	$id =~ s/\'/prime/g;
	$id =~ s/\s/_/g;
	$id =~ s/\-/_/g;
	$id = "path_".$id;
	unless ($pathway = PFBP::Pathway->new(id=>$id)) {
	    print ERR ";ERROR: could not create pathway $id\n";
	}
	$pathway->set_attribute("description", $description);
	$pathway->push_attribute("names", $name);
    }
    close PATH_DESC;
  PFBP::Pathway->index_ids();
  PFBP::Pathway->index_names();

    ### read pathway entities
    open PATH_ENT, $in_file{pathway_molecule} || die ";Error: cannot read pathway molecule file $in_file{pathway_molecule}\n";
    $header = <PATH_ENT>; ### skip header line
    my $line_count = 1;
    while (<PATH_ENT>) {
	$line_count++;
	&MySplit;
	my $path_id = $fields[0];
	my $entity = $fields[1];
	warn "\t$path_id\t$entity\n";
	unless ($pathway = PFBP::Pathway->get_object($path_id)) {
	    print ERR ";ERROR: file $in_file{pathway_molecule} line $line_count: pathway $path_id has not been found in pathway description file\n";
	    next;
	}
	$pathway->push_attribute("entities",$entity);
    }
    close PATH_ENT;
}

sub oldReadPathways {
    open PATHEL, "> $PathwayElementFile" || die  "Error: cannot write pathway element file $PathwayElementFile\n";

    open PATH, "$PathwayFile" || die "Error: cannot read pathway file $PathwayFile\n";
    while (<PATH>) {
	&MySplit;
	### pathway name
	$PathwayName = $fields[0];
	$PathwayName = "undef" if ($PathwayName eq "");

	### pathway accession number
	$PathwayAC = "$PathwayName";
	$PathwayAC =~ s/^\'//g;
	$PathwayAC =~ s/\'$//g;
	$PathwayAC =~ s/\'/prime/g;
	$PathwayAC =~ s/\s/_/g;
	$PathwayAC =~ s/\-/_/g;
	$PathwayAC = "path_".$PathwayAC;

	$pathway{$PathwayAC}->{ac} = $PathwayAC;
	$pathway{$PathwayAC}->{name} = $PathwayName;

	### entities
	$EntityName = $fields[1];
	$lc_name = lc($EntityName);
	if ($AC{$lc_name} ne "") {
	    $EntityAC = $AC{$lc_name};
	    $Type = $entity{$EntityAC}->{type};
	} else {
	    print  ERR ";Error in $PathwayFile: $name had not been previously declared in $in_file{molecules}\n";
	}
	&CreatePathwayEntity($PathwayAC,$EntityAC,$Type,$EntityName);
    }
    close PATH;

    ### convert associations into pathway elements
    foreach $AssocAC (keys %association) {
	$AssocAC = $association{$AssocAC}->{ac};
	$AssocType = $association{$AssocAC}->{type};
	$Label = $association{$AssocAC}->{name};
	$source = $association{$AssocAC}->{source};
	$target = $association{$AssocAC}->{target};
	foreach $PathwayAC (keys %pathway) {
	    next if ($PathwayAC eq "path_undef");
	    if (($member{$source}{$PathwayAC}) && ($member{$target}{$PathwayAC})) {
		$PwelAC = &GetNextPwelAC;
		$FromAC = $member{$source}{$PathwayAC};
		$ToAC = $member{$target}{$PathwayAC};
		print PATHEL "pathway_element($PwelAC,$PathwayAC,$Step,association,$AssocType,$AssocAC,$FromAC,$ToAC,$Label,_,_).\n";
	    }
	}
    }
    close PATHEL;
    ### split pathway_element by pathway
    foreach $PathwayAC (keys %pathway) {
	next if ($PathwayAC eq "path_undef");
	$command = "grep $PathwayAC $PathwayElementFile > ${PathwayAC}.pl";
	system $command;
    }
}

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
    print REPORT "$entity_nb Molecules\n";
    print REPORT "\tTypes of Molecules\n";
    foreach $type (sort keys %moleculeType) {
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
	if ($ARGV[$a] eq "-v") {
	    $verbose = 1;
	} elsif ($ARGV[$a] eq "-vv") {
	    $verbose = 1;
	    $hyperverbose = 1;
	} elsif ($ARGV[$a] eq "-test") {
	    $test = 1;
#    } elsif ($ARGV[$a] eq "-o") {
#      $a++;
#      $output_file = $ARGV[$a];
	} elsif ($ARGV[$a] eq "-format") {
	    $a++;
	    $out_format = $ARGV[$a];
	} elsif ($ARGV[$a] eq "-return") {
	    $a++;
	    @out_fields = split ",", $ARGV[$a];
	    
	} elsif (($ARGV[$a] eq "-h") ||
		 ($ARGV[$a] eq "-help")) {
	    &PrintHelp;
	    exit(0);
	}
    }
}


### print entities and activities in the tabular format for diagram description
sub PrintGraph {
    ### usage : &PrintGraph($graph_file, @activities);
    ### collects all inputs and outputs of these activities
    ### and prints the graphs in .tdf format
    my ($graph_file, @activities) = @_;
    my $title = "Roche T cell signalling data";
    my $xsize = 600;
    my $ysize = 600;
    my $arc_count = 0;
    my %linked_entities = ();
    srand (time);
    
    warn "; printing graph\t$graph_file\n" if ($verbose);
    
    open GRAPH, ">$graph_file" 
	|| die "Error: cannot write $out_file{graph}\n";
    
    print GRAPH "ftype\ttfd\n";
    print GRAPH "title\t$title\n";
    print GRAPH "xsize\t$xsize\n";
    print GRAPH "ysize\t$ysize\n";
    print GRAPH "creator\troche_parser\n";
    print GRAPH "cdate\t", `date`;

    ### collect inputs and outputs for each activity
    foreach my $activity (@activities) {
	foreach $input ($activity->get_attribute("inputs")) {
	    $linked_entities{$input}++;
	}
	foreach $output ($activity->get_attribute("outputs")) {
	    $linked_entities{$output}++;
	}
    }

    ### print nodes corresponding to these entities
    foreach $id (keys %linked_entities) {
	$entity = PFBP::BiochemicalEntity->get_object($id);
	$xpos = int(rand $xsize);
	$ypos = int(rand $ysize);
	$label = $entity->get_name;
	$type = $entity->get_attribute("type");

	&PrintNode;
    }
    
    ### print a node for each activity
    foreach my $activity (@activities) {
	$id = $activity->get_attribute("id");
	$xpos = int(rand $xsize);
	$ypos = int(rand $ysize);
	$type = $activity->get_attribute("type");
	$label = $type."_".$id;
	
	&PrintNode;

	### print arc between the activity node and each input
	foreach $input ($activity->get_attribute("inputs")) {
	    $arc_id = "arc_".$arc_count++;
	    $from = $input;
	    $to = $id; #### activity id
	    $label = $type = "input";
	    &PrintArc;
	}

	### print arc between the activity node and each output
	foreach $output ($activity->get_attribute("outputs")) {
	    $arc_id = "arc_".$arc_count++;
	    $from = $id; #### activity id
	    $to = $output; 
	    $label = $type = "output";
	    &PrintArc;
	}
    }

    close GRAPH;
}

sub PrintArc {
    print GRAPH "arc";
    print GRAPH "\t", $arc_id;
    print GRAPH "\t", $from;
    print GRAPH "\t", $to;
    print GRAPH "\t", $label;
    print GRAPH "\t", $type;
    print GRAPH "\n";
}

sub PrintNode {
    print GRAPH "node";
    print GRAPH "\t", $id;
    print GRAPH "\t", $xpos;
    print GRAPH "\t", $ypos;
    print GRAPH "\t", $label;
    print GRAPH "\t", $type;
    print GRAPH "\n";
}


sub LinkEntities {
    ### usage : @activites = &LinkEntities(@entities);
    ### given a set of entities, return all activities 
    ### that have at one of these as input and one as output
    ### this is a subgraph extraction with a maximum link length of 1 arc
    my @entities = @_;  ### a list of entity IDs
    my @link_activities = ();
    my %seed = ();

    ### put entity ids in a hash
    foreach my $entity (@entities) {
	$seed{uc($entity)} = 1;
    }

    #### collect activities
    foreach my $activity (PFBP::BiochemicalActivity->get_objects()) {
	my $in = 0;
	my $out = 0;
	foreach $input ($activity->get_attribute("inputs")) {
	    $in++ if ($seed{uc($input)});
	}
	foreach $output ($activity->get_attribute("outputs")) {
	    $out++ if ($seed{uc($output)});
	}
	if (($in) and ($out)) {
	    push @link_activities, $activity;
	}
    }

    return @link_activities;
}

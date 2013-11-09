#!/usr/bin/perl
############################################################
#
# $Id: parse_regulation.pl,v 1.10 2011/02/17 05:07:46 rsat Exp $
#
# Time-stamp: <2003-07-10 11:52:58 jvanheld>
#
############################################################
### parse_regulation.plt
### type parse_regulation.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "lib/load_classes.pl";
require "config.pl";
require "lib/util.pl";
require "lib/parsing_util.pl";

################################################################
#
# Main program
#
package main;
{
    ################################################################
    ### default parameters

    #### relational DBMS parameters
    $schema="annotator";
    $user="annotator";
    $password="annotator";

    #### other parameters
    $verbose = 0;
    $debug = 0;
    $clean = 0;

    ### files to parse

#    $in_file{regulation} = "/win/amaze/amaze_team/gaurab/regulation_2001_09_07.txt";
#    $in_file{regulation} = "/win/amaze/amaze_team/georges_cohen/excel_files_georges/regulation/concatenated_regulations.tab";
    if ($hostname eq "brol") {
	$in_file{regulation} = "/win/amaze/amaze_team/georges_cohen/excel_files_georges/regulation/all_regulation_corrected.txt";
    } else {
	$in_file{regulation} = "/rubens/dsk2/dgonze/annotation_regulation/all_regulations.tab";
    }

    $dir{genes} = "$parsed_data/kegg_genes/20021127";
    $dir{polypeptides} = "$parsed_data/swissprot/20021127";


    #### default output directory
    $export_subdir = "regulation";
    $dir{output} = "$parsed_data/$export_subdir/$delivery_date";

#    $dir{output} = "$parsed_data/regulation/$delivery_date";
    $dir{delivery} = "/win/amaze/amaze_programs/amaze_oracle_data";

    $out_format = "obj";


    ################################################################
    #### read command-line arguments
    &ReadArguments();

    #### output directories
    &CheckOutputDir();

#      unless (-d $dir{output}) {
#  	warn "Creating output dir $dir{output}\n"  if ($verbose >= 1);
#  	`mkdir -p $dir{output}`;
#  	die "Error: cannot create directory $dir\n" 
#  	    unless (-d $dir{output});
#      }
#      chdir $dir{output};
    $outfile{regulation} = "$dir{output}/regulation.obj";
    $outfile{stats} = "$dir{output}/regulation.stats.txt";
    $outfile{errors} = "$dir{output}/regulation.errors.txt";
    $outfile{mirror} = "$dir{output}/regulation.mirror.txt";

    ################################################################
    #### check some parameters

    #### check the existence of gene index files
    unless (-d $dir{genes}) {
	die "ERROR: gene directory  $dir{genes} does not exist\n";
    }
    $in_file{gene_names} =  "gunzip -c $dir{genes}/Gene_names.tab.gz | grep -v '^--' | grep -v 'H.sapiens' | ";
    $in_file{genes} = "gunzip -c $dir{genes}/Gene.tab.gz |";
    $in_file{genes} .= " grep -v '^--' | cut -f 1,2  |";
    $in_file{genes} .= " grep -v 'H.sapiens' |";
    $in_file{genes} .= " perl -pe 's|S.cerevisiae|Saccharomyces cerevisiae|' |";
    $in_file{genes} .= " perl -pe 's|E.coli|Escherichia coli|' |";

    #### check existence of polypeptide index files
    unless (-d $dir{polypeptides}) {
	die "ERROR: polypeptide directory  $dir{polypeptides} does not exist\n";
    }
    $in_file{polypeptide_names} = "gunzip -c $dir{polypeptides}/Polypeptide_names.tab.gz  | grep -v '^--' |";
    $in_file{polypeptides} = "gunzip -c $dir{polypeptides}/Polypeptide_organisms.tab.gz  | grep -v '^--'|";


    #### remove all files from the output directory
    if ($clean) {
	warn "; Cleaning output directory $dir{output}\n" if ($verbose >=1);
	system "\\rm -f $dir{output}/*.tab.gz  $dir{output}/*.txt.gz $dir{output}/*.obj.gz" ;
	system "\\rm -f $dir{output}/*.tab  $dir{output}/*.txt $dir{output}/*.obj" ;
    }

    open ERR, ">$outfile{errors}" || die "Error: cannot write error report file $$outfile{errors}\n";
    open MIRROR, ">$outfile{mirror}" || die "Error: cannot write mirror file $$outfile{mirror}\n";

    ### testing mode
    if ($test) {
	warn ";TEST\n" if ($verbose >= 1);
	### fast partial parsing for debugging
	$in_file{regulation} .= " head -10 |";
    }

    ### default verbose message
    &DefaultVerbose() if ($verbose >= 1);

    #### class factories
    $controlOfControls = classes::ClassFactory->new_class(object_type=>"classes::ControlOfControl",
						       prefix=>"act_");
    $transcriptionalRegulations = classes::ClassFactory->new_class(object_type=>"classes::TranscriptRegul",
								prefix=>"trr_");

    $indirectInteractions = classes::ClassFactory->new_class(object_type=>"classes::IndirectInteraction",
							  prefix=>"iin_");

    #### output fields

    @out_fields = qw(
		     pathway
		     source
		     inputType
		     input
		     
		     controlType

		     controlledFrom
		     controlledType
		     controlledTo
		     

		     is_positive
		     description

#		     type_of_inhibition
#		     strength

		     gene_name
		     gene_id
		     factor_name
		     factor_id
		     file
		     pathwayIDs 
		     coenzyme 
		     cofactor 
		     relative_concentration_of_input

		     remark 
		     pubmedIDs 
		     );
    if ($debug) {
	push @out_fields, qw(
			     input_ids
			     output_id
			     org_if_not_Ecoli 
			     org_in_addition_to_Ecoli 
			     remark_pubmedIDs
			     );
	
    }
    
    $transcriptionalRegulations->set_out_fields(@out_fields) unless ($export{all}) ;

    $inductions = classes::ClassFactory->new_class(object_type=>"classes::Induction",
						prefix=>"ind_");

    #### indexes
    &LoadIndexes();

    #### parse reactions
    &ParseRegulation($in_file{regulation});
#    &IdentifyInputsOutputs();

    ### print result
    $transcriptionalRegulations->dump_tables();
    #$controlOfControls->dump_tables();
    #$inductions->dump_tables();


    push @classes, ("classes::TranscriptRegul");
    #push @classes, ("classes::ControlOfControl");
    #push @classes, ("classes::Induction");
    &ExportClasses($outfile{regulation}, $out_format, @classes)  if ($export{obj});


    #### export SQL scripts to load the data in a relational database
    foreach $class_holder ($transcriptionalRegulations) {
	$class_holder->generate_sql(schema=>$schema, 
				    user=>$user,
				    password=>$password,
				    dir=>"$dir{output}/sql_scripts", 
				    prefix=>"");
    }

    ### print some stats after parsing
    &PrintStats($outfile{stats}, @classes);

    ### report execution time
    if ($verbose >= 1) {
	$done_time = `date +%Y-%m-%d.%H%M%S`;
	warn ";\n";
	warn "; job started $start_time";
	warn "; job done    $done_time";
    }

    close ERR;
    close MIRROR;

    &deliver() if ($deliver);

    &CompressParsedData();

    exit(0);


}

################################################################
#
# SUBROUTINES
#


################################################################
### print the help message 
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
        parse_regulation.pl

DESCRIPTION

	Parse Regulation excel sheets (this is a temporary format for
	annotation before we have real annotation tools).

AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  

SYNOPSIS	
	parse_regulation.pl [-h] [-help] [-test] [-v #] 
		[-i infile] [-o outfile]

OPTIONS	
$generic_option_message
	-d	debug
		export additional info for debugging (e.g. indexes
	-genes  gene_dir
	        parsed genes directory. 
	        This directory should contain the result of the script
	        parse_genes.pl. Parsed genes are used for identifying
	        the  target genes.
	-pp	polypeptide_dir
	        parsed polypeptides directory. 
	        This directory should contain the result of the script
	        parse_swissprot.pl. Parsed polypeptides are used for
	        identifying the transcription factors.
	-deliver 
		copy the relevant files in the delivery directory,
		for loading in the aMAZE database. 
	-i	input file
		If ommited, STDIN is used
		This allows to insert the program within a unix pipe
	-o	output file
		If ommited, STDOUT is used. 
		This allows to insert the program within a unix pipe
EndHelp
  close HELP;
}

  


################################################################
### read arguments from the command line
sub ReadArguments {
    my $a = 0;
    while ($a <= $#ARGV) {
#    for $a (0..$#ARGV) {

	&ReadGenericOptions($a);

	    ### debug
	if ($ARGV[$a] eq "-d") {
	    $main::debug = 1;

	    ### deliver
	} elsif ($ARGV[$a] eq "-deliver") {
	    $main::deliver = 1;

	    ### gene directory
	} elsif ($ARGV[$a] eq "-genes") {
	    $a++;
	    $main::dir{genes} = $ARGV[$a];

	    ### gene directory
	} elsif ($ARGV[$a] eq "-pp") {
	    $a++;
	    $main::dir{polypeptides} = $ARGV[$a];

	    ### export all fields
	} elsif ($ARGV[$a] eq "-all") {
	    $main::export{all} = 1;

	    ### input file
	} elsif ($ARGV[$a] eq "-i") {
	    $a++;
	    $in_file{regulation} = $ARGV[$a];
	    
	    ### output file
	} elsif ($ARGV[$a] eq "-o") {
	    $a++;
	    $outfile{regulation} = $ARGV[$a];
	}

	$a++;
    }
}

################################################################
#
# Parse the tab-delimited file containing the description of the
# regulatory interactions
#
sub ParseRegulation {
    my ($input_file) = @_;
    my $filename = `basename $input_file`;

    warn "; Parsing file $input_file\n" if ($verbose >= 1);

    chomp $filename;
    
    my %class_holder;
    $class_holder{activates} = $controlOfControls;
    $class_holder{inhibits} = $controlOfControls;
    $class_holder{'down-regulates'} = $transcriptionalRegulations;
    $class_holder{'up-regulates'} = $transcriptionalRegulations;
    $class_holder{'down-regulates'} = $transcriptionalRegulations;
    $class_holder{'up-regulates'} = $transcriptionalRegulations;
    $class_holder{induces} = $inductions;
    
    my %is_positive;
    $is_positive{activates} = "1";
    $is_positive{inhibits} = "0";
    $is_positive{'up-regulates'} = "1";
    $is_positive{'down-regulates'} = "0";
    $is_positive{induces} = "1";


#    foreach my $c (0..$#cols) {
#	$col{$cols[$c]} = $c;
#    }
    my $l = 0;
    
    ### open the input stream
    if ($input_file) {
	open REG, $input_file 
	    || die "Error: cannot read $input_file\n";
	$in = REG;
    } else {
	$in = STDIN;
    } 

    ################################################################
    #### read header line and parse it to know which column contains
    #### which information
    my $header = <$in>;
    chomp($header);

    my @cols = split "\t", $header;
    for my $c (0..$#cols) {
	$cols[$c] =~ s/ +/ /g;
	$cols[$c] =~ s/ $//g;
	$cols[$c] =~ s/^ //g;
	$cols[$c] =~ s/ /_/g;
    }


    my @single_value = qw( 
			   pathway
			   row_nb
			   source
			   inputType
			   input
			   controlType
			   type_of_inhibition
			   strength
			   controlledType
			   controlledFrom
			   controlledTo
			   org_not_Ecoli
			   );

    my %single_value = ();
    foreach my $key (@single_value) {
	$single_value{$key} = 1;
    }

    foreach my $key (@cols) {
	push @multi_value, $key unless $single_value{$key};
    }


    #### report column contents
    if ($verbose >= 2) {
	warn "; Column contents\n";
	for my $c (1..$#cols) {
	    warn join ("\t", ";", $c, $cols[$c]), "\n";
	}
    }

    @mirror_fields = ('pathway',
		      'pubmedIDs',
		      'source',
		      'inputType',
		      'input',
		      'enzyme_or_protein_controlled',
		      'controlType',
		      'controlledType',
		      'type_of_inhibition',
		      'controlledFrom',
		      'controlledTo',
		      'remark',
		      'remark2',
		      'remark3',
		      'remark4',
		      'remark_4',
		      'method',
		      'parsing_errors');
    print MIRROR join( "\t", @mirror_fields), "\n";
    
    
    #### read the interactions
    while (my $current_line = <$in>) {
	$l++;
	next if ($current_line =~ /^;/);
	chomp $current_line;
	@fields = split "\t", $current_line;

	@current_errors = ();

	%field = ();

	foreach $c (0..$#cols) {
	    $key = $cols[$c];
	    $value = $fields[$c];
	    $value =~ s/^\"//;
	    $value =~ s/\"$//;
	    $field{$key} = $value;
	    warn join ("\t", 
		       $filename,
		       $l,
		       $c, 
		       $key, 
		       $value,
		       $field{$key}
		       ), "\n" if ($verbose >=4);
	}


	#### merge the remarks (Georges annotated them in 4 columns)
	if ($field{remark2} =~ /\S/) {
	    warn join "\t", $l, "merging remark2\n" if ($verbose >= 3);
	    $field{remark} .= "|" if $field{remark};
	    $field{remark} .= $field{remark2};
	}
	if ($field{remark3} =~ /\S/) {
	    warn join "\t", $l, "merging remark3\n" if ($verbose >= 3);
	    $field{remark} .= "|" if $field{remark};
	    $field{remark} .= $field{remark3};
	}
	if ($field{remark4} =~ /\S/) {	
	    warn join "\t", $l, "merging remark4\n" if ($verbose >= 3);
	    $field{remark} .= "|" if $field{remark};
	    $field{remark} .= $field{remark4};
	}
	if ($field{'remark 4'} =~ /\S/) {
	    warn join "\t", $l, "merging remark 4\n" if ($verbose >= 3);
	    $field{remark} .= "|" if $field{remark};
	    $field{remark} .= $field{'remark 4'};
	}
	


	#### check the controlType
	unless (defined($class_holder{$field{controlType}})) {
	    push @current_errors, "unknown control type";
	    &ErrorMessage(join ("\t", 
				$filename,
				$l,
				"unknown control type",
				$field{controlType}), "\n");
	    &PrintMirrorLine();
	    next;
	}
	
	#### select the class for the new interaction
	my $class_holder = $class_holder{$field{controlType}};
	
	if (($class_holder == $transcriptionalRegulations) && 
	    (lc($field{inputType}) ne "protein")) {
	    $class_holder = $indirectInteractions;
	    push @current_errors, "transcriptional regulation with inputType different from protein";
	    &ErrorMessage(join ("\t", 
				$filename,
				$l,
				"transcriptional regulation with inputType different from protein",
				$field{controlledFrom}), "\n");
	    &PrintMirrorLine();
	    next;
	}


	#### when there are multiple outputs, create one regulation per output
	my @control_outputs = split /\|/, $field{controlledFrom};

	unless ($#control_outputs >= 0) {
	    push @current_errors, "the field controlledFrom is empty";
	    &ErrorMessage(join ("\t", 
				$filename,
				$l,
				"the field controlledFrom is empty",
				$field{controlledFrom}), "\n");
	    &PrintMirrorLine();
	    next;
	}


	warn ($field{controlledFrom}, "\t",
	      $#control_outputs + 1, "\t",
	      join ("\t", @control_outputs), "\n" )
	    if ($verbose >=4);
	
	foreach my $control_output (@control_outputs) {
	    $field{controlledFrom} = $control_output;
	    my $current_control = $class_holder->new_object();
	    
	    #### assign single-value attributes
	    foreach my $key (@single_value) {
		$current_control->set_attribute($key, $field{$key});
	    }
	    
	    #### assign multi-value attributes
	    foreach my $key  (@multi_value) {
		my @values = split '\|', $field{$key};
		foreach my $value (@values) {
		    if ($value =~ /\S/) {
			$current_control->new_attribute_value($key, $value);
		    }
		}
	    }
	    
	    #### change the names of attributes to have an uniform source
	    $current_control->set_attribute("biblio_source", $current_control->get_attribute("source"));
	    $current_control->force_attribute("source", $filename); #### for error report
	    $current_control->set_attribute("line", $l);        #### for error report
	    $current_control->set_attribute("is_positive", $is_positive{$field{controlType}});

	    #### bibliographic references
	    my @pubmedIDs = split '\|', $field{pubmedIDs};
	    $current_control->empty_array_attribute("pubmedIDs");
	    foreach my $id (@pubmedIDs) {
		$id = &trim($id);
		if ($id =~ /^\d+$/) {
		    $current_control->new_attribute_value("pubmedIDs", $id);
		} else {
		    $current_control->new_attribute_value("biblioref", $id);
		}
	    }

	    #### report the content of the object
	    warn join ("\t", 
		       $current_control->get_attribute("id"),
		       $current_control->get_attribute("controlType"),
		       ), "\n"
		if ($verbose >= 2);

	    #### check inputs/outputs for transcriptional regulations
	    if ($class_holder == $transcriptionalRegulations) {
		&IdentifyInputsOutputs($current_control);
	    }
	}

	&PrintMirrorLine();
    }
    
    close $in if ($input_file);
}


sub PrintMirrorLine {
    $field{'parsing_errors'} = join (";", @current_errors);
    my @mirror_line = ();
    foreach my $f (@mirror_fields) {
	push @mirror_line, $field{$f}; 
    }
    print MIRROR join ("\t", @mirror_line), "\n";
}

################################################################
#
# Load index files
#
sub LoadIndexes {
    #### load indexes
    warn ("; Loading index files\n") 
	if ($verbose >= 1);


    #### gene names
    $index{name_gene} = classes::Index->new();
    warn ("; ", &AlphaDate(), "\tgene name index ...\n") 
	if ($verbose >= 1);
    $index{name_gene}->load($in_file{gene_names}, 1, 0, reverse=>1);

    #### gene organism
    $index{gene_organism} = classes::Index->new();
#    $in_file{genes} .= " perl -pe 's|H.sapiens|Homo sapiens|' |";
    warn ("; ", &AlphaDate(), "\tgene organism index ...\n") 
	if ($verbose >= 1);
    $index{gene_organism}->load($in_file{genes}, 0, 0);


    #### polypeptide names from swissprot
    $index{name_polypeptide} = classes::Index->new();
    warn ("; ", &AlphaDate(), "\tpolypeptide name index ...\n") 
	if ($verbose >= 1);
    $index{name_polypeptide}->load($in_file{polypeptide_names}, 1, 0, reverse=>1);

    #### polypeptide organisms from swissprot
    $index{polypeptide_organism} = classes::Index->new();
    warn ("; ", &AlphaDate(), "\tpolypeptide organism index ...\n") 
	if ($verbose >= 1);
    $index{polypeptide_organism}->load($in_file{polypeptides}, 0, 1);

    #### report index sizes
    if ($verbose >=1) {
	warn ("; Index sizes\n");
	foreach $key (keys %index) {
	    warn (";\t", 
		  "\t", $index{$key}->get_size(), 
		  "\t", $key,
		  "\n");
	}
    }

    #### export the indexes for checking their contents
    if ($debug) {
	foreach $key (keys %index) {
	    $index{$key}->export("tab", "$dir{output}/index_${key}.tab");
	}
    }
}




################################################################
#
# Identify target genes (in KEGG genes) and transcription factors
# (with Swissprot)
#
sub IdentifyInputsOutputs {
#    foreach my $trreg ($transcriptionalRegulations->get_objects()) {
    
    my ($trreg) = @_;

    warn ";\tIdentifying inputs and outputs for transcriptional regulation", $trreg->get_attribute("id"), "\n" if ($verbose >= 3);

    ################################################################
    #### matching factors
    my $factor_name = $trreg->get_attribute("input");

    $trreg->set_attribute("factor_id", $null);
    $trreg->set_attribute("factor_name", $factor_name);

    my @matching_factors = $index{name_polypeptide}->get_values(&standardize($factor_name));

    if (($#mathing_factors < 0) && ($factor_name =~ /p$/)){
	my $truncated_factor_name = $factor_name;
	$truncated_factor_name =~ s/p$//;
	@matching_factors = $index{name_polypeptide}->get_values(&standardize($truncated_factor_name));
    }



    
    if ($#matching_factors < 0) {
	push @current_errors, "cannot identify polypeptide";
	&ErrorMessage(join ("\t", 	
			    $trreg->get_attribute("source"),
			    $trreg->get_attribute("line"),
			    "Cannot identify polypeptide", 
			    $factor_name), "\n");
    } else {
	if ($#matching_factors > 0) {
	    push @current_errors, join (" ", "Multiple matching factors", $factor_name, @matching_factors);
	    &ErrorMessage (join ("\t", 	
				 $trreg->get_attribute("source"),
				 $trreg->get_attribute("line"),
				 "Multiple matching factors", $factor_name, @matching_factors), "\n");
	}
	$trreg->set_attribute("factor_id", $matching_factors[0]);
	
	foreach my $factor_id (@matching_factors) {
	    warn ("adding polypeptide $factor_id as input for $trreg\n") 
		if ($verbose >=3);
	    $trreg->new_attribute_value("input_ids", $factor_id);
	}
    }

    ################################################################
    #### matching genes
    my $gene_name = $trreg->get_attribute("controlledFrom");
    $trreg->set_attribute("gene_id", $null);
    $trreg->set_attribute("gene_name", $gene_name);
    my @matching_genes = $index{name_gene}->get_values(&standardize($gene_name));

    warn join( "\t", $gene_name, &standardize($gene_name), join( "|", @matching_genes)), "\n" if ($verbose >= 3);

    if ($#matching_genes < 0) {
	push @current_errors, "cannot identify gene";
	&ErrorMessage(join ("\t",
			    $trreg->get_attribute("source"),
			    $trreg->get_attribute("line"),
			    "Cannot identify gene",
			    $gene_name), "\n");
    } else {
	if ($#matching_genes > 0) {
	    push @current_errors, join (" ", "Multiple matching genes", $gene_name, @matching_genes);
	    &ErrorMessage (join ("\t", 	
				 $trreg->get_attribute("source"),
				 $trreg->get_attribute("line"),
				 "Multiple matching genes", $gene_name, @matching_genes), "\n");
	}
	$trreg->set_attribute("gene_id", $matching_genes[0]);
	
	foreach my $gene_id (@matching_genes) {
	    warn ("adding polypeptide $gene_id as output for $trreg\n") 
		if ($verbose >=3);
	    $trreg->new_attribute_value("output_ids", $gene_id);
	}
    }

    ################################################################
    #### description
    my $description = $factor_name;
    $description .= " ".$trreg->get_attribute("controlType");
    $description .= " ".$gene_name;
    $trreg->set_attribute("description", $description);


# }
}

#!/usr/bin/perl
############################################################
#
# $Id: parse_regulation.pl,v 1.5 2002/11/25 19:03:19 jvanheld Exp $
#
# Time-stamp: <2002-11-25 13:03:00 jvanheld>
#
############################################################
### parse_regulation.plt
### type parse_regulation.pl -h for info

#use strict;
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "PFBP_classes.pl";
require "PFBP_config.pl";
require "PFBP_util.pl";
require "PFBP_parsing_util.pl";

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
    $in_file{regulation} = "/win/amaze/amaze_team/georges_cohen/excel_files_georges/regulation/concatenated_regulations.tab";


    $dir{genes} = "$parsed_data/kegg_genes/20020822";
    $dir{polypeptides} = "$parsed_data/swissprot/20021016";


    $dir{output} = "$parsed_data/regulation/$delivery_date";
    $dir{delivery} = "/win/amaze/amaze_programs/amaze_oracle_data";
    unless (-d $dir{output}) {
	warn "Creating output dir $dir{output}\n"  if ($verbose >= 1);
	`mkdir -p $dir{output}`;
	die "Error: cannot create directory $dir\n" 
	    unless (-d $dir{output});
    }
    chdir $dir{output};
    $out_file{regulation} = "$dir{output}/regulation.obj";
    $out_file{stats} = "$dir{output}/regulation.stats.txt";
    $out_file{errors} = "$dir{output}/regulation.errors.txt";

    $out_format = "obj";


    ################################################################
    #### read command-line arguments
    &ReadArguments();


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

    open ERR, ">$out_file{errors}" || die "Error: cannot write error report fle $$out_file{errors}\n";

    ### testing mode
    if ($test) {
	warn ";TEST\n" if ($verbose >= 1);
	### fast partial parsing for debugging
	$in_file{regulation} .= " head -10 |";
    }

    ### default verbose message
    &DefaultVerbose() if ($verbose >= 1);

    #### class factories
    $controlOfControls = PFBP::ClassFactory->new_class(object_type=>"PFBP::ControlOfControl",
						       prefix=>"act_");
    $transcriptionalRegulations = PFBP::ClassFactory->new_class(object_type=>"PFBP::TranscriptionalRegulation",
								prefix=>"trr_");

    $indirectInteractions = PFBP::ClassFactory->new_class(object_type=>"PFBP::IndirectInteraction",
							  prefix=>"iin_");

    #### output fields

    @out_fields = qw(
		     source
		     inputType
		     factor_name
		     factor_id
		     controlType
		     gene_name
		     gene_id
		     is_positive
		     description
		     pathway

		     pubmedIDs 
		     );
    if ($debug) {
	push @out_fields, qw(
			     file
			     input
			     input_ids
			     output_id
			     controlledFrom
			     controlledType
			     controlledTo
			     strength
			     type_of_inhibition
			     coenzyme 
			     cofactor 
			     org_if_not_Ecoli 
			     org_in_addition_to_Ecoli 
			     pathwayIDs 
			     relative_concentration_of_input
			     remark 
			     remark_pubmedIDs);
	
    }
    
    $transcriptionalRegulations->set_out_fields(@out_fields) unless ($export{all}) ;

    $inductions = PFBP::ClassFactory->new_class(object_type=>"PFBP::Induction",
						prefix=>"ind_");

    #### indexes
#    &LoadIndexes();

    #### parse reactions
    &ParseRegulation($in_file{regulation});
    &IdentifyInputsOutputs();

    ### print result
    $transcriptionalRegulations->dump_tables();
    #$controlOfControls->dump_tables();
    #$inductions->dump_tables();


    push @classes, ("PFBP::TranscriptionalRegulation");
    #push @classes, ("PFBP::ControlOfControl");
    #push @classes, ("PFBP::Induction");
    &ExportClasses($out_file{regulation}, $out_format, @classes)  if ($export{obj});


    #### export SQL scripts to load the data in a relational database
    foreach $class_holder ($transcriptionalRegulations) {
	$class_holder->generate_sql(schema=>$schema, 
				    user=>$user,
				    password=>$password,
				    dir=>"$dir{output}/sql_scripts", 
				    prefix=>"");
    }

    ### print some stats after parsing
    &PrintStats($out_file{stats}, @classes);

    ### report execution time
    if ($verbose >= 1) {
	$done_time = `date +%Y-%m-%d.%H%M%S`;
	warn ";\n";
	warn "; job started $start_time";
	warn "; job done    $done_time";
    }

    close ERR;

    &deliver() if ($deliver);

#    system "gzip -f $dir{output}/*.tab $dir{output}/*.txt";
#    system "gzip -f $dir{output}/*.obj" if ($export{obj});




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

	### warn level
	if (($ARGV[$a] eq "-v" ) && 
	    ($ARGV[$a+1] =~ /^\d+$/)){
	    $main::verbose = $ARGV[$a+1];
	    $a++;

	    ### fast test
	} elsif ($ARGV[$a] eq "-test") {
	    $test = 1;

	    ### debug
	} elsif ($ARGV[$a] eq "-d") {
	    $main::debug = 1;

	    ### clean
	} elsif ($ARGV[$a] eq "-clean") {
	    $main::clean = 1;

	    ### clean
	} elsif ($ARGV[$a] eq "-deliver") {
	    $main::deliver = 1;

	    ### output file
 	} elsif ($ARGV[$a] eq "-obj") {
	    $a++;
	    $main::export{obj} = 1;
	    
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
	    $out_file{regulation} = $ARGV[$a];

	    ### help
	} elsif (($ARGV[$a] eq "-h") ||
		 ($ARGV[$a] eq "-help")) {
	    &PrintHelp;
	    exit(0);

	} else {
	    warn "WARNING: unknown option $ARGV[$a] is ignored\n";
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
    $header = <$in>;
    chomp($header);
    #### fix some inconsistencies between Gaurab's and Georges' headers
    $header =~ s/controlled from/controlledFrom/;
    $header =~ s/controlled to/controlledTo/;
    $header =~ s/PubmedID/pubmedID/;
    $header =~ s/type of inhibition/type_of_inhibition/;

    
    my @cols = split "\t", $header;


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


#      my @cols = ();
#      my %col = ();
#      push @cols, "pathwayIDs";
#      push @cols, "pubmedIDs";
#      push @cols, "source";
#      push @cols, "inputType";
#      push @cols, "input";
#      push @cols, "controlType";
#      push @cols, "controlledType";
#      push @cols, "controlledFrom";
#      push @cols, "controlledTo";
#      push @cols, "ligand1";
#      push @cols, "ligand2";
#      push @cols, "org_not_Ecoli";
    

    #### read the interactions
    while (<$in>) {
	$l++;
	next if (/^;/);
	chomp;
	@fields = split "\t";
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

	unless (defined($class_holder{$field{controlType}})) {
	    &ErrorMessage(join ("\t", 
				$filename,
				$l,
				"unknown control type",
				$field{controlType}), "\n");
	    next;
	}
	
	my $class_holder;
	if (lc($field{inputType}) eq "compound") {
	    $class_holder = $indirectInteractions;
	} else {
	    $class_holder = $class_holder{$field{controlType}};
	}

	#### when there are multiple outputs, create one regulation per output
	my @control_outputs = split /\|/, $field{controlledFrom};

	warn join "\t", "HELLO", "output",  $field{controlledFrom}, "\n";

	unless ($#control_outputs >= 1) {
	    &ErrorMessage(join ("\t", 
				$filename,
				$l,
				"the field controlledFrom is empty",
				$field{controlledFrom}), "\n");
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
		    $current_control->new_attribute_value($key, $value);
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
	}
    }
    
    close $in if ($input_file);
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
    $index{name_gene} = PFBP::Index->new();
    warn ("; ", &AlphaDate(), "\tgene name index ...\n") 
	if ($verbose >= 1);
    $index{name_gene}->load($in_file{gene_names}, 1, 0, reverse=>1);

    #### gene organism
    $index{gene_organism} = PFBP::Index->new();
#    $in_file{genes} .= " perl -pe 's|H.sapiens|Homo sapiens|' |";
    warn ("; ", &AlphaDate(), "\tgene organism index ...\n") 
	if ($verbose >= 1);
    $index{gene_organism}->load($in_file{genes}, 0, 0);


    #### polypeptide names from swissprot
    $index{name_polypeptide} = PFBP::Index->new();
    warn ("; ", &AlphaDate(), "\tpolypeptide name index ...\n") 
	if ($verbose >= 1);
    $index{name_polypeptide}->load($in_file{polypeptide_names}, 1, 0, reverse=>1);

    #### polypeptide organisms from swissprot
    $index{polypeptide_organism} = PFBP::Index->new();
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
# Check whether factor-coding genes and target genes have an ID in
# KEGG
#
sub IdentifyInputsOutputs {
    foreach my $trreg ($transcriptionalRegulations->get_objects()) {
	
	#### matching factors
	my $factor_name = $trreg->get_attribute("input");
	$trreg->set_attribute("factor_id", $null);
	$trreg->set_attribute("factor_name", $factor_name);
	my @matching_factors = $index{name_polypeptide}->get_values(&standardize($factor_name));
	
	if ($#matching_factors < 0) {
	    &ErrorMessage(join ("\t", 	
				$trreg->get_attribute("source"),
				$trreg->get_attribute("line"),
				"Cannot identify polypp", 
				$factor_name), "\n");
	} else {
	    if ($#matching_factors > 0) {
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

	#### matching genes
	my $gene_name = $trreg->get_attribute("controlledFrom");
	$trreg->set_attribute("gene_id", $null);
	$trreg->set_attribute("gene_name", $gene_name);
	my @matching_genes = $index{name_gene}->get_values(&standardize($gene_name));

	if ($#matching_genes < 0) {
	    &ErrorMessage(join ("\t",
				$trreg->get_attribute("source"),
				$trreg->get_attribute("line"),
				"Cannot identify gene",
				$gene_name), "\n");
	} else {
	    if ($#matching_genes > 0) {
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

	#### description
	my $description = $factor_name;
	$description .= " ".$trreg->get_attribute("controlType");
	$description .= " ".$gene_name;
	$trreg->set_attribute("description", $description);
    }
}

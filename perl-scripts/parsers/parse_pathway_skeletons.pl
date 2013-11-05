#!/usr/bin/perl
############################################################
#
# $Id: parse_pathway_skeletons.pl,v 1.16 2011/02/17 05:07:46 rsat Exp $
#
# Time-stamp: <2003-07-10 11:52:59 jvanheld>
#
############################################################

if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`"); ### add the program's directory to the lib path
}
require "config.pl";
require "lib/load_classes.pl";
require "lib/util.pl";
require "lib/parsing_util.pl";
require "lib/loading_util.pl";

package main;
{
    ################################################################
    #### default parameters

    #### relational DBMS parameters
    $schema="annotator";
    $user="annotator";
    $password="annotator";

    #### other parameters
    $fast_index = 1;
    $debug = 0;
    $clean = 0;
    $mirror = 0;
    
    $reverse = 0;
    $forward = 1;
    
#    $dir{amaze_export} = "/win/amaze/amaze_data/exported_2001_0719";
    $dir{kegg_parsed} = "$parsed_data/kegg_ligand/20021028";
#    $dir{georges_export} = "/win/amaze/amaze_team/for_georges";
#    $dir{skeletons} = "/win/amaze/amaze_team/gaurab/excel_tables/Pathway_ECno/Amino_acids/Amino_acid_biosyn";
#    $dir{skeletons} = "/win/amaze/amaze_team/gaurab/2001_07/Pathway_ECno/";
#    $dir{skeletons} = "/win/amaze/amaze_team/gaurab/excel_tables/Pathway_ECno/";

    $in_file{georges_identified_compounds} = "/win/amaze/amaze_team/georges_cohen/excel_files_georges/unidentified_compounds/identified_compounds_2001_1007.tab";


    $export_subdir = "pathway_skeletons";

    $dir{output} = "$parsed_data/$export_subdir/$delivery_date";

    $out_format = "obj";

    #### class factory
    $processes = classes::ClassFactory->new_class(object_type=>"classes::Process",
					       prefix=>"proc_");
    $leaves = classes::ClassFactory->new_class(object_type=>"classes::ProcessLeaf",
					    prefix=>"leaf_");


    push @classes, ("classes::Process");
    push @classes, ("classes::ProcessLeaf");

    &ReadArguments();

    #### output directories
    &CheckOutputDir();

    if ($mirror) {
	$dir{mirror} = $dir{output}."/mirror";
	unless (-d $dir{mirror}) {
	    warn "; Creating mirror dir $dir{mirror}\n"  if ($verbose >=1);
	    `mkdir -p $dir{mirror}`;
	    die "Error: cannot create mirror directory $dir{mirror}\n" unless (-d $dir{mirror});
	}
    }
    
    $dir{seed_files} = "$dir{output}/amaze_pathway_seeds";
#    $dir{seed_files} = "/win/amaze/amaze_programs/amaze_graph_analysis/amaze_pathway_seeds";

    $outfile{process_skeleton} = "$dir{output}/process_skeleton.obj";
    $outfile{stats} = "$dir{output}/process_skeleton.stats.txt";
    $outfile{errors} = "$dir{output}/process_skeleton.errors.txt";
    
    #### open error stream
    
    open ERR, ">$outfile{errors}" || die "Error: cannot write error report fle $$outfile{errors}\n";

    unless (-e $outfile{errors}) {
	die "Error: cannot open error file $outfile{errors} with write access";
    }

    my $files_to_parse = `find $dir{skeletons} -follow -name "*.txt" -o -name "*.tab" | xargs`;
    $files_to_parse =~ s|$dir{skeletons}/*||g;
    @files_to_parse = split " ", $files_to_parse;

    unless (-d $dir{skeletons}) {
	die "Error: input directory $dir{skeletons} does not exist\n";
    }

    #### quick test 
    if ($test) {
	@files_to_parse = qw (
			      Escherichia_coli/TCA_glyoxylate_bypass.txt 
			      Saccharomyces_cerevisiae/Isoleucine_and_valine_biosynth.txt
			      );
    }

    #### fields to export
    my @process_out_fields = qw( 
				 id
				 source
				 organism
				 description
				 status

				 names
				 nodes
				 arcs
				 
				 reactions
				 ECs
				 genes
				 substrates
				 products
				 );
    
    
    if ($debug) {
	push @process_out_fields, qw(
				     complete
				     steps
				     matching_reactions
				     problems
				     no_match_reactions
				     duplicate_reactions
				     multi_match_reactions
				     unidentified_compounds
				     missing_direction
				     parsed_file				 
				     );
    }
    $processes->set_out_fields(@process_out_fields);
    $processes->set_attribute_header("nodes", join ("\t", "nodes","step"));
    $processes->set_attribute_header("ECs", join ("\t", "ec","step"));
    $processes->set_attribute_header("genes", join ("\t", "genes","step"));
    $processes->set_attribute_header("products", join ("\t", "products","step"));
    $processes->set_attribute_header("substrates", join ("\t", "substrates","step"));
    $processes->set_attribute_header("arcs", join ("\t", "source", "target") );
    $processes->set_attribute_header("reactions", join ("\t", "reactions","step"));
    
    $leaves->set_out_fields(qw( 
				id
				interaction
				is_direct
				description
				));
    
    ### default verbose message
    if ($verbose >=1) {
	&DefaultVerbose();
	warn "; Files to parse\n\t", join (";\n\t", @files_to_parse), "\n";
    }


    &LoadIndexes();
#    &ReadReactions();
    &ReadProcesses();
    
    $processes->dump_tables();
    $leaves->dump_tables();
    foreach $class_holder ($processes, $leaves) {
	$class_holder->generate_sql(schema=>$schema, 
				    user=>$user,
				    password=>$password,
				    dir=>"$dir{output}/sql_scripts", 
				    prefix=>"");
    }
    &ExportClasses($outfile{process_skeleton}, $out_format, @classes) if ($export{obj});
    &PrintStats($outfile{stats}, @classes);

    &GenerateECSeeds() if ($seeds);
    close ERR;


    ### report execution time
    if ($verbose >=1) {
	$done_time = &AlphaDate;
	warn (";\n",
	      "; job started $start_time",
	      "; job done    $done_time\n")
	}
    
    &CompressParsedData();
    
    exit(0);
}


################################################################
#### export lists of EC numbers for the pathway builder
sub GenerateECSeeds {
    warn "; Generating EC seeds for pathway reconstruction\n"
	if ($verbose >=1);

    #### create directory for storing seed files
    unless (-d $dir{seed_files}) {
	warn "; Creating seed_files dir $dir{seed_files}\n"  if ($verbose >=1);
	`mkdir -p  $dir{seed_files}`;
	die "Error: cannot create directory $dir{seed_files}\n" unless (-d $dir{seed_files});
    }
    
    #### export one seed file per pathway
    foreach $process ($processes->get_objects()) {
	my $organism_prefix = $process->get_attribute("organism");
	$organism_prefix =~ s/\s/_/g;

	#### check the existence of the directory for exporting seed files
	my $export_dir = join "/", $dir{seed_files}, $organism_prefix;
	unless (-d $export_dir) {
	    `mkdir -p  $export_dir`;
	    unless (-d $export_dir) {
		die "Cannot create seed directory $export_dir";
	    }
	}

#	my $file = $process->get_attribute("parsed_file");
#	my $filename = `basename $file .txt`;
#	$filename =~ s/\.tab$//;
#	chomp $filename;
	my @names =  $process->get_attribute("names");
	my $filename = $names[0];
	die join ("\t", "Pathway has no name ", 
		  $process->get_attribute("parsed_file")),"\n" 
		      unless $filename;

	$filename =~ s/\s/_/g;
	my $export_file = join "/", $export_dir, $filename;

	#### export EC number as seeds for pathway reconstruction
#	$export_file .= ".".$process->get_attribute("complete");
	open ECS, ">$export_file.ecs";
	foreach my $ec_ref ($process->get_attribute(ECs)) {
	    my ($ec, $step) = @{$ec_ref};
	    print ECS "$ec\n";
###	    print ECS "${ec}*\n"; #### for acccepting the reverse reaction
	}
	close ECS;
	
	#### export reaction IDs for validation
	#### only export reaction for pathways where all 
	#### reactions are identified in an unequivocal way.

#	next unless ($process->get_attribute("status") eq "OK");
	open REACTIONS, ">$export_file.reactions";
	foreach my $reaction_ref ($process->get_attribute(reactions)) {
	    my ($reaction, $step) = @{$reaction_ref};
	    print REACTIONS "$reaction\n";
	}
	close ECS;
    }
    
}

################################################################
### print the help message 
sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
	parse_process_skeleton.pl

DESCRIPTION
	Parse process reactions from the excel sheets.
	This is a temporary fix until we have real annotation tools.

AUTHOR
	Jacques van Helden (jvanheld\@ucmb.ulb.ac.be)  
	Created		2001/07/20
	
OPTIONS
$generic_option_message
	-d	debug
		add attributes to the process to reflect the 
		parsed data, which ar not required for the actual 
		aMAZE process
	-mirror export a copy of the data files with additional
		columns to indicate whether each equation was
		identified or not.
	-skel  skel_dir
	       skeleton directory
	       Directory with the input files. See below for a
	       description of directory organization and file formats.
	-kegg  kegg_dir
	       parsed kegg directory. 
	       This directory should contain the result of the script
	       parse_kegg.pl. KEGG parsed files are used for
	       identifying compounds and reactions involved in each
	       pathway.
	-seeds
	       Export pathway seeds. Pathway seeds are lists of
	       reactions, which can be used by the pathway builder to
	       draw a pathway diagram.

FILES AND DIRECTORIES

	The skeleton directory should contain one sudirectory per
	organism (the name of the subdirectory is the name of
	organism, wich spaces replaced by underscores).

	Within each organism subdirectory, there is one file per
	pathway.

INPUT FORMAT
        
	Each pathway skeleton comes in a separate tab-delimited file.

	Column  contents
	----------------
	1	step
	2	EC number
	3	gene
	4	from_comp
	5	to_comp
	6	equation of the reaction

	(in George Cohen\'s files, there are some additional columns,
	which are currently ignored)
	7	Name of enzyme
	8	Coenzyme
	9	Remarks
	10	Chromosome
	11	Systematic deletion
	12	pubmedID
EndHelp
  close HELP;
}

################################################################
### read arguments from the command line
sub ReadArguments {
    for my $a (0..$#ARGV) {
	
	&ReadGenericOptions($a);

	################################################################
	#### specific options

	### debug
	if ($ARGV[$a] eq "-d") {
	    $main::debug = 1;

	    ### mirror
	} elsif ($ARGV[$a] eq "-mirror") {
	    $main::mirror = 1;

	    ### skeleton directory
	} elsif ($ARGV[$a] eq "-skel") {
	    $main::dir{skeletons} = $ARGV[$a+1];

	    ### skeleton directory
	} elsif ($ARGV[$a] eq "-kegg") {
	    $main::dir{kegg_parsed} = $ARGV[$a+1];

	    ### ec seeds
	} elsif ($ARGV[$a] eq "-seeds") {
	    $main::seeds = 1;
	}
    }
}

################################################################
#### load indexes for reference reactions and compounds (these should
#### have been previously parsed from KEG)
sub LoadIndexes {
    warn ("Loading index files\n") if ($verbose >=1);

    unless (-d $dir{kegg_parsed}) {
	die "Error: cannot find KEGG directory\t$dir{kegg_parsed}\n";
    }

    #### compound names
    warn ("\tindexing compound names ...\n") if ($verbose >=1);
    $index{compound_name} = classes::Index->new();
    $index{compound_name}->load("gunzip -c $dir{kegg_parsed}/Compound_names.tab.gz | grep -v '^--' |", 0, 1);
    $index{name_compound} = $index{compound_name}->reverse();
    
    #### reactions <-> EC numbers
    warn ("\tindexing reactions <-> EC numbers ...\n") if ($verbose >=1);
    $index{reaction_ec} = classes::Index->new();
    $index{reaction_ec}->load("gunzip -c $dir{kegg_parsed}/Reaction_ecs.tab.gz | grep -v '^--' |", 0, 0);
    $index{ec_reaction} = $index{reaction_ec}->reverse();
    
    #### reactions <-> substrates
    warn ("\tindexing reactions <-> substrates ...\n") if ($verbose >=1);
    $index{reaction_substrate} = classes::Index->new();
    $index{reaction_substrate}->load("gunzip -c $dir{kegg_parsed}/Reactant.tab.gz | awk -F\"\t\" '\$2 == \"substrate\" {print \$3\"\t\"\$4}' |", 0, 0);
    $index{substrate_reaction} = $index{reaction_substrate}->reverse();

    #### reactions <-> products
    warn ("\tindexing reactions <-> products ...\n") if ($verbose >=1);
    $index{reaction_product} = classes::Index->new();
    $index{reaction_product}->load("gunzip -c $dir{kegg_parsed}/Reactant.tab.gz | awk -F\"\t\" '\$2 == \"product\" {print \$3\"\t\"\$4}' |", 0, 0);
    $index{product_reaction} = $index{reaction_product}->reverse();

    #### reactions <-> equations with compound IDs
    warn ("\tindexing reactions <-> equations ...\n") if ($verbose >=1);
    $index{reaction_equation} = classes::Index->new();
    $index{reaction_equation}->load("gunzip -c $dir{kegg_parsed}/Reaction.tab.gz | grep -v '^--' | cut -f 1,3 |", 0 ,1);
    $index{equation_reaction} = $index{reaction_equation}->reverse();

    #### reactions <-> equations with compound names
    warn ("\tindexing reactions <-> equations (names) ...\n") if ($verbose >=1);
    $index{reaction_equation_names} = classes::Index->new();
    $index{reaction_equation_names}->load("gunzip -c $dir{kegg_parsed}/Reaction.tab.gz | grep -v '^--' | cut -f 1,4 |", 0 ,1);
    $index{equation_reaction_names} = $index{reaction_equation}->reverse();

    #### amaze_id <-> kegg_id
#    warn ("\tindexing reaction amaze_id <-> kegg_id ...\n") if ($verbose >=1);
#    $index{amaze_kegg} = classes::Index->new();
#    $index{amaze_kegg}->load("$dir{amaze_export}/amaze_kegg_reactions.tab", 0 , 0);
#    $index{kegg_amaze} = $index{amaze_kegg}->reverse();

    #### patches on compound names
    warn ("\tindexing wrong_name <-> kegg_name ...\n") if ($verbose >=1);
    $index{wrong_kegg} = classes::Index->new();
    $index{wrong_kegg}->load($in_file{georges_identified_compounds}, 1, 1);

    

    #### KEGG reaction ID <-> aMAZE reaction ID
#    warn ("\tamaze_equation index ...\n") if ($verbose >=1);
#    $index{amaze_equation}->load($dir{amaze_export}."/reaction_equation.tab", 0, 1);


    #### report index sizes
    if ($verbose >=1) {
	warn ("Index sizes\n");
	foreach my $key (keys %index) { 
	    warn ("\t", $key, "\t",$index{$key}->get_size(),"\n");
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
#### read process files
sub ReadProcesses {
    #### columns in the input file
    my %col = (
	       step=>0,
	       ec=>1,
	       gene=>2,
	       substrates=>3,
	       products=>4,
	       equation=>5,
	       kegg_id => 6,
	       direction => 7
	       );

    foreach my $file (@files_to_parse) {
	my $process_name = `basename $file .txt`;
	$process_name =~ s/\.tab$//;
	chomp $process_name;
	$process_name =~ s/_/ /g;
	$process_name =~ s/biosyn$/biosynthesis/;
	unless ($process_name =~ /^.[A-Z]/) {
	    $process_name = ucfirst($process_name);
	}

	my $organism =  `dirname $file`;
	chomp $organism;
	$organism =~ s/_/ /g;
	$organism =~ s|/||g;

	#### initialize a new process
	$process = $processes->new_object();
	$process->set_attribute("source","amaze");
	$process->set_attribute("complete",1);
	$process->set_attribute("organism",$organism);
	$process->set_attribute("description","$organism - $process_name");
	$process->set_attribute("parsed_file",$file);
	$process->set_attribute("unidentified_compounds", 0) ;
	$process->set_attribute("matching_reactions", 0) ;
	$process->set_attribute("multi_match_reactions",0);
	$process->set_attribute("missing_direction",0);
	$process->set_attribute("no_match_reactions",0);
	$process->set_attribute("duplicate_reactions",0);
	$process->new_attribute_value("names",$process_name);
	if ($organism =~ /coli/) {
	    $process->set_attribute("annotator", "Gaurab");
	} elsif ($organism =~ /cerevisiae/) {
	    $process->set_attribute("annotator", "Georges Cohen");
	}

	my $product_index = classes::Index->new();
	my $substrate_index = classes::Index->new();

	if ($mirror) {
	    #### mirror file to report the matching reactions for each pathway step
	    my $dir = join "/", $dir{mirror}, $organism;
	    $dir =~ s/\s/_/g;
	    unless (-d $dir) {
		unless (-d $dir) {
		    warn "; Creating mirror dir $dir\n"  if ($verbose >=1);
		    `mkdir -p $dir`;
		    die "Error: cannot create mirror directory $dir\n" unless (-d $dir);
		}
		
	    }

	    #### create a control file for viewing the mirror file as a form in emacs
	    &PrintControlFile($dir, $process_name) ;

	    #### open mirror file
	    my $mirror_file = "$dir/${process_name}.tab";
	    $mirror_file =~ s/ /_/g;
	    warn "; Mirror file\t$mirror_file\n" if ($verbose >=1);
	    open MIRROR, ">$mirror_file"
		|| die "Error: cannot open mirror file $dir{skeletons}/$mirror_file\n";
	    print MIRROR join ("\t", "; step", "EC", "gene", "substrate", "product", "reaction", "reaction_id", "direction", "status", "remark"), "\n";

	}

	#### opening the data file
        warn "; Parsing file\t$dir{skeletons}/$file\n" if ($verbose >=1);
        die "Error: file $dir{skeletons}/$file does not exist\n"
	    unless (-e "$dir{skeletons}/$file");
        open FILE, "$dir{skeletons}/$file"
	    || die "Error: cannot read file $dir{skeletons}/$file\n";

#	my $header = <FILE>;
	my $l = 0;
	my $step = 0;
	while (my $line = <FILE>) {
	    $l++;
	    chomp $line;
            $line =~ s/ +\t/\t/g;
            $line =~ s/\t +/\t/g;
	    next unless ($line =~ /\S/);
	    next if ($line =~ /^;/);
	    next if (($l==1) & ($line =~ /^reaction/i));
	    $line =~ s/\r//;
	    my @fields = split "\t", $line;
	    warn "line $l\t$line\n" if ($verbose >=3);
	    next if ($line =~ /^;/);
	    #### step number
#	    $step = $fields[$col{step}];
	    $step++;

	    my @matching_reactions = ();
	    my %is_direct = ();
#	    my $provided_amaze_id = "";
	    my $provided_equation = "";
	    my $reverse_equation = "";
	    my $provided_direction = &trim($fields[$col{direction}]);

	    ################################################################
	    #### fix a problem with previoous annotation, where AmazeIDs were used to specify the reaciotn
	    my $provided_reaction = &trim($fields[$col{equation}]);
	    if ($provided_reaction =~ /amazereaction/i) {
		$provided_reaction = "";
		$fields[$col{equation}] = "";
	    }



	    ################################################################
	    #### Identify the reaction

	    #### reaction provided in the form of a KEGG ID
	    my $provided_kegg_id = &trim($fields[$col{kegg_id}]); 
	    
	    #### fix a problem with some remarks which are in this column where reaction IDs are supposed to be found
	    unless ($provided_kegg_id =~ /^A{0,1}R0/) {
		$provided_kegg_id = "";
		$fields[$col{kegg_id}] = "";
		warn ("provided_kegg_id\tin column ", 
		      $col{kegg_id} + 1, 
		      "\t", $provided_kegg_id, "\n") if ($verbose >=2);
	    }


	    #### KEGG id sometimes provided in the equation column
	    if (!($provided_kegg_id) && ($provided_reaction =~ /^A{0,1}R0/)) {
		$provided_kegg_id = $provided_reaction;
		warn ("provided_kegg_id\tin column ", 
		      $col{equation} + 1, 
		      "\t", $provided_kegg_id, "\n") if ($verbose >=2);
	    }

	    if ($provided_kegg_id) {
		warn ("provided_kegg_id", 
		      "\t", $provided_kegg_id, "\n") if ($verbose >= 2);
		@provided_reactions = split '\|', $provided_kegg_id;
		if ($#provided_reactions > 0) {
		    warn join ("\t", $file, $l, "Multiple reactions are provided", $provided_kegg_id), "\n" 
			if ($verbose >=2);
		}

		foreach my $reaction_id (@provided_reactions) {
		    if ($reaction_id=~ /^A{0,1}R0/) {
#			warn "HELLO\tprovided reaction\t$reaction_id\n" if ($verbose >= 0);
			push @matching_reactions, $reaction_id;
		    } else {
			&ErrorMessage (join ("\t", $file, $l, "Invalid KEGG ID for a reaction", $provided_kegg_id), "\n");
		    }		
		}

		#### reaction described by its equation
	    } elsif (($provided_reaction =~ /=/) || 
		     ($provided_reaction =~ /->/) || 
		     ($provided_reaction =~ /<-/)) {
		($provided_equation, $reverse_equation) = &CleanAndReverseEquation($provided_reaction);
		warn join ("\t", $file, $l, "provided equation", $provided_equation), "\n" 
		    if ($verbose >=2);

	    } 
	    #### initialize the reaction direction
	    $is_direct{$reaction} = $null;
	    

	    ################################################################
	    #### Parse  and check EC number
	    my $ec = &trim($fields[$col{ec}]);
	    if ($ec) {
		warn join ("\t", $file, $l, "EC number", $ec), "\n" 
		    if ($verbose >=3);
		if ($ec =~ /^[\d\-\.]+$/) {
		    $process->push_expanded_attribute("ECs", $ec, $step); 
		} else {
		    &ErrorMessage(join ("\t", 
					$file, 
					$l, 
					"invalid_EC_number", 
					$ec), 
				  "\n");
		}
	    }

	    ################################################################
	    #### parse gene names
	    my $gene_string = &trim($fields[$col{gene}]);
	    if ($gene_string) {
		my @genes = split '\|', $gene_string;
		foreach $gene (@genes) {
		    warn join ("\t", $file, $l, "gene", $gene), "\n" 
			if ($verbose >=3);
		    $process->push_expanded_attribute("genes", $gene, $step);
		}
	    }

	    ################################################################
	    #### main substrates
	    my $substrate_string = &unquote($fields[$col{substrates}]);
	    my @substrates = ();
	    my @substrate_names = split '\|', $substrate_string;
	    for my $i (0..$#substrate_names) {
		my $substrate_name = &standardize($substrate_names[$i]);
		if ($substrate_name eq "") {
		    &ErrorMessage(join ("\t", $file, $l, "empty_substrate", $substrate_string), "\n");
		} elsif ($index{name_compound}->contains($substrate_name)) {
		    push @substrates, $index{name_compound}->get_first_value($substrate_name);
#		    $process->new_attribute_value("substrates", $substrate_name); #if ($debug);
		    $process->push_expanded_attribute("substrates", $substrate_name, $step); #if ($debug);
		} elsif ($index{wrong_kegg}->contains($substrate_name)) {
		    my $substrate_name = $index{wrong_kegg}->get_first_value($substrate_name);
		    push @substrates, $index{name_compound}->get_first_value($substrate_name);
#		    $process->new_attribute_value("substrates", $substrate_name); #if ($debug);
		    $process->push_expanded_attribute("substrates", $substrate_name, $step); #if ($debug);
		} else {
		    $process->force_attribute("unidentified_compounds", 
					      $process->get_attribute("unidentified_compounds") + 1) ; #if ($debug);
		    &ErrorMessage(join ("\t", $file, $l, "unidentified_compound", $substrate_name), "\n");
		}
	    }
	    if ($verbose >=3) {
		warn "substrate_names\t", join ("|", @substrate_names), "\n" ;
		warn "substrates\t", join ("|", @substrates), "\n" ;
	    }
	    
	    ################################################################
	    #### main products
	    my $product_string = &unquote($fields[$col{products}]);
	    my @products = ();
	    my @product_names = split '\|', $product_string;
	    for my $i (0..$#product_names) {
		my $product_name = &standardize($product_names[$i]);
		if ($product_name eq "") {
		    &ErrorMessage(join ("\t", $file, $l, "empty_product", $product_string), "\n");
		} elsif ($index{name_compound}->contains($product_name)) {
		    push @products, $index{name_compound}->get_first_value($product_name);
		    $process->push_expanded_attribute("products", $product_name, $step); 
		} elsif ($index{wrong_kegg}->contains($product_name)) {
		    my $product_name = $index{wrong_kegg}->get_first_value($product_name);
		    push @products, $index{name_compound}->get_first_value($product_name);
		    $process->push_expanded_attribute("products", $product_name, $step); 
		} else {
		    $process->force_attribute("unidentified_compounds", 
					      $process->get_attribute("unidentified_compounds") + 1) ;
		    &ErrorMessage(join("\t", $file, $l, "unidentified_compound", $product_name), "\n");
		}
	    }
	    if ($verbose >=3) {
		warn "product_names\t", join ("|", @product_names), "\n" ;
		warn "products\t", join ("|", @products), "\n" ;
	    }


	    ################################################################
	    #### match forward reactions
	    my @matching_substrates_fwd = ();
	    if ($substrates[0]) {
		@matching_substrates_fwd = $index{substrate_reaction}->get_values($substrates[0]);
		foreach my $i (1..$#substrates) {
		    my @matching_next_substrate = $index{substrate_reaction}->get_values($substrates[$i]);
		    @matching_substrates_fwd = &intersection(\@matching_substrates_fwd, \@matching_next_substrate);
		}
		warn "matching_substrates_fwd\t", join ("\|", @matching_substrates_fwd), "\n" if ($verbose >=3); 
	    }
	    my @matching_products_fwd = ();
	    if ($products[0]) {
		@matching_products_fwd = $index{product_reaction}->get_values($products[0]);
		foreach my $i (1..$#products) {
		    my @matching_next_product = $index{product_reaction}->get_values($products[$i]);
		    @matching_products_fwd = &intersection(\@matching_products_fwd, \@matching_next_product);
		}
		warn "matching_products_fwd\t", join ("\|", @matching_products_fwd), "\n" if ($verbose >=3); 
	    }
	    
	    my @matching_reactions_fwd = &intersection(\@matching_substrates_fwd, \@matching_products_fwd);
	    foreach my $reaction (@matching_substrates_fwd, @matching_products_fwd) {
		$is_direct{$reaction} = $forward;
	    }
	    warn "matching_reactions_fwd\t", join ("\|", @matching_reactions_fwd), "\n" if ($verbose >=3); 
	    
	    ################################################################
	    #### match reverse reactions
	    my @matching_substrates_rev = ();
	    if ($substrates[0]) {
		@matching_substrates_rev = $index{product_reaction}->get_values($substrates[0]);
		foreach my $i (1..$#substrates) {
		    my @matching_next_substrate = $index{product_reaction}->get_values($substrates[$i]);
		    @matching_substrates_rev = &intersection(\@matching_substrates_rev, \@matching_next_substrate);
		}
		warn "matching_substrates_rev\t", join ("\|", @matching_substrates_rev), "\n" if ($verbose >=3); 
	    }
	    my @matching_products_rev = ();
	    if ($products[0]) {
		@matching_products_rev = $index{substrate_reaction}->get_values($products[0]);
		foreach my $i (1..$#products) {
		    my @matching_next_product = $index{substrate_reaction}->get_values($products[$i]);
		    @matching_products_rev = &intersection(\@matching_products_rev, \@matching_next_product);
		}
		warn "matching_products_rev\t", join ("\|", @matching_products_rev), "\n" if ($verbose >=3); 
	    }
	    
	    my @matching_reactions_rev = &intersection(\@matching_substrates_rev, \@matching_products_rev);
	    foreach my $reaction (@matching_substrates_rev, @matching_products_rev) {
		$is_direct{$reaction} = $reverse;
	    }
	    warn "matching_reactions_rev\t", join ("\|", @matching_reactions_rev), "\n" if ($verbose >=3); 
	    
	    unless ($provided_kegg_id) {
		#### matching reactions in both directions
		push @matching_reactions, @matching_reactions_fwd;
		push @matching_reactions, @matching_reactions_rev;

		#### restrict the list by imposing constraint on  EC number
		if (($#matching_reactions > 0) && ($ec)) {
		    my @matching_ecs = $index{ec_reaction}->get_values($ec);
		    warn "matching_ecs\t", join ("\|", @matching_ecs), "\n" if ($verbose >=3); 
		    @matching_reactions = &intersection(\@matching_reactions, \@matching_ecs);
		}
		
	    }

	    #### in case of doubt, use provided equation
	    if ($#matching_reactions != 0) {
		warn join ("\t", "Trying other columns", 
			   $provided_kegg_id || $provided_equation ), "\n" if ($verbose >=3);
		
		#### reaction equation provided in the reaction column
		if ($provided_equation =~ /\S/) {
		    warn join ("\t", "provided_equation", $provided_equation), "\n" if ($verbose >=3);
		    warn join ("\t", "reverse_equation", $reverse_equation), "\n" if ($verbose >=3);
		    if ($index{equation_reaction}->contains($provided_equation)) {
			my $kegg_reaction_id = $index{equation_reaction}->get_first_value($provided_equation);
			@matching_reactions = ($kegg_reaction_id);
			warn join ("\t", $file, $l, "identified_equation", $kegg_reaction_id, $provided_equation), "\n" if ($verbose >=3);
		    } elsif ($index{equation_reaction}->contains($reverse_equation)) {
			my $kegg_reaction_id = $index{equation_reaction}->get_first_value($reverse_equation);
			$is_direct{$kegg_reaction_id} = $reverse;
			@matching_reactions = ($kegg_reaction_id);
			warn join ("\t", $file, $l, "identified_reverse_equation", $kegg_reaction_id, $reverse_equation), "\n" if ($verbose >=3);
		    } else {
			&ErrorMessage(join ("\t", $file, $l, "unidentified_equation", $provided_equation), "\n");
			warn join ("\t", $file, $l, "unidentified_equation", $provided_equation), "\n" if ($verbose >=3);
		    }
		}

	    }

	    ################################################################
	    #### determine the direction of the reaction (forward or reverse) in each step
	    foreach my $reaction (@matching_reactions) {
		#### use provided direction
		if (($provided_direction =~ /\S/)  &&
		    (($provided_direction eq "0") || ($provided_direction eq "1"))) {
		    $is_direct{$reaction} = $provided_direction;
		}

		#### if nothing is provided, try to deduce from substrates and products
#		unless (defined($is_direct{$reaction})) {
		if ($is_direct{$reaction} eq $null) {
		    my $kegg_id = "";
		    
		    ################################################################
		    #### try to match any kegg substrate or product with the provided substrates/products
		    if ($reaction =~ /^A{0,1}R0/) {
			
			my @kegg_substrates = $index{reaction_substrate}->get_values($kegg_id);
			my @kegg_products = $index{reaction_product}->get_values($kegg_id);
			
			warn join ("\t", "kegg_reactants", join ("|", @kegg_substrates), join ("|", @kegg_products)), "\n" if ($verbose >=3);
			warn join ("\t", 
				   "provided_reactants", 
				   join ("|", @substrates), 
				   join ("|", @substrate_names), 
				   join ("|", @products),
				   join ("|", @product_names)), "\n" if ($verbose >=3);
			
			my @fwd_matches = ();
			push @fwd_matches, &intersection(\@kegg_substrates, \@substrates);
			push @fwd_matches, &intersection(\@kegg_products, \@products);
			
			my @rev_matches = ();
			push @rev_matches, &intersection(\@kegg_substrates, \@products);
			push @rev_matches, &intersection(\@kegg_products, \@substrates);
			
			if ($#fwd_matches >=0) {
			    $is_direct{$reaction} = $forward;
			    warn join ("\t", "Identified direction for reaction", $reaction, $kegg_id, $is_direct{$reaction}), "\n";
			    next;
			} elsif ($#rev_matches >=0) {
			    $is_direct{$reaction} = $reverse;
			    warn join ("\t", "Identified direction for reaction", $reaction, $kegg_id, $is_direct{$reaction}), "\n";
			    next;
			    
			} else {
			    $is_direct{$reaction} = $null;
			    $process->force_attribute("complete", 0);
			    warn join ("\t", "Cannot identify direction for reaction", $reaction, $kegg_id, $is_direct{$reaction}), "\n" if ($verbose >=3);
			    warn join ("\t", 
				       "substr",
				       join ("|", @substrates),
				       "prod",
				       join ("|", @products),
				       "equat",
				       $index{reaction_equation}->get_first_value($kegg_id)
				       ), "\n" if ($verbose >=3);
			    &ErrorMessage (join ("\t", $file, $l, "unidentified_direction", $reaction, $kegg_id), "\n");
			}
		    }
		    
		    $process->force_attribute("missing_direction", $process->get_attribute("missing_direction") + 1) ; #if ($debug);
		    $is_direct{$reaction} = $null;
		    $process->force_attribute("complete", 0);
		    
		    warn join ("\t", $reaction, $is_direct{$reaction}), "\n" if ($verbose >=3);
		    &ErrorMessage(join ( "\t", $file, $l, "arbitrarily_defining_reaction_direction", $reaction), "\n");
		}
	    }

	    ################################################################
	    #### report matching reactions
	    if ($verbose >=3) {
		my $msg = join ("\t", 
				"matching_reactions", 
				$#matching_reactions + 1, 
				join ("\|", @matching_reactions));
		if ($provided_amaz_id || $provided_equation) {
		    $msg .= "\t".join ("\t",, 
				       "provided", 
				       $provided_amaz_id || $provided_equation);
		}
		warn $msg, "\n";
	    } 


	    $mirror_line = join ("\t", @fields[0..5]);

	    ################################################################
	    #### report multiple matches
	    if ($#matching_reactions > 0) {
		
		#### several matching reactions
		$process->force_attribute("multi_match_reactions", $process->get_attribute("multi_match_reactions") + 1) ; #if ($debug);
		$process->force_attribute("complete", 0);
		print MIRROR join ("\t", $mirror_line, "", "",  "MULTI", join ("\|", @matching_reactions)), "\n" if ($mirror);
		
		
		#### report error
		&ErrorMessage (join ("\t",
				     "several matching reactions", 
				     $file,
				     $line,
				     join("\|", @matching_reactions ), 
				     join("\|", @equations ), 
				     ),
			       "\n");
		my $r = 1;
		&ErrorMessage ($file, 
			       "\t", $l, 
			       "\t", "multiple_matches", 
			       "\t", "0",
			       "\tdata\t", $line, 
			       "\n");
		foreach my $reaction (@matching_reactions)  {
		    
		    &ErrorMessage (join ("\t", 
					 $file, 
					 $l, 
					 "multiple_matches", 
					 $r++,
					 $reaction,
					 $index{reaction_equation_names}->get_values($reaction)
					 ), "\n");
		}
	    } elsif ($#matching_reactions == -1) {
		
		#### no reaction matchs the line
		$process->force_attribute("no_match_reactions", $process->get_attribute("no_match_reactions") + 1) ; #if ($debug);
		$process->force_attribute("complete", 0);
		print MIRROR join ("\t", $mirror_line, "", "",  "NO MATCHING REACTION"), "\n" if ($mirror);
		
		#### report error
		&ErrorMessage (join ("\t",
				     $file,
				     $l,
				     "no_matching_reaction",
				     $line
				     ),
			       "\n");
	    } else {
		#### If the reaction has been identified in an
		#### unequivocal way, create a process step

		$process->force_attribute("matching_reactions", 
					  $process->get_attribute("matching_reactions") + 1) ; #if ($debug);
		
		my $reaction = $matching_reactions[$r];
		print MIRROR join ("\t", $mirror_line, $reaction, $is_direct{$reaction}, "OK", 
				   $index{reaction_equation_names}->get_values($reaction)), "\n" if ($mirror);
		my $reaction = $matching_reactions[0];
		my $kegg_id = "ERROR";
		if ($reaction =~ /^A{0,1}R0/) {
		    $kegg_id = $reaction;
		} else {
		    &ErrorMessage(join ("\t", $file, $l, "invalid_reaction_identifier", $reaction), "\n");		
		}
		
		#### check if the direction is known
		if ($is_direct{$reaction} eq $null) {
		    &ErrorMessage(join ("\t", $file, $l, "cannot identify the direction for reaction ", $reaction), "\n");		
#		    next;
		}
		
		$leaf = $leaves->new_object();
		$leaf->set_attribute("interaction", $kegg_id);
		$leaf->set_attribute("is_direct", $is_direct{$reaction});
		if ($index{reaction_equation_names}->contains($kegg_id)) {
		    $leaf->set_attribute("description", $index{reaction_equation_names}->get_first_value($kegg_id));
		} else {
		    $leaf->set_attribute("description", $kegg_id);
		}
		$leaf->set_attribute("step", $step);
		$process->push_expanded_attribute("nodes",$leaf->get_id(), $step);
		$process->push_expanded_attribute("reactions", $reaction, $step); #if ($debug);
		warn ("OK", 
		      "\t",$reaction, 
		      "\t", $index{reaction_equation}->get_values($reaction), 
		      "\n")
		    if ($verbose >= 2);
		
		#### index product names for connection with next reactions
		foreach my $p (@products) {
		    $product_index->add_value($p, $leaf->get_attribute("id"));
		}
		#### index substrate names for connection with next reactions
		foreach my $s (@substrates) {
		    $substrate_index->add_value($s, $leaf->get_attribute("id"));
		}
	    }
	}


	################################################################
	#### create process arcs
	@indexed_products = $product_index->get_keys();
	foreach my $compound (@indexed_products) {
	    warn "Connecting reactions with compound $compound.\n"
		if ($verbose >= 2);
	    my @p_steps = $product_index->get_values($compound);
	    my @s_steps = $substrate_index->get_values($compound);
	    foreach $p (@p_steps) {
		foreach $s (@s_steps) {
		    warn join ("\t", 
			       "arc from step", $p,
			       "to step", $s, 
			       "with compound", $compound, 
			       $index{compound_name}->get_first_value($compound)), "\n"
				   if ($verbose >=3);
		    $process->push_expanded_attribute("arcs", $p, $s);
		}
	    }
	}

	################################################################
	#### check that there is no duplicate reaction
	my @reactions = sort($process->get_attribute("reactions"));
	warn ("; Reactions\n",
	      ";\t", 
	      join ("\n;\t", @reactions), 
	      "\n") if ($verbose >= 3);
	for my $i (0..$#reactions -1) {
	    if ($reactions[$i] eq $reactions[$i+1]) {
		&ErrorMessage (join ("\t", $file, "Duplicate reaction", $reactions[$i]), "\n");
		$process->force_attribute("duplicate_reactions", $process->get_attribute("duplicate_reactions") + 1) ; 
		$process->force_attribute("complete", 0);
	    }
	}

	################################################################
	#### report the number of steps in the process
	$process->set_attribute("steps", $step) ;
	
	################################################################
	#### count the number of problematic reactions
	$process->set_attribute("problems", 
				$process->get_attribute("multi_match_reactions") 
				+ $process->get_attribute("no_match_reactions")
				+ $process->get_attribute("duplicate_reactions")
				);
	if ($process->get_attribute("complete") == 1) {
	    $process->set_attribute("status", "OK");
	} else {
	    $process->set_attribute("status", "ERR");
	}

	close FILE;
	close MIRROR if ($mirror);
    }   
}



################################################################
#### print a control file for displaying each mirror file as a form
#### with emacs
sub PrintControlFile {
    my ($dir, $process_name) = @_;
    my $control_file = "$dir/${process_name}.fctl";
    $control_file =~ s/\s/_/g;
    my $mirror_file_short = $process_name.".tab";
    $mirror_file_short =~ s/\s+/_/g;
    open FCTL, ">$control_file";
    print FCTL<<EndForm;
(setq forms-file \"$mirror_file_short\")
(setq forms-number-of-fields 17)
;; (setq forms-read-only t)                 ; to make sure
(setq forms-field-sep \"\t\")
;;(setq forms-multi-line nil)               ; Don\'t allow multi-line fields.
(setq forms-read-only nil)
(setq forms-format-list (list
        \"\# 1-step        \# \" 1
        \"\\n\# 2-EC          \# \" 2
        \"\\n\# 3-gene        \# \" 3
        \"\\n\# 4-substrate   \# \" 4
        \"\\n\# 5-product     \# \" 5
        \"\\n\# 6-reaction    \# \" 6
        \"\\n\# 7-reaction id \# \" 7
        \"\\n\# 8-direction   \# \" 8
        \"\\n\# 9-status      \# \" 9
        \"\\n\# 10-remark     \# \" 10
        \"\\n\# 11            \# \" 11
        \"\\n\# 12            \# \" 12
        \"\\n\# 13            \# \" 13
        \"\\n\# 14            \# \" 14
        \"\\n\# 15            \# \" 15
        \"\\n\# 16            \# \" 16
        \"\\n\# 17            \# \" 17))
EndForm
    close FCTL;
}


################################################################
#### clean and reverse provided equation for matching it against KEGG index
sub CleanAndReverseEquation {
    my ($provided_reaction) = @_;

    my $direct_equation = $provided_reaction;
    my $reverse_equation = "";

    #### clean the provided equation
    $direct_equation =~ s/\=/\<\=\>/; #### fix cases where the < and > were forgotten
    $direct_equation =~ s/\<\<\=\>\>/\<\=\>/; #### fix cases where the < and > were forgotten
    $direct_equation =~ s/\<\=\>/ \<\=\> /;
    $direct_equation =~ s/\-\>/ \-\> /;
    $direct_equation =~ s/\<\-/ \<\- /;
    $direct_equation =~ s/ +/ /g;
    my $provided_equation = &standardize($direct_equation);
    warn ("provided_equation\t", $provided_equation, "\n") if ($verbose >=3);

    #### calculate the reverse reaction
    my $separator = "";
    my @reaction_sides = ();
    if (($direct_equation =~ /(\<\=\>)/) || 
	($direct_equation =~ /(\-\>)/) || 
	($direct_equation =~ /(\<\-)/)) {
	$separator = $1;
	@reaction_sides = split $separator, $direct_equation;
    }
    warn join ("\t", "Reaction separator", $separator), "\n" if ($verbose >=4);
    warn join ("\t", "Reaction sides", @reaction_sides), "\n" if ($verbose >=4);
    if ($#reaction_sides == 1) {
	$reverse_equation = &trim($reaction_sides[1]);
	$reverse_equation .= " ".$separator." ";
	$reverse_equation .= &trim($reaction_sides[0]);
	$reverse_equation = &standardize($reverse_equation);
	warn ("reverse_equation\t", $reverse_equation, "\n") if ($verbose >=3);
    } else {
	&ErrorMessage (join ("\t", $file, $l, "Cannot reverse equation", $provided_equation), "\n");
    }

    return ($provided_equation, $reverse_equation);
}

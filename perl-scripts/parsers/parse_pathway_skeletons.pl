#!/usr/bin/perl
############################################################
#
# $Id: parse_pathway_skeletons.pl,v 1.6 2002/10/25 17:56:05 jvanheld Exp $
#
# Time-stamp: <2002-10-25 11:54:57 jvanheld>
#
############################################################

require "PFBP_config.pl";
require "PFBP_classes.pl";
require "PFBP_util.pl";
require "PFBP_parsing_util.pl";
require "PFBP_loading_util.pl";

package main;
{
    #### default parameters

    $schema="annotator";
    $user="annotator";
    $password="annotator";

    $fast_index = 1;
    $debug = 0;
    $clean = 0;
    $mirror = 0;
    
    $reverse = 0;
    $forward = 1;
    
#    $dir{amaze_export} = "/win/amaze/amaze_data/exported_2001_0719";
    $dir{skeletons} = "/win/amaze/amaze_data/pathway_skeletons/";
#    $dir{georges_export} = "/win/amaze/amaze_team/for_georges";
#    $dir{skeletons} = "/win/amaze/amaze_team/gaurab/excel_tables/Pathway_ECno/Amino_acids/Amino_acid_biosyn";
#    $dir{skeletons} = "/win/amaze/amaze_team/gaurab/2001_07/Pathway_ECno/";
#    $dir{skeletons} = "/win/amaze/amaze_team/gaurab/excel_tables/Pathway_ECno/";
    $dir{kegg_parsed} = "$parsed_data/kegg_ligand/20020709";

    $in_file{georges_identified_compounds} = "/win/amaze/amaze_team/georges_cohen/excel_files_georges/unidentified_compounds/identified_compounds_2001_1007.tab";


    my $files_to_parse = `find $dir{skeletons} -name "*.txt" | xargs`;
    $files_to_parse =~ s/$dir{skeletons}//g;
    @files_to_parse = split " ", $files_to_parse;

    $dir{output} = "$parsed_data/pathway_skeletons/$delivery_date";

    $dir{seed_files} = "$dir{output}/amaze_pathway_seeds";
#    $dir{seed_files} = "/win/amaze/amaze_programs/amaze_graph_analysis/amaze_pathway_seeds";

    $out_file{process_skeleton} = "$dir{output}/process_skeleton.obj";
    $out_file{stats} = "$dir{output}/process_skeleton.stats.txt";
    $out_file{errors} = "$dir{output}/process_skeleton.errors.txt";
    
    $out_format = "obj";

    #### class factory
    $processes = PFBP::ClassFactory->new_class(object_type=>"PFBP::Process",
					       prefix=>"proc_");
    $leaves = PFBP::ClassFactory->new_class(object_type=>"PFBP::ProcessLeaf",
						    prefix=>"leaf_");


    push @classes, ("PFBP::Process");
    push @classes, ("PFBP::ProcessLeaf");

    &ReadArguments;

    #### output directories
    &CheckOutputDir();
    if ($mirror) {
	$dir{mirror} = $dir{output}."/mirror";
	unless (-d $dir{mirror}) {
	    warn "Creating mirror dir $dir{mirror}\n"  if ($verbose >=1);
	    `mkdir -p $dir{mirror}`;
	    die "Error: cannot create mirror directory $dir{mirror}\n" unless (-d $dir{mirror});
	}
    }

    #### error stream
    open ERR, ">$out_file{errors}" || die "Error: cannot write error report fle $$out_file{errors}\n";
    
    #### test 
    if ($test) {
#	@files_to_parse = qw ( Saccharomyces_cerevisiae/Isoleucine_and_valine_biosynth.txt );
	@files_to_parse = qw ( Escherichia_coli/TCA_glyoxylate_bypass.txt );
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
				));
    
    ### default verbose message
    if ($verbose >=1) {
	&DefaultVerbose;
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
    &ExportClasses($out_file{process_skeleton}, $out_format, @classes) if ($export{obj});
    &PrintStats($out_file{stats}, @classes);

    &GenerateECSeeds() if ($seeds);
    close ERR;


    ### report execution time
    if ($verbose >=1) {
	$done_time = &AlphaDate;
	warn (";\n",
	      "; job started $start_time",
	      "; job done    $done_time\n")
	}

    #### compress result files
    system "gzip -f $dir{output}/*.tab $dir{output}/*.txt";
    system "gzip -f $dir{output}/*.obj" if ($export{obj});

    exit(0);
}

################################################################
#### export data in a special format for pathway reconstruction
################################################################
sub GenerateECSeeds {
    print STDERR "; Generating EC seeds for pathway reconstruction\n"
	if ($verbose >=1);
    unless (-d $dir{seed_files}) {
	warn "Creating seed_files dir $dir{seed_files}\n"  if ($verbose >=1);
	`mkdir -p  $dir{seed_files}`;
	die "Error: cannot create directory $dir{seed_files}\n" unless (-d $dir{seed_files});
    }
    foreach $process ($processes->get_objects()) {
	my $file = $process->get_attribute("parsed_file");
	my $organism_prefix = $process->get_attribute("organism");
	$organism_prefix =~ s/\s/_/g;
	my $basename = `basename $file .txt`;
	chomp $basename;

	#### export EC number as seeds for pathway reconstruction
	my $export_file = $dir{seed_files};
	$export_file .= "/".$organism_prefix."_".$basename;
#	$export_file .= ".".$process->get_attribute("complete");
	open ECS, ">$export_file.ecs";
	foreach $ec ($process->get_attribute(ECs)) {
	    print ECS "$ec\n";
	    print ECS "${ec}*\n"; #### for acccepting the reverse reaction
	}
	close ECS;
	
	#### export reaction IDs for validation
	#### only export reaction for pathways where all 
	#### reactions are identified in an unequivocal way.
	next unless ($process->get_attribute("status") eq "OK");
	open REACTIONS, ">$export_file.reactions";
	foreach $reaction ($process->get_attribute(reactions)) {
	    print REACTIONS "$reaction\n";
	}
	close ECS;
    }
    
}

### print the help message 
### when the program is called with the -h or -help option
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
	-h	detailed help
	-help	short list of options
	-test	fast parsing of partial data, for debugging
	-v #	warn level
		Warn level 1 corresponds to a restricted verbose
		Warn level 2 reports all polypeptide instantiations
		Warn level 3 reports failing get_attribute()
	-obj	export data in object format (.obj file)
		which are human-readable
	-d	debug
		add attributes to the process to reflect the 
		parsed data, which ar not required for the actual 
		aMAZE process
	-mirror export a copy of the data files with additional
		columns to indicate whether each equation was
		identified or not.
	-clean remove all files from the outpu directory before
	       parsing.
	-skel  skeleton directory
	       Directory with the input files. See below for a
	       description of directory organization and file formats.

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
	1	reaction
	2	EC number
	3	gene
	4	from_comp
	5	to_comp
	6	Reaction

	(in George Cohen's files, there are some additional columns,
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

### read arguments from the command line
sub ReadArguments {
    for my $a (0..$#ARGV) {
	
	### warn level
	if (($ARGV[$a] eq "-v" ) && 
	    ($ARGV[$a+1] =~ /^\d+$/)){
	    $main::verbose = $ARGV[$a+1];
	    
	    ### test run
	} elsif ($ARGV[$a] eq "-test") {
	    $main::test = 1;
	    
	    ### output file
 	} elsif ($ARGV[$a] eq "-obj") {
	    $a++;
	    $main::export{obj} = 1;
	    
	    ### debug
	} elsif ($ARGV[$a] eq "-d") {
	    $main::debug = 1;

	    ### mirror
	} elsif ($ARGV[$a] eq "-mirror") {
	    $main::mirror = 1;

	    ### clean
	} elsif ($ARGV[$a] eq "-clean") {
	    $main::clean = 1;

	    ### skeleton directory
	} elsif ($ARGV[$a] eq "-skel") {
	    $main::dir{skeletons} = $ARGV[$a+1];

	    ### ec seeds
	} elsif ($ARGV[$a] eq "-seeds") {
	    $main::seeds = 1;

	    ### help
	} elsif (($ARGV[$a] eq "-h") ||
		 ($ARGV[$a] eq "-help")) {
	    &PrintHelp;
	    exit(0);
	}
    }
}

sub LoadIndexes {
    #### load indexes
    warn ("Loading index files\n") if ($verbose >=1);

    #### compound names
    warn ("\tindexing compound names ...\n") if ($verbose >=1);
    $index{compound_name} = PFBP::Index->new();
    $index{compound_name}->load("gunzip -c $dir{kegg_parsed}/Compound_names.tab.gz | grep -v '^--' |", 0, 1);
    $index{name_compound} = $index{compound_name}->reverse();
    
    #### reactions <-> EC numbers
    warn ("\tindexing reactions <-> EC numbers ...\n") if ($verbose >=1);
    $index{reaction_ec} = PFBP::Index->new();
    $index{reaction_ec}->load("gunzip -c $dir{kegg_parsed}/Reaction_ecs.tab.gz | grep -v '^--' |", 0, 0);
    $index{ec_reaction} = $index{reaction_ec}->reverse();
    
    #### reactions <-> substrates
    warn ("\tindexing reactions <-> substrates ...\n") if ($verbose >=1);
    $index{reaction_substrate} = PFBP::Index->new();
    $index{reaction_substrate}->load("gunzip -c $dir{kegg_parsed}/Reactant.tab.gz | awk -F\"\t\" '\$2 == \"substrate\" {print \$3\"\t\"\$4}' |", 0, 0);
    $index{substrate_reaction} = $index{reaction_substrate}->reverse();

    #### reactions <-> products
    warn ("\tindexing reactions <-> products ...\n") if ($verbose >=1);
    $index{reaction_product} = PFBP::Index->new();
    $index{reaction_product}->load("gunzip -c $dir{kegg_parsed}/Reactant.tab.gz | awk -F\"\t\" '\$2 == \"product\" {print \$3\"\t\"\$4}' |", 0, 0);
    $index{product_reaction} = $index{reaction_product}->reverse();

    #### reactions <-> equations
    warn ("\tindexing reactions <-> equations ...\n") if ($verbose >=1);
    $index{reaction_equation} = PFBP::Index->new();
    $index{reaction_equation}->load("gunzip -c $dir{kegg_parsed}/Reaction.tab.gz | grep -v '^--' | cut -f 1,3 |", 0 ,1);
    $index{equation_reaction} = $index{reaction_equation}->reverse();

    #### amaze_id <-> kegg_id
    warn ("\tindexing reaction amaze_id <-> kegg_id ...\n") if ($verbose >=1);
#    $index{amaze_kegg} = PFBP::Index->new();
#    $index{amaze_kegg}->load("$dir{amaze_export}/amaze_kegg_reactions.tab", 0 , 0);
#    $index{kegg_amaze} = $index{amaze_kegg}->reverse();

    #### patches on compound names
    warn ("\tindexing wrong_name <-> kegg_name ...\n") if ($verbose >=1);
    $index{wrong_kegg} = PFBP::Index->new();
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

#### read process files
sub ReadProcesses {
    #### columns in the input file
    my %col = (
	       step=>0,
	       ec=>1,
	       gene=>2,
	       substrates=>3,
	       products=>4,
	       amaze_id=>5
	       );
    $col{kegg_id} = 5;
    $col{equation} = 5;

    foreach my $file (@files_to_parse) {
	my $process_name = `basename $file .txt`;
	chomp $process_name;
	$process_name =~ s/_/ /g;
	$process_name =~ s/biosyn$/biosynthesis/;
	unless ($process_name =~ /^.[A-Z]/) {
	    $process_name = ucfirst($process_name);
	}

	my $organism =  `dirname $file`;
	chomp $organism;
	$organism =~ s/_/ /g;

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

	my $product_index = PFBP::Index->new();
	my $substrate_index = PFBP::Index->new();

	if ($mirror) {
	    #### mirror file to report the matching reactions for each pathway step
	    my $mirror_file = "$dir{mirror}/mirror_${organism}_${process_name}.tab";
	    $mirror_file =~ s/ /_/g;
	    warn "; Mirror file\t$mirror_file\n" if ($verbose >=1);
	    open MIRROR, ">$mirror_file"
		|| die "Error: cannot open mirror file $dir{skeletons}/$mirror_file\n";
	}

	#### opening the data file
	warn "; Parsing file\t$file\n" if ($verbose >=1);
	die "Error: file $dir{skeletons}/$file does not exist\n"
	    unless (-e "$dir{skeletons}/$file");
	open FILE, "$dir{skeletons}/$file"
	    || die "Error: cannot read file $dir{skeletons}/$file\n";

	my $header = <FILE>;
	my $l = 0;
	my $step = 0;
	while (my $line = <FILE>) {
	    $l++;
	    chomp $line;
	    next unless ($line =~ /\S/);
	    $line =~ s/\r//;
	    my @fields = split "\t", $line;
	    warn "line\t$line\n" if ($verbose >=3);

	    next if ($line =~ /^;/);
	    #### step number
#	    $step = $fields[$col{step}];
	    $step++;

	    my @matching_reactions = ();
	    my %is_direct = ();
	    my $provided_reaction = &trim($fields[$col{amaze_id}]);
	    my $provided_amaze_id = "";
	    my $provided_kegg_id = ""; 
	    my $provided_equation = "";
	    my $reverse_equation = "";
	    my $multi_allowed = 0;


		#### reaction provided in the form of an aMAZE ID
	    if ($provided_reaction =~ /aMAZEReaction/) {
		$provided_amaze_id = $provided_reaction;
		warn ("provided_amaze_id\t", $provided_amaze_id, "\n") if ($verbose >=3);

		#### reaction provided in the form of a KEGG ID
	    } elsif ($provided_reaction =~ /^R0/) {
		$provided_kegg_id = $provided_reaction;
		warn ("provided_kegg_id\t", $provided_kegg_id, "\n") if ($verbose >=3);

		#### reaction described by its equation
	    } elsif (($provided_reaction =~ /=/) || 
		     ($provided_reaction =~ /->/) || 
		     ($provided_reaction =~ /<-/)) {

		#### clean the provided equation
		$provided_reaction =~ s/\=/\<\=\>/; #### fix cases where the < and > were forgotten
		$provided_reaction =~ s/\<\<\=\>\>/\<\=\>/;#### fix cases where the < and > were forgotten
		$provided_reaction =~ s/\<\=\>/ \<\=\> /;
		$provided_reaction =~ s/\-\>/ \-\> /;
		$provided_reaction =~ s/\<\-/ \<\- /;
		$provided_reaction =~ s/ +/ /g;
		$provided_equation = &standardize($provided_reaction);
		warn ("provided_equation\t", $provided_equation, "\n") if ($verbose >=3);

		#### calculate the reverse reaction
		my $separator = "";
		my @reaction_sides = ();
		if (($provided_reaction =~ /(\<\=\>)/) || 
		    ($provided_reaction =~ /(\-\>)/) || 
		    ($provided_reaction =~ /(\<\-)/)) {
		    $separator = $1;
		    @reaction_sides = split $separator, $provided_reaction;
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
	    } 

	    #### EC number
	    my $ec = &trim($fields[$col{ec}]);
	    if ($ec) {
		if ($ec =~ /^[\d\-\.]+$/) {
#		    $process->new_attribute_value("ECs", $ec); #if ($debug);
		    $process->push_expanded_attribute("ECs", $ec, $step); #if ($debug);
		} else {
		    &ErrorMessage(join ("\t", 
					$file, 
					$l, 
					"invalid_EC_number", 
					$ec), 
				  "\n");
		}
#	    } else {
#		&ErrorMessage(join("\t", $file, $l, "no_EC_number"), "\n");
	    }
	    
	    #### parse gene names
	    #	    if ($debug) {
	    my $gene_string = &trim($fields[$col{gene}]);
	    if ($gene_string) {
		my @genes = split '\|', $gene_string;
		foreach $gene (@genes) {
		    $process->push_expanded_attribute("genes", $gene, $step);
		}
	    }
	    
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
#		    $process->new_attribute_value("products", $product_name); #if ($debug);
		    $process->push_expanded_attribute("products", $product_name, $step); #if ($debug);
		} elsif ($index{wrong_kegg}->contains($product_name)) {
		    my $product_name = $index{wrong_kegg}->get_first_value($product_name);
		    push @products, $index{name_compound}->get_first_value($product_name);
#		    $process->new_attribute_value("products", $product_name); #if ($debug);
		    $process->push_expanded_attribute("products", $product_name, $step); #if ($debug);
		} else {
		    $process->force_attribute("unidentified_compounds", 
					    $process->get_attribute("unidentified_compounds") + 1) ; #if ($debug);
		    &ErrorMessage(join("\t", $file, $l, "unidentified_compound", $product_name), "\n");
		}
	    }
	    if ($verbose >=3) {
		warn "product_names\t", join ("|", @product_names), "\n" ;
		warn "products\t", join ("|", @products), "\n" ;
	    }

#  	    unless ($substrates[0]) {
#  		&ErrorMessage("$file\t$l\t", "row without substrate\t", "$line\n"); 
#  		$process->force_attribute("unidentified_compounds", $process->get_attribute("unidentified_compounds") + 1) ; #if ($debug);
#  #		$process->force_attribute("complete", 0);
#  	    }

#  	    unless ($products[0]) {
#  		&ErrorMessage("$file\t$l\t", "row without product\t", "$line\n"); 
#  		$process->force_attribute("unidentified_compounds", $process->get_attribute("unidentified_compounds") + 1) ; #if ($debug);	
#  #		$process->force_attribute("complete", 0);
#  	    }

	    
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

	    #### matching reactions in both directions
	    push @matching_reactions, @matching_reactions_fwd;
	    push @matching_reactions, @matching_reactions_rev;

	    #### restrict the list by imposing constraint on  EC number
	    if (($#matching_reactions > 0) && ($ec)) {
		my @matching_ecs = $index{ec_reaction}->get_values($ec);
		warn "matching_ecs\t", join ("\|", @matching_ecs), "\n" if ($verbose >=3); 
		@matching_reactions = &intersection(\@matching_reactions, \@matching_ecs);
	    }
	    
	    

	    #### in case of doubt, use provided reaction column
	    if ($#matching_reactions != 0) {
		warn join ("\t", "Trying other columns", $provided_amaze_id || $provided_kegg_id || $provided_equation ), "\n" if ($verbose >=3);

		#### aMAZE ID provided directly in the file
		if ($provided_amaze_id =~ /aMAZEReaction/) {
		    @matching_reactions = split '\|', $provided_amaze_id;
		    $multi_allowed = 1;
		    warn join ("\t", $file, $l, "provided_amaze_id", $provided_amaze_id), "\n" if ($verbose >=3);
		    
		} elsif ($provided_kegg_id =~ /R0/) {
		    @matching_reactions = split '\|', $provided_kegg_id;
		    $multi_allowed = 1;

		    warn join ("\t", $file, $l, "provided_kegg_id", $provided_kegg_id), "\n" if ($verbose >=3);
		}

		#### reaction equation provided by Georges Cohen
		#warn join ("\t", "provided_equation", $provided_equation), "\n" if ($verbose >=3);
		#warn join ("\t", "reverse_equation", $reverse_equation), "\n" if ($verbose >=3);
		if ($provided_equation =~ /\S/) {
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

		#### checking reaction direction
		foreach my $reaction (@matching_reactions) {
		    unless (defined($is_direct{$reaction})) {
			my $kegg_id = "";

			#### try to identify the kegg reaction ID
			if ($reaction =~ /R0/) {
			    $kegg_id = $reaction;
#			} elsif ($reaction =~ /aMAZE/) {
#			    #### convert amaze_id into kegg_id
#			    $kegg_id = $index{amaze_kegg}->get_first_value($reaction);
#			    if (defined($is_direct{$kegg_id})) {
#				$is_direct{$reaction} = $is_direct{$kegg_id};
#				warn join ("\t", 
#					   $reaction, 
#					   $is_direct{$reaction}, 
#					   "identified the direction\t$reaction\t$kegg_id\n"), 
#				"\n" if ($verbose >=3);
#				next;
#'			    }
			}
			
			#### try to match any kegg substrate or product with the provided substrates/products
			if ($kegg_id) {
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
	    }

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


	    ################################################################
	    #### report multiple matches
	    if (!($multi_allowed) && ($#matching_reactions > 0)) {
		#### several matching reactions
		$process->force_attribute("multi_match_reactions", $process->get_attribute("multi_match_reactions") + 1) ; #if ($debug);
		$process->force_attribute("complete", 0);
		print MIRROR join ("\t", $line, "MULTI", join ("\|", @matching_reactions)), "\n" if ($mirror);

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
					 $index{reaction_equation}->get_values($reaction)
					 ), "\n");
		}
	    }

	    if ($#matching_reactions == -1) {
		#### no reaction matchs the line
		$process->force_attribute("no_match_reactions", $process->get_attribute("no_match_reactions") + 1) ; #if ($debug);
		$process->force_attribute("complete", 0);
		print MIRROR $line, "\t", "NO MATCHING REACTION\n" if ($mirror);
		

		#### report error
		&ErrorMessage (join ("\t",
				     $file,
				     $l,
				     "no_matching_reaction",
				     $line
				     ),
			       "\n");
	    } else {
		#### the reaction could b identified in an unequivocal way
		#### create a process step

		$process->force_attribute("matching_reactions", 
					$process->get_attribute("matching_reactions") + 1) ; #if ($debug);
		
		$max_reactions = 0; ### TEMPORARY: do not allow more than one reaction per step, to avoid layout problems. 

		foreach my $r (0..$max_reactions) {
		    my $reaction = $matching_reactions[$r];
		    print MIRROR join ("\t", $line, "OK", $reaction, 
				       $index{reaction_equation}->get_values($reaction)), "\n" if ($mirror);
		    my $reaction = $matching_reactions[0];
		    my $kegg_id = "ERROR";
		    if ($reaction =~ /aMAZEReaction/) {
#			if ( $index{amaze_kegg}->contains($reaction)) {
#			    $kegg_id = $index{amaze_kegg}->get_first_value($reaction);
#			} else {
			    &ErrorMessage(join ("\t", $file, $l, "cannot_map_amaze_id_to_kegg_id", $reaction), "\n");
#			}
		    } elsif ($reaction =~ /R0/) {
			$kegg_id = $reaction;
		    } else {
			&ErrorMessage(join ("\t", $file, $l, "invalid_reaction_identifier", $reaction), "\n");		
		    }
		    
		    #### check if the direction is known
		    if ($is_direct{$reaction} eq $null) {
			next;
		    }

		    $leaf = $leaves->new_object();
		    $leaf->set_attribute("interaction", $kegg_id);
		    $leaf->set_attribute("is_direct", $is_direct{$reaction});
		    $leaf->set_attribute("step", $step);
#		    $process->new_attribute_value("nodes",$leaf->get_id());
		    $process->push_expanded_attribute("nodes",$leaf->get_id(), $step);
#		    $process->new_attribute_value("reactions", $reaction); #if ($debug);
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
	}

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

	#### report the number of steps in the process
	$process->set_attribute("steps", $step) ;
	
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


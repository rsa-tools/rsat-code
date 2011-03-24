################################################################
## Library for functions shared by footprint detection programs
## - footprint-discovery
## - footprint-scan
## - get-orthologs

require "RSA2.cgi.lib";		## For sortable HTML tables

## Options for the &doit() command;
local $dry = 0;	  ## Do not run the command, just echo them as warning
local $batch = 0;		## Run the processes on a PC cluster
local $die_on_error = 1;
local $job_prefix = "footprint_disco";
local $cmd;


local $dist_thr = 55;	    ## Distance threshold for operon inference

%task = ();
$main::all_genes = 0;	## Analyze all the genes of the query organism
$main::create_index = 1; ## The index is now (from 2010/11/29) always created. Still to be checked in footprint-scan.

%main::infile = ();		## input files
%main::outfile = ();		## output files
%main::dir = ();		## output directories
$dir{output_root} = "footprints"; ## Default root output directory. Can be changed with the optiojn -o


################################################################
## The procedure &one_command() has been moved to RSAT::util
sub one_command {
  &RSAT::util::one_command(@_);
}

################################################################
## Read a list of organisms from a file, and select those supported on
## the RSAT server (issue a warning for non-supported organisms).
sub ReadOrganismsFromFile {
  ($main::orglist) = &OpenInputFile($main::orglist_file);
  my @organisms = ();
  while (<$main::orglist>) {
    chomp();
    s/\r/\n/g;		  ## Suppress Windows-specific carriage return
    next if /^;/;		## Comment line
    next if /^\#/;		## Header line
    next if /^\--/;		## SQL comment line
    next unless /\S/;		## Empty line
    my ($org_from_list) = split /\s/;
    $org_from_list = &trim($org_from_list); ## Remove leading and trailing spaces
    push @main::org_from_lists, $org_from_list;	## We need a list to report the results in the same order as the org file.
    $main::org_from_lists{$org_from_list} = 1; ## The hash table is convenient for checking which organism is part of the list
  }
  close $main::orglist if ($main::orglist_file);
  if (scalar(@org_from_lists) > 0) {
    &RSAT::message::Info(scalar(@org_from_lists), "organisms specified in file", $main::orglist_file) if ($main::verbose >= 1);
  } else {
    &RSAT::error::FatalError("The organism file does not contain any valid organism name.", $main::orglist_file);
  }

  ## Select only the supported organisms in the input list
  foreach my $org (@org_from_lists) {
    if ($supported_organism{$org}) {
      push @organisms, $org;
    } else {
      &RSAT::message::Warning("Organism", $org, "is not supported on this RSAT server");
    }
  }

  &RSAT::error::FatalError("There is no supported organism in the organism list", $main::orglist_file,
			   "\nUse the RSAT command supported-organisms to obtain the list of supported organisms."
			  )  if (scalar(@organisms) == 0);
  return(@organisms);
}

################################################################
## Select reference organisms, either from a taxon (if the option
## -taxon has been used) or from a user-specified list (option
## -org_list).
sub SelectReferenceOrganisms {

  my @ref_organisms = ();

  ## Check reference taxon or org_list
  &RSAT::error::FatalError("You should select a taxon of interest or provide a list of organisms")
    unless ($main::taxon || $main::orglist_file);

  &RSAT::error::FatalError("-taxon and -org_list option are mutually exclusive")
    if (($main::taxon) && ($main::orglist_file));

  ################################################################
  ## Read the list reference organisms from a user-specified file
  if ($main::orglist_file) {
    @ref_organisms = &ReadOrganismsFromFile($main::orglist_file);
  } else {
    #check if there is is at least one organism for that Taxon
    @ref_organisms = &CheckTaxon($taxon) if $main::taxon;

    ## Check option -no_purge and -org_list
    if ($main::no_purge) {
      &RSAT::error::FatalError("-no_purge option can only be used if -org_list option is active (not with -taxon).");
    }
  }

  return @ref_organisms;
}

################################################################
## Check all parameters required for footprint analysis
## (footprint-discovery or footprint-scan).
sub CheckFootprintParameters {

  ################################################################
  ## If all tasks are requested or if no task is defined, execute all
  ## tasks.
  if ((scalar(keys(%task)) == 0) || ($task{all})) {
    foreach my $task (keys %supported_task) {
      #      &RSAT::message::Debug("Auto adding task", $task);
      $task{$task} = 1;
    }
  }
  &RSAT::message::Info("Footprint analysis tasks: ", join (",", keys %task)) if ($main::verbose >= 2);

  ##############################################################
  ## If orthologs_list is specified omit the tasks that are not longer necesary
  ## Omit: orthologs
  $task{orthologs} = 0 if $main::orthologs_list_file ;

  ################################################################
  ## Check reference organisms (taxon or org list)
  my @ref_organisms = &SelectReferenceOrganisms() unless  $main::orthologs_list_file ;

  ################################################################
  ## Check query organism
  &RSAT::error::FatalError("You must specify a query organism (option -org)")
    unless ($organism_name);
  $organism = new RSAT::organism();
  $organism->check_name($organism_name);
  $organism->set_attribute("name", $organism_name);

  ################################################################
  ## Read query genes from input file
  if ($all_genes) {
    if (defined($supported_organism{$organism_name})) {
      $infile{genes} = $supported_organism{$organism_name}->{'data'};
      $infile{genes} .= "/genome/cds.tab";
    } else {
      &RSAT::error::FatalError("Organism", $organism_name, "is not supported");
    }
  }
  if ($infile{genes}) {
    my ($in) = &OpenInputFile($infile{genes});
    while (<$in>) {
      next if (/^--/);		## Skip mysql-type comment lines
      next if (/^;/);		## Skip comment lines
      next if (/^#/);		## Skip header lines
      next unless (/\S/);	## Skip empty lines
      chomp;
      my @fields = split /\s+/;
      my $query = shift @fields;
      push @query_genes, $query;
    }
    close $in;
  }
  if ( $main::orthologs_list_file) { 
      my $genes=`grep -v "#" $main::orthologs_list_file | grep -v ";" | cut -f1 | sort -u`;
      chomp($genes);
      my @fields = split (/\s+/,$genes);
      push @query_genes, @fields;
  }
  
  &RSAT::message::Info("Genes analyzed taken from orthologs_list file",join(" ", @query_genes)) if ( ($main::verbose >= 0) && ($main::orthologs_list_file) );
  
  ################################################################
  ## Check query genes
  if (scalar(@query_genes) ==0) {
    &RSAT::error::FatalError("You must specify at least one query gene (options -q , -i or -orthologs_list)");
  } else {
    &RSAT::message::Info("Number of query genes", scalar(@query_genes)) if ($main::verbose >= 2);
  }

  ################################################################
  ## Get maximum number of genes if limited
  if ($max_genes) {
    if (scalar(@query_genes)>$max_genes) {
      &RSAT::message::Warning("The analysis has been limited to the first ", $max_genes,"genes");
      @query_genes= splice(@query_genes,0,$max_genes);
      &RSAT::message::Warning(join("\t","Query genes:",@query_genes)) if ($main::verbose >= 2);
    }
  }

}

################################################################
## Check dependency between a task and a list of files
sub CheckDependency {
  my ($task, @files_types) = @_;

  ## in batch mode, the check should not be done since the previous
  ## tasks have not yet been ran
  return(0) if ($batch);

  foreach my $type (@files_types) {
    my $file = $outfile{$type};
    if ((-e $file) || (-e $file.".gz")) {
      &RSAT::message::Info("Checked existence of ", $type, "file required for task", $task, "file", $file)
	if ($main::verbose >= 3);
      return(1);
    } else {
      &RSAT::error::FatalError("Missing", $type, "file required for task", $task, "file", $file)
    }
  }
}

################################################################
## Define a query prefix that will be used as title for the feature
## map and as prefix for the automatic specificationof output files
## (will be used by &getOutfilePrefix()).
sub GetQueryPrefix {
  my $query_prefix;
  if (scalar(@current_query_genes) == 1) {
    $query_prefix = $current_query_genes[0];
  } elsif (scalar(@current_query_genes) <= 10) {
    $query_prefix = join "_", @current_query_genes;
    #  } elsif ($outfile{prefix}) {
    #    $query_prefix = `basename $outfile{prefix}`;
    #    chomp($query_prefix);
  } elsif ($infile{genes}) {
    $query_prefix = `basename $infile{genes} .tab`;
  } else {
    $query_prefix = "multiple_genes";
  }
  &RSAT::message::Info("&GetQueryPrefix()", join("","'", $query_prefix, "'")) if ($main::verbose >= 5);

  $dir{query_prefix} = $query_prefix; ## Index the query subdir, which will also serve for the main index
  return ($query_prefix);
}


################################################################
## Define the result sub-directories and the prefix for all output
## files for a given query (separated queries are dealt successively
## by this function).
##
## The query-specific prefix ($main::outfile{prefix}) is computed
## automatically on the basis of the taxon, query organism, background
## model and query gene(s).
##
## If several genes are analyzed separately (option -sep_genes), the
## result is saved in a separate directory for each gene.
##
## This routine creates the query-specific sub-directory if required.
sub GetOutfilePrefix {

  ## Create the query-specific sub-directory
  my $query_prefix = &GetQueryPrefix();
  $main::used_orgs="";
  if ($taxon){
      $used_orgs=$taxon;
  }
  elsif ($main::orglist_file){
      $used_orgs="org_list";
  }
  elsif ($main::orthologs_list_file){
      $used_orgs="orthologs_list";
  }
 
  $dir{output_per_query} = join("/",$main::dir{output_root}, $used_orgs, $organism_name, $query_prefix);
  &RSAT::util::CheckOutDir($dir{output_per_query});

  ## Compute a query-specific file prefix including the main parameters
  my $outfile_prefix = $query_prefix;
  $outfile_prefix .= "_";
  $outfile_prefix .= join ("_", $organism_name,$used_orgs );

  ## We don't want the bg model in the query prefix, because it is only a parameter for the dyads file (not for the sequences)
  #  if ($bg_model) {
  #    $outfile_prefix .= "_".$bg_model;
  #  }
  if ($infer_operons) {
    $outfile_prefix .= "_operons";
  }
  $outfile{prefix} = join("/", $dir{output_per_query}, $outfile_prefix);
  &RSAT::message::Info("Automatic definition of the output prefix", $outfile{prefix}) if ($main::verbose >= 4);

  &RSAT::message::Info("&GetOutfilePrefix()", 
  		       "query_prefix", $outfile_prefix,
  		       "dir{output_per_query}", $dir{output_per_query},
  		       "outfile{prefix}", $outfile{prefix},
  		      ) if ($main::verbose >= 5);

  return ($outfile_prefix, $query_prefix);
}


################################################################
## Initialize a query-specific output directory.
##
## The root output directory is either defined by the user (option -o)
## or set to the default value "footprints".
##
##
## The output root directory can either be specified by the user, or takes the default value "footprints".
##
##
sub InitQueryOutput {

  my ($outfile_prefix, $query_prefix) = &GetOutfilePrefix();

  ## Create output directory if required
  $dir{output} = `dirname $outfile{prefix}`;
  chomp($dir{output});
  &RSAT::util::CheckOutDir($dir{output});

  ## Open output stream for the log file
  $outfile{log} = $outfile{prefix}."_log.txt";
  $out = &OpenOutputFile($outfile{log});

  ## File for storing the list of query gene names
  $outfile{genes} = $outfile{prefix}."_query_genes.tab";
  #$outfile{genes_info} = $outfile{prefix}."_query_genes_info.tab";
  $genes = &OpenOutputFile($outfile{genes});

  ## Specify other file names
  $outfile{orthologs} = $outfile{prefix}."_ortho_bbh.tab"; ## orthologs of the query gene(s)
  $outfile{query_seq} = $outfile{prefix}."_query_seq.fasta"; ## Upstream sequence(s) of the query gene(s)

  ## type of promoter: either use directly the promoter of each query
  ## gene (ortho) or of the predicted operon leader gene
  if ($infer_operons) {
    $outfile{leader_qgenes} = $outfile{prefix}."_leader_query_genes.tab"; ## Leader genes of the operons containing query genes
    $promoter = "leaders";
  } else {
    $promoter = "ortho";
  }
  $outfile{bbh} = $outfile{prefix}."_".$promoter."_bbh.tab"; ##
  $outfile{seq} = $outfile{prefix}."_".$promoter."_seq.fasta"; 
  $outfile{purged_notclean} = $outfile{prefix}."_".$promoter."_seq_purged_notclean.fasta" unless $main::no_purge;
  $outfile{purged} = $outfile{prefix}."_".$promoter."_seq_purged.fasta" unless $main::no_purge;

  return($outfile_prefix, $query_prefix);
}



################################################################
## Read the options that are common to footprint detection programs.
sub ReadFootprintOptions {

      ## Verbosity

=pod

=item B<-v #>

Level of verbosity (detail in the warning messages during execution)

=cut
  if ($arg eq "-v") {
    if (&IsNatural($arguments[0])) {
      $main::verbose = shift(@arguments);
    } else {
      $main::verbose = 1;
    }

    ## Help message
=pod

=item B<-h>

Display full help message.

=cut
  } elsif ($arg eq "-h") {
    &PrintHelp();

	    ## List of options
=pod

=item B<-help>

Same as -h

=cut
  } elsif ($arg eq "-help") {
    &PrintOptions();

    ## Organism
=pod

=item	B<-org query_organism>

Query organism, to which the query genes belong.

=cut
  } elsif ($arg eq "-org") {
    $main::organism_name = shift(@arguments);

    ## Taxon
=pod

=item	B<-taxon reference_taxon>

Reference taxon, in which orthologous genes have to be collected.

Alternatively, reference organisms can be specified with the option
-org_list.

=cut
  } elsif ($arg eq "-taxon") {
    $main::taxon = shift(@arguments);
    &RSAT::error::FatalError("Options -taxon and -org_list are mutually incompatible") if ($main::orglist_file);

## Organism List
=pod

=item	B<-org_list organisms_list_file>

This option gives the posibility to analyse a user-specified set of
reference organisms rather than a full taxon.

File format: the first word of each line is used as organism ID. Any
subsequent text is ignored. The comment char is ";".

This option is incompatible with the option "-taxon".

=cut
  } elsif ($arg eq "-org_list") {
    $main::orglist_file = shift(@arguments);
    &RSAT::error::FatalError("Options -taxon, -org_list and -orthologs_list are mutually incompatible") if ($main::taxon || $main::orthologs_list_file);

## No Purge
=pod

=item	B<-no_purge>

This option can only be used combined with the -org_list option, this
gives the posibility to analyse a given set of sequences managing
sequence redundancy using a list of "no redundant" organisms.

The file format is one organisms per line, the comment char is ";"

=cut
  } elsif ($arg eq "-no_purge") {
    $main::no_purge = 1;

     ## Orthologs List
=pod

=item	B<-orthologs_list file>

This option gives the posibility to analyse a user-specified set of
orthologs for specific reference organisms instead of using the BBH set of orthologs provided by RSAT.

The query genes included here will be the ones analyzed by the program.

File format: Tab delimited file with three columns. 

  ID of the query gene (in the query organism)
  ID of the reference gene 
  ID of the reference organism

Further columns will be ignored.
The comment char is ";".

This option is incompatible with the option "-taxon", and "-bg_model taxfreq" option.

=cut
  } elsif ($arg eq "-orthologs_list") {
    $main::orthologs_list_file = shift(@arguments);
    &RSAT::error::FatalError("Options -taxon, org_list and -orthologs_list are mutually incompatible") if ($main::taxon || $main::orglist_file);

    ## Query gene
=pod

=item	B<-q query>

Query gene.

This option can be used iteratively on the command line to specify
multiple genes.

=cut
  } elsif ($arg eq "-q") {
    push @main::query_genes, shift(@arguments);


=pod

=item B<-genes>

Specify a file containing a list of genes. Multiple genes can also be
specified by using iteratively the option -q.

=cut
  } elsif ($arg eq "-genes") {
    $main::infile{genes} = shift(@arguments);

=pod

=item B<-all_genes>

Automatically analyze all the genes of a query genome, and store each
result in a separate folder (the folder name is defined
automatically).

=cut
  } elsif ($arg eq "-all_genes") {
    $main::all_genes = 1;
    $main::sep_genes = 1;

=pod

=item B<-max_genes>

Maximal number of genes to analyze.

=cut
  } elsif ($arg eq "-max_genes") {
    $main::max_genes = shift(@arguments);


=pod

=item B<-skip #>

Skip the first # genes (useful for quick testing and for resuming
interrupted tasks).

=cut

  } elsif ($arg eq "-skip") {
    $main::skip = shift(@arguments);
    &RSAT::error::FatalError("Invalid number with option -skip\t$skip") unless &IsNatural($skip);

=pod

=item B<-last #>

Stop after having treated the first # genes (useful for quick
testing).

=cut

  } elsif ($arg eq "-last") {
    $main::last = shift(@arguments);
    &RSAT::error::FatalError("Invalid number with option -last\t$last") unless &IsNatural($last);

	    ## Output prefix
=pod

=item	B<-o output_root_dir>

Main output directory. The results will be dispatched in
sub-directories, defined according to the taxon, query organism and
query gene name(s).

If the main output dir is not specified, the program automatically
sets it to "footprints".

=cut
  } elsif ($arg eq "-o") {
#    $main::outfile{prefix} = shift(@arguments);
    $main::dir{output_root} =  shift(@arguments);

=pod

=item B<-batch>

Generate one command per query gene, and post it on the queue of a PC
cluster.

=cut
  } elsif ($arg eq "-batch") {
    $main::batch = 1;


    ## Dry run
=pod

=item B<-dry>

Dry run: print the commands but do not execute them.

=cut
  } elsif ($arg eq "-dry") {
    $main::dry = 1;;

    ## Don't die on error
=pod

=item B<-nodie>

Do not die in case a sub-program returns an error.  

The option -nodie allows you to circumvent problems with specific
sub-tasks, but this is not recommended because the results may be
incomplete.

=cut

} elsif ($arg eq "-nodie") {
  $main::die_on_error = 0;

=pod

=item B<-sep_genes>

Search footprints for each query gene separately. The results are
stored in a separate folder for each gene. The folder name is defined
automatically.

By default, when several query genes are specified, the program
collects orthologs and analyzes their promoters altogether. The option
I<-sep> allows to automatize the detection of footprint in a set of
genes that will be treated separately.

=cut
  } elsif ($arg eq "-sep_genes") {
    $main::sep_genes = 1;

	  ## infer operons
=pod

=item B<-infer_operons>

Infer operons in order to retrieve the promoters of the predicted
operon leader genes rather than those located immediately upstream of
the orthologs. This method uses a threshold on the intergenic distance.

=cut

	} elsif ($arg eq "-infer_operons") {
	  $main::infer_operons = 1;


=pod

=item B<-task>

Specify a subset of tasks to be executed.

By default, the program runs all necessary tasks. However, in some
cases, it can be useful to select one or several tasks to be executed
separately. For instance, after having collected all the promoter
sequences of ortholog genes, one might desire to run the pattern
detection with various parameter values without having to retrieve the
same sequences each time.

Beware: task selection requires expertise, because most tasks depends
on the prior execution of some other tasks in the workflow. Selecting
tasks before their prerequisite tasks have been completed will provoke
fatal errors.

B<Supported tasks:>

=over

=item For all footprint programs (I<footprint-discovery>,
I<footprint-scan>).

=over

=item I<all>

Run all supported tasks. If no task is specified, all tasks are
performed.

=item I<operons>

Infer operons (using I<infer-operons>. This option should be used only for
Bacteria.

=item I<query_seq>

Retrieve upstream sequence of the query genes (using I<retrieve-seq>).

=item I<orthologs>

Identify theorthologs of the query genes in the selected taxon (using
I<get-orthologs>).

=item I<ortho_seq>

Retrieve upstream sequences of the orthologs (using
I<retrieve-seq-multigenome>).

=item I<purge>

Purge upstream sequences of the orthologs (using I<purge-seq>).


=item I<index>

Generate an HTML index with links to the result files. This option is
used for the web interface, but can also be convenient to index
results, especially when several genes or taxa are analyzed (options
-genes, -all_genes, -all_taxa).

With the option -sep_genes, one index is generated for each gene
separately. An index summarizing the results for all genes can be
generated using the option -synthesis.

=item I<synthesis>

(still to be implemented)

Generate a HTML table with links to the individual result files. The
table contains one row per query gene, one column by output type
(sequences, dyads, maps, ...).

=back

=item Tasks specific to I<footprint-discovery>

=over

=item I<filter_dyads>

Detect all dyads present with at elast one occurrence in the upstream
sequence of the query gene (using I<dyad-analysis>). Those dyads will
be used as filter if the option I<-filter> has been specifed.

=item I<dyads>

Detect significantly over-represented in upstream sequences of
orhtologs (using I<dyad-analysis>).

=item I<map>

Draw feature maps showing the location of over-represented dyads in
upstream sequences of promoters (using I<feature-map>).

=item I<index>

Generate an index file for each gene separately. The index file is in
the gene-specific directory, it is complementary to the general index
file generated with the task "synthesis".

=back

=item Tasks specific to I<footprint-scan>

=over

=item I<occ_sig>

Compute the significance of number of matrix hit occurrences as a
function of the weight score (I<using matrix-scan> and
I<matrix-scan-quick>).

=item I<occ_sig_graph>

Generate graphs showing the distributions of occurrences and their
significances, as a function of the weight score (using >XYgraph>).

=item scan

Scan upstream sequences to detect hits above a given threshold (using
I<matrix-scan>).

=item I<map>

Draw the feature map of the hits (using I<feature-mp>).

=back

=back

=cut
    } elsif ($arg eq "-task") {
      my @requested_tasks = split ",", shift (@arguments);
      foreach my $task (@requested_tasks) {
	next unless $task;
	if ($supported_task{$task}) {
	  $task{$task} = 1;
	} else {
	  &RSAT::error::FatalError("Task '$task' is not supported. \n\tSupported: $supported_tasks");
	}
      }

=pod

=item B<-map_format>

Format for the feature map.

Supported: any format supported by the program feature-map.

=cut
    } elsif ($arg eq "-map_format") {
      $main::map_format = shift(@arguments);


    ## Create HTML Index for each gene separately
=pod

=item B<-index>

Deprecated, replaced by the task "index".

=cut
  } elsif ($arg eq "-index") {
    $main::create_index = 1;
    &RSAT::message::Warning("Option -index is deprecated, indexes are now always created.");

    ## Create a tab-delimited file and a HTML Index for all the results
=pod

=item B<-synthesis>

This option generate synthetic tables (in tab-delimited text and html)
for all the results. It should be combined with the option
-sep_genes. The synthetic tables contain one row per gene, and one
column per parameter. They summarize the results (maximal
significance, top-ranking motifs) and give pointers to the separate
result files.

=cut
  } elsif ($arg eq "-synthesis") {
    $main::synthesis = 1;


  } else {
    return(0); ## No option was found
  }

  return(1);

=pod

=back

=cut

}


################################################################
## Open file for the HTML index
sub OpenIndex {
  my ($program_name) = @_;
  $outfile{index} = $outfile{prefix}."_index.html";
  my $outfile_prefix = $dir{query_prefix};
  $index_list{$outfile_prefix} = $outfile{index};
  $index = &OpenOutputFile($outfile{index});
  print $index "<html>\n";
  $html_title = $outfile_prefix;
  $html_title .= " ".$taxon." " if ($taxon);
  $html_title .= " ".$organism_name if ($organism_name);
  $html_title .= " ".$bg_model if ($bg_model);
  $html_title .= "Adaptive BG model, window size ".$window_size if ($window_size);
  print $index "<head><title>", $html_title , "</title></head>\n" ;
  print $index "<body>\n";
  print $index "<hr size=4 color='#000088'>";

  print $index "<h1 align=center>",$program_name, " result</h1>"  ;
  $html_title2 = "<i>".$organism_name."</i>";
  $html_title2 .= " ".$bg_model unless ($program_name eq "footprint-scan");
  $html_title2 .= " ".$query_prefix;
  $html_title2 .= "; ".$taxon if $taxon;
  print $index "<h2 align=center>",$html_title2 , "</h2>\n";
  print $index "<hr size=2 color='#000088'>";
  print $index "<blockquote>";
  print $index "<table cellspacing=0 cellpadding=3 border=0>\n";
  &IndexOneFile("log", $outfile{log});
  &IndexOneFile("input", $infile{genes}) if (($infile{genes}) && !($main::sep_genes));
}


################################################################
## Compute a prefix for the "main" files (index for all queries, list
## of dyad files, ...)
sub MainPrefix {
  $outfile{main_prefix} = join("_", ($taxon||"org_list"),$organism_name,
			       "bg", $bg_model);
  return $outfile{main_prefix};
}

################################################################
## Compute the name of the main index file. This function is called by
## the perl script footprint-discovery, but also y the CGI script
## footprint-discovery.cgi in the Web interface. The goal of the
## function is to ensure consistency between the index file names.
sub MainIndexFileName {
  my $main_prefix = &MainPrefix();
  $outfile{main_index_file} = $main_prefix."_result_index.html";
  return($outfile{main_index_file});
}

################################################################
## Main index. This is a HTML table with links to the query-specific
## results: one row per query, one column per output type.
sub OpenMainIndex {
  &RSAT::util::CheckOutDir($dir{output_root});

  $outfile{main_index} = $dir{output_root}."/".&MainIndexFileName();

  &RSAT::message::Info("Main index file", $outfile{main_index}) if ($main::verbose >= 1);

  my $main_index = &OpenOutputFile($outfile{main_index});
  print $main_index "<html>\n";
  $html_title = "footprint-discovery";
  $html_title .= " ".$taxon if ($taxon);
  $html_title .= " ".$organism_name if ($organism_name);
  $html_title .= " ".$bg_model if ($bg_model);
  print $main_index "<head><title>", $html_title , "</title></head>\n" ;
  print $main_index &sorttable_script();
  print $main_index "<body>\n";
  print $main_index "<h1>", $html_title, "</h1\n";
  print $main_index "<p><b>Command:</b> footprint-discovery";
  $arguments =  &RSAT::util::hide_RSAT_path(&PrintArguments());
  print $main_index $arguments;
  print $main_index "</p>\n";

  ## Open the index table
  print $main_index "<p><table class='sortable' border='0' cellpadding='3' cellspacing='0'>\n";
#  print $main_index "<table cellspacing=0 cellpadding=3 border=0>\n";
  print $main_index "<tr>\n";
  print $main_index "<th>","Query nb","</th>\n";
  print $main_index "<th>","Query","</th>\n";
  print $main_index "<th>","Top dyad","</th>\n";
  print $main_index "<th>","Max sig","</th>\n";
  print $main_index "<th>","Nb dyads","</th>\n";
  print $main_index "</tr>\n";

  return ($main_index);
}
################################################################
## Main index for footprint-scan. This is a HTML table with links to the query-specific
## results: one row per query, one column per output type.
sub OpenSynthesisFilesScan {
  &RSAT::util::CheckOutDir($dir{output_root});
  $outfile{synthesis} = join( "/",$dir{output_root}, $used_orgs, $organism_name,"result_synthesis");
  $outfile{synthesis_tab} = $outfile{synthesis}.".tab";
  $outfile{synthesis_html} = $outfile{synthesis}.".html";

  my $synthesis_index = &OpenOutputFile($outfile{synthesis_html} ); 
  my $synthesis_table = &OpenOutputFile($outfile{synthesis_tab});

  return ($synthesis_index, $synthesis_table)  ;
}
sub HeaderSynthesisFilesScan {
    ($synthesis_index, $synthesis_table, $occ_one_gene_file) =@_;
    $occ_sig_header = `grep '^#' $occ_one_gene_file`;
    chomp $occ_sig_header;
    $occ_sig_header =~ s/occ_sig_rank/gene_rank/;
    $occ_sig_header =~ s/^#//;
    print $synthesis_table join("\t", "#gene", $occ_sig_header,"name","descr","upstr_neighb_name","file"), "\n";

    print $synthesis_index "<html>\n";
    $html_title = "footprint-scan";
    $html_title .= " ".$taxon if ($taxon);
    $html_title .= " ".$organism_name if ($organism_name);
    $html_title .= " ".$bg_model if ($bg_model);
    print $synthesis_index "<head><title>", $html_title , "</title>";
    print $synthesis_index &sorttable_script();
    $header .= "<style type='text/css'>\n";
    $header .= `cat $ENV{RSAT}/perl-scripts/lib/results.css`;
    $header .=  "</style>\n";
    print $synthesis_index "</head>\n";
    print $synthesis_index "<body>\n";
    print $synthesis_index "<h1>", $html_title, "</h1\n";
    print $synthesis_index "<p><b>Command:</b> footprint-scan";
    print $synthesis_index &PrintArguments();
    print $synthesis_index "</p>\n";

    ## Open the index table
    print $synthesis_index "<p><table class='sortable' border='0' cellpadding='3' cellspacing='0'>\n";
    print $synthesis_index "<tr>\n";
     my $html_sep="<\/th>\n<th>";
    $occ_sig_header=~s/\t/$html_sep/g;
    my $header_index=join("</th>\n<th>", "<th>gene", $occ_sig_header, "OCC_Sig", "SigPlot" ,"FreqPlot","Map", "name","descr","upstr_neighb_name</th>");
    $header_index .= "\n"; 
    print $synthesis_index $header_index;
    print $synthesis_index "</tr>\n";
}

################################################################
## Add one file to the index file
sub IndexOneFile {
  my ($name, $file, %args) = @_;
  my $short_file = `basename $file`;
  print $index "<tr valign=top>\n";
  print $index "<td>", $name, "<td>",  &LinkOneFile($outfile{index}, $file, $short_file), "</td>\n";

  if ($args{image}) {
    print $index "</tr><tr><td colspan=\"2\">(Click on image below)</td></tr><tr><td colspan=\"2\"><a href=".$short_file."><img width=\"100%\" src=".$short_file."></a></td>\n";
  }
  print $index ("</tr>\n\n");
}


################################################################
## Return a HTML link from an index file to a target file
sub LinkOneFile {
  my ($from_file, $to_file, $text) = @_;
  my $path = "";
  my $link = "";
  $text = $text || $to_file;
  if (-e ($to_file)) {
    $path = &RSAT::util::RelativePath($from_file, $to_file);
    $link = join ("",  "<a href='",$path, "'>", $text, "</a>");
  } else {
    $link = join ("",  "<font color='red'>", $text, "</font>");
  }
  return $link;
}


################################################################
## Predict operon leader genes of the query genes
sub InferQueryOperons {
  &RSAT::message::TimeWarn("Get leaders of query genes (d<=".$dist_thr."bp)", $outfile{leader_qgenes}) if ($main::verbose >= 2);
  &CheckDependency("operons", "genes");
  my $cmd = "$SCRIPTS/get-leader-multigenome ";
  $cmd .= " -i ".$outfile{genes};
  $cmd .= " -o ".$outfile{leader_qgenes};
  $cmd .= " -uth interg_dist ".$dist_thr;
  &one_command($cmd) if ($task{operons});
  #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
  &IndexOneFile("leader query genes", $outfile{leader_qgenes});
}


################################################################
## Retrieve promoters of the query organism
sub RetrieveQueryPromoters {
  if ($task{query_seq}) {
    &RSAT::message::TimeWarn("Retrieving promoter sequences for query genes", $outfile{query_seq}) if ($main::verbose >= 2);
    my $cmd = "$SCRIPTS/retrieve-seq ";
    $cmd .= " -org ".$organism_name;
    if ($infer_operons) {
      &CheckDependency("query_seq", "leader_qgenes");
      $cmd .= " -i ".$outfile{leader_qgenes};
    } else {
      &CheckDependency("query_seq", "genes");
      $cmd .= " -i ".$outfile{genes};
    }
    $cmd .= " -noorf";
    $cmd .= " -feattype misc_RNA,cds,trna,rrna";
    $cmd .= " -o ".$outfile{query_seq};
    &RSAT::message::Info("Retrieve seq command", $cmd) if ($main::verbose >= 5);
    &one_command($cmd);
  }
  &IndexOneFile("query sequence", $outfile{query_seq});
}

################################################################
## Detect all dyads in promoters of query genes for dyad filtering
sub ComputeFilterDyads {
  if ($task{filter_dyads}) {
    &RSAT::message::TimeWarn("Computing filter dyads", $outfile{filter_dyads}) if ($main::verbose >= 2);
    &CheckDependency("filter", "query_seq");
    my $cmd = "$SCRIPTS/dyad-analysis -v 1 -return occ -lth occ 1";
    $cmd .= " -i ".$outfile{query_seq};
    $cmd .= " -l 3 -sp 0-20";
    $cmd .= " ".$strands;
    $cmd .= " ".$noov;
    $cmd .= " -o ".$outfile{filter_dyads};
    &one_command($cmd);
    #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
  }
  &IndexOneFile("filter dyads", $outfile{filter_dyads});
}

################################################################
## Detect all matrix hits  in promoters of query genes for scan filtering
sub ComputeFilterScan {
    my ($matrix_format2, @matrix_files2)=@_;
    $main::skip_gene=0;
    &RSAT::message::TimeWarn("Computing filter hits", $outfile{filter_scan}) if ($main::verbose >= 2);
    &CheckDependency("filter", "query_seq");
    my $cmd = "$SCRIPTS/matrix-scan -v 1 -return sites,pval,rank -uth pval " .$main::filter_pval;
    $cmd .= " -uth rank 1 ";
    $cmd .= " -i ".$outfile{query_seq};
    foreach my $file (@matrix_files2) {
	&RSAT::error::FatalError("Matrix file $file does not exist.Matrix file is mandatory.")  unless (-e $file ) ;
	$cmd .= " -m ".$file;
    }
    $cmd .= " -bgfile ".$main::filter_bgfile ; 
    $cmd .= " -matrix_format ".$matrix_format2;
    $cmd .= " ".$strands;
    $cmd .= " -quick ";
    #$cmd .= " ".$noov;
    $cmd .= " -o ".$outfile{filter_scan};
    &one_command($cmd) if ($task{filter_scan});
    #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
    &IndexOneFile("filter scan", $outfile{filter_scan});

    my $filterg=`grep -v ";" $outfile{filter_scan}  | grep -v "#" | cut -f1`;
    chomp($filterg);
    $outfile{genes}= $outfile{prefix}."_filter_query_genes.tab";
    $filter_genes = &OpenOutputFile($outfile{genes});
    print  $filter_genes $filterg;
    close $filter_genes if ($outfile{genes} && $main::filter);
    &RSAT::message::Info("Filter genes  ", $outfile{genes} ) if ($main::verbose >= 1);  
    &IndexOneFile("filter genes", $outfile{genes});
    $main::skip_gene=1 unless ($filterg=~/\w/);
}

################################################################
## Identify ortholog genes
sub GetOrthologs {
  if ($task{orthologs}) {
    &RSAT::message::TimeWarn("Getting orthologs", $outfile{orthologs}) if ($main::verbose >= 2);
    &CheckDependency("orthologs", "genes");
    my $cmd = "$SCRIPTS/get-orthologs";
    $cmd .= " -i ".$outfile{genes};
    $cmd .= " -org ".$organism_name;
    if ($main::tf_ortho_file) {
      $cmd .= " -org_list ". $main::tf_ortho_file ;
      &RSAT::message::Info ("Getting orthologs using option -org_list ", $main::tf_ortho_file) if ($main::verbose >= 2);
    } elsif ($main::orglist_file) {
      $cmd .= " -org_list ". $main::orglist_file ;
      &RSAT::message::Info ("Getting orthologs using option -org_list ", $main::orglist_file ) if ($main::verbose >= 2);
    } else {
      $cmd .= " -taxon ".$taxon ;
    }
    $cmd .= " -return query_name,query_organism -return ident";
    $cmd .= " -uth rank 1";	## BBH criterion
    $cmd .= " -lth ali_len 50";
    $cmd .= " -uth e_value 1e-05";
    $cmd .= " -return e_value";
    $cmd .= " -only_blast";	## only use genome having blast files
    $cmd .= " -o ".$outfile{orthologs};
    &one_command($cmd);
    #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
  }
  elsif ($main::orthologs_list_file){
      #here introduce addaptation for the file so it can be used on next steps
      my $genes=`cut -f1 $outfile{genes} |xargs`;
      chomp($genes);
      $genes=~s/\s+/\|/g;
      my $orthologs=`egrep "$genes" $main::orthologs_list_file | perl -ane 'print "\$F[1]\t\$F[2]\t\$F[0]\t$organism_name\n"' `;
      my $out= &OpenOutputFile($outfile{orthologs});
      print $out $orthologs; 
      close $out;
      &RSAT::message::Info("Orthologs for gene(s)",$genes, "specified by the user can be found in " , $outfile{orthologs}) if ($main::verbose >= 0);
      #die"HELLO";
  }
  &IndexOneFile("orthologs", $outfile{orthologs});
}

################################################################
## Predict operon leader genes for the orthologous genes
sub InferOrthoOperons {
  if  ($task{operons}) {
    &RSAT::message::TimeWarn("Get leaders of query genes (d<=".$dist_thr."bp)", $outfile{bbh}) if ($main::verbose >= 2);
    &CheckDependency("operons", "orthologs");
    my $cmd = "$SCRIPTS/get-leader-multigenome ";
    $cmd .= " -i ".$outfile{orthologs};
    $cmd .= " -o ".$outfile{bbh};
    $cmd .= " -uth interg_dist ".$dist_thr;
    &one_command($cmd) ;
  }
  #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
  &IndexOneFile("leader genes", $outfile{bbh});
}


################################################################
## Retrieve sequences from orthologs
sub RetrieveOrthoSeq {
  if ($task{ortho_seq}) {
    &RSAT::message::TimeWarn("Retrieving promoter sequences of orthologs", $outfile{seq}) if ($main::verbose >= 2);
    &CheckDependency("ortho_seq", "bbh");
    my $cmd = "$SCRIPTS/retrieve-seq-multigenome -ids_only";
    $cmd .= " -i ".$outfile{bbh};
    $cmd .= " -noorf";
    $cmd .= " -feattype CDS,mRNA,tRNA,rRNA,scRNA,misc_RNA" ;
    $cmd .= " -o ".$outfile{seq};
    &one_command($cmd);
    #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
  }
  &IndexOneFile("$promoter sequences", $outfile{seq});
}



################################################################
## Purge sequences
sub PurgeOrthoSeq {
  if ($task{purge}) {
    &RSAT::message::TimeWarn("Purging promoter sequences of orthologs", $outfile{purged}) if ($main::verbose >= 2);
    &CheckDependency("purge", "seq");
    my $cmd = "$SCRIPTS/purge-sequence";
    $cmd .= " -nodie" if ($main::die_on_error == 0);
    $cmd .= " -i ".$outfile{seq};
    $cmd .= " -ml 30 -mis 0 -mask_short 30";
    $cmd .= " ".$strands;
    $cmd .= " -o ".$outfile{purged_notclean};
    &one_command($cmd);

    $cmd = "$SCRIPTS/convert-seq";
    $cmd .= " -i ".$outfile{purged_notclean}; ;
    $cmd .= " -mask non-dna ";
    $cmd .= " -from fasta ";
    $cmd .= " -to fasta ";
    $cmd .= " -dna ";
    $cmd .= " -o ". $outfile{purged} ;
    &one_command($cmd);
  }
  &IndexOneFile("purged sequences", $outfile{purged});
#  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
}



################################################################
## Detect over-representation of matching occurrecnes (hits) of the motif
sub OccurrenceSig {
  if ($task{occ_sig}) {
    &RSAT::message::TimeWarn("Testing over-representation of hits", $outfile{occ_sig}) if ($main::verbose >= 2);
    if ($main::no_purge){
      &CheckDependency("occ_sig", "seq");
    }
    else {
      &CheckDependency("occ_sig", "purged");
    }
    my $cmd = "matrix-scan";
    $cmd .= $ms_parameters;
    $cmd .= " -return distrib,occ_proba,rank -sort_distrib";
    #  $cmd .= " -lth inv_cum 1 -lth occ_sig 0 -uth occ_sig_rank 1";
    if ($main::no_purge){
      $cmd .= " -i ".$outfile{seq};
    }else{
      $cmd .= " -i ".$outfile{purged};
    }
    $cmd .= " -o ".$outfile{occ_sig};
    $cmd .= $occ_sig_opt;
    &one_command($cmd);
  }
  &IndexOneFile("occ sig", $outfile{occ_sig});
}

################################################################
## Draw a graph with occurrence significance profiles
sub OccurrenceSigGraph {
  if ($task{occ_sig_graph}) {
    &RSAT::message::TimeWarn("Graph with over-representation of hits", $outfile{occ_freq_graph}) if ($main::verbose >= 2);
    &CheckDependency("occ_sig_graph", "occ_sig");


    ## Occ frequency graph
    my $cmd = "sort -n -k 2 $outfile{occ_sig} | XYgraph";
    $cmd .= " -xcol 2 -xleg1 'Weight score' -xsize 800 -xgstep1 1 -xgstep2 0.5";
    $cmd .= " -ycol 5,8 -yleg1 'Hit numbers' -ylog 10";
    $cmd .= " -title 'matrix ".$matrix_suffix." ; gene ".$current_gene."'";
    $cmd .= " -lines -legend ";
    $cmd .= " -format ".$plot_format;
    $cmd .= " -o ".$outfile{occ_freq_graph};

    &one_command($cmd);

    ## Occ significance graph
    $cmd = "sort -n -k 2 $outfile{occ_sig} | XYgraph";
    $cmd .= " -xcol 2 -xleg1 'Weight score' -xsize 800  -xgstep1 1 -xgstep2 0.5";
    $cmd .= " -ycol 11 -yleg1 'Binomial significance of hits'";
    $cmd .= " -title 'matrix ".$matrix_suffix." ; gene ".$current_gene."'";
    $cmd .= " -lines -legend ";
    $cmd .= " -format ".$plot_format;
    $cmd .= " -o ".$outfile{occ_sig_graph};

    ##  info lines
    if ($main::draw_info_lines) {
      $cmd .= " -hline violet 0 "; #line showing the significance zero line.
      $cmd .= " -vline violet 0 " ; # line showing the score zero value 

      #Calculate the positive score wirh maximal significance
      my $tab="\t";
      my $new_line="\n"; 
      my $top_sig_cmd = "grep -v '^;' $outfile{occ_sig} | grep -v '^#' | perl -ane 'if(\$F[1]>0) {print join(\"\t\",\@F).\"\n\" }' | head -n 1";
      my $top_sig_pos_score_row = `$top_sig_cmd`;
      if ($top_sig_pos_score_row) {
	my @fields = split "\t", $top_sig_pos_score_row;
	my $sig_max = $fields[1];
	$cmd .= " -vline red ". $sig_max;
      }
    }

    ## options added  from comand line
    $cmd .= " ".$occ_sig_graph_opt;
    &one_command($cmd);
  }
  &IndexOneFile("occ freq graph", $outfile{occ_freq_graph});
  &IndexOneFile("occ sig graph", $outfile{occ_sig_graph});

}

################################################################
## Collect info for the synthetic table
sub GetTopSig {
  &CheckDependency("synthesis", "occ_sig");
  my $top_cmd = "grep -v '^;' $outfile{occ_sig} | grep -v '^#' | awk '\$12==1'";
  my $top_sig_row = `$top_cmd`;
  if ($top_sig_row) {
    my @fields = split "\t", $top_sig_row;

    $top_sig{$current_gene} = $fields[10];;
    $top_score{$current_gene} = $fields[1];
    $top_sig_row{$current_gene} = join ("\t", @fields[0..10]);
#    &RSAT::message::Debug("Top sig", $current_gene, $top_sig{$current_gene}, "score", $top_score{$current_gene}) if ($main::verbose >= 5);

    ## Index occ sig files for the synthetic table
    $occ_sig_file{$current_gene} = $outfile{occ_sig} if (-e $outfile{occ_sig} );
    $occ_freq_graph_file{$current_gene} = $outfile{occ_freq_graph} if (-e $outfile{occ_freq_graph});
    $occ_sig_graph_file{$current_gene} = $outfile{occ_sig_graph} if (-e $outfile{occ_sig_graph});
    $gene_list{$current_gene} = $outfile{genes} if (-e $outfile{genes} )
  }

  ## Index scan files for the synthetic table
  $map_file{$current_gene} = $outfile{map} if (-e $outfile{map});
}


################################################################
## Scan sequences with PSSM to locate sites
sub OrthoScan {
  if ($task{scan}) {
    &RSAT::message::TimeWarn("Scannig sequences to detect sites", $outfile{sites}) if ($main::verbose >= 2);
    &CheckDependency("scan", "seq");
    my $cmd = "matrix-scan";
    $cmd .= $ms_parameters;
    $cmd .= " -return limits,sites,rank";
    $cmd .= " -i ".$outfile{seq};
    $cmd .= " -o ".$outfile{sites};
    $cmd .= " ".$scan_opt;
    &one_command($cmd);
  }
  &IndexOneFile("sites", $outfile{sites});
}

################################################################
## Draw a feature map with the detected sites
sub OrthoMap {
  if ($task{map}) {
    &RSAT::message::TimeWarn("Drawing feature map", $outfile{map}) if ($main::verbose >= 2);
    &CheckDependency("map", "sites");
    my $cmd = "feature-map -i ".$outfile{sites};
    $cmd .= " -scalebar -legend";
    $cmd .= " -xsize 800 -scorethick -minscore 0";
    if ($taxon) { 
      $cmd .= " -title 'matrix hits in ".$taxon." promoters'";
    }
    else {
      $cmd .= " -title 'matrix hits in promoters'";
    }
    $cmd .= " -format ".$map_format;
    $cmd .= " -o ".$outfile{map};
    $cmd .= " ".$map_opt;
    &one_command($cmd);
  }
  &IndexOneFile("map", $outfile{map});
}


1;

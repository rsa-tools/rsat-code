################################################################
## Library for functions shared by footprint detection programs
## - footprint-discovery
## - footprint-scan

%supported_task =  (operons=>1,
		    query_seq=>1,
		    orthologs=>1,
		    ortho_seq=>1,
		    purge=>1,
		    all=>1,
		    );
%task = ();
$supported_tasks = join (",", keys %supported_task);
local $all_genes = 0;         ## Analyze all the genes of the query organism

################################################################
## Treat one command, by either executing it, or concatenating it for
## further batch processing
sub one_command {
  my ($cmd) = @_;
  if ($main::batch) {
    if ($main::batch_cmd =~/\S/) {
      $main::batch_cmd .= " ; $cmd";
    } else {
      $main::batch_cmd = "$cmd";
    }
  } else {
    print $out ("; ", &AlphaDate(), "\n", $cmd, "\n\n") if ($main::verbose >= 1);
    &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
  }
}


################################################################
## Check all parameters required for footprint analysis (discovery or
## scanning).
sub CheckFootprintParameters {

  ## If all tasks are requested or if no task is defined, execute all
  ## tasks.
  if ((scalar(keys(%task)) == 0) || ($task{all})){
    foreach my $task (keys %supported_task) {
#      &RSAT::message::Debug("Auto adding task", $task);
      $task{$task} = 1;
    }
  }
  &RSAT::message::Info("Tasks: ", join (";", keys %task)) if ($main::verbose >= 0);

  ## Check taxon
  &RSAT::error::FatalError("You must specify a taxon (option -taxon)")
    unless ($taxon);
  my @organisms = &CheckTaxon($taxon);

  ## Check organism
  &RSAT::error::FatalError("You must specify a organism (option -org)")
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
      next if (/^--/); ## Skip mysql-type comment lines
      next if (/^;/); ## Skip comment lines
      next if (/^#/); ## Skip header lines
      next unless (/\S/); ## Skip empty lines
      chomp;
      my @fields = split /\s+/;
      my $query = shift @fields;
      push @query_genes, $query;
    }
    close $in;
  }


  ################################################################
  ## Check query genes
  if (scalar(@query_genes) ==0) {
    &RSAT::error::FatalError("You must specify at least one query gene (options -q or -i)");
  } else {
    &RSAT::message::Info("Number of query genes", scalar(@query_genes)) if ($main::verbose >= 1);
  }

  ## Get maximum number of genes if limited
  if ($max_genes){
    if (scalar(@query_genes)>$max_genes){
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
    if ((-e $file) || (-e $file.".gz")){
      &RSAT::message::Info("Checked existence of ", $type, "file required for task", $task, "file", $file)
	if ($main::verbose >= 3);
      return(1);
    } else {
      &RSAT::error::FatalError("Missing", $type, "file required for task", $task, "file", $file)
    }
  }
}

################################################################
## Define a query prefix for the title of the feature map and for
## automatic output file specification
sub GetQueryPrefix {
  my $query_prefix;
  if (scalar(@current_query_genes) == 1) {
    $query_prefix = $current_query_genes[0];
  } elsif (scalar(@current_query_genes) <= 10) {
    $query_prefix = join "_", @current_query_genes;
  } elsif ($outfile{prefix}) {
    $query_prefix = `basename $outfile{prefix}`;
    chomp($query_prefix);
  } elsif ($infile{genes}) {
    $query_prefix = `basename $infile{genes} .tab`;
  }
  &RSAT::message::Info("Query prefix", $query_prefix) if ($main::verbose >= 3);
  return ($query_prefix);
}


################################################################
## Output prefix is mandatory
## If not specified by the user, define it automatically
## If genes are to be analyzed separately, also define output prefix automatically
sub GetOutfilePrefix {
  if (($sep_genes)
      || ($outfile{prefix} eq "")) {
    if ($query_prefix) {
      if ($bg_model) {
	$outfile{prefix} = join( "/", "footprints", $taxon, $organism_name, $query_prefix, $bg_model,
				 join ("_", $query_prefix, $organism_name, $taxon, $bg_model));
      } else {
	$outfile{prefix} = join( "/", "footprints", $taxon, $organism_name, $query_prefix,
				 join ("_", $query_prefix, $organism_name, $taxon));
      }
      &RSAT::message::Info("Automatic definition of the output prefix", $outfile{prefix}) if ($main::verbose >= 3);
    } else {
      &RSAT::error::FatalError("You must define a prefix for the output files with the option -o");
    }
  }
  return ($outfile{prefix});
}


################################################################
## Initialize output directory + output files
sub InitOutput {

  ## Create output directory if required
  $dir{output} = `dirname $outfile{prefix}`;
  chomp($dir{output});
  &RSAT::util::CheckOutDir($dir{output});

  ## Open output stream for the log file
  $outfile{log} = $outfile{prefix}."_log.txt";
  $out = &OpenOutputFile($outfile{log});

  ## File for storing the list of query gene names
  $outfile{genes} = $outfile{prefix}."_query_genes.tab";
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
  $outfile{purged} = $outfile{prefix}."_".$promoter."_seq_purged.fasta";
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

Display full help message

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

=item	B<-o output_prefix>

Prefix for the output files.

If the prefix is not specified, the program can guess a default
prefix, but this is working only if there is a single query gene or
query file.

=cut
  } elsif ($arg eq "-o") {
    $main::outfile{prefix} = shift(@arguments);

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

=cut
  } elsif ($arg eq "-taxon") {
    $main::taxon = shift(@arguments);


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
cses, it can be useful to select one or several tasks to be executed
separately. For instance, after having collected all the promoter
sequences of ortholog genes, one might desire to run the pattern
detection with various parameter values without having to retrieve the
same sequences each time.

Beware: task selection requires expertise, because most tasks depends
on the prior execution of some other tasks in the workflow. Selecting
tasks before their prerequisite tasks have been completed will provoke
fatal errors.

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

=item B<map_format>

Format for the feature map. 
Supported: all formats supported by the program feature-map

=cut
    } elsif ($arg eq "-map_format") {
      $main::map_format = shift(@arguments);


    ## Create HTML Index
=pod

=item B<-index>

Generate an HTML index with links to the result files. This option is
used for the web interface, but can also be convenient to index
results, especially when several genes or taxa are analyzed (options
-genes, -all_genes, -all_taxa).

=cut
  } elsif ($arg eq "-index") {
    $main::create_index = 1;
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
  $outfile{index} = $outfile{prefix}."_index.html";
  $index_list{$query_prefix} = $outfile{index};
  $index = &OpenOutputFile($outfile{index});
  print $index "<html>\n";
  print $index "<head><title>", join (" ", $query_prefix, $taxon, $organism_name, $bg_model), "</title></head>\n";
  print $index "<body>\n";
  print $index "<hr size=4 color='#000088'>";
  print $index "<h1 align=center>Footprint discovery result</h1>";
  print $index "<h2 align=center>", join (" ", $query_prefix, $taxon, "<i>".$organism_name."</i>", $bg_model), "</h2>\n";
  print $index "<hr size=2 color='#000088'>";
  print $index "<blockquote>";
  print $index "<table cellspacing=0 cellpadding=3 border=0>\n";
  &IndexOneFile("log", $outfile{log});
  &IndexOneFile("input", $infile{genes}) if (($infile{genes}) && !($sep_genes));
}

################################################################
## Predict operon leader genes of the query genes
sub InferQueryOperons {
  &RSAT::message::TimeWarn("Get leaders of query genes (d<=".$dist_thr."bp)", $outfile{leader_qgenes}) if ($main::verbose >= 1);
  &CheckDependency("operons", "genes");
  my $cmd = "$SCRIPTS/get-leader-multigenome ";
  $cmd .= " -i ".$outfile{genes};
  $cmd .= " -o ".$outfile{leader_qgenes};
  $cmd .= " -uth interg_dist ".$dist_thr;
  &one_command($cmd) if ($task{operons});
  #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
  &IndexOneFile("leader query genes", $outfile{leader_qgenes}) if ($create_index);
}


################################################################
## Retrieve promoters of the query organism
sub RetrieveQueryPromoters {
  &RSAT::message::TimeWarn("Retrieving promoter sequences for query genes", $outfile{query_seq}) if ($main::verbose >= 1);
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
  $cmd .= " -o ".$outfile{query_seq};
  &one_command($cmd) if ($task{query_seq});
  #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
  &IndexOneFile("query sequence", $outfile{query_seq}) if ($create_index);
}

################################################################
## Detect all dyads in promoters of query genes for dyad filtering
sub ComputeFilterDyads {
  &RSAT::message::TimeWarn("Computing filter dyads", $outfile{filter_dyads}) if ($main::verbose >= 1);
  &CheckDependency("filter", "query_seq");
  my $cmd = "$SCRIPTS/dyad-analysis -v 1 -return occ -lth occ 1";
  $cmd .= " -i ".$outfile{query_seq};
  $cmd .= " -l 3 -sp 0-20";
  $cmd .= " ".$strands;
  $cmd .= " ".$noov;
  $cmd .= " -o ".$outfile{filter_dyads};
  &one_command($cmd) if ($task{filter_dyads});
  #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
  &IndexOneFile("filter dyads", $outfile{filter_dyads}) if ($create_index);
}

################################################################
## Identify ortholog genes
sub GetOrthologs {
  &IndexOneFile("orthologs", $outfile{orthologs}) if ($create_index);
  return(0) unless ($task{orthologs});
  &RSAT::message::TimeWarn("Getting orthologs", $outfile{orthologs}) if ($main::verbose >= 1);
  &CheckDependency("orthologs", "genes");
  my $cmd = "$SCRIPTS/get-orthologs";
  $cmd .= " -i ".$outfile{genes};
  $cmd .= " -org ".$organism_name;
  $cmd .= " -taxon ".$taxon;
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

################################################################
## Predict operon leader genes for the orthologous genes
sub InferOrthoOperons {
  &IndexOneFile("leader genes", $outfile{bbh}) if ($create_index);
  return(0) unless ($task{operons});
  &RSAT::message::TimeWarn("Get leaders of query genes (d<=".$dist_thr."bp)", $outfile{bbh}) if ($main::verbose >= 1);
  &CheckDependency("operons", "orthologs");
  my $cmd = "$SCRIPTS/get-leader-multigenome ";
  $cmd .= " -i ".$outfile{orthologs};
  $cmd .= " -o ".$outfile{bbh};
  $cmd .= " -uth interg_dist ".$dist_thr;
  &one_command($cmd);
  #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
}


################################################################
## Retrieve sequences from orthologs
sub RetrieveOrthoSeq {
  &IndexOneFile("$promoter sequences", $outfile{seq}) if ($create_index);
  return(0) unless ($task{ortho_seq});
  &RSAT::message::TimeWarn("Retrieving promoter sequences of orthologs", $outfile{seq}) if ($main::verbose >= 1);
  &CheckDependency("ortho_seq", "bbh");
  my $cmd = "$SCRIPTS/retrieve-seq-multigenome -ids_only";
  $cmd .= " -i ".$outfile{bbh};
  $cmd .= " -noorf";
  $cmd .= " -o ".$outfile{seq};
  &one_command($cmd);
  #  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
}



################################################################
## Purge sequences
sub PurgeOrthoSeq {
  &IndexOneFile("purged sequences", $outfile{purged}) if ($create_index);
  return(0) unless ($task{purge});
  &RSAT::message::TimeWarn("Purging promoter sequences of orthologs", $outfile{purged}) if ($main::verbose >= 1);
  &CheckDependency("purge", "seq");
  my $cmd = "$SCRIPTS/purge-sequence";
  $cmd .= " -i ".$outfile{seq};
  $cmd .= " -ml 30 -mis 0 -mask_short 30";
  $cmd .= " ".$strands;
  $cmd .= " -o ".$outfile{purged};
  &one_command($cmd);
#  print $out "\n; ", &AlphaDate(), "\n", $cmd, "\n\n"; &doit($cmd, $dry, $die_on_error, $main::verbose, $batch, $job_prefix);
}



################################################################
## Detect over-representation of matching occurrecnes (hits) of the motif
sub OccurrenceSig {
  &IndexOneFile("occ sig", $outfile{occ_sig}) if ($create_index);
  return(0) unless ($task{occ_sig});
  &RSAT::message::TimeWarn("Testing over-representation of hits", $outfile{occ_sig}) if ($main::verbose >= 1);
  &CheckDependency("occ_sig", "purged");
  my $cmd = "matrix-scan";
  $cmd .= $ms_parameters;
  $cmd .= " -return distrib,occ_proba,rank -sort_distrib";
  #  $cmd .= " -lth inv_cum 1 -lth occ_sig 0 -uth occ_sig_rank 1";
  $cmd .= " -i ".$outfile{purged};
  $cmd .= " -o ".$outfile{occ_sig};
  $cmd .= $occ_sig_opt;
  &one_command($cmd);
}

################################################################
## Draw a graph with occurrence significance profiles
sub OccurrenceSigGraph {
  &IndexOneFile("occ freq graph", $outfile{occ_freq_graph}) if ($create_index);
  &IndexOneFile("occ sig graph", $outfile{occ_sig_graph}) if ($create_index);

  return(0) unless ($task{occ_sig_graph});

  &RSAT::message::TimeWarn("Graph with over-representation of hits", $outfile{occ_freq_graph}) if ($main::verbose >= 1);
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
  my $cmd = "sort -n -k 2 $outfile{occ_sig} | XYgraph";
  $cmd .= " -xcol 2 -xleg1 'Weight score' -xsize 800  -xgstep1 1 -xgstep2 0.5";
  $cmd .= " -ycol 11 -yleg1 'Binomial significance of hits'";
  $cmd .= " -title 'matrix ".$matrix_suffix." ; gene ".$current_gene."'";
  $cmd .= " -lines -legend ";
  $cmd .= " -format ".$plot_format;
  $cmd .= " -o ".$outfile{occ_sig_graph};
  $cmd .= " ".$occ_sig_graph_opt;
  &one_command($cmd);
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
    $occ_sig_file{$current_gene} = $outfile{occ_sig} if (-e $outfile{occ_sig});
    $occ_freq_graph_file{$current_gene} = $outfile{occ_freq_graph} if (-e $outfile{occ_freq_graph});
    $occ_sig_graph_file{$current_gene} = $outfile{occ_sig_graph} if (-e $outfile{occ_sig_graph});
  }

  ## Index scan files for the synthetic table
  $map_file{$current_gene} = $outfile{map} if (-e $outfile{map});
}


################################################################
## Scan sequences with PSSM to locate sites
sub OrthoScan {
  &IndexOneFile("sites", $outfile{sites}) if ($create_index);
  return(0) unless ($task{scan});
  &RSAT::message::TimeWarn("Scannig sequences to detect sites", $outfile{sites}) if ($main::verbose >= 1);
  &CheckDependency("scan", "seq");
  my $scan_cmd = "matrix-scan";
  $scan_cmd .= $ms_parameters;
  $scan_cmd .= " -return limits,sites,rank";
  $scan_cmd .= " -i ".$outfile{seq};
  $scan_cmd .= " -o ".$outfile{sites};
  $scan_cmd .= " ".$scan_opt;
  &one_command($scan_cmd);
}

################################################################
## Draw a feature map with the detected sites
sub OrthoMap {
  &IndexOneFile("map", $outfile{map}) if ($create_index);
  return(0) unless ($task{map});
  &RSAT::message::TimeWarn("Drawing feature map", $outfile{map}) if ($main::verbose >= 1);
  &CheckDependency("map", "sites");
  my $cmd = "feature-map -i ".$outfile{sites};
  $cmd .= " -scalebar -legend";
  $cmd .= " -xsize 800 -scorethick";
  $cmd .= " -title 'matrix hits in ".$taxon." promoters'";
  $cmd .= " -format ".$map_format;
  $cmd .= " -o ".$outfile{map};
  &one_command($cmd);
}


1;

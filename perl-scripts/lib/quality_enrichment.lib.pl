################################################################
## Library for functions shared by motif enrichment detection programs
## - matrix-quality
## - matrix-enrichment

require "RSA2.cgi.lib";		## For sortable HTML tables
use List::MoreUtils qw(uniq);

################################################################
## Export the matrix in tab-delimited format. This will be used
## for permuting the matrix.
sub ExportTabMatrix {
  my ($matrix) = @_;

  &RSAT::message::TimeWarn("Exporting matrix in tab-delimited format",  $main::outfile{matrix_tab})
    if ($main::verbose >= 2);

  my $verbose_bk = $verbose;
  $verbose = 0;
  $matrix_handle = &OpenOutputFile($main::outfile{matrix_tab});
  print $matrix_handle $matrix->toString(sep=>"\t",
					 type=>"counts",
					 format=>"tab",
					 pipe=>"", ## We suppress the pipe for permute-table
					);
  close $matrix_handle;
  $verbose = $verbose_bk;
}




################################################################
## Export the matrix in transfac format. This will be used
## for permuting the matrix.
sub ExportTransfacMatrix {
  my ($matrix) = @_;

  &RSAT::message::TimeWarn("Exporting matrix in transfac format",  $outfile{matrix_transfac})
    if ($main::verbose >= 2);

  my $verbose_bk = $verbose;
  $verbose = 0;
  $matrix_handle = &OpenOutputFile($main::outfile{matrix_transfac});
  print $matrix_handle $matrix->toString(type=>"counts",
					 format=>"transfac",
					);
  close $matrix_handle;
  $verbose = $verbose_bk;
}


################################################################
## Export the matrix in tab-delimited format with additional
## information + the logos.
sub ExportMatrixInfo {
  my ($matrix) = @_;

  ## Compute information (logos, consensus, parameters)
  &RSAT::message::TimeWarn("Exporting matrix information",  $outfile{matrix_info})
    if ($main::verbose >= 2);
  my $cmd = $SCRIPTS."/convert-matrix -v 1";
  $cmd .= " -from transfac -i ".$main::outfile{matrix_transfac};
  $cmd .= " -to tab -o ".$outfile{matrix_info};
  $cmd .= " -bgfile ".$outfile{bg_file_inclusive};
  $cmd .= " -bg_format inclusive";
  $cmd .= " -return counts,frequencies,weights,info,parameters,sites,logo";
  $cmd .= " -logo_format ".$logo_formats;
  $cmd .= " -logo_opt '-e -M -t ".$matrix_name." ' ";
  $cmd .= " -logo_file ". $outfile{matrix_logo};
  &doit($cmd, $dry, $die_on_error, $verbose, 0, $job_prefix);

  
}


################################################################
## Export the matrix in tab-delimited format. This will be used
## for permuting the matrix.
sub PermuteMatrixColumns {
  ## Define the permutation number for each sequence type and the max permutation number
    foreach my $seq_type (@local_seq_types) {
    unless (defined($perm_nb{$seq_type})) {
      $perm_nb{$seq_type} = 0;
      }
  }
  local $perm_nb_max = &RSAT::stats::checked_max(0, values %perm_nb);

  return if ($perm_nb_max == 0);

  &RSAT::message::TimeWarn("Permuting matrix columns", $perm_nb_max, "permutations")
    if ($main::verbose >= 2);


  ## Define the names of the column-permuted matrices (required for the index)
  for my $i (1..$perm_nb_max) {
    $outfile{'matrix_perm_col_'.$i} = $matrix_prefix{$matrix_name}."_matrix_perm_col_".$i.".tab";
  }

  ## Define file names for sequence type-specific permuted matrices
  ## (each sequence type can have its particular number of
  ## permutations)
  print $main::out "; Sequence sets (name, permutations, file)";
  foreach my $seq_type (@main::local_seq_types) {
  &RSAT::message::Debug("Defining file names for column-permuted matrices",
			  "seq_type=".$seq_type,
			  "perm_nb=".$perm_nb{$seq_type},
			 ) if ($main::verbose >= 10);

    #print $main::out join("\t", ";", $seq_type, $perm_nb{$seq_type}, $seqfile{$seq_type}), "\n";
    $outfile{'perm_col_matrices_'.$seq_type.'_'.$perm_nb{$seq_type}.'perm'} = $matrix_prefix{$matrix_name}."_".$seq_type."_matrix_perm_col_all_".$perm_nb{$seq_type}.".tab";
    push @files_to_index, 'perm_col_matrices_'.$seq_type.'_'.$perm_nb{$seq_type}.'perm' if ($perm_nb{$seq_type} > 0);
  }

  if ($task{permute}) {

    ## Remove previous version of the column-permuted matrix files
    ## before appending the new permuted columns
      foreach my $seq_type (@main::local_seq_types) {
	  
      #    $outfile{'perm_col_matrices_'.$seq_type.'_'.$perm_nb{$seq_type}.'perm'} = $matrix_prefix{$matrix_name}."_".$seq_type."_matrix_perm_col_all_".$perm_nb{$seq_type}.".tab";
      $init_matrix_cmd = "rm -f ".$outfile{'perm_col_matrices_'.$seq_type.'_'.$perm_nb{$seq_type}.'perm'};
      &doit($init_matrix_cmd, $dry, $die_on_error, $verbose, 0, $job_prefix);
    }

    ## Generate the column-permuted matrices
    for my $i (1..$perm_nb_max) {

      ## Perform one column permutation
#      $outfile{'matrix_perm_col_'.$i} = $matrix_prefix{$matrix_name}."_matrix_perm_col_".$i.".tab";
      my $permute_matrix_cmd = $SCRIPTS."/permute-table -rownames -entire_col";
      $permute_matrix_cmd .= " -i ".$outfile{matrix_tab};
      $permute_matrix_cmd .= " -o ".$outfile{'matrix_perm_col_'.$i};

      ## Append the column-permuted matrix to the permuted matrices
      ## for each sequence set (the number of required permutation can
      ## vary between sequence sets)
      foreach my $seq_type (sort keys %seqfile) {
	  if (defined($perm_nb{$seq_type}) && ($i <= $perm_nb{$seq_type})) {
	      #print "seq_type ".$seq_type ." +++ perm_nb ".$perm_nb{$seq_type}." \n" ;
	      #print "+ ".$main::outfile{'perm_col_matrices_'.$seq_type.'_'.$perm_nb{$seq_type}.'perm'}." +\n";
	      #die "boom";
	      $permute_matrix_cmd .= "; cat ".$outfile{'matrix_perm_col_'.$i}." >> ".$outfile{'perm_col_matrices_'.$seq_type.'_'.$perm_nb{$seq_type}.'perm'};
	      
	}
      }

      ## Execute the command
      &doit($permute_matrix_cmd, $dry, $die_on_error, $verbose, 0, $job_prefix);
    }
  }
}



################################################################
## Calculate score distribution
sub CalcTheorScoreDistribution {
    my ($matrix_tab_file,  $out_file) = @_;
    
    &RSAT::message::TimeWarn("Calculating theoretical distribution for matrix", $matrix_tab_file)
	if ($main::verbose >= 2);
    
    my $matrix_distrib_cmd = $SCRIPTS."/matrix-distrib";
    $matrix_distrib_cmd .= " -v 1";
    $matrix_distrib_cmd .= " -m ".$matrix_tab_file;
    $matrix_distrib_cmd .= " -matrix_format tab";
    $matrix_distrib_cmd .= " -pseudo ".$main::pseudo_counts;
    $matrix_distrib_cmd .= " -bgfile ".$outfile{bg_file_inclusive};
    $matrix_distrib_cmd .= " -bg_format inclusive";
    $matrix_distrib_cmd .= " -bg_pseudo ".$main::bg_pseudo if (defined($main::bg_pseudo));
    $matrix_distrib_cmd .= " -decimals ".$decimals;
    $matrix_distrib_cmd .= " -o ".$out_file;
    
    ## Execute the command
    &RSAT::message::TimeWarn("Matrix-distrib command: ", $matrix_distrib_cmd)
	if ($main::verbose >= 2);
    &doit($matrix_distrib_cmd, $dry, $die_on_error, $verbose, $batch, $job_prefix);
}


################
## Draw the NWD comparison plot : by motif or by sequence

sub Draw_NWD{
    my ($prefix,@nwd_files)= @_;
    
    
    my $ycols = "";
    if ((scalar(@nwd_files)>=2)){
	$ycols = join ",", 2..(scalar(@nwd_files)+1);
    }else{
	$ycols = 2;
    }
    
    my $nwd_input_files =" -i ". join (" -i " ,@nwd_files)." ";
    my $ic_column= 1;
    my $sc_column= 4;
    my $nwd_outfile_prefix= $prefix ;
    my $nwd_compare_scores_file=$nwd_outfile_prefix."_compare-scores.tab";
    
    my $compare_nwd_cmd =  $SCRIPTS."/compare-scores " ;
    $compare_nwd_cmd .= " ". $nwd_input_files ." " ;
    $compare_nwd_cmd .= " -ic ". $ic_column . " " ;
    $compare_nwd_cmd .= " -sc ". $sc_column . " ";
    $compare_nwd_cmd .= " -numeric  -basename ";
    $compare_nwd_cmd .= " -o ".  $nwd_compare_scores_file . " ";
    &doit($compare_nwd_cmd, $dry, $die_on_error, $verbose, 0, $job_prefix);
    &RSAT::message::Info("Merging NWD files with compare-scores ", $compare_nwd_cmd ) if ($main::verbose >= 5);
    
    
    my $nwd_xygrpah_file=$nwd_outfile_prefix."_compare-scores.";
    
    foreach my $image_format (@image_formats) {
	my $XYgraph_nwd_cmd = $SCRIPTS."/XYgraph ".$main::rplot_option." " ;
	$XYgraph_nwd_cmd .= " -i ".  $nwd_compare_scores_file ." ";
	$XYgraph_nwd_cmd .= " -format ". $image_format  ." " ;
	$XYgraph_nwd_cmd .= " -xcol 1  -ycol ".$ycols;
	$XYgraph_nwd_cmd .= " -lines -xlog " ;
	$XYgraph_nwd_cmd .= " -yleg1 'NWD' -xleg1 'log10(Pvalue)' ";
	$XYgraph_nwd_cmd .= " -legend -pointsize 0 ";
	$XYgraph_nwd_cmd .= " -o ". $nwd_xygrpah_file.$image_format." ";
	&doit($XYgraph_nwd_cmd, $dry, $die_on_error, $verbose, 0, $job_prefix);
	&RSAT::message::Info(" Drawing NWD  curves ",  $XYgraph_nwd_cmd) if ($main::verbose >= 5);   
	
    }
    return($nwd_compare_scores_file, $nwd_xygrpah_file )
	
	
	
}



################
## Subroutine to draw occurrence significance plot comparing distribitions across files: by motif or by sequence
sub Draw_OCC{
    my ($prefix,@occ_files)= @_;
    
    my $occ_input_files =" -i ". join (" -i " ,@occ_files)." ";
    my $ic_column= 7;
    my $sc_column= 11;
    my $occ_outfile_prefix= $prefix ;
    my $occ_compare_scores_file=$occ_outfile_prefix."_compare-scores.tab";
    
    my $compare_occ_cmd =  $SCRIPTS."/compare-scores " ;
    $compare_occ_cmd .= " ". $occ_input_files ." " ;
    $compare_occ_cmd .= " -ic ". $ic_column . " " ;
    $compare_occ_cmd .= " -sc ". $sc_column . " ";
    $compare_occ_cmd .= " -numeric  -basename ";
    $compare_occ_cmd .= " -o ".  $occ_compare_scores_file . " ";
    &doit($compare_occ_cmd, $dry, $die_on_error, $verbose, 0, $job_prefix);
    &RSAT::message::Info("Merging OCC files with compare-scores ", $compare_occ_cmd ) if ($main::verbose >= 10);
    
    
    my $occ_xygrpah_file=$occ_outfile_prefix."_compare-scores.";
    
    my $ycols = "";
    if ((scalar(@occ_files)>=2)){
	$ycols = join ",", 2..(scalar(@occ_files)+1);
    }else{
	$ycols = 2;
    }   
    

    #print join ("++",@main::image_formats );
    #die "BOOM";
    
    
    foreach my $image_format (@main::image_formats) {
	my $XYgraph_occ_cmd = $SCRIPTS."/XYgraph ".$main::rplot_option." " ;
	$XYgraph_occ_cmd .= " -i ".  $occ_compare_scores_file ." ";
	$XYgraph_occ_cmd .= " -format ". $image_format  ." " ;
	$XYgraph_occ_cmd .= " -xcol 1  -ycol ".$ycols;
	$XYgraph_occ_cmd .= " -lines -xlog 10 " ;
	$XYgraph_occ_cmd .= " -yleg1 'Binomial significance of hit number (OCC)' -xleg1 'log10(Score Pvalue)' ";
	$XYgraph_occ_cmd .= " -legend -pointsize 0 ";
	$XYgraph_occ_cmd .= " -hline red 100 -hline violet 0 -ysize 400 -force_lines "; ## Parammeters suggested by Morgane      
	$XYgraph_occ_cmd .= " -o ". $occ_xygrpah_file.$image_format." ";
	&doit($XYgraph_occ_cmd, $dry, $die_on_error, $verbose, 0, $job_prefix);
	&RSAT::message::Info(" Drawing OCC  curves ",  $XYgraph_occ_cmd) if ($main::verbose >= 5);   
    }
    return($occ_compare_scores_file, $occ_xygrpah_file )
	
	
	
}


################################################################
## Compute the score distribution in one sequence set
sub CalcSequenceDistrib {
  ## Arguments are local, because they are needed in sub-routines
  local ($sequence_file, $matrix_file, $matrix_format, $seq_type, $index, @args) = @_;

  ## Only tab format is supported
  if ($matrix_format ne "tab") {
    &RSAT::error::FatalError("&CaclSequenceDistrib() only supports tab format in quick scan mode", $seq_type, $matrix_format, $matrix_file);
  }

  ## Define the output file for the current sequence type
  $outfile{'empirical_distrib_'.$seq_type} = $matrix_prefix{$matrix_name}."_scan_".$seq_type."_score_distrib.tab";

  ## Add the file to the list for the comparison of distributions
  if ($index) {
    push @files_to_index, 'empirical_distrib_'.$seq_type;
    &RSAT::message::Debug("Adding file to index ",  'empirical_distrib_'.$seq_type, $outfile{'empirical_distrib_'.$seq_type}  ) if ($main::verbose >= 10);
    push @distrib_files, $outfile{'empirical_distrib_'.$seq_type}; $file_nb{$seq_type} = scalar(@distrib_files);
  }


  local $matrix_scan_cmd = "";
  if (($quick) &&
      !($scanopt{$seq_type}) ## Scanning options may be incompatible with matrix-scan-quick -> if specified, we pass the command to matrix-scan
     ) {
    $matrix_scan_cmd = $quick_scan_cmd;
  } else {
    $matrix_scan_cmd = $SCRIPTS."/matrix-scan -v ".$main::verbose;
    #    $matrix_scan_cmd .= " -quick"; ## Run in quick mode if possible
    #    $matrix_scan_cmd .= " -m ".$matrix_file;
    #    $matrix_scan_cmd .= " -top_matrices 1";
    #    $matrix_scan_cmd .= " -matrix_format ".$matrix_format;
    $matrix_scan_cmd .= " -matrix_format tab"; ## We use tab as matrix format for compatibiliy with matrix-scan-quick
    $matrix_scan_cmd .= " -bg_format inclusive"; ## We use inclusive as bg format for compatibiliy with matrix-scan-quick
  }
  $matrix_scan_cmd .= " -i ".$sequence_file;
  $matrix_scan_cmd .= " -m ".$matrix_file;
  $matrix_scan_cmd .= " -pseudo ".$main::pseudo_counts;
  $matrix_scan_cmd .= " -decimals ".$decimals;
  $matrix_scan_cmd .= " -bgfile ".$outfile{bg_file_inclusive};
  $matrix_scan_cmd .= join(" ", "", @args);

  ## Sequence type-Specific options
  &RSAT::message::TimeWarn("\tScanning options for ".$seq_type,  $scanopt{$seq_type})
    if ((defined($main::scanopt{$seq_type})) && ($main::verbose >= 2));
  if (defined($main::scanopt{$seq_type})) {
    $matrix_scan_cmd .= " ".$main::scanopt{$seq_type};
  }

  if ($scanopt{$seq_type}) {
    ## Scanning options may be ignored by the option -return distrib
    ## -> if specified, we detect sites and use classfeq do
    ## determine the distirbution of weight scores
    &AddSequenceDistribOptions_classfreq();
  } else {
    &AddSequenceDistribOptions_direct();
  }

  $matrix_scan_cmd .= " > ".$outfile{'empirical_distrib_'.$seq_type};

  &RSAT::message::Info("Scanning to compute distribution", $matrix_scan_cmd) if ($main::verbose >= 2);

  ## Print the complete command in the log file
  print $main::out "\n; ", &AlphaDate(), "\tComputing score distribution\n";
  printf $main::out ";\t%-22s\t%s\n", "Sequence type", $seq_type;
  printf $main::out ";\t%-22s\t%s\n", "Sequence file", $sequence_file;
  if (defined($main::scanopt{$seq_type})) {
    printf $main::out ";\t%-22s\t%s\n", "Type-specific options", $scanopt{$seq_type};
  }
  printf $main::out "; %s\n%s\n", "Command:", $matrix_scan_cmd;
  print $main::out "\n";

  ## Execute the command
  if ($task{scan}) {
    &RSAT::message::TimeWarn("Computing observed distribution",
			     "\nseq_type=".$seq_type,
			     "\nmatrix_file=".$matrix_file,
			     "\nseq_file=".$sequence_file,
			     "\nout_file=".$outfile{'empirical_distrib_'.$seq_type},
			    )
      if ($main::verbose >= 2);

    &doit($matrix_scan_cmd, $dry, $die_on_error, $verbose, $batch, $job_prefix);
  }
  return($outfile{'empirical_distrib_'.$seq_type});
}



################################################################
## Options to compute the empirical distribution using matrix-scan
## -return distrib (direct computation).
sub AddSequenceDistribOptions_direct {
  $matrix_scan_cmd .= " -return distrib ";
}


################################################################
## Compare the score distribution files
sub CompareDistrib {
  my ($score_column, @distrib_files) = @_;

  $outfile{distrib_compa} = $matrix_prefix{$matrix_name}."_score_distrib_compa";

  if ($task{compare}) {
    &RSAT::message::TimeWarn("Comparing score distributions",  $outfile{distrib_compa})
      if ($main::verbose >= 5);
    &RSAT::message::Info("\n", "distrib_files", @distrib_files)
      if ($main::verbose >= 2);

    ################################################################
    ## Compare the distributions
    my $distrib_compa_cmd = $SCRIPTS."/compare-scores ";
    $distrib_compa_cmd .= " -numeric";
    $distrib_compa_cmd .= " -sc1 4"; # score column for the theoretical distribution
    $distrib_compa_cmd .= " -sc ".$score_column; # score column for the observed distributions
    $distrib_compa_cmd .= " -suppress ".$matrix_prefix{$matrix_name}."_scan_";
    $distrib_compa_cmd .= " -suppress ".$matrix_prefix{$matrix_name}."_";
    $distrib_compa_cmd .= " -suppress _score_distrib.tab ";
    #$distrib_compa_cmd .= " -suppress ".$dir{output}." ";
    #$distrib_compa_cmd .= " -suppress ".$matrix_name." ";
    $distrib_compa_cmd .= " -o ".$outfile{distrib_compa}.".tab";
    $distrib_compa_cmd .= " -files ";
    $distrib_compa_cmd .= join(" ", $outfile{'matrix_theoretical_distrib'}, @distrib_files);
    &doit($distrib_compa_cmd, $dry, $die_on_error, $verbose, $batch, $job_prefix);
    &RSAT::message::TimeWarn("Comparing distribution scores from different files ", $distrib_compa_cmd ) if ($main::verbose >= 5);
    #<STDIN>;
    
  }
  if ($task{graphs}) {
      &RSAT::message::TimeWarn("Generating comparison graphs") if ($main::verbose >= 2);

    ## Generate the graphs for each image format
    foreach my $image_format (@image_formats) {

      ## General options for all the graphs below
      my $all_graph_options = " -i ".$outfile{distrib_compa}.".tab";
      $all_graph_options .= " -format ".$image_format." -lines -pointsize 0";
      $all_graph_options .= " ".$graph_options;

      ## Alternative options for the large graphs and for the icons, respectively
      my $large_graph_options = " -title1 '".$matrix_name."'";
#      $large_graph_options .= " -title2 ".$matrix_prefix{$matrix_name};
      $large_graph_options .= " -legend ";
      $large_graph_options .= " -xsize 800 -ysize 400 ";
      $large_graph_options .= " -xleg1 'matrix score' ";
      $large_graph_options .= " -yleg1 'dCDF (log scale)' ";

      my $icon_options;


      ################################################################
      ## Draw a graph with all the decreasing cumulative distributions
      my $XYgraph_cmd = $SCRIPTS."/XYgraph ".$main::rplot_option." ".$all_graph_options;

      my $ycols = join ",", 2..(scalar(@distrib_files)+2);
      
      $XYgraph_cmd .= " -xcol 1 -ycol ".$ycols;
      $XYgraph_cmd .= " -ymin 0  -ymax 1 ";
      $XYgraph_cmd .= " -xgstep1 5 -xgstep2 1 -ygstep1 0.1 -ygstep2 0.02";
      $XYgraph_cmd .= " -gp 'set size ratio 0.5' ";
      $graph_file_opt = $large_graph_options." ".$distrib_options." -o ".$outfile{distrib_compa}.".".$image_format;
      &doit($XYgraph_cmd.$graph_file_opt, $dry, $die_on_error, $verbose, $batch, $job_prefix);
      &RSAT::message::Info("Distribution comparison graph", $outfile{distrib_compa}.".".$image_format) if ($main::verbose >= 2);
      print $main::out ";\n; XYgraph command\n", $XYgraph_cmd.$graph_file_opt, "\n";

      ## Generate the icon
      unless ($noicon) {
	$icon_options = " -xsize 120 -ysize 120 -o ".$outfile{distrib_compa}."_small.".$image_format;
	&doit($XYgraph_cmd.$icon_options, $dry, $die_on_error, $verbose, $batch, $job_prefix);
      }

      ################################################################
      ## Draw a graph with all the decreasing cumulative distributions
      ## and a logarithmic Y axis
      $XYgraph_cmd = $SCRIPTS."/XYgraph ".$main::rplot_option." ".$all_graph_options;
      $XYgraph_cmd .= " -xcol 1 -ycol ".$ycols;
      $XYgraph_cmd .= " -xgstep1 5 -xgstep2 1";
      $XYgraph_cmd .= " -ymax 1 -ylog 10";
      $XYgraph_cmd .= " -gp 'set size ratio 0.5' ";
      $graph_file_opt = $large_graph_options." ".$distrib_options." -o ".$outfile{distrib_compa}."_logy.".$image_format;
      &doit($XYgraph_cmd.$graph_file_opt, $dry, $die_on_error, $verbose, $batch, $job_prefix);
      &RSAT::message::Info("Distribution comparison graph (log Y)", $outfile{distrib_compa}."_logy.".$image_format)
	if ($main::verbose >= 2);
      print $main::out ";\n; XYgraph command\n", $XYgraph_cmd.$graph_file_opt, "\n";

      ## Generate the icon
      unless ($noicon) {
	$icon_options = " -xsize 120 -ysize 120 -o ".$outfile{distrib_compa}."_logy_small.".$image_format;
	&doit($XYgraph_cmd.$icon_options, $dry, $die_on_error, $verbose, $batch, $job_prefix);
	&RSAT::message::Info("Distribution comparison icon (log Y)", $outfile{distrib_compa}."_logy_small.".$image_format)
	  if ($main::verbose >= 2);
      }

      ################################################################
      ## Draw a ROC curve
      if ($draw_roc){
	  my $ref_column = 2;
	  if ($roc_ref) {
	      if (defined($file_nb{$roc_ref})) {
		  $ref_column = 2 + $file_nb{$roc_ref};
	      } else {
		  if ($roc_ref ne "theor") {
		      &RSAT::message::Warning($roc_ref, "Invalid reference distribution for the ROC curve: should be one of the input sequence types, or 'theor'.");
		      $roc_ref = "Forced to use theoretical";
		  }
	      }
	  }
	  
	  $ycols = join ",", 2..(scalar(@distrib_files)+2);
	  #      $large_graph_options =~ s/-xsize 800/-xsize 400/;
	  $XYgraph_cmd = $SCRIPTS."/XYgraph ".$main::rplot_option." ".$all_graph_options;
	  $XYgraph_cmd .= " -xcol ".$ref_column;
	  $XYgraph_cmd .= " -ycol ".$ycols;
	  $XYgraph_cmd .= " -ygstep1 0.1 -ygstep2 0.02";
	  # $XYgraph_cmd .= " -ymin 0  -ymax 1 ";
	  # $XYgraph_cmd .= " -xmin 0  -xmax 1 ";
	  $XYgraph_cmd .= " -ymax 1 ";
	  $XYgraph_cmd .= " -xmax 1 ";
	  my $roc_file_opt = $large_graph_options.$roc_options." -o ".$outfile{distrib_compa}."_roc.".$image_format;
	  $roc_file_opt .= " -xleg1 'FPR (Reference = ".$roc_ref.")' ";
	  $roc_file_opt .= " -yleg1 'Site Sn + other distributions' ";
	  
	  ################################################################
	  ## Draw a ROC curve with non-logarithmic axes
	  ## Beware: this curve is generally not informative, so I inactivate this drawing.
	  ## In case it would appear useful for some purpose, I would add an option "-ROC_nolog"
	  my $ROC_nolog = 0;
	  if ($ROC_nolog) {
	      &doit($XYgraph_cmd.$roc_file_opt, $dry, $die_on_error, $verbose, $batch, $job_prefix);
	      &RSAT::message::Info("ROC curve graph", $outfile{distrib_compa}."_roc.".$image_format) if ($main::verbose >= 2);
	      print $main::out ";\n; XYgraph command\n", $XYgraph_cmd.$roc_file_opt, "\n";
	      
	      ## Generate the icon for the ROC curve
	      unless ($noicon) {
		  $icon_options = " -xsize 120 -ysize 120 -o ".$outfile{distrib_compa}."_roc_small.".$image_format;
		  &doit($XYgraph_cmd.$icon_options, $dry, $die_on_error, $verbose, $batch, $job_prefix);
		  &RSAT::message::Info("ROC curve icon", $outfile{distrib_compa}."_roc_small.".$image_format) if ($main::verbose >= 2);
	      }
	  }
	  
	  
	  ################################################################
	  ## Draw a ROC curve with xlog This is the relevant way to
	  ## display the ROC curve with pattern matching, because we are
	  ## only interested in the low FPR values (< 10-3), which are not
	  ## visible on the non-log representations.
	  $XYgraph_cmd =~ s/XYgraph/XYgraph -xlog 10/;
	  $roc_file_opt =~ s/_roc/_roc_xlog/;
	  
	  &doit($XYgraph_cmd.$roc_file_opt, $dry, $die_on_error, $verbose, $batch, $job_prefix);
	  &RSAT::message::Info("ROC curve graph (log X)", $outfile{distrib_compa}."_roc_xlog.".$image_format) if ($main::verbose >= 2);
	  print $main::out ";\n; XYgraph command\n", $XYgraph_cmd.$roc_file_opt, "\n";
	  
	  ## Generate the icon for the ROC curve
	  unless ($noicon) {
	      $icon_options = " -xsize 120 -ysize 120 -o ".$outfile{distrib_compa}."_roc_xlog_small.".$image_format;
	      &doit($XYgraph_cmd.$icon_options, $dry, $die_on_error, $verbose, $batch, $job_prefix);
	      &RSAT::message::Info("ROC curve icon (log X)", $outfile{distrib_compa}."_roc_xlog_small.".$image_format) if ($main::verbose >= 2);
	  }
	  
	  ################################################################
	  ## Draw a ROC curve with xylog
	  my $ROC_xylog = 0;
	  if ($ROC_xylog) {
	      $XYgraph_cmd =~ s/XYgraph/XYgraph -ylog 10/;
	      $roc_file_opt =~ s/_roc_xlog/_roc_xylog/;
	      &doit($XYgraph_cmd.$roc_file_opt, $dry, $die_on_error, $verbose, $batch, $job_prefix);
	      &RSAT::message::Info("ROC curve graph (log XY)", $outfile{distrib_compa}."_roc_xylog.".$image_format) if ($main::verbose >= 2);
	      print $main::out ";\n; XYgraph command\n", $XYgraph_cmd.$roc_file_opt, "\n";
	      
	      ## Generate the icon for the ROC curve
	      unless ($noicon) {
		  $icon_options = " -xsize 120 -ysize 120 -o ".$outfile{distrib_compa}."_roc_xylog_small.".$image_format;
		  &doit($XYgraph_cmd.$icon_options, $dry, $die_on_error, $verbose, $batch, $job_prefix);
		  &RSAT::message::Info("ROC curve icon (log XY)", $outfile{distrib_compa}."_roc_xylog_small.".$image_format) if ($main::verbose >= 2);
	      }
	  }
      }
      
      unless ($no_cv) {
	  if ($task{theor_cv} ||$task{compara}) {
	      
	      $outfile{th_distrib_compa} = $matrix_prefix{$matrix_name}."_theoretical_score_distrib_compa";

	      ################################################################
	      ## Compare the theoretical distributions
	      my $distrib_compa_cmd = $SCRIPTS."/compare-scores ";
	      $distrib_compa_cmd .= " -numeric";
	      $distrib_compa_cmd .= " -sc 4";	# score column for the theoretical distribution
	      $distrib_compa_cmd .= " -suppress ".$matrix_prefix{$matrix_name}."_";
	      $distrib_compa_cmd .= " -suppress .tab";
	      $distrib_compa_cmd .= " -o ".$outfile{th_distrib_compa}.".tab";
	      $distrib_compa_cmd .= " -files ";
	      $distrib_compa_cmd .= join(" ", $outfile{'matrix_theoretical_distrib'}, @th_distrib_files);
	      &doit($distrib_compa_cmd, $dry, $die_on_error, $verbose, $batch, $job_prefix);
	      
	      ################################################################
	      ## draw a graph with the theoretical distributions of partial and complete matrix
	      ## General options for all the graphs below
	      $all_graph_options =~ s/$outfile{distrib_compa}/$outfile{th_distrib_compa}/g;
	      
	      ################################################################
	      ## Draw a graph with all the decreasing cumulative distributions
	      my $XYgraph_cmd = $SCRIPTS."/XYgraph ".$main::rplot_option." ".$all_graph_options;
	      my $ycols = join ",", 2..(scalar(@th_distrib_files)+2);
	      $XYgraph_cmd .= " -xcol 1 -ycol ".$ycols;
	      $XYgraph_cmd .= " -ymin 0  -ymax 1 ";
	      $XYgraph_cmd .= " -gp 'set size ratio 0.5' ";
	      $graph_file_opt = $large_graph_options." ".$distrib_options." -o ".$outfile{th_distrib_compa}.".".$image_format;
	      &doit($XYgraph_cmd.$graph_file_opt, $dry, $die_on_error, $verbose, $batch, $job_prefix);
	      print $main::out ";\n; XYgraph command\n", $XYgraph_cmd.$graph_file_opt, "\n";
	      
	      ################################################################
	      ## Draw a graph with all the decreasing cumulative distributions
	      ## and a logarithmic Y axis
	      $XYgraph_cmd = $SCRIPTS."/XYgraph ".$main::rplot_option." ".$all_graph_options;
	      $XYgraph_cmd .= " -xcol 1 -ycol ".$ycols;
	      $XYgraph_cmd .= " -ymax 1 -ylog 10";
	      $XYgraph_cmd .= " -gp 'set size ratio 0.5' ";
	      $graph_file_opt = $large_graph_options." ".$distrib_options." -o ".$outfile{th_distrib_compa}."_logy.".$image_format;
	      &doit($XYgraph_cmd.$graph_file_opt, $dry, $die_on_error, $verbose, $batch, $job_prefix);
	      print $main::out ";\n; XYgraph command\n", $XYgraph_cmd.$graph_file_opt, "\n";
	  }	  
      }
    }
  }
  
}





################################################################
## Calculate NWD
################################################################
## Read matrix quality file and calculate the NWD curve

sub Calculate_NWD {
    my ($m_w,$comp_distrib_file,$nwd_seq_type)= @_;
    $outfile{'distrib_nwd'.$nwd_seq_type} = $matrix_prefix{$matrix_name}."_score_distrib_".$nwd_seq_type."_nwd.tab"; 
    push @files_to_index, 'distrib_nwd'.$nwd_seq_type;
    $nwd_all_files{$matrix_name}{$nwd_seq_type}=$outfile{'distrib_nwd'.$nwd_seq_type};
    &RSAT::message::Debug("Adding file to index ", "distrib_nwd",  $matrix_prefix{$matrix_name}."_score_distrib_".$nwd_seq_type."_nwd.tab") if ($main::verbose >= 10);
    $main::out_nwd = &OpenOutputFile( $outfile{'distrib_nwd'.$nwd_seq_type});
    #die "width " . $m_w. "file " . $comp_distrib_file;
    my $dists;
    my ($dist_file) =  &OpenInputFile($comp_distrib_file) ;

    my $head=1;
    my $case_dist= $nwd_seq_type;
    my $base_dist= "theor";

    &RSAT::message::TimeWarn("Calculating NWD between ", $base_dist ," and ",$case_dist , " distributions ")if ($main::verbose >= 5);   


    my $case_col;
    my $base_col;
    my $j=0;
    my $point=0;
    my %p_val_score_case=();
    my %p_val_score_base=();   

    #print join ("\t",";Pvalue","score_". $base_dist,"score_".$case_dist,"NWD")."\n";

    print $main::out_nwd join ("\t","#Pvalue","score_". $base_dist,"score_".$case_dist,"NWD_".$case_dist."_vs_".$base_dist)."\n";

    my %sort_pval;

    ## Read distribution comparison file
    while (<$dist_file>){
	#print $_ ; <STDIN>;
	next if (/'^;'/);		# skip comment lines
	next if (/'^--'/);	# skip mysql-type comment lines

	if ((/^#/) && ($head)){
	    $head=0;
	    #print "header ".$_."\n" ; <STDIN>;
	    my @head = split/\t+/;
	    ## Identify column with the correct header matching the base and case sequence sets IDs
	    foreach my $i (@head) {
		if ( ($i =~/$base_dist/) && !$base_col ){
		    $base_col= $j ;
		}
		elsif ( ($i =~/$case_dist/) && !$case_col ){
		    $case_col= $j ;	
		}
		$j++;
		last if($case_col && $base_col);
	    }
	    @head2=@head;
	    shift(@head2);

	    # If nine if the columns matched the case name die on error.
	    &RSAT::error::FatalError("Please specify an adequate distribution to calculate the NWD\n",
				     "Select one of the following distributions: ", 
				     join("\t",@head2)) unless $case_col;

	    &RSAT::message::Debug("Case column",$case_dist,"#",$case_col) if ($main::verbose >= 10);
	    &RSAT::message::Debug("Base column",$base_dist,"#",$base_col) if ($main::verbose >=10);
	    next;
	}
	next if (/'^#'/ ) ;		# skip coments once the header has been saved
	#next unless ($case_col);
	@line = split /\t+/ ;

	my $score = $line[0];
	my $case_pval= $line[$case_col] ;
	my $base_pval= $line [$base_col] ;
	#my @scores=($case_score,$base_score);


	if (($case_pval =~ /NUL/)
	    || ($base_pval =~ /NUL/)
	    ){
	    next;
	}else{
	    $round_case_pval=sprintf("%.1e", $case_pval);
	    my $round_base_pval= $base_pval ;
	    $sort_pval{$round_case_pval}=$case_pval;
	    push(@{$p_val_score_case{$round_case_pval}}, $score);
	    push(@{$p_val_score_base{$round_base_pval}}, $score);
	    &RSAT::message::Debug("Line point",$score,$round_case_pval,$round_base_pval) if ($main::verbose >= 10);
	}

    }

    my @pvals_list =  (keys(%p_val_score_case),keys(%p_val_score_base));


    %hashTemp = map { $_ => 1 } @pvals_list; ## Declare and fill a hash with pvalue as keys
    @pvals_list = sort keys %hashTemp; ## Sort the keys of this pval hash and restore them in the pval list array

    my %hash_print;
    foreach my $pval (sort {$b cmp $a}(@pvals_list)){
	next unless $p_val_score_case{$pval}; ## If one of the distributions didn't show the pval skip
	next unless $p_val_score_base{$pval};

	&RSAT::message::Debug("Intersection of score value on Pval ",$pval) if ($main::verbose >= 10);
	#print $pval."\n";<STDIN>;
       	my $NWD="";
	my $case_max_score= &RSAT::stats::max(@{$p_val_score_case{$pval}});
	my $base_max_score= &RSAT::stats::max(@{$p_val_score_base{$pval}});
	$NWD = ($case_max_score - $base_max_score) / $m_w ; ## Case max score for the give pvalue minus
	                                                    ## Base max score for the given pvalue
	                                                    ## Divided by the matrix-width to correct

	$main::key_diferences_results{$pval}{$matrix_name}=$NWD if ($NWD);

	&RSAT::message::Debug("Score diference ",$matrix_name,"Pval", $pval," $case_max_score - $base_max_score  $m_w " ,$main::key_diferences_results{$pval}{$matrix_name}=$NWD) if ($main::verbose >= 10);

	## Fill hash to print the NWD table
	$hash_print{$pval}=join ("\t",$pval,$base_max_score ,$case_max_score,$NWD)."\n";
	    
    }

    ## Sort printing hash and print the content to a text file 
    foreach my $pval ( sort {$sort_pval{$b} <=> $sort_pval{$a}}  keys %sort_pval ){
	next unless $hash_print{$pval};
	print $main::out_nwd $hash_print{$pval} ;
    }
    close $dist_file ;
    &RSAT::message::TimeWarn(" NWD subroutine done")  if ($main::verbose >= 5);  

    return ($outfile{'distrib_nwd'.$nwd_seq_type});
}



################################################################
## Calculate OCC
################################################################
## Read matrix quality file and calculate the OCC curve

sub Calculate_OCC { 
    local ($sequence_file, $matrix_file, $matrix_format, $seq_type, $index, @args) = @_; 
    ## Only tab format is supported
    if ($matrix_format ne "tab") {
	&RSAT::error::FatalError("&Calculate_OCC only supports tab format in quick scan mode", $seq_type, $matrix_format, $matrix_file);
    }

    ## Define the output file for the current sequence type
    $outfile{'occ_proba_'.$seq_type} = $matrix_prefix{$matrix_name}."_occ_proba_".$seq_type.".tab";  

    ## Add the file to the list for the comparison of distributions
    if ($index) {
	push @files_to_index, 'occ_proba_'.$seq_type;
	&RSAT::message::Debug("Adding file to index ", 'occ_proba_'.$seq_type,  $outfile{'occ_proba_'.$seq_type}  ) if ($main::verbose >= 2);
    }
    
    #-decimals 1 -bg_pseudo 0.01 -n score -lth score -5 \
    
    local $matrix_scan_cmd = "";
    local $decimals_occ=0;
    $matrix_scan_cmd = $SCRIPTS."/matrix-scan -v ".$main::verbose;
    $matrix_scan_cmd .= " -quick "; 
    $matrix_scan_cmd .= " -matrix_format tab"; ## We use tab as matrix format for compatibiliy with matrix-scan-quick
    $matrix_scan_cmd .= " -bg_format inclusive"; ## We use inclusive as bg format for compatibiliy with matrix-scan-quick
    
    $matrix_scan_cmd .= " -i ".$sequence_file;
    $matrix_scan_cmd .= " -m ".$matrix_file;
    $matrix_scan_cmd .= " -pseudo ".$main::pseudo_counts; 
    $matrix_scan_cmd .= " -2str ";
    $matrix_scan_cmd .= "  -return distrib  -return occ_proba ";
    $matrix_scan_cmd .= " -decimals ".$decimals_occ;
    $matrix_scan_cmd .= " -bgfile ".$outfile{bg_file_inclusive};
    
    $matrix_scan_cmd .= " > ".  $outfile{'occ_proba_'.$seq_type};

    &RSAT::message::TimeWarn("Before scanning for OCC plots", $matrix_scan_cmd) if ($main::verbose >= 5);

    &doit($matrix_scan_cmd, $dry, $die_on_error, $verbose, 0, $job_prefix) if $task{scan} || $task {compare
};
    &RSAT::message::TimeWarn("Scanning to compute occurence probability for OCC plots", $matrix_scan_cmd) if ($main::verbose >= 5);
     
    return ($outfile{'occ_proba_'.$seq_type});
  
}

###############################################################
## Calculate empirical distributions 
###############################################################

sub CalculateEmpiricalDistributions{
    foreach my $seq_type (@local_seq_types) {
	  &RSAT::message::TimeWarn("Analyzing sequence type", $seq_type, $seqfile{$seq_type}) if ($main::verbose >= 5);
	  &CalcSequenceDistrib($seqfile{$seq_type}, $outfile{matrix_tab}, 'tab', $seq_type, 1,  @matrix_scan_options) ;
	  ## Score sequences with the permuted matrices
	  if (($seqfile{$seq_type}) &&
	      (defined($perm_nb{$seq_type})) &&
	      ($perm_nb{$seq_type} > 0)) {

	      &RSAT::message::Warning("Calculating permuted matrices and their emprical distributions for",$seq_type, "Number of permutations", $perm_nb{$seq_type}) if ($main::verbose >= 5);
	      
	      ## Calculate the separate distributions for each permuted matrix
	      ## (this highlights the variability but the graph is noisy)
	      for my $i (1..$perm_nb{$seq_type}) {
		  $perm_suffix = $seq_type."_perm_col_".$i;
		  if (defined($scanopt{$seq_type})) {
		      $scanopt{$perm_suffix} = $scanopt{$seq_type};
		  }
		  push @perm_distrib_files, &CalcSequenceDistrib($seqfile{$seq_type}, $outfile{'matrix_perm_col_'.$i}, "tab", $perm_suffix, $perm_separate_distrib,  @matrix_scan_options) if $task{scan};
	      }
	      
	      ## Compute the distribution for all the permutation tests
	      
	      ## Define the output file for the regrouped permutation tests
	      my $perm_suffix = $seq_type.'_'.$perm_nb{$seq_type}.'perm';
	      $outfile{$perm_suffix} =  $matrix_prefix{$matrix_name}."_scan_".$perm_suffix."_score_distrib.tab"; push @files_to_index, $perm_suffix;
	      &RSAT::message::Debug("Adding file to index ",  $perm_suffix, $outfile{$perm_suffix}  ) if ($main::verbose >= 10);
	      push @distrib_files, $outfile{$perm_suffix}; $file_nb{$perm_suffix} = scalar(@distrib_files);
	      
	      ## Run compare-scores to compute the dCDF of the mergeed permutation test
	      my $merge_cmd = $SCRIPTS."/compare-scores -v 2 ";
	      $merge_cmd .= " -ic 1 -numeric -sc 2";
	      $merge_cmd .= " -files ";
	      $merge_cmd .= join " ", @perm_distrib_files;
	      my $last_col = scalar(@perm_distrib_files) + 1;
	      $merge_cmd .= " | ".$SCRIPTS."/row-stats -before -col 2-".$last_col;
	      &RSAT::message::Debug("Merging permuted distributions", $merge_cmd) if ($main::verbose >= 3);
	      
	      ## Compute the cumulative and decreasing cumlative
	      ## distributions
	      
	      my @weights = ();
	      my @occ = ();
	      my @cum_occ = ();
	      my %merged_occ = ();
	      my $cum_occ = 0;
	      if ($task{scan}){
		  open MERGE, "$merge_cmd |";
		  while (<MERGE>) {
		      chomp();
		      next if /^;/;
		      next if /^#/;
		      next unless /\S/;
		      my @fields = split /\t/, $_;
		      my $weight = $fields[4];
		      my $occ = $fields[2];
		      $cum_occ += $occ;
		      push @weights, $weight;
		      push @occ, $occ;
		      push @cum_occ, $cum_occ;
		  }
		  close MERGE;
		  my $total_occ = $cum_occ[$#cum_occ];
	      
		  
		  ## Print the merged distribution
		  my $merged_distrib = &OpenOutputFile($outfile{$perm_suffix});
		  print $merged_distrib join ("\t", "#weight", "occ", "cum", "dcum", "dCDF"), "\n";
		  for my $i (0..$#weights) {
		      my $dcum_occ = $total_occ - $cum_occ[$i]+$occ[$i];
		      my $dcdf = $dcum_occ / $total_occ;
		      print $merged_distrib join ("\t",
						  $weights[$i],
						  $occ[$i],
						  $cum_occ[$i],
						  $dcum_occ,
						  sprintf("%7g", $dcdf)
			  ), "\n";
		  }
		  close $merged_distrib;
		  &RSAT::message::TimeWarn("Exported merged distribution", $outfile{$perm_suffix}) if ($main::verbose >= 10);    
	      }
	  }
      }


}

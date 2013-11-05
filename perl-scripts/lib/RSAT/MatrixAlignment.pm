##############################################################
#
# Class MatrixReader
#
package RSAT::MatrixAlignment;

use RSAT::GenericObject;
use RSAT::matrix;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

=pod

=head1 NAME

RSAT::MatrixAlignment

=head1 DESCRIPTION

Class holding methods for matrix alignment and clustering.

=cut



################################################################
## Method definition



=pod

=item B<AlignMatrixPair>

Generate a single matrix summarizing the alignment of two matrices.

The offset should be provided. The program computes the sum or the
mean count for each row (residue) and each column (aligned position).

Usage:
  my $merged_matrix =
    &RSAT::MatrixAlignment::AlignMatrixPair($matrix1, $matrix2, $offset, $strand, %args);

The method can also return the shifted original matrices, which are
convenient for drawing aligned logos.

  my ($merged_matrix, $matrix1_shifted, $matrix2_shifted) =
     &RSAT::MatrixAlignment::AlignMatrixPair($matrix1, $matrix2, $offset, $strand, %args);

Arguments:
  stat=> mean|sum   compute the mean or the sum of the counts, respectively.

=cut

sub AlignMatrixPair {
  my ($matrix1, $matrix2, $offset, $strand, %args) = @_;
  my $stat = $args{stat} || "sum";
  my $rc = "";
  $rc = "_rc" if ($strand eq "R"); ## Suffix added to matrix ID, name and descrption to indicate reverse complement

  my $id1 = $matrix1->get_attribute("id");
  my $desc1 = $id1;

  my $id2 = $matrix2->get_attribute("id");
  my $desc2 = $id2.$rc;

  my $name1 = $matrix1->get_attribute("name");
  if (($name1) && ($name1 ne $id1)) {
    $desc1 .= " (".$name1.")";
  }
  my $name2 = $matrix2->get_attribute("name");
  $name2 .= $rc;
  if (($name2) && ($name2 ne $id2.$rc)) {
    $desc2 .= " (".$name2.")";
  }

  &RSAT::message::TimeWarn("Pairwise matrix alignment", $desc1, $desc2, "offset=".$offset,  "strand=".$strand)
    if ($main::verbose >= 3);

  my @counts1 = $matrix1->getMatrix();
  my @counts2;
  if ($strand eq "R") {
    @counts2 = $matrix2->getCountRC();
  } else {
    @counts2 = $matrix2->getMatrix();
  }
  my @alphabet = $matrix1->getAlphabet();

  my $ncol1 = $matrix1->get_attribute("ncol");
  my $ncol2 = $matrix2->get_attribute("ncol");
  my $nrow = $matrix1->get_attribute("nrow");

  ## Compute aligned matrix positions
  my $end1 = &RSAT::stats::min($ncol1, $ncol2+$offset);
  my $start1 = &RSAT::stats::max(1, $offset+1, $end1-$ncol2+1);
  my $w = $end1-$start1+1;
  my $ncol = $ncol1 + $ncol2 - $w;
  my $start2 = &RSAT::stats::max(1, 1- $offset);
  my $end2 = &RSAT::stats::min($start2+$w-1, $ncol2);

  ## Create the shifted matrix 1.  This matrix is not necessary for
  ## computing the aligned matrix, but convenient for visualization
  ## (aligned matrices, aligned logos).
  my $shifted_matrix1 = new RSAT::matrix();
  $shifted_matrix1->force_attribute("ncol", $ncol);
  $shifted_matrix1->setAlphabet_lc(@alphabet);
#  $shifted_matrix1->force_attribute("nrow", scalar(@alphabet));
  my @shifted_counts1 = ();

  ## Create the shifted matrix 2.  This matrix is not necessary for
  ## computing the aligned matrix, but convenient for visualization
  ## (aligned matrices, aligned logos).
  my $shifted_matrix2 = new RSAT::matrix();
  $shifted_matrix2->force_attribute("ncol", $ncol);
  $shifted_matrix2->setAlphabet_lc(@alphabet);
#  $shifted_matrix2->force_attribute("nrow", scalar(@alphabet));
  my @shifted_counts2 = ();

  ## Create the alignment matrix
  my $aligned_matrix = new RSAT::matrix();
  $aligned_matrix->force_attribute("id", $id1."_".$id2.$rc."_".$stat),
  $aligned_matrix->force_attribute("ncol", $ncol);
  $aligned_matrix->setAlphabet_lc(@alphabet);
#  $aligned_matrix->force_attribute("nrow", scalar(@alphabet));
  my @ali_counts = ();

  ## Compute the cell values of the aligned matrix
  my $count1;
  my $count2;
  my $shift1 = 0;
  my $shift2 = 0;
  if ($offset >= 0) {
    $shift1 = 0;
    $shift2 = $offset;
  } else {
    $shift1 = -$offset;
    $shift2 = 0;
  }
  $shifted_matrix1->set_parameter("shift",$shift1);
  $shifted_matrix1->force_attribute("id", $id1."_shift".$shift1);
  $shifted_matrix1->force_attribute("identifier", $id1."_shift".$shift1);
  $shifted_matrix2->set_parameter("shift",$shift2);
  $shifted_matrix2->force_attribute("id", $id2.$rc."_shift".$shift2);
  $shifted_matrix2->force_attribute("identifier", $id2.$rc."_shift".$shift2);
  for my $c (1..$ncol) {
    for my $r (1..$nrow) {
      $c1 = $c - $shift1;
      $c2 = $c - $shift2;

      ## Counts in the first matrix
      if (($c1 < 1) || ($c1 > $ncol1)) {
	$count1 = 0;
      } else {
	$count1 = $counts1[$c1-1][$r-1] || 0;
      }
      $shifted_counts1[$c-1][$r-1] = $count1;

      ## Counts in the second matrix
      if (($c2 < 1) || ($c2 > $ncol2)) {
	$count2 = 0;
      } else {
	$count2 = $counts2[$c2-1][$r-1] || 0;
      }
      $shifted_counts2[$c-1][$r-1] = $count2;

      my $sum = $count1 + $count2;
      my $count;
      if ($stat eq "sum") {
	$count = $sum;
      } elsif ($stat eq "mean") {
	$count = $sum/2;
      } else {
	&RSAT::error::FatalError("&AlignMatrixPair", "Invalid value for the stat argment. Supported: mean|sum.");
      }
      $ali_counts[$c-1][$r-1] = $count;

#       &RSAT::message::Debug("alignment", 
# 			    $id1,
# 			    $id2.$rc,
# 			    "offset=".$offset, 
# 			    "r=".$r,
# 			    "c=".$c."/".$ncol,
# 			    "c1=".$c1."/".$ncol1,
# 			    "c2=".$c2."/".$ncol1,
# 			    "count1=".$shifted_counts1[$c][$r],
# 			    "count2=".$shifted_counts2[$c][$r],
# 			    "sum=".$sum,
# 			    "count=".$count,
# 			   ) if ($main::verbose >= 10);

    }
#    die "HELLO" if ($c > 3);
  }

  ## Set the counts of the shifted matrix 1
  $shifted_matrix1->setMatrix($nrow, $ncol, @shifted_counts1);
  $shifted_matrix1->calcConsensus();

  ## Set the counts of the shifted matrix 2
  $shifted_matrix2->setMatrix($nrow, $ncol, @shifted_counts2);
  $shifted_matrix2->calcConsensus();

  ## Set the counts of the aligned matrix
  $aligned_matrix->setMatrix ($nrow, $ncol, @ali_counts);
  $aligned_matrix->calcConsensus();
  return($aligned_matrix, $shifted_matrix1, $shifted_matrix2);
}




=pod

=item B<AlignMatricesOneToN>

Generate a single matrix summarizing the alignment of several
matrices.

Usage:

  my ($merged_matrix, @shifted_matrices) = 
    &RSAT::MatrixAlignment::AlignMatricesOneToNOneToN($matrix1, 
                                                  \@matching_matrices,
                                                  \@matching_offsets,
                                                  \@matching_strands,
                                                  \@matching_scores, %args);


Arguments:
  stat=> mean|sum   compute the mean or the sum of the counts, respectively.

The method returns the merged matrix and the list of shifted matrices,
which are convenient for drawing aligned logos.

=cut

sub AlignMatricesOneToN {
  my ($matrix1, $matching_compa_ids_ref,  $matching_matrices_ref, $matching_offsets_ref, $matching_strands_ref, $matching_scores_ref, %args) = @_;
  my @shifted_matrices = ();
  my @matching_compa_ids = @{$matching_compa_ids_ref};
  my @matching_matrices = @{$matching_matrices_ref};
  my @matching_offsets = @{$matching_offsets_ref};
  my @matching_strands = @{$matching_strands_ref};
  my @matching_scores = @{$matching_scores_ref};
  my $stat = $args{stat} || "sum";

  ## Parameters for the first matrix (the one on which the other ones
  ## will be aligned)
  my $id1 = $matrix1->get_attribute("id");
  my $name1 = $matrix1->get_attribute("name") || $id1;
  my @counts1 = $matrix1->getMatrix();
  my @alphabet = $matrix1->getAlphabet();
  my $ncol1 = $matrix1->get_attribute("ncol");
  my $nrow = $matrix1->get_attribute("nrow");

  ## Create the shifted matrix 1.  This matrix is not necessary for
  ## computing the aligned matrix, but convenient for visualization
  ## (aligned matrices, aligned logos).
  my $shifted_matrix1 = new RSAT::matrix();
  my $max_offset = &RSAT::stats::max(@matching_offsets);
  my $min_offset = &RSAT::stats::min(@matching_offsets);

  my @ncol2 = ();
  my @start1 = ();
  my @end1 = ();
  my @start2 = ();
  my @end2 = ();
  my @w = (); ## Number of aligned columns
  my @ncol = (); ## Total number of columns per pairwise alignment
  foreach my $m (0..$#matching_matrices) {
    my $matrix2 = $matching_matrices[$m];
    my $offset = $matching_offsets[$m];
    my $ncol2 = $matrix2->get_attribute("ncol"); push @ncol2, $ncol2;
    my $end1 = &RSAT::stats::min($ncol1, $ncol2+$offset); push @end1, $end1;
    my $start1 = &RSAT::stats::max(1, $offset+1, $end1-$ncol2+1); push @start1, $start1;
    my $w = $end1-$start1+1; push @w, $w;
    my $ncol = $ncol1 + $ncol2 - $w; push @ncol, $ncol;
    my $start2 = &RSAT::stats::max(1, 1- $offset); push @start2, $start2;
    my $end2 = &RSAT::stats::min($start2+$w-1, $ncol2); push @end2, $end2;
  }
  my $ncol = &RSAT::stats::max(@ncol);


  my $shift1 = 0;
  if ($min_offset < 0) {
    $shift1 = -$min_offset;
  }
  $shifted_matrix1->force_attribute("ncol", $ncol);
  $shifted_matrix1->setAlphabet_lc(@alphabet);
  $shifted_matrix1->set_parameter("shift",$shift1);
  $shifted_matrix1->force_attribute("id", $id1."_shift".$shift1);
  $shifted_matrix1->force_attribute("name", $name1);
  $shifted_matrix1->force_attribute("identifier", $id1."_shift".$shift1);

  ## Generate the shifted matrix 1
  my @shifted_counts1 = ();
  my $count1;
  for my $c (1..$ncol) {
    for my $r (1..$nrow) {
      $c1 = $c - $shift1;

      ## Counts in the first matrix
      if (($c1 < 1) || ($c1 > $ncol1)) {
	$count1 = 0;
      } else {
	$count1 = $counts1[$c1-1][$r-1] || 0;
      }
      $shifted_counts1[$c-1][$r-1] = $count1;
    }
  }
  ## Set the counts of the shifted matrix 1
  $shifted_matrix1->setMatrix($nrow, $ncol, @shifted_counts1);
  $shifted_matrix1->calcConsensus();
  push @shifted_matrices, $shifted_matrix1;

  &RSAT::message::Debug("&RSAT::MatrixAlignment::AlignMatricesOneToN()",
			$id1,
			"ncol1=".$ncol1,
			"max_offset=".$max_offset,
			"min_offset=".$min_offset,
			"ncol=".$ncol,
			"shift1=".$shift1,
		       ) if ($main::verbose >= 4);

  ## Description of the first matrix (reference for the 1-to-n alignment)
  my $desc1 = $id1;
  if (($name1) && ($name1 ne $id1)) {
    $desc1 .= " (".$name1.")";
  }
  $desc1 .= "; m=0 (reference)";
  $desc1 .= "; ncol1=".$ncol1;
  $desc1 .= "; shift=".$shift1;
  $desc1 .= "; ncol=".$ncol;
  $desc1 .= "; ".$shifted_matrix1->get_attribute("consensus.IUPAC");
  $shifted_matrix1->set_parameter("shift", $shift1);
  $shifted_matrix1->set_parameter("description", $desc1);
#  &RSAT::message::TimeWarn("Aligning matrices 1-to-n", $desc1, "n=".scalar(@matching_matrices)) if ($main::verbose >= 10);

#  &RSAT::message::Debug("Shifted matrix 0","\n".$shifted_matrix1->toString(format=>'tab', type=>"counts")) if ($main::verbose >= 10);

  my $n2 = scalar(@matching_matrices);
  foreach my $m (0..$#matching_matrices) {
    my $compa_id = $matching_compa_ids[$m];
    my $matrix2 = $matching_matrices[$m];
    my $offset = $matching_offsets[$m];
    my $strand = $matching_strands[$m];
    my $rc = "";
    $rc = "_rc" if ($strand eq "R"); ## Suffix added to matrix ID, name and descrption to indicate reverse complement

    my $score = $matching_scores[$m];
    my $id2 = $matrix2->get_attribute("id");
    my $name2 = $matrix2->get_attribute("name") || $id2;
    $name2 .= $rc;


    ## Get the direct or reverse complementary matrix, depending on the alignment strand
    my @counts2;
    if ($strand eq "R") {
      @counts2 = $matrix2->getCountRC();
    } else {
      @counts2 = $matrix2->getMatrix();
    }

    ## Get the pre-computed positions for the alignment
    my $ncol2 = $ncol2[$m];
    my $start1 = $start1[$m];
    my $end1 = $end1[$m];
    my $start2 = $start2[$m];
    my $end2 = $end2[$m];
    my $w = $w[$m];

    ## Compute the cell values of the shifted matrix
    my $shift2 = $shift1 + $offset;

    ## Create the shifted matrix 2.  This matrix is not necessary for
    ## computing the aligned matrix, but convenient for visualization
    ## (aligned matrices, aligned logos).
    my $shifted_matrix2 = new RSAT::matrix();
    $shifted_matrix2->set_attribute("compa_id", $compa_id);
    $shifted_matrix2->force_attribute("ncol", $ncol);
    $shifted_matrix2->setAlphabet_lc(@alphabet);
    #  $shifted_matrix2->force_attribute("nrow", scalar(@alphabet));
    $shifted_matrix2->set_parameter("shift",$shift2);
    $shifted_matrix2->force_attribute("id", $id2.$rc."_shift".$shift2);
    $shifted_matrix2->force_attribute("name", $name2);
    $shifted_matrix2->force_attribute("identifier", $id2.$rc."_shift".$shift2);
    my @shifted_counts2 = ();

#     ## Create the alignment matrix
#     my $aligned_matrix = new RSAT::matrix();
#     $aligned_matrix->force_attribute("id", $id1."_".$id2.$rc."_".$stat);
#     $aligned_matrix->force_attribute("ncol", $ncol);
#     $aligned_matrix->setAlphabet_lc(@alphabet);
#     #  $aligned_matrix->force_attribute("nrow", scalar(@alphabet));
#     my @ali_counts = ();

    ## Description of the shifted matrix (aligned on the ref matrix)
    my $desc2 = $id1." versus ".$id2.$rc;
    if (($name2) && ($name2 ne $id2.$rc)) {
      $desc2 .= " (".$name2.")";
    }
    $desc2 .= "; m=".($m+1)."/".$n2;
    $desc2 .= "; ncol2=".$ncol2;
    $desc2 .= "; w=".$w;
    $desc2 .= "; offset=".$offset;
    $desc2 .= "; strand=".$strand;
    $desc2 .= "; shift=".$shift2;
    $desc2 .= sprintf("; score=%7g", $score);
#    &RSAT::message::TimeWarn("Aligning matrices 1-to-n", $desc2) if ($main::verbose >= 10);

    for my $c (1..$ncol) {
      for my $r (1..$nrow) {
	$c2 = $c - $shift2;

	## Counts in the second matrix
	my $count2;
	if (($c2 < 1) || ($c2 > $ncol2)) {
	  $count2 = 0;
	} else {
	  $count2 = $counts2[$c2-1][$r-1] || 0;
	}
	$shifted_counts2[$c-1][$r-1] = $count2;

# 	my$count1 = $shifted_counts1[$c-1][$r-1];
# 	my $sum = $count1 + $count2;
# 	my $count;
# 	if ($stat eq "sum") {
# 	  $count = $sum;
# 	} elsif ($stat eq "mean") {
# 	  $count = $sum/2;
# 	} else {
# 	  &RSAT::error::FatalError("&AlignMatricesOneToN", "Invalid value for the stat argment. Supported: mean|sum.");
# 	}
# 	$ali_counts[$c-1][$r-1] = $count;
# 	#       &RSAT::message::Debug("alignment", 
# 	# 			    $id1,
# 	# 			    $id2.$rc,
# 	# 			    "offset=".$offset, 
# 	# 			    "r=".$r,
# 	# 			    "c=".$c."/".$ncol,
# 	# 			    "c1=".$c1."/".$ncol1,
# 	# 			    "c2=".$c2."/".$ncol1,
# 	# 			    "count1=".$shifted_counts1[$c][$r],
# 	# 			    "count2=".$shifted_counts2[$c][$r],
# 	# 			    "sum=".$sum,
# 	# 			    "count=".$count,
# 	# 			   ) if ($main::verbose >= 10);

      }
    }

    ## Set the counts of the shifted matrix 2
    $shifted_matrix2->setMatrix($nrow, $ncol, @shifted_counts2);
    $shifted_matrix2->calcConsensus();
    $desc2 .= "; ".$shifted_matrix2->get_attribute("consensus.IUPAC");
    $shifted_matrix2->set_parameter("description", $desc2);
    push @shifted_matrices, $shifted_matrix2;

    &RSAT::message::Debug("Shifted matrix ".($m+1),"\n".$shifted_matrix2->toString(format=>'tab', type=>"counts")) if ($main::verbose >= 10);


    ## Set the counts of the aligned matrix
#    $aligned_matrix->setMatrix ($nrow, $ncol, @ali_counts);
#    $aligned_matrix->calcConsensus();
  }

  return($aligned_matrix, @shifted_matrices);
}

return 1;

__END__


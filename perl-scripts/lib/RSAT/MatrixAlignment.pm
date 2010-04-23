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

=pod

=item B<AlignMatrices>

Generate a single matrix summarizing the alignment of two matrices.

The offset should be provided. The program computes the sum or the
mean count for each row (residue) and each column (aligned position).

Usage:
  my $aligned_matrix = &RSAT::MatrixAlignment::AlignMatrices($matrix1, $matrix2, $offset, %args);

The method can also return the shifted original matrices, which are
convenient for drawing aligned logos.

  my ($aligned_matrix, $shifted_matrix1, $shifted_matrix2) =
         &RSAT::MatrixAlignment::AlignMatrices($matrix1, $matrix2, $offset, %args);

Arguments:
  stat=> mean|sum   compute the mean or the sum of the counts, respectively.

=cut 

sub AlignMatrices {
  my ($matrix1, $matrix2, $offset, %args) = @_;
  my $stat = $args{stat} || "sum";

  my $desc1 = $id1 = $matrix1->get_attribute("id");
  my $desc2 = $id2 = $matrix2->get_attribute("id");

  my $name1 = $matrix1->get_attribute("name");
  if (($name1) && ($name1 ne $id1)) {
    $desc1 .= " (".$name1.")";
  }
  my $name2 = $matrix2->get_attribute("name");
  if (($name2) && ($name2 ne $id2)) {
    $desc2 .= " (".$name2.")";
  }

  &RSAT::message::TimeWarn("Aligning matrices", $desc1, $desc2, "offset=".$offset)
    if ($main::verbose >= 2);

  my @counts1 = $matrix1->getMatrix();
  my @counts2 = $matrix2->getMatrix();
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
  ## computing the alignmed matrix, but convenient for visualization
  ## (aligned matrices, aligned logos).
  my $shifted_matrix1 = new RSAT::matrix();
  $shifted_matrix1->force_attribute("ncol", $ncol);
  $shifted_matrix1->setAlphabet_lc(@alphabet);
  $shifted_matrix1->set_attribute("nrow", scalar(@alphabet));
  my @shifted_counts1 = ();

  ## Create the shifted matrix 2.  This matrix is not necessary for
  ## computing the alignmed matrix, but convenient for visualization
  ## (aligned matrices, aligned logos).
  my $shifted_matrix2 = new RSAT::matrix();
  $shifted_matrix2->force_attribute("ncol", $ncol);
  $shifted_matrix2->setAlphabet_lc(@alphabet);
  $shifted_matrix2->set_attribute("nrow", scalar(@alphabet));
  my @shifted_counts2 = ();

  ## Create the alignment matrix
  my $aligned_matrix = new RSAT::matrix();
  $aligned_matrix->force_attribute("id", $id1."_".$id2."_".$stat),
  $aligned_matrix->force_attribute("ncol", $ncol);
  $aligned_matrix->setAlphabet_lc(@alphabet);
  $aligned_matrix->set_attribute("nrow", scalar(@alphabet));
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
  $shifted_matrix2->force_attribute("id", $id2."_shift".$shift2);
  $shifted_matrix2->force_attribute("identifier", $id2."_shift".$shift2);
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
	&RSAT::error::FatalError("&AlignMatrices", "Invalid value for the stat argment. Supported: mean|sum.");
      }
      $ali_counts[$c-1][$r-1] = $count;

#       &RSAT::message::Debug("alignment", 
# 			    $id1,
# 			    $id2,
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

return 1;

__END__


###############################################################
#
# Class Operon
#
package RSAT::Operon;
$_count = 0;
$_prefix = "operon";

use RSAT::GenericObject;
use RSAT::error;
use RSAT::GenomeFeature;

@ISA = qw( RSAT::GenericObject );

### class attributes

=pod

=head1 NAME

    RSAT::Operon

=head1 DESCRIPTION

Operon class. 

=cut

################################################################
=pod

    =item new()

    Create a new Operon.

=cut
sub new {
  my ($class, %args) = @_;
  my $self = bless {
		    %args,
		   }, $class;
  $self->init();
  return $self;
}



################################################################
=pod

=item B<add_gene()>

Add a gene object to the operon.

=cut

sub add_gene {
  my ($self, $gene) = @_;    
  &RSAT::message::Info("Operon", 
		       $self->get_attribute("id"),
		       "Adding gene",
		       $gene->get_attribute("id"),
		       $gene->get_attribute("name"),
		      ) if ($main::verbose >= 5);

  if (defined($self->{"gene_hash"}->{$gene})) {
    &RSAT::message::Warning("Gene", $gene->get_attribute("id"), "is already member of the operon", $self->get_attribute("id") ) if ($main::verbose >= 5);
  } else {
    $self->add_hash_attribute("gene_hash", $gene, 1);
    $self->push_attribute("genes", $gene);
  }
}

# ################################################################
# =pod
#
# =item B<get_genes()>
#
# Get the list of member genes.
#
# =cut
#
# sub get_genes {
#   my ($self) = @_;    
#   return  ($self->get_attribute("genes"));
# }




################################################################
=pod 

=item B<get_strand>

Return the strand of the operon. If required, calculate this strand from the member genes.

=cut

sub get_strand {
  my ($self) = @_;
  unless ($self->get_attribute("strand")) {
    $self->check_strand();
  }
  return $self->get_attribute("strand");
}

################################################################
=pod

=item B<check_strand>

Check strand consistency: all genes of an operon must by definition be on the same strand.
Assign the common strand of all genes to the operon. 

Return 1 if all genes have the same strand, 0 otherwise.

=cut

sub check_strand {
  my ($self) = @_;
  my @genes = $self->get_attribute("genes");
  my $operon_strand;
  my $is_consistent = 1;
  foreach my $gene (@genes) {
    my $gene_strand = $gene->get_attribute("strand");
    if ($operon_strand) {
      unless ($gene_strand eq $operon_strand) {
	$is_consistent = 0;
	&RSAT::message::Warning(join("\t", "Operon", $self->getattribute("id"), $operon_strand,
				     "Strand inconsistency", $gene->get_attribute("id"), $gene_strand,
				    )); 
	$operon_strand = "ERROR";
	last;
      }
    } else {
      $operon_strand = $gene_strand;
    }
    $self->force_attribute("strand", $operon_strand);
  }
  return($is_consistent);
}



################################################################
=pod 

=item B<get_sorted_genes>

Return sorted genes of the operon. If required, sort the genes with
$self->sort_genes().

=cut

sub get_sorted_genes {
  my ($self) = @_;
  unless (defined($self->{sorted_genes})) {
    $self->sort_genes();
  }
  return $self->get_attribute("sorted_genes");
}

################################################################
=pod 

=item B<sort_genes>  

Sort genes according to their relative chromosomal position, and to
the strand of the operon.

Genes are always sorted according to the relative positions of their
start codons.  Thus, if the operon is on the direct strand (D), genes
are sorted by increasing values of left limits. If the operon is on
the reverse strand (R), genes are sorted by decreasing values of right
limits.

=cut

sub sort_genes {
  my ($self) = @_;
  my @genes = $self->get_attribute("genes");
  my $operon_strand = $self->get_strand();

  ## Sort genes
  my @sorted_genes = ();
  if (scalar(@genes) > 0) {
    if ($operon_strand eq "D") {
      @sorted_genes = sort {$a->get_attribute("left") <=> $b->get_attribute("left") } @genes;
    } else {
      @sorted_genes = sort {$b->get_attribute("right") <=> $a->get_attribute("right") } @genes;
    }
  }
  $self->set_array_attribute("sorted_genes", @sorted_genes);

  ## Collect sorted gene IDs
  my @sorted_gene_ids = ();
  foreach my $gene (@sorted_genes) {
    push @sorted_gene_ids, $gene->get_attribute("id");
  }
  $self->set_array_attribute("sorted_gene_ids", @sorted_gene_ids);

  ## Collect sorted gene names
  my @sorted_gene_names = ();
  foreach my $gene (@sorted_genes) {
    push @sorted_gene_names, $gene->get_attribute("name");
  }
  $self->set_array_attribute("sorted_gene_names", @sorted_gene_names);

  ## Set the operon name
  $self->force_attribute("name", join ("-", @sorted_gene_names));

}



################################################################
=pod

=item B<set_leader>

Set the identity of the operon leader gene. This specifies three
attributes: "leader", "leader_id" and "leader_name".

=cut
sub set_leader {
  my ($self, $feature) = @_;
  $self->force_attribute("leader", $feature);
  $self->force_attribute("leader_id", $feature->get_attribute("id"));
  $self->force_attribute("leader_name", $feature->get_attribute("name"));
}

################################################################
=pod

=item B<set_trailer>

Set the identity of the operon trailer gene. This specifies three
attributes: "trailer", "trailer_id" and "trailer_name".

=cut
sub set_trailer {
  my ($self, $feature) = @_;
  $self->force_attribute("trailer", $feature);
  $self->force_attribute("trailer_id", $feature->get_attribute("id"));
  $self->force_attribute("trailer_name", $feature->get_attribute("name"));
}

################################################################

=pod

=item <get_gene_nb>

Return the number of genes in the operon.

=cut

sub gene_nb {
  my ($self) = @_;
  my @genes = $self->get_attribute("genes");
  return (scalar(@genes));
}

return 1;

__END__


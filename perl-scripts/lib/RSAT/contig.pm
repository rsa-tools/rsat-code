###############################################################
#
# Manipulation of contigs
#
package RSAT::contig;

use RSAT::GenericObject;
@ISA = qw( RSAT::GenericObject );

=pod

=head1 NAME

    RSAT::contig

=head1 DESCRIPTION

Object for manipulating contigs. 

=cut


################################################################
=pod

=item new()

Create a new contig.

=cut
sub new {
    my ($class, %args) = @_;
    my $contig = bless {
	}, $class;
    return $contig;
}


################################################################
=pod

=item init

Initalize the contig.

=cut
sub init() {
    my ($self) = @_;
    $self->set_genes();
}

################################################################
=pod

=item get_prefix

Get the prefix.

=cut
sub get_prefix {
    my ($class) = @_;
    return $_prefix;
}

################################################################
=pod

=item get_genes

Get the list of genes annotated for this contig.

=cut
sub get_genes {
    my ($self) = @_;
    return @{$self->{genes}};
}

################################################################
=pod

=item set_genes

Specify the list of genes annotated for this contig.

=cut
sub set_genes {
    my ($self, @genes) = @_;
    @{$self->{genes}} = @genes;
}

################################################################
=pod

=item add_gene

Add a gene to the list. 

=cut
sub add_gene {
    my ($self, $gene) = @_;
    push @{$self->{genes}}, $gene;
}

################################################################
=pod

=item count_genes

Return the number of genes on this contig.

=cut
sub count_genes {
    my ($self) = @_;
    return $#{$self->{genes}} + 1;
}


################################################################
=pod

=item get_organism

Return the organism (perl object of class RSAT::organism) to which
this contig belongs.

=cut
sub get_organism {
    my ($self) = @_;
    return $self->{organism};
}

################################################################
=pod

=item set_organism

Specify the organism (perl object of class RSAT::organism) to which
this contig belongs.

=cut
sub set_organism {
    my ($self, $new_organism) = @_;
    $self->{organism} = $new_organism;
}


################################################################
=pod

=item get_sequence

Return the sequence (string). This is a wrapper, which passes the
arguments to the embedded RSAT::Sequence object.

=cut
sub get_sequence {
  my ($self, @args) = @_;
  my $sequence_object = $self->get_sequence_object();
  return $sequence_object->get_sequence(@args);
}

################################################################
=pod

=item get_length

Return the length of the contig. This is a wrapper, which passes the
request to the embedded RSAT::Sequence object.

=cut
sub get_length {
    my ($self) = @_;
    my $sequence_object = $self->get_sequence_object();
    return $sequence_object->get_length();
}

################################################################
=pod

=item get_sequence_object

Return the sequence object (RSAT::sequence) of this contig. 

=cut
sub get_sequence_object {
    my ($self) = @_;
    return $self->get_attribute("sequence");
}


################################################################
=pod

=item set_sequence_object

Specify the sequence object (RSAT::sequence) of this contig. 

=cut
sub set_sequence_object {
    my ($self, $new_sequence) = @_;
    $self->{sequence} = $new_sequence;
}


return 1;


__END__

=pod

=back


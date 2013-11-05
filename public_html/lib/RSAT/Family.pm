###############################################################
#
# Class Family
#
package RSAT::Family;

use RSAT::GenericObject;
use RSAT::error;
@ISA = qw( RSAT::GenericObject );

### class attributes
$_count = 0;
$_prefix = "ctg_";
@_objects = ();
%_name_index = ();
%_id_index = ();
%_attribute_count = ();
%_attribute_cardinality = (id=>"SCALAR",
			   name=>"SCALAR",
			   organism=>"SCALAR",
			   size=>"SCALAR"
			  );

=pod

=head1 NAME

    RSAT::Family

=head1 DESCRIPTION

Class used to store a family (cluster) of genes.

=cut



################################################################

=pod

=item get_members()

Return the list of members

=cut

sub get_members {
    my ($self) = @_;
    return @{$self->{members}};
}

################################################################

=pod

=item get_size()

return the number of members

=cut
sub get_size {
    my ($self) = @_;
    my $size = scalar(@{$self->{members}});
    return $size;
}

################################################################

=pod

=item new_member($new_member, $allow_duplicates)

Add a member to the family, if it has not yet been inserted.

By default, a member fcan be inserted only once in a family.  If the argument
"allow_duplicates" is set to 1, this check is disabled (the list of members
can contain several times the same entry).

Usage: $family->new_member($member);

A score can be assigned to the new member with the argument score=>$score

Usage: $family->new_member($member,score=>$score);

A hash table with the scores can be otbained with the command

%scores = $family->get_attribute("scores");

=cut
sub new_member {
  my ($self, $new_member, $allow_dup, %args) = @_;
  if ($allow_duplicates) {
    $self->push_attribute("members", $new_member);
  } else {
    if ($self->is_member($new_member)) {
      &RSAT::message::Warning(join("\t", "Family", $self->get_attribute("name"), "skipped duplicate member", $new_member)) if ($main::verbose >= 3);
    } else {
      &RSAT::message::Info(join("\t", "Family", $self->get_attribute("name"), "Adding member", $new_member)) if ($main::verbose >= 4);
      $self->add_hash_attribute("member_index", $new_member, 1);
      $self->push_attribute("members", $new_member);
    }
  }
  if (defined($args{score})) {
    $self->add_hash_attribute("scores", $new_member, $args{score});
    #	&RSAT::message::Debug("Family", $self->get_attribute("name"), "member", $new_member, "score", $args{score}) if ($main::verbose >= 10);
  }
}

################################################################

=pod

=item set_members(@members)

Set the member list for the class.

It is more efficient to set the list in one shot than to call
iteratively the method &new_member().

=cut

sub set_members {
  my ($self, @members) = @_;
  $self->set_array_attribute("members", @members);
}


################################################################

=pod

=item is_member($member)

Check whether an element is already member of the family.

=cut
sub is_member {
    my ($self, $member) = @_;
    my %member_index = $self->get_attribute("member_index");
    if ($member_index{$member}) {
	return 1;
    } else {
	return 0;
    }
}

return 1;

__END__


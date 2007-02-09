###############################################################
#
# Class go
#
package RSAT::go;

use RSAT::GenericObject;
use RSAT::error;
use RSAT::util;
@ISA = qw( RSAT::GenericObject );

### class attributes
@ISA = qw( RSAT::GenericObject );
=pod

=head1 NAME

    RSAT::go

=head1 DESCRIPTION

    Implementation of the gene ontology complete directed acyclic graph. 
    go class. This class allows to look for the ancesters and children of the each gene ontology  annotation.
    
=cut


################################################################

=pod

=item B<read_from_obo>

Read the graph from a tab-delimited text file.


 Title    : read_from_obo
 Usage    : $graph->read_from_obo($input_file)
 Function : Read the graph from a obo text file.
 Returns  : void

=cut



sub read_from_obo {
  my ($self, $inputfile) = @_;
  &RSAT::message::TimeWarn("Loading GO from obo file", $inputfile) if ($main::verbose >= 2);
  ($main::in) = &RSAT::util::OpenInputFile($inputfile); 

  my @ids;
  $ids[0] = "-1";
  my $idcpt = 0;
  my $name;
  my $namespace;
  my $def;
  my @isa = ("-1");
  my @parenttable;
  my @childrentable;
  my %indextable;
  my %functionname;
  my %name_id;
  my %namespace;
  my %definition;
  while (my $ligne = <$main::in>) {

    chomp($ligne);
    if ($ligne eq "[Term]" && $ids[0] ne "-1") {
      my $parentnb = 0;
      foreach my $id (@ids) {
        my $parentcpt;
        if (defined($indextable{$id})) {
          $parentcpt = $indextable{$id};
        } else {
          $parentcpt = $idcpt;
          $indextable{$id} = $parentcpt;
          $idcpt++;
        }
        foreach my $parent (@isa) {
          $functionname{$id} = $name;
          $name_id{$name} = $id;
          $definition{$id} = $def;
          $namespace{$id} = $namespace;
          ### Add the parents in @parenttable 
          $parenttable[$parentcpt][$parentnb] = $parent;
          $parentnb++;
          ### Add the children in @childrentable 
          my $childrencpt;
          if (defined($indextable{$parent})) {
            $childrencpt = $indextable{$parent};
          } else {
            $childrencpt = $idcpt;
            $indextable{$parent} = $childrencpt;
            $idcpt++;
          }
          push @{$childrentable[$childrencpt]}, $id;
          
        }
      }
      @ids = ();
      $name = "";
      $namespace = "";
      $def = "";
      @isa = ("-1");
    }
    if ($ligne =~ /^id\: /) {
      $ligne =~ s/id\: //;
      $ids[0] = $ligne;
    } elsif ($ligne =~ /^alt_id\: /) {
      $ligne =~ s/alt_id\: //;
      push(@ids, $ligne);    
    } elsif ($ligne =~ /namespace\: /) {
      $ligne =~ s/namespace\: //;
      $namespace = $ligne;  
    } elsif ($ligne =~ /name\: /) {
      $ligne =~ s/name\: //;
      $name = $ligne;
    } elsif ($ligne =~ /is_a\: /) {
      my $parent = substr($ligne,6,10);
      if ($isa[0] eq "-1") {
        @isa = ();
      }
      push(@isa, $parent);
    } elsif ($ligne =~ /def\: /) {
      $ligne =~ s/def\: //;
      $def = $ligne;
    }
  }
  $self->set_array_attribute("children", @childrentable);
  $self->set_array_attribute("parents", @parenttable);
  $self->set_hash_attribute("indices", %indextable);
  $self->set_hash_attribute("functionname", %functionname);
  $self->set_hash_attribute("namespace", %namespace);
  $self->set_hash_attribute("definition", %definition);
  $self->set_hash_attribute("name_id", %name_id);
}
    
################################################################

=pod

=item B<get_parents>

Returns the GOID of the parents of a specific goid


 Title    : get_parents
 Usage    : $go->get_parents(go_id)
 Returns  : the list of all parents classes of the specified class

=cut

sub get_parents {
  my ($self, $goid) = @_;
  my %indextable = $self->get_attribute("indices");
  my @parenttable = $self->get_attribute("parents");
  my @rep;
  $rep[0] = $goid;
  for (my $i = 0; $i < scalar(@rep); $i++) {
    my $index = $indextable{$rep[$i]};
    my $parents = $parenttable[$index];
    for my $j (0 .. $#{$parents}) {
      if ($parents->[$j] ne -1) {
        push(@rep, $parents->[$j]);
      }
    }
  }
  my %seen = ();
  my @uniq = ();
  foreach my $item (@rep) {
    unless ($seen{$item}) {
      # if we get here, we have not seen it before
      $seen{$item} = 1;
      push(@uniq, $item);
    }
  }
  my $repRef = \@uniq;
  return($repRef);
}
################################################################

=pod

=item B<get_children>

Returns the GOID of the children of a specific goid


 Title    : get_children
 Usage    : $go->get_children(go_id)
 Returns  : the list of all children classes of the specified class

=cut

sub get_children {
  my ($self, $goid) = @_;
  my %indextable = $self->get_attribute("indices");
  my @childrentable = $self->get_attribute("children");
  my @rep;
  $rep[0] = $goid;
  for (my $i = 0; $i < scalar(@rep); $i++) {
    my $index = $indextable{$rep[$i]};
    my $children = $childrentable[$index];
    for my $j (0 .. $#{$children}) {
      if ($children->[$j] ne -1) {
        push(@rep, $children->[$j]);
      }
    }
  }
  my %seen = ();
  my @uniq = ();
  foreach my $item (@rep) {
    unless ($seen{$item}) {
      # if we get here, we have not seen it before
      $seen{$item} = 1;
      push(@uniq, $item);
    }
  }
  my $repRef = \@uniq;
  return($repRef);
}


################################################################

=pod

=item B<get_parents>

Read the graph from a tab-delimited text file.


 Title    : exists
 Usage    : $go->exists(go_id)
 Returns  : if this goid is annotated in the go

=cut

sub exists {
  my ($self, $arg) = @_;
  my %indextable = $self->get_attribute("indices");
  my %name_id = $self->get_attribute("name_id");
  my $index = $indextable{$arg};
  my $id = $name_id{$arg};
  my $repRef = 0;
  if (defined($index) || defined($id)) {
    $repRef = 1;
  }
  return($repRef);
}
################################################################

=pod

=item B<get_namespace>

 Title    : get_namespace
 Usage    : $go->print_go_class(goid)
 Returns  : a string indicating the namespace of the specified go_class

=cut

sub get_namespace {
  my ($self, $goid) = @_;
  my %namespace = $self->get_attribute("namespace");
  my $rep = $namespace{$goid};
  if (!defined($rep)) {
    &RSAT::message::Warning("Unable to find the name space of class $goid\n") if ($main::verbose >= 2);
  }
  return($rep);
}

################################################################

=pod

=item B<get_goid>


 Title    : get_goid
 Usage    : $go->print_go_class(name)
 Returns  : a string indicating the namespace of the specified go_class

=cut

sub get_goid {
  my ($self, $name) = @_;
  my %namespace = $self->get_attribute("name_id");
  my $rep = $namespace{$name};
  if (!defined($rep)) {
    &RSAT::message::Warning("Unable to find the id of class $name\n") if ($main::verbose >= 2);
  }
  return($rep);
}
################################################################

=pod

=item B<get_name>

Read the graph from a tab-delimited text file.


 Title    : get_name
 Usage    : $go->print_go_class(goid)
 Returns  : a string indicating the namespace of the specified go_class

=cut

sub get_name {
  my ($self, $goid) = @_;
  my %functionname = $self->get_attribute("functionname");
  my $rep = $functionname{$goid};
  if (!defined($rep)) {
    &RSAT::message::Warning("Unable to find the name of class $goid\n") if ($main::verbose >= 2);
  }
  return($rep);
}


################################################################

=pod

=item B<get_definition>

Read the graph from a tab-delimited text file.


 Title    : get_definition
 Usage    : $go->print_go_class(goid)
 Returns  : a string indicating the definition of the specified go_class

=cut

sub get_definition {
  my ($self, $goid) = @_;
  my %definition = $self->get_attribute("definition");
  my $rep = $definition{$goid};
  if (!defined($rep)) {
    &RSAT::message::Warning("Unable to find the description of class $goid\n") if ($main::verbose >= 2);
  }
  return($rep);
}




return 1;

__END__


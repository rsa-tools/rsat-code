###############################################################
#
# Class to handle organisms supported in RSAT
#

package RSAT::OrganismManager;

use RSAT::GenericObject;
use RSAT::error;
use RSAT::message;
use RSAT::SequenceOnDisk;
use RSAT::GenomeFeature;
use RSAT::Index;
use RSAT::stats;
use RSAT::organism;
use Storable qw(nstore retrieve);

@ISA = qw( RSAT::GenericObject );

################################################################
## Class variables

## Fields for supported organisms
our @supported_org_fields = qw(ID name data last_update taxonomy up_from up_to genome seq_format);

## Name of the table containing the list of supported organisms
our $organism_table_name = "supported_organisms.tab";

=pod

=head1 NAME

    RSAT::OrganismManager

=head1 DESCRIPTION

Object for handling organisms supported in RSAT.

=cut



################################################################
=pod

=item B<supported_org_fields()>

Return the list of organism fields (parameters used to describe each
organism).

=cut
sub supported_org_fields {
  return (@supported_org_fields);
}

################################################################
=pod

=item B<load_supported_organisms>

Load the list of supported organisms with their parameters (taxon,
directory, upstream length, ...).

This command can be called several times in order to load several
lists of supported organisms stored in separate tables.

=cut
sub load_supported_organisms {
  my ($organism_table) = @_;
  $organism_table = $organism_table || $ENV{RSAT}."/data/".$organism_table_name;

  unless (-e $organism_table) {
      &RSAT::message::Warning("The tabular file with the list of supported organism cannot be read");
      &RSAT::message::Warning("Missing file",  $organism_table);
      return();
  }
  my ($table_handle) = &RSAT::util::OpenInputFile($organism_table);
  
  my @fields = ();
  my $linecpt = 0;
  while (my $line = <$table_handle>) {
    $linecpt++;
    next if $line =~ /^;/;
    chomp $line;
    if ($line =~ /^#/) { # Load header
      $line =~ s/#//;
      @fields = split /\t/, $line;
    } else {
      @linecp = split /\t/, $line;
      &RSAT::error::FatalError("Number of fields in the header does not correspond to the number of fields at line $linecpt\n") if (scalar (@linecp) != scalar (@fields));
      for (my $i = 1; $i < scalar @fields; $i++) {
        $main::supported_organism{$linecp[0]}->{$fields[$i]} = $linecp[$i];
      }
    }
  }
}

################################################################
=pod

=item B<export_supported_organisms>

Export the list of supported organisms with their parameters (taxon,
directory, upstream length, ...).

=cut
sub export_supported_organisms {
  my ($organism_table, @fields) = @_;
  $organism_table = $organism_table || $ENV{RSAT}."/data/supported_organisms.tab";
  my ($table_handle) = &RSAT::util::OpenOutputFile($organism_table);
  print $table_handle &supported_organism_table("header", @fields);
  &RSAT::message::Info("Exporting table with supported organisms", $organism_table) if ($main::verbose >= 1);
}

################################################################
=pod

=item B<supported_organism_table>

Return a string with a tab-delimited table containing one row per
organism and one column per field (feature of an organism).

Usage:

 &RSAT::OrganismManager::supported_organism_table($header, @fields)

=cut
sub supported_organism_table {
  my ($header, @fields) = @_;
#  &RSAT::message::Debug("&RSAT::OrganismManager::supported_organism_table()", "fields", join( ";", @fields)) if ($main::verbose >= 3);

  ## Default fields
  if (scalar(@fields) == 0) {
#    @fields = qw(name data last_update features genome seq_format taxonomy synonyms up_from up_to);
    @fields = qw(ID name data last_update taxonomy up_from up_to genome seq_format);
  }

  ## Check if the requested fields are supported
  my %supported_org_fields;
  foreach my $field (@supported_org_fields) {
    $supported_org_fields{$field} = 1;
  }
  foreach my $field (@fields) {
    unless (defined($supported_org_fields{$field})) {
      &RSAT::error::FatalError($field, "is not a valid field for &RSAT::OrganismManager::supported_organism_table()");
    }
  }

  ## Add the header row
  my $table = "";
  if ($header) {
    $table .= "#";
    $table .= join ("\t", @fields);
    $table .= "\n";
  }

  ## Add fields for each organism
  foreach my $org (sort keys %main::supported_organism) {
    $main::supported_organism{$org}->{'ID'} = $org;
    my @values = ();
    foreach my $field (@fields) {
      push @values, $main::supported_organism{$org}->{$field};
    }
    my $row = join ("\t", @values);
#    &RSAT::message::Debug($org, $row) if ($main::verbose >= 3);
    $table .= $row."\n";
  }
  return($table);
}

################################################################
=pod

=item is_supported

Indicates whether a given organism name is supported on this RSAT
site.  Return value is boolean (1=true, 0=false).

=cut
sub is_supported {
    my ($organism_name) = @_;
    if (defined($main::supported_organism{$organism_name})) {
	return 1;
    } else {
	return 0;
    }
}


################################################################
=pod

=item B<get_supported_organisms>

return a list with the IDs of the supported organisms.

=cut

sub get_supported_organisms {
  return sort keys %main::supported_organism;
}

################################################################
=pod

=item check_name

Check if the specified organism has been installed on this RSAT
site. If not, die.

=cut
sub check_name {
  my ($organism_name) = @_;
  unless ($organism_name) {
    &RSAT::error::FatalError("you should specify an organism name");
  }
  my $supported = &is_supported($organism_name);
  unless ($supported) {
    &RSAT::error::FatalError("Organism $organism_name is not supported.",
			     "Use the command supported-organisms for a list of supported organisms");
  }
}


################################################################
=pod

=item serial_file_name()

Return the name of the file containing the serialized organism.

Usage

 $organism->serial_file_name($imp_pos, $synonyms)

=cut

sub serial_file_name {
  my ($self, $imp_pos, $synonyms) = @_;
  my $serial_dir = $ENV{RSAT}."/public_html/tmp";
  my $serial_file = join ("", $self->get_attribute("name"),
			  "_imp_pos",$imp_pos,
			  "_synonyms",$synonyms,
			  "_",join("_", $self->get_attribute("feature_types")),
			  ".serial");
  return ($serial_dir."/".$serial_file);
}


################################################################
=pod

=item load_and_serialize()

Serialize the organism, i.e. store it as a binary file. The serialized
organism can be reloaded faster than with the flat files.

=cut
sub load_and_serialize {
  use Storable qw(nstore);
  my ($self, $imp_pos, $synonyms) = @_;
  $self->LoadFeatures($annotation_table, $imp_pos);
  $self->LoadSynonyms() if ($synonyms);
  my $serial_file = $self->serial_file_name($imp_pos, $synonyms);
  nstore $self, $serial_file;
  system ("chmod 777 $serial_file");
  &RSAT::message::TimeWarn("Serialized organism", $organism, $serial_file) 
    if ($main::verbose >= 3);
}

################################################################
=pod

=item is_serialized()

Return 1 if there is an up-to-date serialized version of the organism.

Return zero otherwise, i.e. if either the serialized file either does
not exist, or is older than the contig file.

=cut
sub is_serialized {
  my ($self, $imp_pos, $synonyms) = @_;

  ## Get last modifiction date of the contig file
  my $organism_name = $self->get_attribute("name");
  my $ctg_file = join("", $main::supported_organism{$organism_name}->{'data'}, "/genome/","contigs.txt");
  my ($ctg_dev,$ctg_ino,$ctg_mode,$ctg_nlink,$ctg_uid,$ctg_gid,$ctg_rdev,$ctg_size,
      $ctg_atime,$ctg_mtime,$ctg_ctime,$ctg_blksize,$ctg_blocks)
    = stat($ctg_file);

  ## Get last modifiction date of the serialized file
  my $serial_file = $self->serial_file_name($imp_pos, $synonyms);
  my ($serial_dev,$serial_ino,$serial_mode,$serial_nlink,$serial_uid,$serial_gid,$serial_rdev,$serial_size,
      $serial_atime,$serial_mtime,$serial_ctime,$serial_blksize,$serial_blocks)
    = stat($serial_file);

  ## Compare modification dates
  if (-e $serial_file) {
    if ($serial_mtime > $ctg_mtime) {
      &RSAT::message::Info("Serialized file is up-to-date", $serial_file) if ($main::verbose >= 3);
      return (1);
    } else {
      &RSAT::message::Info("Serialized file is obsolete", $serial_file) if ($main::verbose >= 3);
      return (0);
    }
  } else {
      &RSAT::message::Info("Serialized file does not exist", $serial_file) if ($main::verbose >= 3);
      return (0);
  }
}

return 1;


__END__



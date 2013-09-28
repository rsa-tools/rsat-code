###############################################################
#
# Class to handle organisms supporte in RSAT
#

package RSAT::OrganismManager;

use RSAT::util;
use RSAT::GenericObject;
use RSAT::error;
use RSAT::message;
use RSAT::SequenceOnDisk;
use RSAT::GenomeFeature;
use RSAT::Index;
use RSAT::stats;
use RSAT::organism;
use Storable qw(nstore retrieve);
use RSAT::Tree;
use RSAT::TreeNode;

@ISA = qw( RSAT::GenericObject );


################################################################
## Class variables

## Fields for supported organisms
our @supported_org_fields = qw(ID
			       name
			       data
			       last_update
			       taxonomy
			       up_from
			       up_to
			       genome
			       seq_format
			       source
			       nb
			       selected_taxon
			       taxon_depth
			       org_per_taxon
			      );

## Name of the table containing the list of supported organisms
our $organism_table_name = "supported_organisms.tab";

## Null value for undefined fields
our $null = "<NA>";


## Taxonomic tree (should be loaded only once)
our $tree;
our %count_per_taxon = ();

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
  $organism_table = $organism_table || $ENV{RSAT}."/public_html/data/".$organism_table_name;

  unless (-e $organism_table) {
      &RSAT::message::Warning("The tabular file with the list of supported organism cannot be read");
      &RSAT::message::Warning("Missing file",  $organism_table);
      return();
  }
  my ($table_handle) = &RSAT::util::OpenInputFile($organism_table);
  my @fields = ();
  my @values = ();
  my $l = 0;

  while (my $line = <$table_handle>) {
    $l++;
    next if $line =~ /^;/;
    chomp $line;
    if ($line =~ /^#/) { # Load header
      $line =~ s/#//;
      @fields = split /\t/, $line;
    } else {
#      $line =~ s|\$ENV\{RSAT\}|BOUM|g;
#      my $tmp = $ENV{RSAT};
      #      &RSAT::error::FatalError($tmp, $ENV{RSAT});
      #      $line =~ s|\$ENV\{RSAT\}|$tmp|g;
      #      $line =~ s|\/+|\/|g;
      @values = split /\t/, $line;
      &RSAT::error::FatalError("&RSAT::OrganismManager::load_supported_organisms()\n",
			       "Number of fields in the header does not correspond to the number of fields\n",
			       "file", $organism_table, "line", $l, "\n") if (scalar (@values) != scalar (@fields));
      for (my $i = 1; $i < scalar @fields; $i++) {
	my $field = $fields[$i];
	my $value = $values[$i];
	$value =~ s|\$ENV\{RSAT\}|$ENV{RSAT}|;
#	if ($value =~ /^\$ENV\{RSAT\}/) {
#	  $value = $ENV{RSAT}; #"/".$'; ## '
#	}
       $main::supported_organism{$values[0]}->{$field} = $value;
      }
    }
  }
}

# ################################################################
# ## Update one organism in the tab-delimited file
# sub UpdateConfigTab {
#   my ($org, %args) = @_;

#   ## Organism name
#   my $name = "";
#   if ($args{name}) {
#     $name = $args{name};
#   } else {
#     $name = $org;
#     $name =~ s/\_/ /;
#   }
#   $supported_organism{$organism_short_name}->{'name'} = $args{name} || $organism_full_name;

#   ## Data directory
#   $supported_organism{$organism_short_name}->{'data'} = $args{data} || $ENV{RSAT}."/public_html/data/genomes/".$org;

#   ## Update date
#   my $now = `date '+%Y/%m/%d %H:%M:%S'`;
#   $supported_organism{$organism_short_name}->{'last_update'} = $args{last_update} ||  $install_date;

#   $supported_organism{$organism_short_name}->{'source'} = $args{source} || $source;
#   ## OLIVIER SAND SHOULD CHEXK IF THIS RESTRICTION FOR ensembl IS STILL VALID
#   unless ($main::source eq 'ensembl') {
#     $supported_organism{$organism_short_name}->{'features'} = $args{features} || $outfile{features};
#     $supported_organism{$organism_short_name}->{'genome'} = $args{genome} || $outfile{genome};
#     $supported_organism{$organism_short_name}->{'seq_format'} = $args{seq_format} || "filelist";
#   }
#   $supported_organism{$organism_short_name}->{'taxonomy'} = $args{taxonomy} || $taxonomy;
#   if (defined($outfile{synonyms})) {
#     $supported_organism{$organism_short_name}->{'synonyms'} = $args{synonyms} || $outfile{synonyms};
#   }
#   $supported_organism{$organism_short_name}->{'up_to'} = $args{up_to} || $up_to;
#   $supported_organism{$organism_short_name}->{'up_from'} = $args{up_from} || $up_from;

#   if ($main::verbose >= 0) {
#     &RSAT::message::Debug("");
#   }
# #  &RSAT::message::Debug("new_org_config", $new_org_config) if ($main::verbose >= 0);

#   ## Export the updated table of supported organisms
#   &export_supported_organisms($config_table);

# }


################################################################

=pod

=item B<export_supported_organisms>

Export the list of supported organisms with their parameters (taxon,
directory, upstream length, ...).

Usage:

=over

=item Table name set automatically

  &RSAT::OrganismManager::export_supported_organisms();

=item Specify a custom table name

  &RSAT::OrganismManager::export_supported_organisms(file=>$my_organism_table_file);

=item Save a backup copy of the previous table

  &RSAT::OrganismManager::export_supported_organisms(backup=>1);


=cut
sub export_supported_organisms {
  my (%args) = @_;
  my $organism_table = $args{file} || $ENV{RSAT}."/public_html/data/supported_organisms.tab";
  my $store_backup = $args{backup} || 0;

  ## Save a copy of the previous table
  if ($store_backup) {
    $organism_table_bk =  &RSAT::util::make_temp_file($ENV{RSAT}."/public_html/data/","supported_organisms_backup", 1,0);
    &RSAT::message::TimeWarn("Storing backup copy of previous organism table", $organism_table_bk) if ($main::verbose >= 2);
    my $cmd = "rsync -uptL ".$organism_table." ".$organism_table_bk;
    $cmd .= "; gzip $organism_table_bk";
    system($cmd);
  }

  ## Write the table to a temporary file. This avoids problems if the
  ## process is interrupted during the table writing.
  $organism_table_tmp =  &RSAT::util::make_temp_file($ENV{RSAT}."/public_html/data/","supported_organisms_tmp", 1,0);
  &RSAT::message::TimeWarn("Storing updated organism table to temporary file", $organism_table_tmp) if ($main::verbose >= 2);
  my ($table_handle) = &RSAT::util::OpenOutputFile($organism_table_tmp);

  print $table_handle &supported_organism_table("header", 1, @fields);
  &RSAT::message::Warning("Make sure that the file RSA.config does not load the old format file",$ENV{RSAT}."/public_html/data/supported_organisms.pl") if ($main::verbose >= 3);

  ## Rename the updated table to make it effective
  system("mv ".$organism_table_tmp." ".$organism_table);
  &RSAT::message::Info("Exported supported organisms", $organism_table) if ($main::verbose >= 1);
}

################################################################

=pod

=item B<supported_organism_table>

Return a string with a tab-delimited table containing one row per
organism and one column per field (feature of an organism).

Usage:

  &RSAT::OrganismManager::supported_organism_table($header, $relative_path, $taxon, @fields)

Arguments:

=over

=item header

When the argument is not null, the first row of the table is a header
indicating the column contents.

=item relative_path

When non null, the data paths are given relative to the $RSAT
directory. Otherwise, absolute paths are returned.

=item taxon

Restrict only return organisms belonging to a given taxon.

=back

=cut
sub supported_organism_table {
  my ($header,$relative_path, $taxon, $depth, @fields) = @_;
  my $table = "";

  &RSAT::message::Debug("&RSAT::OrganismManager::supported_organism_table()", "taxon: ".$taxon, "fields", join( ";", @fields)) 
    if ($main::verbose >= 3);

  ## Default fields
  if (scalar(@fields) == 0) {
#    @fields = qw(name data last_update features genome seq_format taxonomy synonyms up_from up_to);
#    @fields = qw(ID name data last_update taxonomy up_from up_to genome seq_format);
    @fields = @supported_org_fields;
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
  if ($header) {
    $table .= "#";
    $table .= join ("\t", @fields);
    $table .= "\n";
  }

  ## Select the organisms
  my @selected_organisms = ();
  if ($taxon) {
    @selected_organisms = &GetOrganismsForTaxon($taxon, $depth);
  } else {
    @selected_organisms = sort keys %main::supported_organism;
    if ($depth != 0) {
      @selected_organisms = &OneOrgPerTaxonomicDepth($depth, @selected_organisms);
    }
  }


  ## Add fields for each organism
  my $n = 0;
  foreach my $org (@selected_organisms) {
    $n++;
    $main::supported_organism{$org}->{'ID'} = $org;
    $main::supported_organism{$org}->{'nb'} = $n; ## Reset nb in case of filtering
    my @values = ();
    foreach my $field (@fields) {
      my $value = $null;
      if (defined($main::supported_organism{$org}->{$field})) {
	$value = $main::supported_organism{$org}->{$field};
	if ($relative_path) {
	  $value =~ s|$ENV{RSAT}|\$ENV\{RSAT\}\/|;
	  $value =~ s|\/+|\/|g;
	}
#      } elsif ($field eq "nb") {
#	$value = $n;
      } else {
	$value = $null;
	&RSAT::message::Warning("Field", $field, "has no value for organism", $org);
      }
      push @values, $value;
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
#### Check if an organism is supported on the current installation,
#### and open streams to read its genome sequence.
####
#### Usage
#### -----
#### Automatic selection of genome and feature file : 
####    &CheckOrganism($organism_name); 
####
#### Manual specification of input files :
#### &CheckOrganism($organism_name, 
####                $annotation_table, 
####                $input_sequence_format);
sub CheckOrganism {
    my ($organism_name, $annotation_table, $input_sequence_file, $input_sequence_format) = @_;
    my $organism_object = new RSAT::organism();
    $organism_object->OpenContigs($organism_name, 
				  $annotation_table, 
				  $input_sequence_file, 
				  $input_sequence_format);
    return $organism_object;
}


################################################################
## Collect all organisms belonging to a given taxon
sub GetOrganismsForTaxon {
  my ($taxon, $depth) = @_;
  my @organisms = ();
  unless ($tree) {
    $tree = new RSAT::Tree();
  }
  $tree->LoadSupportedTaxonomy("Organisms", \%main::supported_organism);
  my $node = $tree->get_node_by_id($taxon);
  if ($node){
    @organisms = $node->get_leaves_names();

    ## Cut the taxonomic tree by selecting only one organism for each
    ## taxon at a given depth of the taxonomic tree.
    if ($depth != 0) {
      @organisms = &OneOrgPerTaxonomicDepth($depth, @organisms);
    }
  }else{
    &RSAT::error::FatalError("Taxon $taxon is not supported\n");
  }
  return @organisms;
}

################################################################
## Select one organism per taxonomic depth.
##
## Positive depth values: from root (level 1) to leaves.
##
## Negative depth values: from leaves (level -1) to root. Beware:
## level 0 is the species -> level -1 is the last taxonomy level
sub OneOrgPerTaxonomicDepth {
  my ($depth, @organisms) = @_;
  my @selected_organisms = ();
  foreach my $org (@organisms) {
    my $taxonomy = $main::supported_organism{$org}->{'taxonomy'};
    my @taxonomy = split (/;\s*/, $taxonomy);
    my $max_depth = scalar(@taxonomy);
    my $selected_depth;
    if ($depth < 0) {
      $selected_depth = $max_depth + $depth;
    } else {
      $selected_depth = &RSAT::stats::min($depth, $max_depth);
    }
    my $selected_taxon = $taxonomy[$selected_depth -1];
    $main::supported_organism{$org}->{'taxon_depth'} = $selected_depth;
    $main::supported_organism{$org}->{'selected_taxon'} = $selected_taxon;
    unless (defined($count_per_taxon{$selected_taxon})) {
      push @selected_organisms, $org;
    }
    $count_per_taxon{$selected_taxon}++;
    &RSAT::message::Debug("&SelectOneOrgPerTaxon()",
			  $org,
			  $taxonomy,
			  "depth=".$depth,
			  "max_depth=".$max_depth, $selected_taxon, $count_per_taxon{$selected_taxon})
      if ($main::verbose >= 5);
  }

  foreach my $org (@selected_organisms) {
    $main::supported_organism{$org}->{'org_per_taxon'} = $count_per_taxon{$main::supported_organism{$org}->{'selected_taxon'}};
  }
  return (@selected_organisms);
}

################################################################
## Check if there is at least one supported organism belonging to the
## specified taxon on this RSAT server.
sub CheckTaxon {
  my ($query_taxon, $warn_only) = @_;

  ## Load the taxonomic tree of organisms supported in RSAT
  unless ($tree) {
    $tree = new RSAT::Tree();
  }

  &RSAT::message::Info("CheckTaxon()", $query_taxon, $tree) if ($main::verbose >= 0);

  unless ($tree->get_attribute("is_loaded")) {
    my @supported_organisms = sort keys (%supported_organism);
    $tree->LoadSupportedTaxonomy("Organisms", \%main::supported_organism, 1);
  }

  my @organisms = $tree->get_node_descendents_names($query_taxon, "DFS", "leaf");

  &RSAT::message::Info("CheckTaxon", $query_taxon, scalar(@organisms), "organisms") if ($main::verbose >= 4);
  if (scalar(@organisms) > 0) {
    return @organisms;
  } else {
    my $message = join("\t", "There is no supported organism for taxon", $query_taxon, "\nType supported-organisms -format full for a list of supported taxa");
    if ($warn_only) {
      &RSAT::message::Warning($message);
    } else {
      &RSAT::error::FatalError($message);
    }
  }
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
    &RSAT::error::FatalError("You should specify an organism name");
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



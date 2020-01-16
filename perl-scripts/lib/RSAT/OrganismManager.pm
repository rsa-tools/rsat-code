###############################################################
##
## Class to handle organisms supported in this instance of RSAT.
##


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
			       taxid
			       source
			       last_update
			       nb
			       seq_format
			       up_from
			       up_to
			       taxonomy
			       data
			       genome
			       genome_assembly
			       genome_version
			       download_date
			       variant_available
			       variant_source
			       path_to_variant_files
                   blast_available
			      );

#			       selected_taxon
#			       taxon_depth
#			       org_per_taxon

## Name of the table containing the list of supported organisms
our $organism_table_name = "supported_organisms.tab";
#our $organism_table_name= "supported_organisms2.tab"; #2019 version with added columns for consistency between organisms

## Flag indicating whether the organism table has been loaded
our $organism_table_loaded = 0;

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
        if ($main::verbose >= 2) {
            &RSAT::message::Warning("The tabular file with the list of supported organism cannot be read");
            &RSAT::message::Warning("Missing file",  $organism_table);
        }
        return();
    }
    my ($table_handle) = &RSAT::util::OpenInputFile($organism_table);
    my @fields = ();
    my $l = 0;
    
    while (my $line = <$table_handle>) {
        my @values = ();
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
                #	my $value = &RSAT::util::hide_RSAT_path($values[$i]);
                $main::supported_organism{$values[0]}->{$field} = $value;
            }
        }
    }

    ## Activate the "loaded" flag
    $organism_table_loaded = 1;
}


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
  my $organism_table = $args{file} || $ENV{RSAT}."/public_html/data/".$organism_table_name;
  my $store_backup = $args{backup} || 0;
  my $source = $args{source} || "";
  my $taxon = $args{taxon} || "";
  my $depth = $args{depth} || 0;

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
  &RSAT::message::TimeWarn("Storing updated organism table to temporary file", $organism_table_tmp) if ($main::verbose >= 3);
  my ($table_handle) = &RSAT::util::OpenOutputFile($organism_table_tmp);

  print $table_handle &supported_organism_table("header", 1, $source, $taxon, $group, $depth, 0, 0, @fields);

  ## Rename the updated table to make it effective
  system("mv ".$organism_table_tmp." ".$organism_table);
  &RSAT::message::Info("Exported supported organisms", $organism_table) if ($main::verbose >= 2);
}

=pod

=item B<UniquePerTaxon>

Select a single organism per species or per genus (the other
taxonomical levels are not supported because we do not dispose of the
information about taxonomic depth). 

Usage:

  @organisms_filtered = &RSAT::OrganismManager::UniquePerTaxon($taxon, @organisms);

Arguments:

=over

=item taxon

Taxonomical level at which the filtering must be done. Supported:
genus, species.

=item organisms

List of organisms

=back

=cut

sub UniquePerTaxon {
  my ($taxon, @organisms) = @_;
  unless (($taxon eq "species") || ($taxon eq "genus")) {
    &RSAT::error::FatalError("&RSAT::OrganismManager::UniquePerTaxon()", $taxon, "is not a valid taxon. Supported: genus, species");
  }
  &RSAT::message::Info("Filtering organisms per", $taxon) if ($main::verbose >= 4);
  my @filtered_organisms = ();
  my %orgs_per_genus = ();
  my %orgs_per_species = ();
  foreach my $org (@organisms) {
    my ($genus, $species) = split("_", $org);
    $species = $genus."_".$species;
    $orgs_per_genus{$genus}++;
    $orgs_per_species{$species}++;
    next if (($taxon eq "genus") && ($orgs_per_genus{$genus} > 1));
    next if (($taxon eq "species") && ($orgs_per_species{$species} > 1));
#      &RSAT::message::Debug($org, $genus, $species, $orgs_per_species{$species}, $orgs_per_genus{$genus}) if ($main::verbose >= 10);
    push @filtered_organisms, $org;
  }
  return(@filtered_organisms);
}

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
  my ($header,$relative_path, $source, $taxon, $group, $depth, $with_variant, $with_blast, @fields) = @_;
  my $table = "";

  ## Check arguments
  &RSAT::error::FatalError("&RSAT::OrganismManager::supported_organism_table()",
			   "arguments 'taxon' and 'group' are mutually exclusive")
      if (($taxon) && ($group));
  
  &RSAT::message::Debug("&RSAT::OrganismManager::supported_organism_table()", 
			"taxon: ".$taxon, 
			"group: ".$group, 
			"depth: ".$depth, 
			"fields", join( ";", @fields)) 
    if ($main::verbose >= 5);

  ## Default fields
  if (scalar(@fields) == 0) {
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
  } elsif (($group) && ($group ne "None")) {
    @selected_organisms = &GetOrganismsForGroup($group, $depth);
    push @selected_organisms, &get_demo_organisms(@group_specificity);

    ## Filter out duplicates between selected organisms and demo organisms
    my %selected = ();
    foreach my $org (@selected_organisms) {
      $selected{$org} = 1;
    }
    @selected_organisms = sort(keys(%selected));
    
  } else {
    @selected_organisms = sort keys %main::supported_organism;
    if ($depth != 0) {
      @selected_organisms = &OneOrgPerTaxonomicDepth($depth, @selected_organisms);
    }
  }

  &RSAT::message::Debug("Group", $group, "Selected organisms: ", scalar(@selected_organisms)) if ($main::verbose >= 5);

  ## Select unique organisms per genus or species if required
  @selected_organisms = &RSAT::OrganismManager::UniquePerTaxon("species", @selected_organisms) if ($main::unique_species);
  @selected_organisms = &RSAT::OrganismManager::UniquePerTaxon("genus", @selected_organisms) if ($main::unique_genus);

  ## Add fields for each organism
  my $n = 0;
  foreach my $org (@selected_organisms) {
    $n++;

    ## Filter by source
    if ($source) {
      next unless lc($main::supported_organism{$org}->{'source'}) eq lc($source);
    }

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
      } else {
	$value = $null;
	&RSAT::message::Warning("Field", $field, "has no value for organism", $org) if ($main::verbose >= 3);
      }
      push @values, $value;
    }
      if(($with_variant == 0 && $with_blast == 0) || ($with_variant == 1 && $main::supported_organism{$org}->{'variant_available'} eq "1") || ($with_blast == 1 && $main::supported_organism{$org}->{'blast_available'} eq "1")){
          my $row = join ("\t", @values);
          #    &RSAT::message::Debug($org, $row) if ($main::verbose >= 3);
          $table .= $row."\n";
      }
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

    &load_supported_organisms() unless ($organism_table_loaded);

    if (defined($main::supported_organism{$organism_name})) {
	return 1;
    } else {
	return 0;
    }
}


=pod

=item B<get_supported_organisms>

Return a list with all supported organisms on this RSAT instance.

=cut
sub get_supported_organisms {
  return sort keys %main::supported_organism;
}

=pod

=item B<get_supported_organisms_web>

Return a list with all supported organisms on the web server for this
RSAT instance. This can differ from the complete list of organisms
supported on the command line, because some servers are dedicated to
specific taxonomic groups.

The taxonomic range supported on a web interface is defined by the
property "group_specificity" in the file $RSAT/RSAT_config.props.

=cut
sub get_supported_organisms_web {
  &RSAT::message::Info("Getting supported organisms on web site", $RSAT{rsat_www}) if ($main::verbose >= 3);

  my  @selected_organisms = ();

  ## Multiple groups can be specified
  my @group_specificity = split(/,/, $ENV{group_specificity});

  if (scalar(@group_specificity) >= 1) {
    ## Collect organisms for each group specificity
    foreach my $group_specificity (@group_specificity) {
      $group_specificity = ucfirst(lc($group_specificity));
      push @selected_organisms , &RSAT::OrganismManager::GetOrganismsForGroup($group_specificity);
      if (scalar(@selected_organisms) < 1) {
	if (lc($ENV{group_specificity}) eq "none") {
	  &RSAT::error::FatalError("No organism supported on this server");
	} else {
	  &RSAT::error::FatalError("No organism supported on this server for group", $ENV{group_specificity});
	}
      }
    }
    push @selected_organisms, &get_demo_organisms(@group_specificity);

  } else {
    @selected_organisms = &RSAT::OrganismManager::get_supported_organisms();
  }


  @selected_organisms = &RSAT::util::sort_unique(@selected_organisms);
  return (@selected_organisms);
}


######### For menu
sub get_supported_organisms_formenu {
    my  @selected_organisms = ();
    
    ## Multiple groups can be specified
    my @group_specificity = split(/,/, $ENV{group_specificity});
    
    if (scalar(@group_specificity) >= 1) {
        ## Collect organisms for each group specificity
        foreach my $group_specificity (@group_specificity) {
            $group_specificity = ucfirst(lc($group_specificity));
            push @selected_organisms , &RSAT::OrganismManager::GetOrganismsForGroup($group_specificity);
            if (scalar(@selected_organisms) < 1) {
                return ();
            }
        }
        push @selected_organisms, &get_demo_organisms(@group_specificity);
        
    } else {
        @selected_organisms = &RSAT::OrganismManager::get_supported_organisms();
    }
    
    
    @selected_organisms = &RSAT::util::sort_unique(@selected_organisms);
    return (@selected_organisms);
}

################################################################
## Return the list of organisms required for Web demos
sub get_demo_organisms {
    my (@group_specificity) = @_;
    foreach my $group (@group_specificity) {
	$group_specificity{$group} = 1;
    }
    my @selected_organisms = ();
    unless ($group_specificity{"Fungi"}) {
	if (&is_supported("Saccharomyces_cerevisiae")) {
	    ## old nomenclature, probably outdated, I should define an alias for the model yeast strain
	    push @selected_organisms, "Saccharomyces_cerevisiae";
	} else {
	    push @selected_organisms, &GetOrganismsForTaxon("Saccharomyces");
	}
    }
    unless ($group_specificity{"Metazoa"}) {
	if (&is_supported("Drosophila_melanogaster")) {
	    push @selected_organisms, "Drosophila_melanogaster"; ## Required for matrix-scan demo. SHOULD WE IMPOSE THIS GENOME JUST FOR THAT, OR HAVE GROUP-SPEFICIC DEMOS ?
	}
    }
    unless (($group_specificity{"Bacteria"}) || ($group_specificity{"Prokaryotes"})) {
	if (&is_supported("Escherichia_coli_K_12_substr__MG1655_uid57779")) {
	    push @selected_organisms, "Escherichia_coli_K_12_substr__MG1655_uid57779";
	}elsif (&is_supported("Escherichia_coli_GCF_000005845.2_ASM584v2")){
	    push @selected_organisms, "Escherichia_coli_GCF_000005845.2_ASM584v2";
	} else {
	    push @selected_organisms, &GetOrganismsForTaxon("Escherichia");
	}
	
    }
    return (@selected_organisms);
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
  my ($taxon, $depth, $die_if_noorg) = @_;
  &RSAT::message::Info("Collecting organisms for taxon", $taxon, "depth: ".$depth) if ($main::verbose >= 4);

  my @organisms = ();

  ## Load the taxonomy of the organisms supported on this RSAT
  ## instance.
  unless ($tree) {
    $tree = new RSAT::Tree();
  }
  $tree->LoadSupportedTaxonomy("Organisms", \%main::supported_organism);

  ## Identify the tree node corresponding to the query taxon
  my $node = $tree->get_node_by_id($taxon);


  if ($node){

    ## Get all organisms belonging to the query taxon, i.e. the leaves
    ## descending from the selected tree node.
    @organisms = $node->get_leaves_names();
    &RSAT::message::Debug("GetOrganismsForTaxon()", $taxon, scalar(@organisms)) if ($main::verbose >= 5);

    &RSAT::message::Debug("node:", $node->get_attribute("id"), "nb organisms:", scalar(@organisms)) if ($main::verbose >= 5);
#  die("HELLO\n");

    ## If depth argument has been specified, cut the taxonomic tree by
    ## selecting only one organism for each taxon at a given depth of
    ## the taxonomic tree.
    if (defined($depth) && ($depth != 0)) {
      @organisms = &OneOrgPerTaxonomicDepth($depth, @organisms);
    }
  } else {
    $message = join ("\t", "Taxon", $taxon, "is not supported on server", $ENV{rsat_site});
    if ($die_if_noorg) {
      &RSAT::error::FatalError($message);
    } #else {
      #&RSAT::message::TimeWarn($message);
      #}
  }
  
  ## Select unique organisms per genus or species if required
  @organisms = &RSAT::OrganismManager::UniquePerTaxon("species", @organisms) if ($main::unique_species);
  @organisms = &RSAT::OrganismManager::UniquePerTaxon("genus", @organisms) if ($main::unique_genus);

  @organisms = &RSAT::util::sort_unique(@organisms);
  &RSAT::message::Info("Collected",scalar(@organisms),"organisms for taxon", $taxon) if ($main::verbose >= 3);
  return(@organisms);
}


################################################################
## Collect all organisms belonging to a given group (in the sense of
## EnsemblGenomes).
##
## Supported groups: Fungi, Prokaryotes, Bacteria, Archaea, Protists,
## Metazoa, Plants.
##
## Note that some of the "groups" correspond to a specific taxon
## defined by its systematic name (e.g. Metazoa, Fungi) or by its
## common name (Plants, Prokaryotes), whilst others are defined
## according ot the common usage (e.g. Protists) but do not properly
## speaking correspond to a taxonomic group. These non-taxonomic
## groups are converted as follows:
##
## - "Protists" is converted to 
##   "(Eukaryota NOT (Metazoa OR Fungi OR Viridiplantae)) OR EnsemblProtists"
## - "Plants" is converted to Viridiplantae
## - "Prokaryotes" is converted to "Bacteria OR Archaea" 
sub GetOrganismsForGroup {
  my ($group_specificities) = @_;

  &RSAT::message::Debug("GetOrganismsForGroup", $group_specificities) if ($main::verbose >= 4);
  my @selected_organisms = ();
  my @specific_taxa = ();

  
  my @supported_groups  = qw(Fungi
                             Prokaryotes
                             Bacteria
                             Archaea
                             Metazoa
                             Protists
                             Fungi
                             Plants
                             Teaching
                             None
                             );
  my $supported_groups = (join ", ", @supported_groups);

  foreach my $group_specificity (split(",", $group_specificities)) {
    ## Convert "groups" to corresponding taxa
    &RSAT::message::Debug("Converting groups to taxa", $group_specificity) if ($main::verbose >= 4);
    if ($group_specificity eq "Fungi") {
      push @specific_taxa, "Fungi";
    } elsif ($group_specificity eq "Prokaryotes") {
      push @specific_taxa, "Bacteria", "Archaea";
    } elsif ($group_specificity eq "Bacteria") {
      push @specific_taxa, "Bacteria";
    } elsif ($group_specificity eq "Archaea") {
      push @specific_taxa, "Archaea";
    } elsif ($group_specificity eq "Metazoa") {
      push @specific_taxa, "Metazoa";
    } elsif ($group_specificity eq "Protists") {
      ## Very tricky way to specify "Protists", which is not a
      ## taxonomic group: I select all organisms that are neither
      ## Metazoa nor Fungi
      my %non_protist = ();
      push @selected_organisms, &GetOrganismsForTaxon("EnsemblProtists");
      my @eukaryotes = &GetOrganismsForTaxon("Eukaryota");
      my @metazoa = &GetOrganismsForTaxon("Metazoa");
      my @fungi = &GetOrganismsForTaxon("Fungi");
      my @plants = &GetOrganismsForTaxon("Viridiplantae");
      my @non_protists = (@metazoa, @fungi, @plants);
      &RSAT::message::Info("Collected", scalar(@eukaryotes),"eukaroytes",
			   "\n", "Discarding non-protist eukaryotes: ", scalar(@non_protists), 
			   "\n", scalar(@metazoa), "Metazoa",
			   "\n", scalar(@fungi), "Fungi",
			   "\n", scalar(@plants), "Plants"
	  ) if ($main::verbose >= 5);
      foreach my $org (@non_protists) {
	$non_protist{$org} = 1;
      }
      foreach my $org (@eukaryotes) {
	if ($non_protist{$org}) {
	  &RSAT::message::Debug("Discarding non-protist", $org) if ($main::verbose >= 10);
	} else {
	  &RSAT::message::Debug("Adding protist", $org) if ($main::verbose >= 10);
	  push (@selected_organisms, $org);
	}
      }

    } elsif ($group_specificity eq "Plants") {
      push @specific_taxa, "Viridiplantae";

    } elsif ($group_specificity eq "Teaching") {

      ## Define a list of model organisms
      my @model_organisms = qw(
                             Arabidopsis_thaliana.TAIR10.29
                             Bacillus_subtilis_168_uid57675
                             Drosophila_melanogaster
                             Escherichia_coli_K_12_substr__MG1655_uid57779
                             Homo_sapiens_GRCh37
                             Homo_sapiens_GRCh38
                             Mus_musculus_GRCm38
                             Saccharomyces_cerevisiae
                            );

      ## Check that each organism is properly instaled.
      foreach my $org (@model_organisms) {
	if ($main::supported_organism{$org}) {
	  push @selected_organisms, $org
	}
      }
    } elsif ($group_specificity eq "None") {
      @selected_organisms = &get_supported_organisms();

    } else {
      &RSAT::error::FatalError($group_specificity, "Invalid group specificity. Supported groups: ", $supported_groups);
    }
  }
  
  ## Add organisms from selected taxa
  foreach my $taxon (@specific_taxa) {
    push @selected_organisms, &GetOrganismsForTaxon($taxon);
  }
  @selected_organisms = &RSAT::util::sort_unique(@selected_organisms);
  return(@selected_organisms);
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
    my @taxonomy = split (/\s*;\s*/, $taxonomy);
#    die "HELLO";
#    &RSAT::message::Debug("&RSAT::OrganismManager::OneOrgPerTaxonomicDepth", $org, "taxonomy", join("::", @taxonomy)) if ($main::verbose >= 2);
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

  &RSAT::message::Info("RSAT::OrganismManager::CheckTaxon()", $query_taxon) if ($main::verbose >= 3);

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
  my ($organism_name, $no_die) = @_;

  unless ($organism_name) {
    if ($no_die) {
      &RSAT::message::Warning("You should specify an organism name");
    } else {
      &RSAT::error::FatalError("You should specify an organism name");
    }
  }
  my $supported = &is_supported($organism_name);
  unless ($supported) {
    my $message = "Organism $organism_name is not supported. Run the command supported-organisms";
    if ($no_die) {
      &RSAT::message::Warning($message);
      return(0);
    } else {
      &RSAT::error::FatalError($message);
    }
  }
  return(1);
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


################################################################
## Name of the expected oligo frequency file, 
## given the organism, oligo length and background model
##
## Usage:
##
## For an organism-specific background model
##    my $background_file = &ExpectedFreqFile($organism_name, $oligo_length,
##				$background,str=>$str,noov=>$overlap,type=>$type,
##                              warn=>0|1,nowarn=>0|1,
##                              taxon=>0);
##
## For a taxon-wide background model
##    my $background_file = &ExpectedFreqFile($taxon_name, $oligo_length,
##				$background,str=>$str,noov=>$overlap,type=>$type,
##                              warn=>0|1,nowarn=>0|1,
##                              taxon=>1);
 sub ExpectedFreqFile {
  my ($org_or_taxon, $oligo_length, $background, %args) = @_;
  &RSAT::message::Info("RSAT::OrganismManager::ExpectedFreqFile()", $org_or_taxon) if ($main::verbose >= 3);
  my $exp_freq_file = "";

  my $type = $args{type} || "oligo";
  my $str = $args{str} || "-2str";
  my $overlap = $args{noov} || "-noov";
  my $overlap_suffix;
  my $taxon = $args{taxon} || 0;
  if ($overlap eq "-noov") {
    $overlap_suffix = "-noov";
  } else { 
    $overlap_suffix = "-ovlp";
  }
  

  #### Check organism or taxon
  my $data_dir = "";
  if ($taxon) {
    &CheckTaxon($org_or_taxon, 1);
    $data_dir = $ENV{RSAT}."/public_html/data";
    $data_dir .= "/taxon_frequencies/".$org_or_taxon."/";
    $data_dir =~ s/\s/_/g;
  } else {
#    &CheckOrganismName($org_or_taxon);
    &check_name($org_or_taxon);
    $data_dir = $main::supported_organism{$org_or_taxon}->{'data'};
    $data_dir .= "/oligo-frequencies/";
  }

  #### check oligo length
  unless ($oligo_length) {
    &RSAT::error::FatalError("You must specify the oligonucleotide length for the background model (expected frequency file)");
  }

  unless (&RSAT::util::IsNatural($oligo_length)) {
    &RSAT::error::FatalError("Invalid oligonucleotide length : must be a natural number");
  }

  #### check background
  my %supported_bg = ('upstream'=>1,
		      'upstream-noorf'=>1,
		      'upstream-rm'=>1,
		      'upstream-noorf-rm'=>1,
		      'protein'=>1,
#		      'intergenic'=>1,
		      'genomic'=>1
		     );

  my @sequence_types = ('upstream_mRNA','upstream_CDS','firstintron','intron','utr');
  my @maskcoding = ('-maskcoding');
#  my @noorf = ('-noorf','-nogene','');
  my @maskrepeats = ('-rm','');

  foreach my $sequence_type(@sequence_types) {
      foreach my $maskcod(@maskcoding) {
	  foreach my $maskrep(@maskrepeats) {
	      my $bg_type = $sequence_type.$maskcod.$maskrep;
		  $supported_bg{$bg_type} = 1;
	  }
      }
  }

  $supported_bg = join ",", sort keys %supported_bg;
  unless ($supported_bg{$background}) {
    &RSAT::error::FatalError($background, "is not a supported background model type. Supported: ", $supported_bg);
  }

  ## Residue type
  my $residue_type;
  if ($background eq 'protein') {
    $residue_type = "pept";
  } else {
    $residue_type = "nt";
  }

  ## Oligo length
  if ($type eq 'oligo') {
    if ($background eq "protein") {
      &RSAT::error::FatalError($oligo_length, "is not a supported length for oligopeptides, should be comprised between 1 and 3")
	unless (($oligo_length >= 1) && ($oligo_length <= 3));
    } else {
      &RSAT::error::FatalError($oligo_length, "is not a supported length for oligonucleotides, should be comprised between 1 and 8")
	unless (($oligo_length >= 1) && ($oligo_length <= 8));
    }
  } elsif ($type eq 'dyad') {
    &RSAT::error::FatalError($oligo_length, "is not a supported length for monads, should be comprised between 1 and 3")
      unless (($oligo_length >= 1) && ($oligo_length <= 3));
  }

  $exp_freq_file = $data_dir;
  if ($type eq "dyad") {
    ## Dyad frequencies
    $exp_freq_file = join("",
			  $data_dir,
			  "dyads_",$oligo_length,"nt_sp0-20_",
			  $background,
			  "_",
			  $org_or_taxon,
			  $overlap_suffix.$str,
			  ".freq");

  } elsif ($background eq "upstreamL") {
    ## ???
    $exp_freq_file .= "_".$org_or_taxon;
    $exp_freq_file .= "_allup500";
    $exp_freq_file .= $oligo_length.$residue_type;
    $exp_freq_file .= $str unless ($background eq "protein");
    $exp_freq_file .= $overlap_suffix;
    $exp_freq_file .= "_freq.tab";

  } else {
    ## Oligonucleotide frequencies
    $exp_freq_file .= $oligo_length.$residue_type;
    $exp_freq_file .= "_".$background;
    $exp_freq_file .= "_".$org_or_taxon;
    $exp_freq_file .= $overlap_suffix;
    $exp_freq_file .= $str unless ($background eq "protein");
    $exp_freq_file .= ".freq";
  }

  ## Check the existence of the exp freq file
  unless ((-e $exp_freq_file) ||
	  (-e $exp_freq_file.".gz")) {
    if ($args{warn}) {
      &RSAT::message::Warning("Cannot find expected frequency file $exp_freq_file");
    } elsif ($args{nowarn}) {
      ## Do nothing
    } else {
      &RSAT::error::FatalError("Cannot find expected frequency file $exp_freq_file");
    }
  }
  &RSAT::message::Info("Expected frequency file", $exp_freq_file) if ($main::verbose >= 4);
  return($exp_freq_file);
}


return 1;


__END__



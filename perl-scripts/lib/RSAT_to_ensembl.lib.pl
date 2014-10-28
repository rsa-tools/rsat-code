#!/usr/bin/perl

package main;

############################################################################
############################################################################
########################## FTP ENSEMBL FONCTION ############################
#  our $ensembl_url = "ftp://ftp.ensemblgenomes.org/pub/protists/current/fasta/";

################################################################
## Define global variables
## 
## Note: this is not very clean. I (JvH) should pass these variables
## in the appropriate methods and check that everything still works
## fine.
our $ensembl_version_safe = $ENV{ensembl_version_safe} || 75;
our $ensemblgenomes_version_safe = $ENV{ensemblgenomes_version_safe} || 22;

################################################################
## Functions
################################################################

################################################################
## Return the safe version of ensembl or ensembl_genomes
## (this version should be defined in RSAT_config.props).
sub get_ensembl_version_safe {
  my ($db) = @_;

  if ($db eq "ensembl") {
    return $ensembl_version_safe;
  }

  elsif ($db eq "ensembl_genomes") {
    return $ensemblgenomes_version_safe;
  }
}

################################################################
## Return the ensembl version to be used.
sub get_ensembl_version {

  ## TEMPORARY: for the time being, the version is only used by
  ##  reference to ensembl, not ensemblgenomes, even if the queried db is
  ##  ensembl_genomes. I (JvH) should cliarify this.  
  my $db = "ensembl";
  # my ($db) = @_;

  my $ensembl_version = $ENV{ensembl_version} || 76;
  my $ensemblgenomes_version = $ENV{ensemblgenomes_version} || 23;  
  ## By preference, return the version defined in the property file
  ## (RSAT_config.props) for this database)
  $prop_name = $db."_version";
  if (defined($ENV{$prop_name})) {
    return ($ENV{$prop_name});
  } else {
    return &get_ensembl_version_latest($db);
  }
}


################################################################
## Get the latest ensembl version for a species from Ensembl ftp site.
##
## THIS IS REALLY TRICKY, I (JvH) should see with Ensembl if their API
## includes a way to get the current release.
sub get_ensembl_version_latest {
  my ($db) = @_;
  my $latest_ensembl_release = 0;
  my $ftp = &Get_ftp($db);

  &RSAT::message::TimeWarn("Getting ensembl version from FTP site", $ftp) 
      if ($main::verbose >= 2);

  my @available_releases = qx{wget -S --spider $ftp 2>&1};

#  &RSAT::message::Debug("Available releases", join ("\n\t", @available_releases)) if ($main::verbose >= 10);

  foreach my $release (@available_releases) {
    next if ($release =~ /current/);
    next unless ($release =~ /^[dl].*[^\.]release/);
    chomp($release);
    my @token = split("release-",$release);
    $latest_ensembl_release = $token[-1] if ($latest_ensembl_release < $token[-1]);
  }

  &RSAT::message::Info("Latest ensembl release", $latest_ensembl_release) if ($main::verbose >= 0);

  return ($latest_ensembl_release);
}

################################################################
## Check that the user-selected version is correct.  This version can
## neither be lower than the safe version, nor higher than the latest
## version.
sub check_ensembl_version {
  my ($db, $ensembl_version) = @_;

  ## TEMPORARY: for the time being, the version is only used by
  ##  reference to ensembl, not ensemblgenomes, even if the queried db is
  ##  ensembl_genomes. I (JvH) should cliarify this.  
  $db = "ensembl";

  my $safe_ensembl_version = &get_ensembl_version_safe($db);
  my $latest_ensembl_version = &get_ensembl_version_latest($db);

  ## Report checking parameters
  if ($main::verbose >= 2) {
    &RSAT::message::Info("Latest ensembl version", $latest_ensembl_version);
    &RSAT::message::Info("Safe ensembl version", $safe_ensembl_version);
    &RSAT::message::Info("Selected ensembl version", $ensembl_version);
  }

  if ($ensembl_version eq "safe") {
    ## Automatic selection of the safe version
    $ensembl_version = $safe_ensembl_version;
    
  } elsif ($ensembl_version eq "latest") {
    ## Automatic selection of the latest version
    $ensembl_version = $latest_ensembl_version;
    
  } else {

    ## Check that selected version is not smaller than safest one (defined in RSAT_config.props)
    &RSAT::error::FatalError($ensembl_version, 
			     "Invalid  version for", $db, 
			     "Minimun safe version: ".$safe_ensembl_version) 
	if ($ensembl_version < $safe_ensembl_version);

    ## Check that selected version is not higher than the latest supported version
    &RSAT::error::FatalError($ensembl_version, 
			     "Invalid  version for", $db, 
			     "Latest version: ".$latest_ensembl_version) 
	if ($ensembl_version > $latest_ensembl_version);
  }

  return(1);
}

################################################################
## Get the Ensembl FTP site, according to the selected DB
sub Get_ftp {
  my ($db) = @_;
  if ($db eq "ensembl") {
    return "ftp://ftp.ensembl.org/pub/";
  } elsif ($db eq "ensembl_genomes") {
    return "ftp://ftp.ensemblgenomes.org/pub/";
  } else {
    &RSAT::error::FatalError($db, "Is not a valid name for Ensembl databases. Supported: ensembl, ensembl_genomes.");
  }
}

################################################################
## Get the URL of the ftp site to download the genome fasta files.
## This is really tricky, but it is the simplest way we found to
## download genome sequences.
sub Get_fasta_ftp {
  my ($db,$ensembl_version) = @_;
  
  if ($db eq "ensembl") {                                                    ## Ensembl
    if ($ensembl_version < 47) {                                             # Version  1 to 46
      return ();
    } else {                                                                 # Version 47 to ??
      return (&Get_ftp($db)."release-".$ensembl_version."/fasta/");
    }
  }

  elsif ($db eq "ensembl_genomes") {                                         ## Ensembl genomes
    my @sites = ("fungi","bacteria","metazoa","plants","protists");
    if ($ensembl_version < 3) {                                              # Version  1 to 3
      return ();
    } else {                                                                 # Version 3 to ??
      my @fasta_ftp = ();

      foreach $site (@sites) {
        my $site_ftp = &Get_ftp($db)."release-".$ensembl_version."/".$site."/fasta/";

        if ($site eq "bacteria") {
          my @available_files = qx{wget -S --spider $site_ftp 2>&1};
          foreach $file (@available_files) {
            next unless ($file =~ /^d/);
            @token = split(" ",$file);
            next if ($token[-1] =~ /^\./);
            push (@fasta_ftp,$site_ftp.$token[-1]."/");
          }
        } else {
    
          push (@fasta_ftp,$site_ftp);
        }
      }
      return @fasta_ftp;
    }
  }
}


################################################################
##
## Get the URL of the peptidic sequences fasta file on Ensembl ftp
## site.
##
## This URL can then be used to download peptidic sequences in a quick
## way.
##
## Beware: the method is only valid for Ensembl version >= 72.
sub Get_pep_fasta_ftp {
  my ($db,$species,$ensembl_version) = @_;
  my @fasta_ftps = &Get_fasta_ftp($db,$ensembl_version);

  if ($db eq "ensembl") {                                                    ## Ensembl
    my $pep_fasta_ftp = $fasta_ftps[0].$species."/pep/";
    my @available_pep_fasta = qx{wget -S --spider $pep_fasta_ftp 2>&1};

    foreach my $fasta (@available_pep_fasta) {
      chomp ($fasta);
      my @token = split(" ",$fasta);
      return $pep_fasta_ftp.$token[-1] if ( $fasta =~ 'all' );
    }
  } elsif ($db eq "ensembl_genomes") {                                         ## Ensembl genomes
    foreach my $fasta_ftp (@fasta_ftps) {
      my @available_species = qx{wget -S --spider $fasta_ftp 2>&1};
      foreach my $spe (@available_species) {
        chomp($spe);
        next unless ($spe =~ /$species/);
        my $pep_fasta_ftp = $fasta_ftp.$species."/pep/";
        my @available_pep_fasta = qx{wget -S --spider $pep_fasta_ftp 2>&1};
        foreach my $fasta (@available_pep_fasta) {
          chomp ($fasta);
          my @token = split(" ",$fasta);
          return $pep_fasta_ftp.$token[-1] if ( $fasta =~ 'all' );
        }
      }
    }
  }
}


## Get ftp variation url
sub Get_variation_ftp {
  my ($db,$ensembl_version) = @_;

  if ($db eq "ensembl") {                                                    ## Ensembl

    if ($ensembl_version < 60) {                                             # Version  1 to 59
      return ();
    } elsif ($ensembl_version < 64) {                                        # Version 60 to 63
      return (&Get_ftp($db)."release-".$ensembl_version."/variation/");
    } else {                                                                 # Version 64 to ??
      return (&Get_ftp($db)."release-".$ensembl_version."/variation/gvf/");
    }
  }

  elsif ($db eq "ensembl_genomes") {                                         ## Ensembl genomes
    my @sites = ("fungi","bacteria","metazoa","plants","protists");

    if ($ensembl_version < 17) {                                             # Version  1 to 16
      return ();
    } else {                                                                 # Version 17 to ??
      my @variation_ftp = ();

      foreach $site (@sites) {
        my $site_ftp = &Get_ftp($db)."release-".$ensembl_version."/".$site."/";

        my @available_files = qx{wget -S --spider $site_ftp 2>&1};
        foreach $file (@available_files) {
          next unless ($file =~ /^d/);
          push (@variation_ftp,$site_ftp."gvf/") if ($file =~ /gvf/);
        }
      }

      return @variation_ftp;
    }
  }
}

## Get ftp variation species url
sub Get_variation_species_ftp {
  my ($db,$species,$ensembl_version) = @_;
  my @variation_ftps = &Get_variation_ftp($db,$ensembl_version);

  if ($db eq "ensembl") {                                                    ## Ensembl

    if ($ensembl_version == 61) {                                            # Version 61
      return $variation_ftps[0].ucfirst($species)."/";
    } else {                                                                 # Version 60, 62 to ??
      return $variation_ftps[0].$species."/";
    }
  }

  elsif ($db eq "ensembl_genomes") {                                         ## Ensembl genomes

                                                                             # Version 17 to ??
    foreach my $variation_ftp (@variation_ftps) {
      my @available_species = qx{wget -S --spider $variation_ftp 2>&1};
      foreach my $spe (@available_species) {
        chomp($spe);
        next unless ($spe =~ /$species/);
        return $variation_ftp.$species."/" if ($spe eq $species);
      }
    }
  }
}


##Get ftp gvf file url
sub Get_gvf_ftp {
  my ($db,$species,$ensembl_version) = @_;

  if ($db eq "ensembl") {                                                    ## Ensembl

    if ($ensembl_version == 60) {                                            # Version 60
      return &Get_variation_species_ftp($db,$species,$ensembl_version).$species.".gvf.gz";

    } else {                                                                 # Version 61 to ??
      return &Get_variation_species_ftp($db,$species,$ensembl_version).ucfirst($species).".gvf.gz";
    }
  }

  elsif ($db eq "ensembl_genomes") {                                         ## Ensembl genomes

                                                                             # Version 17 to ??
    return &Get_variation_species_ftp($db,$species,$ensembl_version).$species.".gvf.gz";
  }
}


################################################################
## Get an the main taxon (bacteria, fungi, metazoa, ...) for each
## species supported in an ansembl database. The result is returned as
## a has table, with species names as keys and taxa as values.
sub Get_species_taxon {
  my ($db, $ensembl_version) = @_;
  my %species_taxon = ();

  my @fasta_url = &Get_fasta_ftp($db,$ensembl_version);

  foreach my $url (@fasta_url) {
    my @token = split('/',$url);
    my $taxon = $token[5];

    my @available_files = qx{wget -S --spider $url 2>&1};
    foreach $file (@available_files) {
      next unless ($file =~ /^d/);
      my @token = split(" ",$file);
      next if ($token[-1] =~ /^\./);
      $species_taxon{$token[-1]} = $taxon;
    }
  }
  return %species_taxon;
}

############################################################################
############################################################################
########################## API Ensembl FCT ################################# 
############################################################################ 

# get API host name
sub Get_host_port {
  my ($db) = @_;

  &RSAT::message::TimeWarn("Getting host port for database", $db) if ($main::verbose >= 2);

  if ($db eq "ensembl") {
    return ('ensembldb.ensembl.org','5306');
  } elsif ($db eq "ensembl_genomes") {
    return ("mysql.ebi.ac.uk","4157");
  }
}

############################################################################
############################################################################
#### Specification of local directories for installing Ensembl on RSAT #####
############################################################################ 

=pod

=item B<Get_species_dir_name()>

Compute the directory name for a given species, assembly version.

Parameters:
    species            mandatory
    assembly_version   mandatory
    ensembl_version    mandatory
=cut

sub Get_species_dir_name {
  my ($species, $assembly_version, $ensemb_version) = @_;
  my $dir_name = ucfirst($species);

  ## Previous directory naming convention, temporarily maintained for
  ## backward compatibility.
  $old_naming = 0;
  if ($old_naming) {
      $dir_name .= "_".$main::db;
      if ($main::db eq "ensembl_genomes") {
	  $dir_name .= "-".$ensembl_version;
      } else {
	  $dir_name .= "-".$ensembl_version;
      }
      $dir_name .= "_".$assembly_version;
  } else {
      ## New directory naming convention (2014-10-28, JvH  AMR)
      $dir_name .= "_".$assembly_version;
      $dir_name .= "_".$main::db.$ensembl_version;
  }

  &RSAT::message::Info("&Get_species_dir_name() result", $dir_name) if ($main::verbose >= 0);
  return($dir_name);
}


=pod

=item B<Get_assembly_vesion()>

Return the genome assembly that corresponds to a specific Ensembl
version for a given species.

=cut
sub Get_assembly_version {
  my ($species,$ensembl_version) = @_;
  $species = ucfirst($species);
  $supported_file = &Get_supported_file();

  if (-f $supported_file ) {
    my ($file) = &OpenInputFile($supported_file);

    while (<$file>) {
	chomp();
	my ($id,$name,$dir) = split("\t");
	my ($species_f,$assembly_version_f,$ensembl_version_f) = split(" ",$name);
	return $assembly_version_f if ($species_f eq $species && $ensembl_version_f eq $ensembl_version);
    }
  }
  return "";
}

=pod

=item B<Get_data_dir()>

Return the main data directory for this RSAT server.

=cut
sub Get_data_dir {
  return $ENV{'RSAT'}."/data";
}

=pod

=item B<Get_genomes_dir()>

Return the directory where genomes are stored on this RSAT server.

=cut 
sub Get_genomes_dir {
    my $data_dir = &Get_data_dir();
    return $data_dir."/genomes";
}

=pod

=item B<Get_species_dir()>

Return the directory where the current species (downloaded from
Ensembl) is installed on this RSAT server.

=cut
sub Get_species_dir {
  my ($species,$assembly_version,$ensembl_version) = @_;
  $species = ucfirst($species);
  $supported_file = &Get_supported_file();

  my %assembly_directory = ();


  my $species_dir = join("/", &Get_genomes_dir(),
			 &Get_species_dir_name($species,$assembly_version,$ensembl_version));
  &RSAT::message::Info("&Get_species_dir() result", $species_dir) if ($main::verbose >= 0);
  return($species_dir);
}

=pod

=item B<Get_assembly_from_ensembl_version()>

Given a user-specified ensembl version, return the corresponding
genome assembly. The information is read from the file
supported_organisms.tab.

=cut
sub Get_assembly_from_ensembl_version {  
  ## Temporarily inactivate the previous options, which were searching
  ## for a compatible genome already installed in the table describing
  ## ensembl supported organisms. We now impose the directory to be
  ## always specified inthe same way. Modif by Jacques van Helden and
  ## Alejandra medina-Rivera, 2014-10-28.
    if (-f $supported_file ) {
	## Open the file containing the list of supported Ensembl species
	my ($file) = &OpenInputFile($supported_file);
	
	foreach (<$file>) {
	    chomp();
	    my ($id,$name,$dir) = split("\t");
	    $dir =~ s|\$ENV\{RSAT\}|$ENV{RSAT}|g;
	    my ($spe,$ass,$ens) = split(" ",$name);
	    
	    if ($ensembl_version && $assembly_version) {
		  ## If the directory has already been defined in the
		## supported organisms file, return it from there
		return $dir if (($spe eq $species) && ($ass eq $assembly_version) && ($ens eq $ensembl_version));
	    } elsif ($ensembl_version) {
		return $dir if (($spe eq $species) && ($ens eq $ensembl_version));
	    } else {
		$assembly_directory{$ens} = $dir if (($spe eq $species) && ($ass eq $assembly_version));
	    }
	}
	foreach (sort{$b<=>$a} (keys(%assembly_directory))) {
	    return $assembly_directory{$_};
	}
    }
}

=pod

=item B<Get_genome_dir()>

Return the directory in which the genome data (sequences + features)
will be installed for a given ensembl species.

=cut
sub Get_genome_dir {
  my ($species, $assembly_version,$ensembl_version) = @_;

  my $genome_dir = &Get_species_dir($species, $assembly_version,$ensembl_version);
  $genome_dir .= "/genome";
  &RSAT::message::Info("&Get_genome_dir() result", $genome_dir) if ($main::verbose >= 0);

  return($genome_dir);
}


=pod

=item B<Get_variation_dir()>

Return the directory in which the variations will be installed for a
given ensembl species.

=cut
sub Get_variation_dir {
  my ($species, $assembly_version,$ensembl_version) = @_;

  my $variation_dir = &Get_species_dir($species, $assembly_version,$ensembl_version);
  $variation_dir .= "/variations";
  &RSAT::message::Info("&Get_variation_dir() result", $variation_dir) if ($main::verbose >= 0);
  return ($variation_dir);
}


=pod

=item B<Get_supported_file()>

Return the path to the tab-delimited file describing all species
installed from Ensembl.

=cut
sub Get_supported_file {
  my ($data_dir) = @_;
  $data_dir = &Get_data_dir unless ($data_dir);
  return $data_dir."/supported_organisms_ensembl.tab";
}

## Contigs.txt
sub Get_contigs_file {
  my ($genome_dir) = @_;
  return $genome_dir."/contigs.txt";
}

## Contig.tab
sub Get_contig_file {
  my ($genome_dir) = @_;
  return $genome_dir."/contig.tab";
}

## Feature.tab
sub Get_feature_file {
  my ($species, $assembly_version,$ensembl_version,$name) = @_;
  $name =~ s/ /_/g;
  $name = lc($name);
  return &Get_genome_dir($species, $assembly_version,$ensembl_version).$name.".tab";
}


############################ Fct get file_chr name

## Get list of sequence file
sub Get_file_seq_name {
  my ($genome_dir) = @_;
  my %chr_file = ();
  my %file_info = ();

  ## Get $accession and seq_id
  my $contig = &Get_contig_file($genome_dir);
  if (-f $contig) {
    my ($file) = &OpenInputFile($contig);
    while (<$file>) {
      next if (/--/);
      chomp();
      my ($chr,$acc) = split("\t");
      $file_info{$acc} = $chr;
    }
    close $file;
  } else {
    &RSAT::error::FatalError("$contig is missing.");
  }

  ##  Get seq_file_name and seq_id
  my $contigs = &Get_contigs_file($genome_dir);
  if (-f $contigs) {
    my ($file) = &OpenInputFile($contigs);
    while (<$file>) {
      chomp();
      my ($file_name,$acc) = split("\t");
      $chr_file{$file_info{$acc}} = $file_name;
    }
    close $file;
  } else {
    &RSAT::error::FatalError("$contigs is missing.");
  }

  return %chr_file;
}

return 1;

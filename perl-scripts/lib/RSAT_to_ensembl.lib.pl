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
our $ensembl_version_safe = $ENV{ensembl_version_safe} || 72;
our $ensemblgenomes_version_safe = $ENV{ensemblgenomes_version_safe} || 19;

## Fields of the table describing the supported organisms obtained
## from Ensembl. These fields are used by several methods
our @supported_header_fields = ("id",
#			 "name", 
			 "species",
			 "assembly_version",
			 "db",
			 "ensembl_version",
			 "update_date",
			 "species_directory",
	);

################################################################
## Functions
################################################################

################################################################
## Return the safe version of ensembl or EnsemblGenomes
## (this version should be defined in RSAT_config.props).
sub get_ensembl_version_safe {
  my ($db) = @_;

  if ($db eq "ensembl") {
    return $ensembl_version_safe;
  } elsif (lc($db) eq "ensemblgenomes") {
    return $ensemblgenomes_version_safe;
  }
}

=pod

=item B<&get_ensembl_version()>

Return the ensembl release to be installed.  The default ensembl
release is specified in the local property file
$RSAT/RSAT_config.props. This default value can be overwritten with
the option -version.

=cut
sub get_ensembl_version {

  ## TEMPORARY: for the time being, the version is only used by
  ##  reference to ensembl, not ensemblgenomes, even if the queried db is
  ##  EnsemblGenomes. I (JvH) should cliarify this.  
  my $db = "ensembl";
  # my ($db) = @_;

  ## By preference, return the version defined in the property file
  ## (RSAT_config.props) for this database). If not defined, get the
  ## latest version available on the Ensembl or EnsemblGenomes server.
  if ($db eq "ensemblgenomes") {
      ## Note: currently ignored, see TEMPORARY above.
      my $ensemblgenomes_version = $ENV{ensemblgenomes_version} || &get_ensembl_version_latest("ensemblgenomes");
      return($ensemblgenomes_version);
  } else {
      $ensembl_version = $ENV{ensembl_version} || &get_ensembl_version_latest("ensembl");
      return ($ensembl_version);
  }

#  $prop_name = $db."_version";
#  if (defined($ENV{$prop_name})) {
#    return ($ENV{$prop_name});
#  } else {
#    return &get_ensembl_version_latest($db);
#  }
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

  &RSAT::message::TimeWarn("Getting latest ensembl version from FTP site", $ftp) 
      if ($main::verbose >= 3);

  my @available_releases = qx{wget -S --spider $ftp 2>&1};

#  &RSAT::message::Debug("Available releases", join ("\n\t", @available_releases)) if ($main::verbose >= 10);

  foreach my $release (@available_releases) {
    next if ($release =~ /current/);
    next unless ($release =~ /^[dl].*[^\.]release/);
    chomp($release);
    my @token = split("release-",$release);
    $latest_ensembl_release = $token[-1] if ($latest_ensembl_release < $token[-1]);
  }

  &RSAT::message::Info("&get_ensembl_version_latest() result", $latest_ensembl_release) if ($main::verbose >= 5);

  return ($latest_ensembl_release);
}

################################################################
## Check that the user-selected version is correct.  This version can
## neither be lower than the safe version, nor higher than the latest
## version.
sub check_ensembl_version {
  my ($db,$ensembl_version) = @_;

  ## TEMPORARY: for the time being, the version is only used by
  ##  reference to ensembl, not ensemblgenomes, even if the queried db is
  ##  EnsemblGenomes. I (JvH) should cliarify this.  
  $db = "ensembl";

  my $safe_ensembl_version = &get_ensembl_version_safe($db);
  my $latest_ensembl_version = &get_ensembl_version_latest($db);

  ## Report checking parameters
  if ($main::verbose >= 2) {
    &RSAT::message::Info("Latest ensembl version", $latest_ensembl_version);
    &RSAT::message::Info("Safe ensembl version", $safe_ensembl_version);
    &RSAT::message::Info("Selected ensembl version",$ensembl_version);
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
  } elsif (lc($db) eq "ensemblgenomes") {
    return "ftp://ftp.ensemblgenomes.org/pub/";
  } else {
    &RSAT::error::FatalError($db, "Is not a valid name for Ensembl databases. Supported: ensembl, ensemblgenomes.");
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

  elsif (lc($db) eq "ensemblgenomes") {                                         ## Ensembl genomes
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
  } elsif (lc($db) eq "ensemblgenomes") {                                         ## Ensembl genomes
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


=pod

Return the URL of the FTP site containing variations for the selected
version of ensembl.

=cut
sub Get_variation_ftp {
  my ($db,$ensembl_version) = @_;
  my %variation_ftp = ();

  ## Ensembl
  if ($db eq "ensembl") {                                                    
    if ($ensembl_version < 60) {                                             # Version  1 to 59
      return ();
    } elsif ($ensembl_version < 64) {                                        # Version 60 to 63
      return (&Get_ftp($db)."release-".$ensembl_version."/variation/");
    } else {                                                                 # Version 64 to ??
      return (&Get_ftp($db)."release-".$ensembl_version."/variation/gvf/");
    }

    ## Variations from EnsemblGenomes are distributed on a different
    ## ftp site, and the version numbers differ
  } elsif (lc($db) eq "ensemblgenomes") {                                         
    my @taxa = ("fungi","bacteria","metazoa","plants","protists");
    if ($ensembl_version < 17) {                                             # Version  1 to 16
      return ();
    } else {                                                                 # Version 17 to ??
      foreach $taxon (@taxa) {
        my $gvf_ftp = &Get_ftp($db).$taxon."/release-".$ensembl_version."/gvf/";
	&RSAT::message::Info("Getting list of GVF files for taxon", $taxon, "FTP", $gvf_ftp) if ($main::verbose >= 2);
        my @available_files = qx{wget -S --spider $gvf_ftp 2>&1};
        foreach my $line (@available_files) {
	  chomp($line);
#	  &RSAT::message::Debug($line) if ($main::verbose >= 10);
          next unless ($line =~ /^drw+.*\s+(\S+\_\S+)\s*$/);
	  my $species = $1;
	  my $species_ftp = $gvf_ftp.$species."/";
	  $variation_ftp{$species} = $species_ftp;
	}
      }
      return(%variation_ftp);
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

  elsif (lc($db) eq "ensemblgenomes") {                                         ## Ensembl genomes

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

  elsif (lc($db) eq "ensemblgenomes") {                                         ## Ensembl genomes

                                                                             # Version 17 to ??
    return &Get_variation_species_ftp($db,$species,$ensembl_version).$species.".gvf.gz";
  }
}


################################################################
## Get an the main taxon (bacteria, fungi, metazoa, ...) for each
## species supported in an ansembl database. The result is returned as
## a has table, with species names as keys and taxa as values.
##
## JvH: THIS IS TRICKY: uses an ftp server. I should rewrite it using
## the Lookup interface or something else, see with Dan Staines.
sub Get_species_taxon {
  my ($db,$ensembl_version) = @_;
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
  &RSAT::message::TimeWarn("Getting host port for database", $db) if ($main::verbose >= 3);

  if ($db eq "ensembl") {
      return ('ensembldb.ensembl.org','5306');
  } elsif (lc($db) eq "ensemblgenomes") {
      return('mysql-eg-publicsql.ebi.ac.uk', '4157');
#      return ("mysql.ebi.ac.uk","4157");
  }
}



=pod

=head1 B<LoadRegistry>

Usage:
    LoadRegistry($registry);

Establish connection to Ensembl and EnsemblGenomes.


BEWARE: the API relies on SQL queries through non-conventional
ports. These ports must be authorized by the local firewall.

List of ports for Ensembl
    http://www.ensembl.org/info/data/mysql.html

List of ports for EnsemblGenomes
    http://ensemblgenomes.org/info/access/mysql

=cut


sub LoadRegistry {
  ($registry, $db, $ensembl_version) = @_;

  &RSAT::message::TimeWarn("Loading registry from", $db, "version ".$ensembl_version) if ($main::verbose >= 2);
# my $ensembl_version = 78;
# Bio::EnsEMBL::Registry->load_registry_from_db(
#            -host => 'mysql-eg-publicsql.ebi.ac.uk',
#            -port => '4157',
#            -user => 'anonymous',
#            -db_version => $ensembl_version,
#    -verbose=>1
#    );

  if ($db eq "ensembl") {
    $registry->load_registry_from_db(
      -host => 'ensembldb.ensembl.org',
      -port => '5306',
      -user => 'anonymous',
      -db_version => $ensembl_version,
#      -verbose=>0
	);
  } elsif ($db eq "ensemblgenomes") {
    $registry->load_registry_from_db(
      -host => 'mysql-eg-publicsql.ebi.ac.uk', 
      -port => '4157',
      -user => 'anonymous',
      -db_version => $ensembl_version,
#      -verbose=>0
	);
  } elsif ($db eq "ensemblall") {
    $registry->load_registry_from_multiple_dbs 
	(
	 {-host => 'mysql-eg-publicsql.ebi.ac.uk',
	  -port => 4157, 
	  -user => 'anonymous',
	  -db_version => $ensembl_version,
#	  -verbose=>0
	 },
	 {-host => 'ensembldb.ensembl.org',
	  -port => 5306,
	  -user    => 'anonymous',
	  -db_version => $ensembl_version,
#	  -verbose=>0
	 }
	);
  } else {
    &RSAT::error::FatalError("Invalid db for ensembl queries. Supported: ensembl, ensemblgenomes, ensemblall.");
  }
  my $nb_species = scalar(@{ $registry->get_all_DBAdaptors(-group => 'core') });
  &RSAT::message::TimeWarn("Loaded registry with", $nb_species, "species") if ($main::verbose >= 3);
}

############################################################################
############################################################################
#### Specification of local directories for installing Ensembl on RSAT #####
############################################################################ 

=pod

=item B<Get_full_species_ID()>

Compute the directory name for a given species and assembly version.

I<Parameters>

=over

=item I<species>

Mandatory argument.

E.g.: Homo_sapiens, Escherichia_coli_str_k_12_substr_mg1655

Species should correspond to species names supported by ensembl, with
underscores to replace the spaces.

=item i<assembly_version>

Ex: 
  GRCh38 (for Homo sapiens), 
   GCA_000005845.2 (for Escherichia_coli_str_k_12_substr_mg1655).

Optional. If not specified, the program attempts to find get
assembly version in the table of organisms previously installed from
Ensembl.

=item I<ensembl_version>

Optional. If not specified, the program uses the default Ensembl
version, which is specified in the file $RSAT/RSAT_config.props.

=back

The full ID will be used to select the organism in the RSAT tools. It
also servers as name for the directory in which the genome has to be
installed.

By default, the full species ID is built by concatenating species and
assembly, separated by an underscore character.

Ex: 

=cut

sub Get_full_species_ID {
  my ($species, $assembly_version, $ensembl_version, $species_suffix) = @_;

  &RSAT::message::Debug("&Get_full_species_ID()", $db, $species, $assembly_version, $ensembl_version, $species_suffix) if ($main::verbose >= 5);

  ## Check that Ensembl version has been provided. If not, take
  ## default one.
  unless ($ensembl_version) {
      $ensembl_version = &get_ensembl_version();
      &RSAT::message::Debug("&Get_full_species_ID() called without ensembl_version argument",
			    "Using default",$ensembl_version) if ($main::verbose >= 5);
  }

  ## Check that the assembly version has been provided. If not, guess
  ## it.
  unless ($assembly_version) {
      &RSAT::message::Debug("&Get_full_species_ID() called without assembly_version argument") if ($main::verbose >= 5);
      $assembly_version = &Get_assembly_version($species,$ensembl_version,$species_suffix);
      &RSAT::message::Debug("Got from &Get_assembly_version()", $assembly_version) if ($main::verbose >= 5);
  }


  ## Full ID convention (2014-10, JvH  AMR)
  ## [Species]_[assembly_version]_[db][ensembl_version]
  my $full_species_id = ucfirst($species);
  $full_species_id .= "_".$assembly_version;
#  $full_species_id .= "_".$main::db.$ensembl_version; ## We prefer to avoid creating one folder for each new release of ensembl
  $full_species_id .= "_".$species_suffix if ($species_suffix);

  &RSAT::message::Info("&Get_full_species_ID() result", $full_species_id) if ($main::verbose >= 5);
  return($full_species_id);
}


=pod

=item B<Get_assembly_version()>

Return the genome assembly that corresponds to a specific Ensembl
version for a given species.

=cut
sub Get_assembly_version {
  my ($species,$ensembl_version,$species_suffix) = @_;
  my $assembly_version = "";
  &RSAT::message::Debug("&Get_assembly_version()", 
			"main::db=".$main::db,
			"species=".$species, 
			"ensembl_version=".$ensembl_version,
			"species_suffix=".$species_suffix,
      ) 
      if ($main::verbose >= 5);
  $supported_file = &Get_supported_file();

  ## Check if the organism is installed in the tab-delimited file of organisms
  if (-f $supported_file ) {
    my ($file) = &OpenInputFile($supported_file);

    my $l=0;
    while (<$file>) {
	$l++;
	next if (/^;/); ## Skip comment lines
	next if (/^#/); ## Skip header line
	next unless (/\S/); ## Skip empty lines
	chomp();
	my (@fields) = split("\t");
	foreach my $field  (@supported_header_fields) {
	  ## Automatically fill attributes corresponding to the column header
	  $var_name = "db_".$field;
	  $$var_name = shift(@fields);
	}
	
	&RSAT::message::Debug("Get_assembly_version", "line=".$l, 
			      "\n\tquery", $species, $main::db,$ensembl_version,
			      "\n\tdb", $db_species, $db_db, $db_ensembl_version,
	    ) if ($main::verbose >= 5);
	if ((lc($species) eq lc($db_species)) 
	    && ($db_db eq $main::db)
	    && ($ensembl_version eq $db_ensembl_version)
	    ) {
	    $assembly_version = $db_assembly_version;
	    &RSAT::message::Info("&Get_assembly_version() result", $assembly_version) if ($main::verbose >= 4);
	    return($assembly_version);
	}
    }
  }

  ## If not found, issue a warning then die
  &RSAT::message::Warning("&Get_assembly_version() could not identify species", $species, 
			  "from", $main::db.$ensembl_version, "in the organism table\n", $supported_file);
  &RSAT::error::FatalError($species, "genome does not seem to be installed in version ".$ensembl_version." of ".$main::db, 
			   "\nTry the following command:\n\t", 
			   "install-ensembl-genome -v 2 -species ".$species);
}

=pod

=item B<Get_data_dir()>

Return the main data directory for this RSAT server.

=cut
sub Get_data_dir {
    my $data_dir = $ENV{'RSAT'}."/public_html/data";
    &RSAT::message::Info("&Get_data_dir() result", $data_dir) if ($main::verbose >= 5);
    return $data_dir;
}

=pod

=item B<Get_genomes_dir()>

Return the directory where genomes are stored on this RSAT server.

=cut 
sub Get_genomes_dir {
    my $data_dir = &Get_data_dir();
    my $genomes_dir = $data_dir."/genomes";
    &RSAT::message::Info("&Get_genomes_dir() result", $genomes_dir) if ($main::verbose >= 5);
    return $genomes_dir;
}

=pod

=item B<Get_species_dir()>

Return the directory where the current species (downloaded from
Ensembl) is installed on this RSAT server.

=cut
sub Get_species_dir {
  my ($species,$assembly_version,$ensembl_version,$species_suffix) = @_;
  &RSAT::message::Debug("&Get_species_dir()", "species=".$species, "assembly_version=".$assembly_version, "ensembl_version=".$ensembl_version) 
      if ($main::verbose >= 5);

  $species = ucfirst($species);
  $supported_file = &Get_supported_file();

  my %assembly_directory = ();

  my $species_dir = &Get_species_dir_from_supported_file($species);

  ## Define species directory based on species name, assembly and ensembl_version
  unless ($species_dir) {
      $species_dir = join("/", &Get_genomes_dir(),
			  &Get_full_species_ID($species,$assembly_version,$ensembl_version,$species_suffix));
  }
  &RSAT::message::Info("&Get_species_dir() result", $species_dir) if ($main::verbose >= 5);
  return($species_dir);
}

=pod

=item B<Get_species_dir_from_supported_file()>

Given a user-specified ensembl version, identify the corresponding
assembly version, and the full species ID (species name + assembly +
ensembl version). The information is read from the table of
ensembl-specific organisms
(${RSAT}/public_html/data/supported_organisms_ensembl.tab).

=cut
sub Get_species_dir_from_supported_file {  
  my ($species) =  @_;
    if (-f $supported_file ) {
	## Open the file containing the list of supported Ensembl species
	my ($file) = &OpenInputFile($supported_file);
	
	while (my $line = <$file>) {
	    chomp($line);
	    my ($id,$name,$dir) = split("\t", $line);

	    ## The full RSAT path should not be writen explicitly in
	    ## the files.
	    if ($dir) {
		$dir =~ s|\$ENV\{RSAT\}|$ENV{RSAT}|g;
	    }

	    if ($name) {
		## Note (JvH, 2014-10-30): the "species name" actually
		## includes the species name (with _ to separate substrain
		## etc), the assembly, and the ensembl version. This is
		## not very clean. We should have a file with the
		## different information types in separated fields.
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
	}

	## ??? THIS SHOULD NOT WORK: the return cannot be included in
	## a loop ! (Note by JvH, 2014-10-30)
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
  my ($species, $assembly_version,$ensembl_version,$species_suffix) = @_;
  &RSAT::message::Debug("&Get_genome_dir()", "species=".$species, "assembly_version=".$assembly_version, "ensembl_version=".$ensembl_version) 
      if ($main::verbose >= 5);

  my $genome_dir = &Get_species_dir($species, $assembly_version,$ensembl_version,$species_suffix);
  $genome_dir .= "/genome";
  &RSAT::message::Info("&Get_genome_dir() result", $genome_dir) if ($main::verbose >= 5);

  return($genome_dir);
}


=pod

=item B<Get_variation_dir()>

Return the directory in which the variations will be installed for a
given ensembl species.

=cut
sub Get_variation_dir {
  my ($species, $assembly_version,$ensembl_version, $species_suffix) = @_;

  my $variation_dir = &Get_species_dir($species, $assembly_version,$ensembl_version,$species_suffix);
  $variation_dir .= "/variations";
  &RSAT::message::Info("&Get_variation_dir() result", $variation_dir) if ($main::verbose >= 5);
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

=pod

=item B<Get_feature_file()>

Returns the full path for a tab-delimited file containing features of
a given type.

Usage: 

 my $feature_file = ($species, $assembly_version,$ensembl_version,$species_suffix,$type);

Where $type is the feature type (e.g. CDS, gene, mrna).

=cut
sub Get_feature_file {
  my ($species, $assembly_version,$ensembl_version,$species_suffix,$type) = @_;
  $type =~ s/ /_/g;
  $type = lc($type);
  my $file = &Get_genome_dir($species, $assembly_version,$ensembl_version,$species_suffix)."/".$type.".tab";
  return ($file);
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
      &RSAT::error::FatalError("Missing contig table", 
			       "\n\t".$contig,
			       "\n\tYou should first install genomic features (download-ensembl-features)",
	  );
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

=pod

    Update the tab-delimited file with the description of supported
    genomes downloaded from Ensembl.

=cut
sub UpdateEnsemblSupported {
    my	$supported_organism_file = &Get_supported_file();
    &RSAT::message::TimeWarn("Updating supported organism file", $supported_organism_file) if ($main::verbose >= 2);
    
    ## Hash table to store previous species description lines
    my %species_description = ();

    ## Find the current species ID
    my $current_species_id = &Get_full_species_ID($species,$assembly_version,$ensembl_version,$species_suffix);
	
    ## Read the list of previously installed organisms if it exists.
    if (-f $supported_organism_file) {
	my ($s_o_file) = &OpenInputFile($supported_organism_file);
	
	## Read the whole file of supported organisms from ensembl,
	## and store species description lines in a hash indexed by
	## full species ID, in order to sort them after having changed
	## the current species fields.
	my $l = 0;
	while (<$s_o_file>) {
	    $l++;
	    next if (/^;/); ## Skip comment lines
	    next if (/^#/); ## Skip header line
	    next unless (/\S/); ## Skip empty lines
	    chomp();
	    my @fields = split("\t");
	    my $full_species_id = $fields[0];
	    $species_description{$full_species_id} = $_;
	}
	close $s_o_file;
    }
    

    ## Build the line for the currently installed species
    my $id = &Get_full_species_ID($species,$assembly_version,$ensembl_version,$species_suffix);
    my $name = $id; $name =~ s/_/ /g;

    my $new_org_config = join ("\t", 
			       $id,
#			       $name, 
			       $species,
			       $assembly_version,
			       $main::db,
			       $ensembl_version , #&get_ensembl_version,
			       &AlphaDate(),
			       &Get_species_dir($species,$assembly_version,$ensembl_version,$species_suffix),
	);

    ## Avoid to expose the full RSAT path
    $new_org_config =~ s|$ENV{RSAT}|\$\{RSAT\}\/|g;
    $new_org_config =~ s|\/\/|/|g;

    ## Index the new species description
    $species_description{$current_species_id} = $new_org_config;
       
    ## Write the updated table of supported organisms from Ensembl
    my $s_o_file = &OpenOutputFile($supported_organism_file);

    ## Print the header with column content
    print $s_o_file "#", join ("\t", @supported_header_fields), "\n";

    ## Print the table of supported organisms
    foreach my $id (sort keys %species_description) {
	print $s_o_file $species_description{$id}, "\n";
    }
#    print $s_o_file join("",@other_species);
    close $s_o_file;
    
    &RSAT::message::Info("Ensembl genome installed in folder", $genome_dir) if ($main::verbose >= 1);
}

return 1;

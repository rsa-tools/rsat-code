#!/usr/bin/env perl

package main;

############################################################################
############################################################################
########################## FTP ENSEMBL FONCTION ############################

################################################################
## Define global variables
##
## Note: this is not very clean. I (JvH) should pass these variables
## in the appropriate methods and check that everything still works
## fine.
our $ensembl_release_safe = $ENV{ensembl_release_safe} || 72;
our $ensemblgenomes_release_safe = $ENV{ensemblgenomes_release_safe} || 19;

## Fields of the table describing the supported organisms obtained
## from Ensembl. These fields are used by several methods
our @supported_header_fields = ("ID",
             "name",    #species
             "taxid",
             "source",  #db
             "last_update",
             "nb",
             "seq_format",
             "up_from",
             "up_to",
             "taxonomy",
             "data",    #species_directory
             "genome",
             "genome_assembly", #assembly
             "genome_version", #ensembl_release
             "download_date", #update_date
             "variant_available",
             "variant_source",
			 "path_to_variant_files",
             "blast_available"
	);

################################################################
## Functions
################################################################

################################################################
## Return the safe release of ensembl or EnsemblGenomes
## (this release should be defined in RSAT_config.props).
sub get_ensembl_release_safe {
  my ($db) = @_;

  if ($db eq "ensembl") {
       &RSAT::message::Info("&get_ensembl_release_safe() result", $ensembl_release_safe) if ($main::verbose >= 5);
    return $ensembl_release_safe;
  } elsif (lc($db) eq "ensemblgenomes") {
       &RSAT::message::Info("&get_ensembl_release_safe() result", $ensemblgenomes_release_safe) if ($main::verbose >= 5);
    return $ensemblgenomes_release_safe;
  }
}

=pod

=item B<&get_ensembl_release()>

Return the ensembl release to be installed.  The default ensembl
release is specified in the local property file
$RSAT/RSAT_config.props. This default value can be overwritten with
the option -release.

=cut
sub get_ensembl_release {

  ## TEMPORARY: for the time being, the release is only used by
  ##  reference to ensembl, not ensemblgenomes, even if the queried db is
  ##  EnsemblGenomes. I (JvH) should clarify this.
  my ($db) = @_;

  ## By preference, return the release defined in the property file
  ## (RSAT_config.props) for this database). If not defined, get the
  ## latest release available on the Ensembl or EnsemblGenomes server.
  if ($db eq "ensemblgenomes") {
      ## Note: currently ignored, see TEMPORARY above.
      my $ensemblgenomes_release = $ENV{ensemblgenomes_release} || &get_ensembl_release_latest("ensemblgenomes");
       &RSAT::message::Info("&get_ensembl_release() result", $ensemblgenomes_release) if ($main::verbose >= 5);
      return($ensemblgenomes_release);
  } else {
      $ensembl_release = $ENV{ensembl_release} || &get_ensembl_release_latest("ensembl");
       &RSAT::message::Info("&get_ensembl_release() result", $ensembl_release) if ($main::verbose >= 5);
      return ($ensembl_release);
  }

#  $prop_name = $db."_release";
#  if (defined($ENV{$prop_name})) {
#    return ($ENV{$prop_name});
#  } else {
#    return &get_ensembl_release_latest($db);
#  }
}


################################################################
## Get the latest ensembl release for a species from Ensembl ftp site.
##
## THIS IS REALLY TRICKY, I (JvH) should see with Ensembl if their API
## includes a way to get the current release.
sub get_ensembl_release_latest {
  my ($db) = @_;
  my $latest_ensembl_release = 0;
  my $ftp = &Get_ftp($db);

  &RSAT::message::TimeWarn("Getting latest ensembl release from FTP site", $ftp)
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

  &RSAT::message::Info("&get_ensembl_release_latest() result", $latest_ensembl_release) if ($main::verbose >= 5);

  return ($latest_ensembl_release);
}

################################################################
## Check that the user-selected release is correct.  This release can
## neither be lower than the safe release, nor higher than the latest
## release.
sub check_ensembl_release {
  my ($db,$ensembl_release) = @_;

  ## TEMPORARY: for the time being, the release is only used by
  ##  reference to ensembl, not ensemblgenomes, even if the queried db is
  ##  EnsemblGenomes. I (JvH) should cliarify this.
    #$db = "ensembl";

  my $safe_ensembl_release = &get_ensembl_release_safe($db);
  my $latest_ensembl_release = &get_ensembl_release_latest($db);

  ## Report checking parameters
  if ($main::verbose >= 2) {
    &RSAT::message::Info("Latest ensembl release", $latest_ensembl_release);
    &RSAT::message::Info("Safe ensembl release", $safe_ensembl_release);
    &RSAT::message::Info("Selected ensembl release",$ensembl_release);
  }

  if ($ensembl_release eq "safe") {
    ## Automatic selection of the safe release
    $ensembl_release = $safe_ensembl_release;

  } elsif ($ensembl_release eq "latest") {
    ## Automatic selection of the latest release
    $ensembl_release = $latest_ensembl_release;

  } else {

    ## Check that selected release is not smaller than safest one (defined in RSAT_config.props)
    &RSAT::error::FatalError($ensembl_release,
			     "Invalid  release for", $db,
			     "Minimun safe release: ".$safe_ensembl_release)
	if ($ensembl_release < $safe_ensembl_release);

    ## Check that selected release is not higher than the latest supported release
    &RSAT::error::FatalError($ensembl_release,
			     "Invalid  release for", $db,
			     "Latest release: ".$latest_ensembl_release)
	if ($ensembl_release > $latest_ensembl_release);
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
  my ($db,$ensembl_release) = @_;

  if ($db eq "ensembl") {                                                    ## Ensembl
    if ($ensembl_release < 47) {                                             # Release  1 to 46
      return ();
    } else {                                                                 # Release 47 to ??
      return (&Get_ftp($db)."release-".$ensembl_release."/fasta/");
    }
  }

  elsif (lc($db) eq "ensemblgenomes") {                                         ## Ensembl genomes
    my @sites = ("fungi","bacteria","metazoa","plants","protists");
    if ($ensembl_release < 3) {                                              # Release  1 to 3
      return ();
    } else {                                                                 # Release 3 to ??
      my @fasta_ftp = ();

      foreach $site (@sites) {
        my $site_ftp = &Get_ftp($db)."release-".$ensembl_release."/".$site."/fasta/";

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
## Beware: the method is only valid for Ensembl release >= 72.
sub Get_pep_fasta_ftp {
  my ($db,$species,$ensembl_release) = @_;
  my @fasta_ftps = &Get_fasta_ftp($db,$ensembl_release);

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

Return the URL of the FTP site containing GVF variations for the selected release of ensembl

Valid URLS as of Jan2020:
ftp://ftp.ensemblgenomes.org/pub/release-46/plants/variation/vcf/triticum_turgidum
ftp://ftp.ensemblgenomes.org/pub/release-46/metazoa/variation/gvf/aedes_aegypti_lvpagwg

=cut
sub Get_variation_ftp {
  my ($db,$ensembl_release) = @_;
  my %variation_ftp = ();

  ## Ensembl
  if ($db eq "ensembl") {
    if ($ensembl_release < 60) {                                             # Release  1 to 59
      return ();
    } elsif ($ensembl_release < 64) {                                        # Release 60 to 63
      return (&Get_ftp($db)."release-".$ensembl_release."/variation/");
    } else {                                                                 # Release 64 to ??
      return (&Get_ftp($db)."release-".$ensembl_release."/variation/gvf/");
    }

    ## Variations from EnsemblGenomes are distributed on a different
    ## ftp site, and the release numbers differ
  } elsif (lc($db) eq "ensemblgenomes") {
    my @divisions = ("fungi","bacteria","metazoa","plants","protists");
    if ($ensembl_release < 17) {                                             # Release  1 to 16
      return ();
    } else {                                                                 # Release 17 to ??
      foreach $taxon (@diviions) {
        my $gvf_ftp = &Get_ftp($db)."/release-".$ensembl_release."/$taxon/variation/gvf/";
        &RSAT::message::Info("Getting list of GVF files for taxon", $taxon, "FTP", $gvf_ftp) if ($main::verbose >= 2);
        my @available_files = qx{wget -S --spider $gvf_ftp 2>&1};
        foreach my $line (@available_files) {
          chomp($line);
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
  my ($db,$species,$ensembl_release) = @_;
  my @variation_ftps = &Get_variation_ftp($db,$ensembl_release);

  if ($db eq "ensembl") {                                                    ## Ensembl

    if ($ensembl_release == 61) {                                            # Release 61
      return $variation_ftps[0].ucfirst($species)."/";
    } else {                                                                 # Release 60, 62 to ??
      return $variation_ftps[0].$species."/";
    }
  }

  elsif (lc($db) eq "ensemblgenomes") {                                         ## Ensembl genomes

                                                                             # Release 17 to ??
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
  my ($db,$species,$ensembl_release) = @_;

  if ($db eq "ensembl") {
      ## Note: the case of GVF file names changed between release 90 and 91.
      ## Until release 90 (included), species name was with first letter in uppercase.
      ## Since release 91, genus and species names are fully in lowercases.
      if ($ensembl_release <= 90) {
	  return &Get_variation_species_ftp($db,$species,$ensembl_release).ucfirst($species).".gvf.gz";
      } else {
	  return &Get_variation_species_ftp($db,$species,$ensembl_release).$species.".gvf.gz";
      }
  } elsif (lc($db) eq "ensemblgenomes") {
      return &Get_variation_species_ftp($db,$species,$ensembl_release).$species.".gvf.gz";
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
  my ($db,$ensembl_release) = @_;
  my %species_taxon = ();

  my @fasta_url = &Get_fasta_ftp($db,$ensembl_release);

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
  ($registry, $db, $ensembl_release) = @_;

  &RSAT::message::TimeWarn("Loading registry from", $db, "release ".$ensembl_release) if ($main::verbose >= 2);
# my $ensembl_release = 78;
# Bio::EnsEMBL::Registry->load_registry_from_db(
#            -host => 'mysql-eg-publicsql.ebi.ac.uk',
#            -port => '4157',
#            -user => 'anonymous',
#            -db_release => $ensembl_release,
#    -verbose=>1
#    );

  if ($db eq "ensembl") {
    $registry->load_registry_from_db(
      -host => 'ensembldb.ensembl.org',
      -port => '5306',
      -user => 'anonymous',
      -db_version => $ensembl_release,
#      -verbose=>0
	);
  } elsif ($db eq "ensemblgenomes") {
    $registry->load_registry_from_db(
      -host => 'mysql-eg-publicsql.ebi.ac.uk',
      -port => '4157',
      -user => 'anonymous',
      -db_version => $ensembl_release,
#      -verbose=>0
	);
  } elsif ($db eq "ensemblall") {
    $registry->load_registry_from_multiple_dbs
	(
	 {-host => 'mysql-eg-publicsql.ebi.ac.uk',
	  -port => 4157,
	  -user => 'anonymous',
	  -db_version => $ensembl_release,
#	  -verbose=>0
	 },
	 {-host => 'ensembldb.ensembl.org',
	  -port => 5306,
	  -user    => 'anonymous',
	  -db_version => $ensembl_release,
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

Compute the directory name for a given species and assembly.

I<Parameters>

=over

=item I<species>

Mandatory argument.

E.g.: Homo_sapiens, Escherichia_coli_str_k_12_substr_mg1655

Species should correspond to species names supported by ensembl, with
underscores to replace the spaces.

=item i<assembly>

Ex:
  GRCh38 (for Homo sapiens),
   GCA_000005845.2 (for Escherichia_coli_str_k_12_substr_mg1655).

Optional. If not specified, the program attempts to find get
assembly in the table of organisms previously installed from
Ensembl.

=item I<ensembl_release>

Optional. If not specified, the program uses the default Ensembl
release, which is specified in the file $RSAT/RSAT_config.props.

=back

The full ID will be used to select the organism in RSAT tools. It
also serves as name for the directory in which the genome has to be
installed.

By default, the full species ID is built by concatenating species and
assembly, separated by an underscore character. Example: 

However, historically other formats have co-existed. For instance, 'install-organisms' 
installs organisms from NCBI and uses strain instead of assembly ([species]_[strain])

Finally, makefiles/ensemblgenomes_FTP_client.mk produces ID such as [species].[assembly].[release]
Example: Oryza_longistaminata.O_longistaminata_v1.0.43

=cut

sub Get_full_species_ID {
  my ($species, $assembly, $ensembl_release, $species_suffix) = @_;

  &RSAT::message::Debug("&Get_full_species_ID()", $db, $species, $assembly, $ensembl_release, $species_suffix) if ($main::verbose >= 5);

  ## Check that Ensembl release has been provided. If not, take
  ## default one.
  unless ($ensembl_release) {
      $ensembl_release = &get_ensembl_release();
      &RSAT::message::Debug("&Get_full_species_ID() called without ensembl_release argument",
			    "Using default",$ensembl_release) if ($main::verbose >= 5);
  }

  ## Check that the assembly has been provided. If not, guess
  ## it.
  unless ($assembly) {
      &RSAT::message::Debug("&Get_full_species_ID() called without assembly argument") if ($main::verbose >= 5);
      $assembly = &Get_assembly($species,$ensembl_release,$species_suffix);
      &RSAT::message::Debug("Got from &Get_assembly()", $assembly) if ($main::verbose >= 5);
  }

  # compose full ID by concatenating bits
  my $full_species_id;
  my $genome_data_dir = Get_data_dir() . '/genomes/';

  # 1) try the makefiles/ensemblgenomes_FTP_client.mk
  $full_species_id = ucfirst($species) .'.'. $assembly .'.'. $ensembl_release; 
  if(-d $genome_data_dir.$full_species_id){
    &RSAT::message::Info("&Get_full_species_ID() result", $full_species_id) if ($main::verbose >= 5);
    return($full_species_id);
  } 
  
  # 2) try the NCBI way
  if($species_suffix){
    $full_species_id = ucfirst($species) .'.'. $species_suffix;
    if(-d $genome_data_dir.$full_species_id){
      &RSAT::message::Info("&Get_full_species_ID() result", $full_species_id) if ($main::verbose >= 5);
      return($full_species_id);
    }
  }

  # 3) finally try the default RSAT ID format
  ## Full ID convention (2014-10, JvH  AMR)
  ## [Species]_[assembly]_[db][ensembl_release]
  $full_species_id = ucfirst($species);
  $full_species_id .= "_".$assembly;
  #$full_species_id .= "_".$main::db.$ensembl_release; ## We prefer to avoid creating one folder for each new release of ensembl
  $full_species_id .= "_".$species_suffix if ($species_suffix);

  &RSAT::message::Info("&Get_full_species_ID() result", $full_species_id) if ($main::verbose >= 5);
  return($full_species_id);
}


=pod

=item B<Get_assembly()>

Return the genome assembly that corresponds to a specific Ensembl
release for a given species.

=cut
sub Get_assembly {
    my ($species,$ensembl_release,$species_suffix) = @_;
    my $assembly = "";
    &RSAT::message::Debug("&Get_assembly()",
    "main::db=".$main::db,
    "species=".$species,
    "ensembl_release=".$ensembl_release,
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
            
            &RSAT::message::Debug("Get_assembly", "line=".$l,
            "\n\tquery", $species, $main::db,$ensembl_release,
            "\n\tdb", $db_species, $db_db, $db_ensembl_release,
            ) if ($main::verbose >= 5);
            
            $db_name =~ s/\s/_/g;
            if ((lc($species) eq lc($db_name))
                && ($db_source eq $main::db)
                && ($ensembl_release eq $db_genome_version)
                ) {
                    $assembly = $db_genome_assembly;
                    &RSAT::message::Info("&Get_assembly() result", $assembly) if ($main::verbose >= 4);
                    return($assembly);
                }
        }
    }
    
    ## If not found, issue a warning then die
    &RSAT::message::Warning("&Get_assembly() could not identify species", $species,
    "from", $main::db.$ensembl_release, "in the organism table\n", $supported_file);
    &RSAT::error::FatalError($species, "genome does not seem to be installed in release ".$ensembl_release." of ".$main::db,
    "\nTry the following command:\n\t",
    "install-ensembl-genome -v 2 -db ".$main::db." -species ".$species);
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
  my ($species,$assembly,$ensembl_release,$species_suffix) = @_;
  &RSAT::message::Debug("&Get_species_dir()", "species=".$species, "assembly=".$assembly, "ensembl_release=".$ensembl_release)
      if ($main::verbose >= 5);

  $species = ucfirst($species);
  $supported_file = &Get_supported_file();

  my %assembly_directory = ();

  my $species_dir = &Get_species_dir_from_supported_file($species);

  ## Define species directory based on species name, assembly and ensembl_release
  unless ($species_dir) {
      $species_dir = join("/", &Get_genomes_dir(),
			  &Get_full_species_ID($species,$assembly,$ensembl_release,$species_suffix));
  }
  &RSAT::message::Info("&Get_species_dir() result", $species_dir) if ($main::verbose >= 5);
  return($species_dir);
}

=pod

=item B<Get_species_dir_from_supported_file()>

Given a user-specified ensembl release, identify the corresponding
assembly, and the full species ID (species name + assembly +
ensembl release). The information is read from the table of
ensembl-specific organisms
(${RSAT}/public_html/data/supported_organisms_ensembl.tab).

=cut
sub Get_species_dir_from_supported_file {
  my ($species) =  @_;

  my %assembly_directory;
  my $supported_file = &Get_supported_file();

  if (-f $supported_file ) {
		## Open the file containing the list of supported Ensembl species
		my ($file) = &OpenInputFile($supported_file);

		while (my $line = <$file>) {
	    	chomp($line);
	    	my ($id,$name,$ass,$db,$ens,$update,$dir) = split("\t", $line);

				# Make sure first letter of species is upper case
				$spe = ucfirst($name);
				$species = ucfirst($species);

   	    	## The full RSAT path should not be writen explicitly in
	    	## the files.
	    	if ($dir) {
					$dir =~ s/\$\{RSAT\}/$ENV{RSAT}/g;
	    	}

	    	if ($spe) {
					 # If release and assembly parameters have been queried, use them for evaluations
					if ($ensembl_release && $assembly) {
						return $dir if (($spe eq $species) && ($ass eq $assembly) && ($ens eq $ensembl_release));

					# If ONLY release parameter has been queried, use it for evaluations
					} elsif ($ensembl_release) {
						return $dir if (($spe eq $species) && ($ens eq $ensembl_release));

					# If ONLY assembly parameter has been queried, use it for evaluations
					# and stote results in a hash in order to select the most recent release
					} else {
						$assembly_directory{$ens} = $dir if (($spe eq $species) && ($ass eq $assembly));
					}
	    	}
		}

		# Sort directories by latest release and then return it
		my @sort_release_index = (sort{$b<=>$a} (keys(%assembly_directory)));
		if(scalar(@sort_release_index)){
            return $assembly_directory{$sort_release_index[0]};
        } else {
            return '';
        }
  }
}

=pod

=item B<Get_genome_dir()>

Return the directory in which the genome data (sequences + features)
will be installed for a given ensembl species.

=cut
sub Get_genome_dir {
  my ($species, $assembly,$ensembl_release,$species_suffix) = @_;
  &RSAT::message::Debug("&Get_genome_dir()", "species=".$species, "assembly=".$assembly, "ensembl_release=".$ensembl_release)
      if ($main::verbose >= 5);

  my $genome_dir = &Get_species_dir($species, $assembly,$ensembl_release,$species_suffix);
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
  my ($species, $assembly,$ensembl_release, $species_suffix) = @_;

  my $variation_dir = &Get_species_dir($species, $assembly,$ensembl_release,$species_suffix);
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
  return $data_dir."/supported_organisms.tab";
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

 my $feature_file = ($species, $assembly,$ensembl_release,$species_suffix,$type);

Where $type is the feature type (e.g. CDS, gene, mrna).

=cut
sub Get_feature_file {
  my ($species, $assembly,$ensembl_release,$species_suffix,$type) = @_;
  $type =~ s/ /_/g;
  $type = lc($type);
  my $file = &Get_genome_dir($species, $assembly,$ensembl_release,$species_suffix)."/".$type.".tab";
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
    my ($species, $assembly,$ensembl_release, $species_suffix) = @_;
    my	$supported_organism_file = &Get_supported_file();
    &RSAT::message::TimeWarn("Updating supported organism file", $supported_organism_file) if ($main::verbose >= 2);
    
    ## Hash table to store previous species description lines
    my %species_description = ();
    
    ## Find the current species ID
    my $current_species_id = &Get_full_species_ID($species,$assembly,$ensembl_release,$species_suffix);
    
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
    
    my $name = $species;
    $name =~ s/_/ /g;
    ## Build the line for the currently installed species
    my $new_org_config = join ("\t",
    $current_species_id,
    ucfirst($name),
    "<NA>",
    $main::db,
    &AlphaDate(),
    "<NA>",
    "<NA>",
    "<NA>",
    "<NA>",
    "<NA>",
    &Get_species_dir($species,$assembly,$ensembl_release,$species_suffix),
    &Get_species_dir($species,$assembly,$ensembl_release,$species_suffix) . "/configs.txt",
    $assembly,
    $ensembl_release , #&get_ensembl_release,
    &AlphaDate(),
    "<NA>",
    "<NA>",
    "<NA>",
    "<NA>"
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

=pod
 
 Update the tab-delimited file with the blasts installed for
 get-orthologs
 
=cut

##########################################
###
sub UpdateBlastSupported {
    my ($species_id) = @_;
    my $supported_organism_file = &Get_supported_file();
    &RSAT::message::TimeWarn("Updating supported organism file", $supported_organism_file) if ($main::verbose >= 2);
    
    ## Hash table to store previous species description lines
    my %species_description = ();
    if (-f $supported_organism_file) {
        my ($s_o_file) = &OpenInputFile($supported_organism_file);
        
        ## Read the whole file of supported organisms from ensembl,
        ## and store species description lines in a hash indexed by
        ## full species ID, in order to sort them after having changed
        ## the current species fields.
        my $l = 0;
        while (<$s_o_file>){
            $l++;
            next if (/^;/); ## Skip comment lines
            next if (/^#/); ## Skip header line
                next unless (/\S/); ## Skip empty lines
            chomp();
            my @fields = split("\t");
            my $current_species_id = $fields[0];
            #check if the species is already in the list
            if($species_id eq $current_species_id){
                # Build the line for the species
                $fields[18] = 1;
                $species_description{$current_species_id} = join("\t", @fields);
            }else{
                $species_description{$current_species_id} = $_;
            }
        }
        close $s_o_file;
    }
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


=pod

    Update the tab-delimited file with the variants installed for
		variation-info

=cut

##########################################
###
sub UpdateVariationsSupported {
    my ($species, $assembly, $source, $release) = @_;
    #my @supported_header_fields = ("id", "species", "assembly", "db", "release", "update_date");
    my $data_dir = &Get_data_dir();
    #my $supported_organism_file = $data_dir."/supported_variation_info.tab";
    my $supported_organism_file = &Get_supported_file();
    &RSAT::message::TimeWarn("Updating supported organism file", $supported_organism_file) if ($main::verbose >= 2);

    ## Hash table to store previous species description lines
    my %species_description = ();

    ## Find the current species ID
    #my $current_species_id = ucfirst($species)."_".$assembly;
    my $current_species_id = Get_full_species_ID($species, $assembly, $release );

    my $already_exist = 0;
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
                #check if the species is already in the list
                if($current_species_id eq $full_species_id){
                    $already_exist = 1;
                    # Build the line for the species
                    if($fields[12] eq "<NA>") { $fields[12] = $assembly; }
                    if($fields[13] eq "<NA>") { $fields[13] = $release; }
                    $fields[14] = &AlphaDate();
                    $fields[15] = 1;
                    $fields[16] = $source;
                    $fields[17] = &Get_variation_dir($species, $assembly,$ensembl_release, $species_suffix);
                    $species_description{$full_species_id} = join("\t", @fields);
                }else{
                    $species_description{$full_species_id} = $_;
                }
			}
			close $s_o_file;
    }


    ## Build the line for the currently installed species
    ##my $id = &Get_full_species_ID($species,$assembly,$ensembl_release,$species_suffix);
    ##my $name = $id; $name =~ s/_/ /g;
    if($already_exist == 0){
        my $name = $species;
        $name =~ s/_/ /g;
        my $new_org_config = join ("\t",
                       $current_species_id,
                             ucfirst($name),
                            "<NA>",
                            $source,
                            &AlphaDate(),
                            "<NA>",
                            "<NA>",
                            "<NA>",
                            "<NA>",
                            "<NA>",
                            "<NA>",
                            "<NA>",
                             $assembly,
                             $release,
                       &AlphaDate(),
                            1,
                            $source,
                            &Get_variation_dir($species, $assembly,$release, $species_suffix),
                            "<NA>"
                        );

        ## Avoid to expose the full RSAT path
        ##$new_org_config =~ s|$ENV{RSAT}|\$\{RSAT\}\/|g;
        ##$new_org_config =~ s|\/\/|/|g;

        ## Index the new species description
        $species_description{$current_species_id} = $new_org_config;
    }
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

    &RSAT::message::Info("Ensembl variations installed for ", $species) if ($main::verbose >= 1);
		return;
}


return 1;

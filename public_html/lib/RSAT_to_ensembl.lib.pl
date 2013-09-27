#!/usr/bin/perl

package main;

############################################################################
############################################################################
########################## FTP ENSEMBL FONCTION ############################
#  our $ensembl_url = "ftp://ftp.ensemblgenomes.org/pub/protists/current/fasta/";

## Define global variables
our $ensembl_version_safe = $ENV{ensembl_version_safe} || 72;

## Functions

# get $ensembl_version_safe
sub Get_ensembl_version_safe {
  my ($db) = @_;

  if ($db eq "ensembl") {
    return $ensembl_version_safe;
  }

  elsif ($db eq "ensembl_genomes") {
   return 18;
  }
}

# get $ensembl_ftp
sub Get_ftp {
  my ($db) = @_;

  if ($db eq "ensembl") {
    return "ftp://ftp.ensembl.org/pub/";
  }

  elsif ($db eq "ensembl_genomes") {
    return "ftp://ftp.ensemblgenomes.org/pub/";
  }
}

## Get ftp fasta url
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


## Get ftp pep fasta species url
#######  Quick way to download features 
#######  Work only for version 72-73

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
  }

  elsif ($db eq "ensembl_genomes") {                                         ## Ensembl genomes
                                                                             
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


## Get the latest ensembl version for a species
sub Get_ensembl_version {
  my ($db) = @_;

  my $ftp = &Get_ftp($db);

  my $current_release = 0;

  &RSAT::message::TimeWarn("Getting ensembl version", $ftp) if ($main::verbose >= 2);
  my @available_release = qx{wget -S --spider $ftp 2>&1};


  foreach (@available_release) {
    next if (/current/);
    next unless (/^[dl].*[^\.]release/);
    chomp();
    my @token = split("release-",$_);
    $current_release = $token[-1] if ($current_release < $token[-1]);
  }

  return $current_release;
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

sub Get_species_dir_name {
  my ($species,$assembly_version,$ensemb_version) = @_;
  $species = ucfirst($species);
  return $species."_ensembl_".$ensembl_version."_".$assembly_version;
}

sub Get_assembly_version {
  my ($data_dir,$species,$ensembl_version) = @_;
  $species = ucfirst($species);
  $supported_file = &Get_supported_file($data_dir);

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

############################ Fct local dir

sub Get_data_dir {
  return $ENV{'RSAT'}."/data/";
}

sub Get_genomes_dir {
  my ($data_dir) = @_;
  return $data_dir."/genomes/";
}

sub Get_species_dir {
  my ($data_dir,$species,$assembly_version,$ensembl_version) = @_;
  $species = ucfirst($species);
  $supported_file = &Get_supported_file($data_dir);
  
  my %assembly_directory = ();

  ## Open the file containing the list of supported Ensembl species
  if (-f $supported_file ) {
    my ($file) = &OpenInputFile($supported_file);

    foreach (<$file>) {
      chomp();
      my ($id,$name,$dir) = split("\t");
      $dir =~ s|\$ENV\{RSAT\}|$ENV{RSAT}|g;
      my ($spe,$ass,$ens) = split(" ",$name);

      if ($ensembl_version && $assembly_version) {
        return $dir if ($spe eq $species && $ass eq $assembly_version && $ens eq $ensembl_version);
      } elsif ($ensembl_version) {
        return $dir if ($spe eq $species && $ens eq $ensembl_version);
      } else {
        $assembly_directory{$ens} = $dir if ($spe eq $species && $ass eq $assembly_version);
      }
    }

    foreach (sort{$b<=>$a} (keys(%assembly_directory))) {
      return $assembly_directory{$_};
    }
  }

  return &Get_genomes_dir($data_dir).&Get_species_dir_name($species,$assembly_version,$ensembl_version)."/";
}


sub Get_genome_dir {
  my ($data_dir,$species, $assembly_version,$ensembl_version) = @_;
 
  return &Get_species_dir($data_dir, $species, $assembly_version,$ensembl_version)."genome/";
}

sub Get_variation_dir {
  my ($data_dir,$species, $assembly_version,$ensembl_version) = @_;
  return &Get_species_dir($data_dir, $species, $assembly_version,$ensembl_version)."variations/";
}

############################ Fct get file

## supported_organims_ensembl.tab
sub Get_supported_file {
  my ($data_dir) = @_;
  return $data_dir."/supported_organisms_ensembl.tab";
}

## Contigs.txt
sub Get_contigs_file {
  my ($genome_dir) = @_;
  return $genome_dir."contigs.txt";
}

## Contig.tab
sub Get_contig_file {
  my ($genome_dir) = @_;
  return $genome_dir."contig.tab";
}

## Feature.tab
sub Get_feature_file {
  my ($data_dir,$species, $assembly_version,$ensembl_version,$name) = @_;
  $name =~ s/ /_/g;
  $name = lc($name);
  return &Get_genome_dir($data_dir,$species, $assembly_version,$ensembl_version).$name.".tab";
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

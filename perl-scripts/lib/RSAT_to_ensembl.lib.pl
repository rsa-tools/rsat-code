#!/usr/bin/perl

package main;

############################################################################
############################################################################
########################## FTP ENSEMBL FONCTION ############################
#  our $ensembl_url = "ftp://ftp.ensemblgenomes.org/pub/protists/current/fasta/";

## Define global variables
our $ensembl_rsync = $ENV{ensembl_rsync} || "rsync://ftp.ensembl.org/ensembl/pub";
our $ensembl_version_safe = $ENV{ensembl_version_safe} || 72;


## Functions

# get $ensembl_ftp
sub Get_ensembl_ftp() {
  my ($db) = @_;
  
  if ($db eq "ensembl") {
    return "ftp://ftp.ensembl.org/pub/";
    
  } 
}

# get $ensembl_version_safe
sub Get_ensembl_version_safe() {
  return $ensembl_version_safe;
}


## Get the latest ensembl version for a species
sub Get_ensembl_version() {
  my ($db) = @_;

  if ($db eq "ensembl") {
  	
  	  my $ftp = &Get_ensembl_ftp($db)."current_fasta/homo_sapiens/dna/";
      my @available_fasta = qx{wget -S --spider $ftp 2>&1};

	  foreach (@available_fasta) {
	  	print $_;
	    next unless (/Homo_sapiens/);
	    $_ =~ s/Homo_sapiens\.//g;
	    my @token = split(".dna",$_);
	    my @token2 = split(/\./,$token[0]);
	    
	    return $token2[-1];
	  }
  }
}


## Get ftp fasta url
sub Get_fasta_ftp() {
  my ($db,$ensembl_version) = @_;
  return &Get_ensembl_ftp($db)."release-".$ensembl_version."/fasta/";     # Version 60 to 72
}

## Get ftp variation url
sub Get_variation_ftp() {
  my ($db,$ensembl_version) = @_;

  if ($ensembl_version < 64) {
    return &Get_ensembl_ftp($db)."release-".$ensembl_version."/variation/";     # Version 60 to 63
  } else {
    return &Get_ensembl_ftp($db)."release-".$ensembl_version."/variation/gvf/";     # Version 64 to 72
  }
}

## Get ftp species url
sub Get_species_ftp() {
  my ($db,$species,$ensembl_version,$type) = @_;


  if ($type eq "fasta") {
    return &Get_fasta_ftp($db,$ensembl_version).$species."/dna/";      # Version 48 to 72
  }

  if ($type eq "variation") {
    if ($ensembl_version == 61) {
      return &Get_variation_ftp($db,$ensembl_version).ucfirst($species)."/";      # Version 61
    } else {
      return &Get_variation_ftp($db,$ensembl_version).$species."/";     # Version 60 and 62 to 72
    }
  }
}


##Get ftp gvf file url
sub Get_gvf_ftp() {
  my ($db,$species,$ensembl_version) = @_;

  if ($ensembl_version == 60) {
    return &Get_species_ftp($db,$species,$ensembl_version,'variation').$species.".gvf.gz";     # Version 60
  } else {
    return &Get_species_ftp($db,$species,$ensembl_version,'variation').ucfirst($species).".gvf.gz";      # Version 61 to 72
  }
}

## Get the genome version for a species
sub Get_assembly_version() {
  my ($db,$species,$ensembl_version) = @_;
  my $species_fasta_ftp = &Get_species_ftp($db,$species,$ensembl_version,'fasta');

  my @available_fasta = qx{wget -S --spider $species_fasta_ftp 2>&1};

  $species = ucfirst($species);
  foreach (@available_fasta) {
    next unless (/$species/);
    my @token = split($species.".");
    my @token2 = split(".dna",$token[-1]);
    my @token3 = split(/\./,$token2[0]);
    return join '.', @token3[0..$#token3-1];
  }
}


############################################################################
############################################################################
########################## API Ensembl FCT ################################# 
############################################################################ 

# get API host name
sub Get_host_port() {
  my ($db) = @_;
  
  if ($db eq "ensembl") {
    return ('ensembldb.protist.org','5306');
  } elsif ($db eq "ensembl_genomes") {
    return ("mysql.ebi.ac.uk","4157");
  }
  
}

sub Get_lastest_ensembl_version_api() {
  my (@db_adaptors,$species) = @_;

  foreach my $db_adaptor (@db_adaptors) {
    if ($db_adaptor->species() eq $species) {
      my $db_connection = $db_adaptor->dbc();
      my @token = split ("_",$db_connection->dbname());
      reutrn $token[-1];
    }
  }
}


############################################################################
############################################################################
#### Specification of local directories for installing Ensembl on RSAT #####
############################################################################ 


############################ Fct get local dir

sub Get_data_dir() {
	return $ENV{'RSAT'}."/data/";
}

sub Get_genomes_dir() {
	my ($data_dir) = @_;
    return $data_dir."genomes/";
}

## Get the local directory for the user-specified species
sub Get_species_dir() {
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

      if ($ensembl_version) {
         return $dir if ($name =~ /$species.*$assembly_version.*$ensembl_version/);
      } else {
         my ($spe,$ass,$ens) = split(" ",$name);
         $assembly_directory{$ens} = $dir if ($name =~ /$species.*$assembly_version/);
      }
    } 

    foreach (sort{$b<=>$a} (keys(%assembly_directory))) {
      return $assembly_directory{$_};
    }
  }

  return &Get_genomes_dir($data_dir).&Get_species_dir_name($species,$assembly_version,$ensembl_version);
}


sub Get_species_dir_name() {
  my ($species,$assembly_version,$ensemb_version) = @_;
  return $species."_ensembl_".$assembly_version."_".$ensembl_version."/";
}


## Genome dir
sub Get_genome_dir() {
  my ($data_dir,$species, $assembly_version,$ensembl_version) = @_;
  return &Get_species_dir($data_dir, $species, $assembly_version,$ensembl_version)."genome/";
}

## Variation dir
sub Get_variation_dir() {
  my ($data_dir,$species, $assembly_version,$ensembl_version) = @_;
  return &Get_species_dir($data_dir, $species, $assembly_version,$ensembl_version)."variations/";
}

############################ Fct get file

## supported_organims_ensembl.tab
sub Get_supported_file() {
  my ($data_dir) = @_;
  return $data_dir."/supported_organisms_ensembl.tab";
}

## Contigs.txt
sub Get_contigs_file() {
  my ($genome_dir) = @_;
  return $genome_dir."contigs.txt";
}

## Contig.tab
sub Get_contig_file() {
  my ($genome_dir) = @_;
  return $genome_dir."contig.tab";
}

## Feature.tab
sub Get_feature_file() {
  my ($species, $assembly_version,$ensembl_version,$name) = @_;
  $name =~ s/ /_/g;
  return &Get_genome_dir($species, $assembly_version,$ensembl_version).$name.".tab";
}


## variation.gvf
sub Get_variation_file() {
  my ($genome_dir) = @_;
  return $genome_dir."contig.tab";
}


############################ Fct get file_chr name

## Get list of sequence file
sub Get_file_seq_name() {
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

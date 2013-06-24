#!/usr/bin/perl

package main;


############################################################################
############################################################################
############################## RSYNC FONCTION ##############################
#  our $ensembl_url = "ftp://ftp.ensemblgenomes.org/pub/protists/current/fasta/";

####### Variables
our $ensembl_rsync = $ENV{ensembl_rsync} || "rsync://ftp.ensembl.org/ensembl/pub";
our $ensembl_version_safe = $ENV{ensembl_version_safe} || 72;


####### Fct

# get $ensembl_rsync
sub Get_ensembl_rsync() {
  return $ensembl_rsync;
}

# get $ensembl_version_safe
sub Get_ensembl_version_safe() {
  return $ensembl_version_safe;
}


## Get the latest ensembl version for a species
sub Get_ensembl_version() {
  my @available_fasta = qx{rsync -navP $ensembl_rsync/current_fasta/homo_sapiens/dna/ "."};

  foreach (@available_fasta) {
    next unless (/Homo_sapiens/);

    $_ =~ s/$species\.//g;
    my @token = split(".dna",$_);
    my @token2 = split(/\./,$token[0]);
    return $token2[-1];
  }
}


## Get rsync fasta url
sub Get_fasta_rsync() {
  my ($ensembl_version) = @_;
  return $ensembl_rsync."/release-".$ensembl_version."/fasta/";     # Version 48 to 72
}

## Get rsync variation url
sub Get_variation_rsync() {
  my ($ensembl_version) = @_;

  if ($ensembl_version < 64) {
    return $ensembl_rsync."/release-".$ensembl_version."/variation/";     # Version 60 to 63

  } else {
    return $ensembl_rsync."/release-".$ensembl_version."/variation/gvf/";     # Version 64 to 72
  }
}

## Get rsync species url
sub Get_species_rsync() {
  my ($url,$species,$type,$ensembl_version) = @_;

  if ($type eq "fasta") {
    return $url.$species."/dna/";      # Version 48 to 72
  }
  
  if ($type eq "variation") {
    
    if ($ensembl_version == 61) {
      return $url.ucfirst($species)."/";      # Version 61
    } else {
      return $url.$species."/";     # Version 60 and 62 to 72
    }
  }
}

##Get rsync gvf file url

sub Get_gvf_rsync() {
  my ($url,$species,$ensembl_version) = @_;

  if ($ensembl_version == 60) {
    return $url.$species.".gvf.gz";     # Version 60
  } else {
    return $url.ucfirst($species).".gvf.gz";      # Version 61 to 72
  }
}


## Get the genome version for a species
sub Get_assembly_version() {
  my ($species_rsync,$species) = @_;
  my @available_fasta = qx{rsync -navP $species_rsync/ "."};

  foreach (@available_fasta) {
    next unless (/$species/);

    $_ =~ s/$species\.//g;
    my @token = split(".dna",$_);
    my @token2 = split(/\./,$token[0]);
    return join '.', @token2[0..$#token2-1];
  }
}



#{} []


############################################################################
############################################################################
#### Specification of local directories for installing Ensembl on RSAT #####
############################################################################

############################ Variables
our $supported_file = $ENV{'RSAT'}."/data/supported_organisms_ensembl.tab";


############################ Fct get local dir

## Get the local directory for the user-specified species
sub Get_species_dir() {
  my ($species,$assembly_version,$ensembl_version) = @_;

  ## Open the file containing the list of supported Ensembl species
  my ($file) = &OpenInputFile($supported_file);

  foreach (<$file>) {
      chomp();
      my ($id,$name,$dir) = split("\t");
      return $dir if ($name =~ /$species.*$assembly_version.*$ensembl_version/);
  }

  return $genomes_dir.&Get_species_dir_name($species,$assembly_version,$ensembl_version);
}


sub Get_species_dir_name() {
  my ($species,$assembly_version,$ensemb_version) = @_;
  return $species."_ensembl_".$assembly_version."_".$ensembl_version;
}


## Genome dir
sub Get_genome_dir() {
  my ($species, $assembly_version,$ensembl_version) = @_;
  return &Get_species_dir($species, $assembly_version,$ensembl_version)."/genome/";
}

## Variation dir
sub Get_variation_dir() {
  my ($species_dir) = @_;
  return $species_dir."/variations/";
}

############################ Fct get file

## supported_organims_ensembl.tab
sub Get_supported_file() {
  return $supported_file;
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

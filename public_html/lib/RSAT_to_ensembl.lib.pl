#!/usr/bin/perl

package main;

############################################################################
############################################################################
############################## RSYNC FONCTION ##############################
#  our $ensembl_url = "ftp://ftp.ensemblgenomes.org/pub/protists/current/fasta/";

## Define global variables
our $ensembl_rsync = $ENV{ensembl_rsync} || "rsync://ftp.ensembl.org/ensembl/pub";
our $ensembl_version_safe = $ENV{ensembl_version_safe} || 72;


## Functions

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
    $_ =~ s/Homo_sapiens\.//g;
    my @token = split(".dna",$_);
    my @token2 = split(/\./,$token[0]);
    
    return $token2[-1];
  }
}


## Get rsync fasta url
sub Get_fasta_rsync() {
  my ($ensembl_version) = @_;
  return $ensembl_rsync."/release-".$ensembl_version."/fasta/";     # Version 60 to 72
}

## Get rsync fasta url
sub Get_feature_rsync() {
  my ($ensembl_version) = @_;
  return $ensembl_rsync."/release-".$ensembl_version."/embl/";     # Version 48 to 72
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
  my ($species,$ensembl_version,$type) = @_;

  if ($type eq "fasta") {
    return &Get_fasta_rsync($ensembl_version).$species."/dna/";      # Version 48 to 72
  }
 
  if ($type eq "feature") {
    return &Get_feature_rsync($ensembl_version).$species."/";      # Version 48 to 72
  } 

  if ($type eq "variation") {
    if ($ensembl_version == 61) {
      return &Get_variation_rsync($ensembl_version).ucfirst($species)."/";      # Version 61
    } else {
      return &Get_variation_rsync($ensembl_version).$species."/";     # Version 60 and 62 to 72
    }
  }
}

## Get rsunc dat file url
sub Get_dat_rsync() {
	  my ($species,$ensembl_version) = @_;
	  my $species_feature_rsync = &Get_species_rsync($species,$ensembl_version,'feature');
	  my @seq_available = ();
	  
   my @available_feature = qx{rsync -navP $species_feature_rsync "."};
	
	 foreach (@available_feature ) {
    my $species_ucf = ucfirst($species);
    next unless (/$species_ucf/);
    next unless (/chrom/ && $ensembl_version >= 68);    ##Remove non chromosomal or contig seq (Version >68 only)
    
    my @info = split(/chrom/);
    next if ($info[1] =~ "_" && $ensembl_version >= 68);    ##Remove non chromosomal or contig seq (Version >68 only)
      
    push(@seq_available, $_);
  }
  return @seq_available;
}


##Get rsync gvf file url
sub Get_gvf_rsync() {
  my ($species,$ensembl_version) = @_;

  if ($ensembl_version == 60) {
    return &Get_species_rsync($species,$ensembl_version,'variation').$species.".gvf.gz";     # Version 60
  } else {
    return &Get_species_rsync($species,$ensembl_version,'variation').ucfirst($species).".gvf.gz";      # Version 61 to 72
  }
}

## Get the genome version for a species
sub Get_assembly_version() {
  my ($species,$ensembl_version) = @_;
  my $species_fasta_rsync = &Get_species_rsync($species,$ensembl_version,'fasta');

  my @available_fasta = qx{rsync -navP $species_fasta_rsync "."};

  $species = ucfirst($species);
  foreach (@available_fasta) {
    next unless (/$species/);

    $_ =~ s/$species\.//g;
    my @token = split(".dna",$_);
    my @token2 = split(/\./,$token[0]);
    return join '.', @token2[0..$#token2-1];
  }
}

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
  $species = ucfirst($species);

  my %assembly_directory = ();

  ## Open the file containing the list of supported Ensembl species
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

  return $genomes_dir.&Get_species_dir_name($species,$assembly_version,$ensembl_version);
}


sub Get_species_dir_name() {
  my ($species,$assembly_version,$ensemb_version) = @_;
  return $species."_ensembl_".$assembly_version."_".$ensembl_version."/";
}


## Genome dir
sub Get_genome_dir() {
  my ($species, $assembly_version,$ensembl_version) = @_;
  return &Get_species_dir($species, $assembly_version,$ensembl_version)."genome/";
}

## Variation dir
sub Get_variation_dir() {
  my ($species, $assembly_version,$ensembl_version) = @_;
  return &Get_species_dir($species, $assembly_version,$ensembl_version)."variations/";
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

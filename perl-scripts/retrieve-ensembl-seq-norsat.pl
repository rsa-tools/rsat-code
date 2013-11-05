#!/usr/bin/perl -w
############################################################
#
# $Id: retrieve-ensembl-seq-norsat.pl,v 1.7 2011/04/10 13:49:58 jvanheld Exp $
#
# Time-stamp
#
############################################################
#use strict;
use DBI();

BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
#	push (@INC, "/home/rsat/rsa-tools/perl-scripts/lib/");
#	push (@INC, "/home/rsat/rsa-tools/");
    }
}

#require "RSA.lib";
#require "RSA.seq.lib";
#require RSAT::util;

## EnsEMBL libraries
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
################################################################
#### main package
package main;
{

  ################################################################
  #### initialise parameters
  local $start_time = &RSAT::util::StartScript();

  local $verbose = 0;
  local $feattype = "mrna";    # other values: Gene, Intron, Exon, CDS and UTR
  local $type = "upstream";
  local $from = -800;
  local $to = -1;
  local $noorf = 0;
  local $nogene = 0;
  local $rm = 0;
  local $all_transcripts = 0;
  local $first_intron = 0;
#  local $utr = "all";    # other values: 5prime and 3prime
  local $non_coding = 0;
  local $all = 0;
  local $query_file;
  local @queries;
  local $left_limit;
  local $right_limit;
  local $strand = 1;
  local $chrom;
  local $ft_file;
  local $ft_file_format = "gft";
  local $mask_coding = 0;

  local $output_file;

  ## Connection to the EnsEMBL MYSQL database
  local $ensembl_host = 'ensembldb.ensembl.org';  # db at EBI (use outside BIGRE)
  local $ensembl_user = "anonymous";
  local $dbname = '';
  local $org = '';
  local $dbversion = '';

  ################################################################
  ## Read arguments
  &ReadArguments();

  ################################################################
  ## Check argument values

  # Either an organism or a database name must be provided
  unless ($org || $dbname) {
    die "; You must provide either an organism name (-org) or a database name (-dbname)\n";
  }

  # Database name has priority over organism if both provided
  if ($org && $dbname) {
    $org = '';
  }

  # A query must be provided
  unless (@queries || $query_file || $all || $ft_file || ($left_limit && $right_limit)) {
    die "; You must either provide an EnsEMBL gene ID (-q), a query file (-i), left and right limits or use the '-all' option\n";
  }

  ### verbose ###
  if ($verbose >= 1) {
    print "; retrieve-ensembl-seq ";
    print "\n";
  }

  ################################################################
  ## If option -org is used, connect to ensembldb to get list of 
  ## databases and pick the latest one corresponding to chosen organism
  if ($org) {
      my $dbh = DBI->connect("DBI:mysql:host=$ensembl_host", "$ensembl_user", "", {'RaiseError' => 1});
      my $sth = $dbh->prepare("SHOW DATABASES");
      $sth->execute();
      while (my $ref = $sth->fetchrow_hashref()) {
	  if ($ref->{'Database'} =~ /($org)_core_\d+/) {
	      $dbname = $ref->{'Database'};
	  }
      }
      $sth->finish();
      $dbh->disconnect();
  } else {  # get organism name from dbname
      $org = $dbname;
      $org =~s/_core_.+//;
  }
  
  ################################################################
  ## Get EnsEMBL db version from db name
  unless ($dbversion) {
      $dbversion = $dbname;
      $dbversion =~ s/($org)_core_//;
      $dbversion =~ s/_.+//;
  }
  
  ################################################################
  ## Open a new connection to EnsEMBL database, but this time we specify the DB name

  my $registry = "Bio::EnsEMBL::Registry";

  $registry->load_registry_from_db(
				   -host => $ensembl_host,
				   -user => $ensembl_user,
				   -db_version => $dbversion,
				   -verbose => "0" );

  local $db = Bio::EnsEMBL::Registry->get_DBAdaptor($org, "core");

#    local $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $ensembl_host, -user => $ensembl_user, -dbname => $dbname);


  ################################################################
  ### Open output stream
  if ($output_file) {
    $fh = 'OUT';
    open $fh, ">".$output_file || die "cannot open file ".$output_file."\n";
  } else {
    $fh = *STDOUT;
  }
  #### print verbose
#  &Verbose() if ($verbose);

  ################################################################
  ## Get query (one gene, several genes or all genes)

  local $slice_adaptor = $db->get_sliceAdaptor();

  # Left and right limits
  if ($left_limit && $right_limit && $chrom) {
    local $chromosome = $slice_adaptor -> fetch_by_region('chromosome', $chrom);
    # Get sequence (repeat masked or not)
    $sequence = &GetSequence($left_limit, $right_limit);
    my $size = $right_limit - $left_limit + 1;

    my $rsat_strand;
    if ($strand == 1) {
      $rsat_strand = "D";
    } else {
      $rsat_strand = "R";
    }

    # Export sequence to file
    print $fh ">$chrom-$left_limit-$right_limit\t$chrom-$left_limit-$right_limit; from 1 to $size; size: $size; location: $chrom $left_limit $right_limit $rsat_strand\n";
    print $fh "$sequence\n";

  # Feature file
  } elsif ($ft_file) {
    open FEAT, $ft_file;
    my $ft_name;
    while ($line = <FEAT>) {
      chomp($line);
      next if (($line =~/^[#|;]/)||($line eq ""));
      if ($ft_file_format eq "ft") {
	($chrom, $ft_type, $ft_id, $strand, $left_limit, $right_limit,@other_comments) = split (/\t/,$line);
      } elsif ($ft_file_format eq 'gft') {
	($ft_id, $ft_type, $ft_name, $chrom, $left_limit, $right_limit, $strand, @other_comments) = split (/\t/,$line);
      }

      # Extract only chromosome number if necessary
      $chrom =~ s/chromosome:[\w\.]*?://;
      $chrom =~ s/:.*//;

      # Tranforms strand in ensembl format
      $strand =~ s/F/1/;
      $strand =~ s/R/-1/;
      $strand =~ s/D/1/;
      $strand =~ s/W/1/;
      $strand =~ s/C/-1/;
      $strand =~ s/>/1/;
      $strand =~ s/</-1/;

      local $chromosome = $slice_adaptor -> fetch_by_region('chromosome', $chrom);
      # Get sequence (repeat masked or not)
      $sequence = &GetSequence($left_limit, $right_limit);
      my $size = $right_limit - $left_limit + 1;

      my $rsat_strand;
      if ($strand == 1) {
	$rsat_strand = "D";
      } elsif ($strand == -1) {
	$rsat_strand = "R";
      }

      # Export sequence to file
      if ($ft_id) {
	print $fh ">$ft_id\t$ft_id; from 1 to $size; size: $size; location: $chrom $left_limit $right_limit $rsat_strand\n";
      } else {
	print $fh ">$chrom-$left_limit-$right_limit\t$chrom-$left_limit-$right_limit; from 1 to $size; size: $size; location: $chrom $left_limit $right_limit $rsat_strand\n";
      }
      print $fh  "$sequence\n";
    }

    # All genes
  } elsif ($all) {
    my @slices = @{$slice_adaptor->fetch_all("chromosome")};
    foreach my $slice (@slices) {
      my @genes = @{$slice->get_all_Genes()};
      foreach my $gene (@genes) {
	&Main($gene);
      }
    }
    # Query
  } else {
#    my $gene_id = "CG40293"; #gene on D strand
#    my $gene_id = "CG18001"; #gene on R strand
#    @genes = @{$gene_adaptor->fetch_all_by_external_name('BRCA2')};

    my $gene_adaptor = $db->get_GeneAdaptor();
    # Input file
    if ($query_file) {
	open IN, $query_file;
	while ($line = <IN>) {
	    my $gene;
	    $line =~s/\t//;
	    chomp($line);
	    if ($line =~ /ENST/) {
		$gene = $gene_adaptor -> fetch_by_transcript_stable_id($line);
	    } elsif ($line =~ /ENSP/) {
		$gene = $gene_adaptor -> fetch_by_translation_stable_id($line);
	    } else {
		my $gene_id = $line;
		$gene = $gene_adaptor -> fetch_by_stable_id($gene_id);
	    }
	    &Main($gene);
	}
    } else {
	foreach my $id (@queries) {
	    my $gene;
	    if ($id =~ /ENST/) {
		$gene = $gene_adaptor -> fetch_by_transcript_stable_id($id);
	    } elsif ($id =~ /ENSP/) {
		$gene = $gene_adaptor -> fetch_by_translation_stable_id($id);
	    } else {
		$gene = $gene_adaptor -> fetch_by_stable_id($id);
	    }
	    &Main($gene);
	}
    }
}


  ################################################################
  ## Report execution time
  my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
  warn $exec_time if ($main::verbose >= 1); ## only report exec time if verbosity is specified

  ################################################################
  ## Close output stream
  close $fh if ($output_file);

  exit(0);
}
################################################################
################### SUBROUTINE DEFINITION ######################
################################################################

################################################################
#### read arguments
sub ReadArguments {
    foreach my $a (0..$#ARGV) {
	### verbose  
	if ($ARGV[$a] eq "-v") {
#	    if (&IsNatural($ARGV[$a+1])) {
		$verbose = $ARGV[$a+1];
#	    } else {
#		$verbose = 1;
#	    }
	    
	    ### detailed help
	} elsif ($ARGV[$a] eq "-h") {
	    &PrintHelp();
	    
	    ### list of options
	} elsif ($ARGV[$a] eq "-help") {
	    &PrintShortHelp();
	    
	    ### output file
	} elsif ($ARGV[$a] eq "-o") {
	    $output_file = $ARGV[$a+1];
	    
	    ### EnsEMBL database server (host)
	} elsif ($ARGV[$a] eq "-ensemblhost") {
	    $ensembl_host = $ARGV[$a+1];
	    
	    ### EnsEMBL database name
	} elsif ($ARGV[$a] eq "-dbname") {
	    $dbname = $ARGV[$a+1];
	    
	    ### EnsEMBL database version
	} elsif ($ARGV[$a] eq "-dbversion") {
	    $dbversion = $ARGV[$a+1];
	    
	    ### organism
	} elsif ($ARGV[$a] eq "-org") {
	    $org = lc($ARGV[$a+1]);
	    
	    ### Limits of region to extract
	} elsif ($ARGV[$a] eq "-from") {
	    $from = $ARGV[$a+1];
	    
	} elsif ($ARGV[$a] eq "-to") {
	    $to = $ARGV[$a+1];
	    
	    ### chromosome name
	} elsif ($ARGV[$a] eq "-chrom") {
	    $chrom = $ARGV[$a+1];
	    
	    ### Left limit
	} elsif ($ARGV[$a] eq "-left") {
	    $left_limit = $ARGV[$a+1];
	    
	    ### Right limit
	} elsif ($ARGV[$a] eq "-right") {
	    $right_limit = $ARGV[$a+1];
	    
	    ### Strand
	} elsif ($ARGV[$a] eq "-strand") {
	    $strand = $ARGV[$a+1];
	    
	    ### Feature file
	} elsif ($ARGV[$a] eq "-ftfile") {
	    $ft_file = $ARGV[$a+1];
	    
	    ### Feature file format
	} elsif ($ARGV[$a] eq "-ftfileformat") {
	    $ft_file_format = $ARGV[$a+1];
	    
	    ### Sequence type
	} elsif ($ARGV[$a] eq "-type") {
	    $type = $ARGV[$a+1];
	    
	    ### Feature type
	} elsif ($ARGV[$a] eq "-feattype") {
	    $feattype = lc($ARGV[$a+1]);
	    
	    ### Noorf
	} elsif ($ARGV[$a] eq "-noorf") {
	    $noorf = 1;
	    
	    ### Mask coding
	} elsif ($ARGV[$a] eq "-maskcoding") {
	    $mask_coding = 1;
	    
	    ### Nogene
	} elsif ($ARGV[$a] eq "-nogene") {
	    $nogene = 1;
	    
	    ### Query
	} elsif ($ARGV[$a] eq "-q") {
	    @queries = (@queries, $ARGV[$a+1]);
	    
	    ### Query file
	} elsif ($ARGV[$a] eq "-i") {
	    $query_file = $ARGV[$a+1];
	    
	    ### All genes
	} elsif ($ARGV[$a] eq "-all") {
	    $all = 1;
	    
	    ### All transcripts
	} elsif ($ARGV[$a] eq "-alltranscripts") {
	    $all_transcripts = 1;
	    
	    ### Repeat mask
	} elsif ($ARGV[$a] eq "-rm") {
	    $rm = 1;
	    
	    ### First intron
	} elsif ($ARGV[$a] eq "-firstintron") {
	    $first_intron = 1;
	    
	    ### UTRs
#    } elsif ($ARGV[$a] eq "-utr") {
#      $utr = $ARGV[$a+1];

	    ### Non coding
	} elsif ($ARGV[$a] eq "-noncoding") {
	    $non_coding = 1;
	}
    }
}
################################################################
#### verbose message
#sub Verbose {
#  print "; retrieve-ensembl-seq.pl ";
#  &PrintArguments();
#}
################################################################
### Get sequence(s) relative to a feature
sub Main {
  my ($gene) = @_;
  my $gene_id = $gene -> stable_id();
  my $gene_name = $gene -> external_name();
  unless ($gene_name) {
    $gene_name = "";
  }

  local $chromosome = $gene -> slice;
  my $chromosome_name = $chromosome -> name();

  my $gene_start = $gene -> start();
  my $gene_end = $gene -> end();
  $strand = $gene -> strand();

  my $description = $gene -> description();

  my $rsat_strand;
  if ($strand == 1) {
    $rsat_strand = "D";
  } else {
    $rsat_strand = "R";
  }

  if ($feattype eq "gene") {
    my ($left, $right) = &GetLimits($gene_id, $gene_start, $gene_end);

    # Get sequence (repeat masked or not)
    $sequence = &GetSequence($left, $right);
    my $size = $new_to - $new_from + 1;

    # Export sequence to file
    print $fh ">$gene_id\t$gene_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand\n";
    print $fh "$sequence\n";

  } else { # feattype = mrna, cds, introns...
    # Get transcripts
    @transcripts = @{$gene->get_all_Transcripts};

    # Initialize transcripts limits
    my $three_primest_id = $transcripts[0] -> display_id();
    my $five_primest_id = $transcripts[0] -> display_id();
    my $three_primest_start;
    my $five_primest_end;
    if ($strand == 1) {
      $three_primest_start = $transcripts[0] -> start();
      $five_primest_end = $transcripts[0] -> end();
    } else {
      $three_primest_start = $transcripts[0] -> end();
      $five_primest_end = $transcripts[0] -> start();
    }

    foreach my $transcript(@transcripts) {
      my $transcript_id = $transcript -> display_id();
      my $transcript_start = $transcript -> start();
      my $transcript_end = $transcript -> end();

      # Find all transcripts minimal limits
      if ($strand == 1 && $transcript_start > $three_primest_start) {
	$three_primest_id = $transcript_id;
	$three_primest_start = $transcript_start;
      } elsif ($strand == -1 && $transcript_end < $three_primest_start) {
	$three_primest_id = $transcript_id;
	$three_primest_start = $transcript_end;
      }
      if ($strand == 1 && $transcript_end < $five_primest_end) {
	$five_primest_id = $transcript_id;
	$five_primest_end = $transcript_end;
      } elsif ($strand == -1 && $transcript_start > $five_primest_end) {
	$five_primest_id = $transcript_id;
	$five_primest_end = $transcript_start;
      }

      if ($feattype eq 'mrna' && $all_transcripts) {
	my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $transcript_start, $transcript_end);

	# Output sequence
	$sequence = &GetSequence($left, $right);
	my $size = $new_to - $new_from + 1;

	# Export sequence to file
	print $fh ">$gene_id-$transcript_id\t$gene_id-$transcript_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand\n";
	print $fh "$sequence\n";
      }

      # Introns
      if ($feattype eq "intron") {
	my $i = 0; # index for introns since they don't have an id in ensembl
	my $intron1;
	my $intron1_id;
	my $start1;
	my $end1;
	if ($strand == 1) {
	  $start1 = $transcript_end;
	} else {
	  $end1 = $transcript_start;
	}

	foreach my $intron (@{$transcript -> get_all_Introns()}) {
	  $i++;
	  my $intron_start = $intron -> start();
	  my $intron_end = $intron -> end();

	  if ($first_intron) {
	    if ($strand == 1) {
	      if ($intron -> start() < $start1) {
		$intron1 = $intron;
		$intron1_id = $transcript_id."-".$i;
		$start1 = $intron -> start();
		$end1 = $intron -> end();
	      }
	    } else {
	      if ($intron -> end() > $end1) {
		$intron1 = $intron;
		$intron1_id = $transcript_id."-".$i;
		$start1 = $intron -> start();
		$end1 = $intron -> end();
	      }
	    }
	  }

	  unless ($first_intron) {
	    $sequence = &GetSequence($intron -> start(), $intron -> end());
	    my $size = ($intron -> end() - $intron -> start()) + 1;

	    # Export sequence to file
	    print $fh ">$gene_id-$transcript_id-$i\t$gene_id-$transcript_id-$i; from 1 to $size; size: $size; location: $chromosome_name $intron_start $intron_end $rsat_strand\n";
	    print $fh "$sequence\n";
	  }
	}

	if ($first_intron) {

	  $sequence = &GetSequence($start1, $end1);
	  my $size = ($end1 - $start1) + 1;

	  # Export sequence to file
	  print $fh ">$gene_id-$intron1_id\t$gene_id-$intron1_id; from 1 to $size; size: $size; location: $chromosome_name $start1 $end1 $rsat_strand\n";
	  print $fh "$sequence\n";
	}
      }

      my $coding_region_start = $transcript->coding_region_start();
      my $coding_region_end = $transcript->coding_region_end();

      # Exons
      if ($feattype eq "exon") {
	foreach my $exon (@{$transcript -> get_all_Exons()}) {
	  my $exon_id = $exon -> stable_id();
	  my $exon_start = $exon -> start();
	  my $exon_end = $exon -> end();

	  if ($non_coding) {
	    if ($coding_region_start > $exon_start && $coding_region_start < $exon_end) {
	      $sequence = &GetSequence($exon -> start(), $coding_region_start - 1);
	      my $non_coding_exon_right = $coding_region_start - 1;
	      my $size = $coding_region_start - $exon -> start();

	      # Export sequence to file
	      print $fh ">$gene_id-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $exon_start $non_coding_exon_right $rsat_strand\n";
	      print $fh "$sequence\n";
	    } elsif ($coding_region_end < $exon_end && $coding_region_end > $exon_start) {
	      $sequence = &GetSequence($coding_region_end + 1, $exon_end);
	      my $non_coding_exon_left = $coding_region_end + 1;
	      my $size = $exon -> end() - $coding_region_end;

	      # Export sequence to file
	      print $fh ">$gene_id-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $non_coding_exon_left $exon_end $rsat_strand\n";
	      print $fh "$sequence\n";
	    } elsif ($coding_region_start > $exon_end) {
	      $sequence = &GetSequence($exon -> start(), $exon -> end());
	      my $size = ($exon -> end() - $exon -> start()) + 1;

	      # Export sequence to file
	      print $fh ">$gene_id-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $exon_start $exon_end $rsat_strand\n";
	      print $fh "$sequence\n";
	    } elsif ($coding_region_end < $exon_start) {
	      $sequence = &GetSequence($exon -> start(), $exon -> end());
	      my $size = ($exon -> end() - $exon -> start()) + 1;

	      # Export sequence to file
	      print $fh ">$gene_id-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $exon_start $exon_end $rsat_strand\n";
	      print $fh "$sequence\n";
	    }
	  } else {
	    $sequence = &GetSequence($exon -> start(), $exon -> end());
	    my $size = ($exon -> end() - $exon -> start()) + 1;

	    # Export sequence to file
	    print $fh ">$gene_id-$exon_id\t$gene_id-$exon_id; from 1 to $size; size: $size; location: $chromosome_name $exon_start $exon_end $rsat_strand\n";
	    print $fh "$sequence\n";
	  }
	}
      }

      # UTR
      if ($feattype eq "utr") {
	my $utr5_start;
	my $utr5_end;
	my $utr3_start;
	my $utr3_end;
	if ($strand == 1) {
	  $utr5_start = $transcript_start;
	  $utr5_end = $coding_region_start - 1;
	  $utr3_start = $coding_region_end + 1;
	  $utr3_end = $transcript_end;
	} else {
	  $utr3_start = $transcript_start;
	  $utr3_end = $coding_region_start - 1;
	  $utr5_start = $coding_region_end + 1;
	  $utr5_end = $transcript_end;
	}
	$utr5_sequence = &GetSequence($utr5_start, $utr5_end);
	$utr3_sequence = &GetSequence($utr3_start, $utr3_end);
	my $utr5_size = $utr5_end - $utr5_start + 1;
	my $utr3_size = $utr3_end - $utr3_start + 1;

	# Export sequence to file
	print $fh ">$gene_id-$transcript_id-5prime_UTR\t$gene_id-$transcript_id-5prime_UTR; from 1 to $utr5_size; size: $utr5_size; location: $chromosome_name $utr5_start $utr5_start $rsat_strand\n";
	print $fh "$utr5_sequence\n";
	print $fh ">$gene_id-$transcript_id-3prime_UTR\t$gene_id-$transcript_id-3prime_UTR; from 1 to $utr3_size; size: $utr3_size; location: $chromosome_name $utr3_start $utr3_start $rsat_strand\n";
	print $fh "$utr3_sequence\n";
      }

      # CDS
      if ($feattype eq 'cds') {
	my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $coding_region_start, $coding_region_end);

	# Output sequence
	$sequence = &GetSequence($left, $right);
	my $size = $new_to - $new_from + 1;
	my $cds_id = $transcript -> translation() -> stable_id();

	# Export sequence to file
	print $fh ">$gene_id-$transcript_id-$cds_id\t$gene_id-$transcript_id-$cds_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand\n";
	print $fh "$sequence\n";
      }
    }

    if ($feattype eq 'mrna' && !$all_transcripts) {

      my $ref_transcript;

      if ($strand == 1) {
	my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $three_primest_start, $five_primest_end);
      } else {
	my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $five_primest_end, $three_primest_start);
      }

      # Output sequence
      $sequence = &GetSequence($left, $right);
      my $size = $new_to - $new_from + 1;
      if ($type eq "upstream") {
	$ref_transcript = $three_primest_id;
      } else {
	$ref_transcript = $five_primest_id;
      }

      # Export sequence to file
      print $fh ">$gene_id-$ref_transcript\t$gene_id-$ref_transcript; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand\n";
      print $fh "$sequence\n";
    }
  }
}
################################################################
#### Calculate left and right limits
sub GetLimits {
  my ($gene_id, $start, $end) = @_;

  if ($type eq "upstream" && $strand == 1) {
    if ($nogene || $noorf) {
      # Get expanded slice and identify closest neighbour
      $expanded_slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome -> seq_region_name(), $start + $from, $end);
      $closest_neighbour_limit = &GetNeighbours($gene_id, $start, $end, $expanded_slice);
    }
#    if (($nogene || $noorf) && $closest_neighbour_limit > $start + $from) {
    if (($nogene || $noorf) && $closest_neighbour_limit > $start + $from && $closest_neighbour_limit < $start) {
      $left = $closest_neighbour_limit + 1;
    }else {
      $left = $start + $from;
    }
    $right = $start + $to;
    # Test if whithin chromosome limits
    if ($left < 1) {
      $left = 1;
    }
    if ($right < 1) {
      $right = 1;
    }
    if ($left > $chromosome->end()) {
      $left = $chromosome->end();
    }
    if ($right > $chromosome->end()) {
      $right = $chromosome->end();
    }
    # Recalculate from and to (in case left and/or right had to be reassigned above)
    $new_from = $left - $start;
    $new_to = $right - $start;
  } elsif ($type eq "upstream" && $strand == -1) {
    $left = $end - $to;
    if ($nogene || $noorf) {
      # Get expanded slice and identify closest neighbour
      $expanded_slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome -> seq_region_name(), $start, $end - $from);
      $closest_neighbour_limit = &GetNeighbours($gene_id, $start, $end, $expanded_slice);

    }
#    if (($nogene || $noorf) && $closest_neighbour_limit < $end - $from) {
    if (($nogene || $noorf) && $closest_neighbour_limit < $end - $from && $closest_neighbour_limit > $end) {
      $right = $closest_neighbour_limit - 1;
    } else {
      $right = $end - $from;
    }
    # Test if whithin chromosome limits
    if ($left < 1) {
      $left = 1;
    }
    if ($right < 1) {
      $right = 1;
    }
    if ($left > $chromosome->end()) {
      $left = $chromosome->end();
    }
    if ($right > $chromosome->end()) {
      $right = $chromosome->end();
    }
    # Recalculate from and to (in case left and/or right had to be reassigned above)
    $new_from = $end - $right;
    $new_to = $end - $left;
  } elsif ($type eq "downstream" && $strand == 1) {
    $left = $end + $from;
    if ($nogene || $noorf) {
      # Get expanded slice and identify closest neighbour
      $expanded_slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome -> seq_region_name(), $start, $end + $to);
      $closest_neighbour_limit = &GetNeighbours($gene_id, $start, $end, $expanded_slice);

    }
#    if (($nogene || $noorf) && $closest_neighbour_limit < $end + $to) {
    if (($nogene || $noorf) && $closest_neighbour_limit < $end + $to && $closest_neighbour_limit > $end) {
      $right = $closest_neighbour_limit - 1;
    } else {
      $right = $end + $to;
    }
    # Test if whithin chromosome limits
    if ($left < 1) {
      $left = 1;
    }
    if ($right < 1) {
      $right = 1;
    }
    if ($left > $chromosome->end()) {
      $left = $chromosome->end();
    }
    if ($right > $chromosome->end()) {
      $right = $chromosome->end();
    }
    # Recalculate from and to (in case left and/or right had to be reassigned above)
    $new_from = $left - $end;
    $new_to = $right - $end;
  } elsif ($type eq "downstream" && $strand == -1) {
    if ($nogene || $noorf) {
      # Get expanded slice and identify closest neighbour
      $expanded_slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome -> seq_region_name(), $start - $to, $end);
      $closest_neighbour_limit = &GetNeighbours($gene_id, $start, $end, $expanded_slice);

    }
#    if (($nogene || $noorf) && $closest_neighbour_limit > $start - $to) {
    if (($nogene || $noorf) && $closest_neighbour_limit > $start - $to && $closest_neighbour_limit < $start) {
      $left = $closest_neighbour_limit + 1;
    } else {
      $left = $start - $to;
    }
    $right = $start - $from;
    # Test if whithin chromosome limits
    if ($left < 1) {
      $left = 1;
    }
    if ($right < 1) {
      $right = 1;
    }
    if ($left > $chromosome->end()) {
      $left = $chromosome->end();
    }
    if ($right > $chromosome->end()) {
      $right = $chromosome->end();
    }
    # Recalculate from and to (in case left and/or right had to be reassigned above)
    $new_from = $start - $right;
    $new_to = $start - $left;
  }
  return ($left, $right, $new_from, $new_to);
}
################################################################
#### Get neighbouring genes
sub GetNeighbours {
  my ($gene_id, $start, $end, $expanded_slice) = @_;

  # Initialize closest neighbour limit
  my $closest_neighbour_limit;
  if ($strand == 1 && $type eq "upstream") {
    $closest_neighbour_limit = $start - abs($from);
  } elsif ($strand == 1 && $type eq "downstream") {
    $closest_neighbour_limit = $end + abs($to);
  } elsif ($strand == -1 && $type eq "upstream") {
    $closest_neighbour_limit = $end + abs($from);
  } elsif ($strand == -1 && $type eq "downstream") {
    $closest_neighbour_limit = $start - abs($to);
  }

  # Get neighbour genes
  my @neighbour_genes = @{$expanded_slice->get_all_Genes()};

  foreach my $neighbour_gene(@neighbour_genes) {
    if ($neighbour_gene->stable_id() eq $gene_id) {    # the query gene itself is in the expanded slice
      next
    }
    my $neighbour_name = $neighbour_gene->external_name();
    unless ($neighbour_name) {
      $neighbour_name = "";
    }

    # This is for debug only
    my $neighbour_strand = $neighbour_gene->seq_region_strand();
    if ($neighbour_strand == 1) {
      $neighbour_strand = "D";
    } else {
      $neighbour_strand = "R";
    }
    my $neighbour_start = $neighbour_gene -> seq_region_start();
    my $neighbour_end = $neighbour_gene -> seq_region_end();

    #Find closest neighbour limit
    if ($nogene) {    # neighbour limits are closest gene limits
      if ($strand == 1 && $type eq "upstream") {
	if ($neighbour_gene->seq_region_end() > $closest_neighbour_limit) {
	  $closest_neighbour_limit = $neighbour_gene->seq_region_end();
	}
      } elsif ($strand == 1 && $type eq "downstream") {
	if ($neighbour_gene->seq_region_start() < $closest_neighbour_limit) {
	  $closest_neighbour_limit = $neighbour_gene->seq_region_start();
	}
      } elsif ($strand == -1 && $type eq "upstream") {
	if ($neighbour_gene->seq_region_start() < $closest_neighbour_limit) {
	  $closest_neighbour_limit = $neighbour_gene->seq_region_start();
	}
      } elsif ($strand == -1 && $type eq "downstream") {
	if ($neighbour_gene->seq_region_end() > $closest_neighbour_limit) {
	  $closest_neighbour_limit = $neighbour_gene->seq_region_end();
	}
      }
    } elsif ($noorf) {    # neighbour limits are closest CDS limits
      my @transcripts = @{$neighbour_gene->get_all_Transcripts};
      foreach my $transcript(@transcripts) {
	my $coding_region_start = $transcript->coding_region_start();
	my $coding_region_end = $transcript->coding_region_end();
	if (($strand == 1 && $type eq "upstream") || ($strand == -1 && $type eq "downstream")) {
	  if ($type eq "upstream") {
	  $coding_region_end = $coding_region_end + $start - abs($from) - 1;
	} else {
	  $coding_region_end = $coding_region_end + $start - abs($to) - 1;
	}
	  if ($coding_region_end > $closest_neighbour_limit) {
	    $closest_neighbour_limit = $coding_region_end;
	  }
	} elsif (($strand == 1 && $type eq "downstream") || ($strand == -1 && $type eq "upstream")) {
	  $coding_region_start = $coding_region_start + $start - 1;
	  if ($coding_region_start < $closest_neighbour_limit) {
	    $closest_neighbour_limit = $coding_region_start;
	  }
	}
      }
    }
  }
  return ($closest_neighbour_limit);
}
################################################################
#### Get sequence
sub GetSequence {
  my ($left, $right) = @_;
  my $sequence;
  unless ($rm) {
    $sequence = $chromosome->subseq($left, $right);
  }
  if ($rm) {
#    my @repeat_features = @{$chromosome->get_all_RepeatFeatures()};
#    foreach my $repeat_feature(@repeat_features) {
#      my $repeat_consensus = $repeat_feature->repeat_consensus;
#      my $repeat_type = $repeat_consensus->repeat_type();
#      my $repeat_class = $repeat_consensus->repeat_class();
#      print $repeat_type, "\t", $repeat_class, "\n";
#    }
    my $mini_slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome -> seq_region_name(), $left, $right);
    my $masked_slice = $mini_slice->get_repeatmasked_seq();
#    my $masked_slice = $mini_slice->get_repeatmasked_seq([],0,{"repeat_class_Dust" => 0});
#    my $masked_slice = $mini_slice->get_repeatmasked_seq(['TRF']);
    $sequence = $masked_slice->seq();
#    my $unmasked_seq = $slice->seq();
#    my $hardmasked_seq = $slice->get_repeatmasked_seq();
#    my $softmasked_seq = $slice->get_repeatmasked_seq(undef, 1);
  }
  if ($mask_coding) {
    my $retrieved_slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome -> seq_region_name(), $left, $right);
    my $transcript_adaptor = $db->get_TranscriptAdaptor();
    my @retrieved_transcripts = @{$transcript_adaptor->fetch_all_by_Slice($retrieved_slice)};
    foreach $retrieved_transcript (@retrieved_transcripts) {
      my @retrieved_exons = @{$retrieved_transcript->get_all_translateable_Exons};
      foreach my $retrieved_exon (@retrieved_exons) {
	my $coding_length = $retrieved_exon->end() - $retrieved_exon->start() + 1;
	my $masking;
	if ($retrieved_exon->seq_region_start() >= $left && $retrieved_exon->seq_region_end() <= $right) {
	  foreach (1..$coding_length) {
	    $masking .= "N";
	  }
	  substr($sequence, $retrieved_exon->start() - 1, $coding_length) = $masking;
	} elsif ($retrieved_exon->seq_region_start() < $left && $retrieved_exon->seq_region_end() >= $left && $retrieved_exon->seq_region_end() <= $right) {
	  $coding_length = $retrieved_exon->end();
	  foreach (1..$coding_length) {
	    $masking .= "N";
	  }
	  substr($sequence, 0, $coding_length) = $masking;
	} elsif ($retrieved_exon->seq_region_start() >= $left && $retrieved_exon->seq_region_start() <= $right && $retrieved_exon->seq_region_end() > $right) {
	  $coding_length = $right - $retrieved_exon->seq_region_start() + 1;
	  foreach (1..$coding_length) {
	    $masking .= "N";
	  }
	  substr($sequence, $retrieved_exon->start() - 1, $coding_length) = $masking;
	} elsif ($retrieved_exon->seq_region_start() < $left && $retrieved_exon->seq_region_end() > $right) {
	  $coding_length = $right - $left + 1;
	  foreach (1..$coding_length) {
	    $masking .= "N";
	  }
	  substr($sequence, 0, $coding_length) = $masking;
	}
      }
    }
  }
#  if ($strand == -1) {
#    $sequence = &ReverseComplement($sequence);
#  }
  return $sequence;
}
################################################################
#### detailed help message
sub PrintHelp {
#    my $supported_organisms = &ListSupportedOrganisms("sizes");
    open(HELP, "| less");
    print HELP<<End_help;
USAGE
	retrieve-ensembl-seq [-type type] -org organism | -dbname database
			[-from] [-to] [-noorf] [-rm] [-maskcoding]
			[-o outpufile] -q query_orf | -i query file | -all

DESCRIPTION
	Returns upstream, downstream, intronic, exonic or UTR  DNA sequences for list of query
	genes.

CATEGORY
	genomics
	sequences

REMARK  Requires local instal of the EnsEMBL Perl API (see http://www.ensembl.org/info/software/api_installation.html)

OPTIONS
	-org organism
	        No caps, underscore between words (eg 'homo_sapiens')

	        If this option is not used, the option -dbname must be used
	         instead.

	        (type 'supported-organism | grep EnsEMBL' to obtain the list of supported
	         organisms)

        -ensemblhost
                address of ensembl database server (default is EBI server)

	-dbname	name of EnsEMBL database
		(alternative to organism)

	-feattype
		Feature type.
		Supported: gene, mrna (default), CDS, intron, exon, utr

	-type	sequence type
		Currently supported sequence types
			upstream (default)
			downstream

	-q query
		The query should be an EnsEMBL gene identifier (eg 'ENSG00000177799').
		Multiple queries can be entered by reiteratively using the -q
		option.

        -i     query file. The first word of each line is taken as a query.
                This option is incompatible with -q.

	-all	return all genomic upstream/downstream regions

	-o	name of the output file

        -from #1 -to #2
                where #1 and #2 are numbers. #2 should be higher than #1.
                limits of the region to extract, relative to feattype start or end
                (=position 0). Use negative values for upstream sequence.
                        example: -from -800 -to -1
                        will extract the 800 bp upstream the feattype start or end.
			 (this is the default

	-noorf	the upstream/downstream sequence can only contain non-coding sequence.
		i.e. the -from values is modified if a predicted orf
		is encountered within its range.
		The weaknesses of using this option are that
		- all predicted orf do not correspond to real orf,
		- there is no a priori reason to exclude a regulatory site
		  which would overlap the upstream coding sequence.

        -nogene the upstream/downstream sequence can only contain non-transcribed sequence.

        -maskcoding
                all coding sequence is replaced by N in the retrieved sequence

	-rm     Use the repeat masked version of the genome.  Attention :
		repeated regions are annotated for some genomes only.

        -alltranscripts
                Get sequences for all transcript of genes.
                Use purge-sequence if you do motif discovery afterwards

        -firstintron
                With feattype intron, get only first intron sequence

        -noncoding
                With feattype exon, get only non-coding (part of) exons

        -chrom  Chromosome name or number (to use with -left and -right)

        -left   Left limit of sequence to retrieve

        -right  Right limit of sequence to retrieve

        -strand Strand of seauence to retrieve when using -left and -right. Values: 1, -1

        -ftfile Feature file

        -ftfileformat
                Feature file format. Supported: ft, gft

End_help
    close HELP;
    exit;
}

################################################################
#### list of options
sub PrintShortHelp {
  open(HELP, "| less");
  print HELP<<End_short_help;
retrieve-seq options
--------------------
-org		organism
-ensemblhost    address of ensembl database server (default is EBI server)
-dbname         name of ensembl db
-feattype	accepted feature types. Supported: gene, mrna, cds, intron, exon, utr
-type		upstream | downstream | orf | random
-q		query
-i              query file
-all		returns all genomic upstream regions
-o		followed by the name of the outputfile.
-from #1 -to #2	limits of the region to extract, relative to feattype start (upstream) or end (downstream)
-noorf		the upstream/downstream sequence can only contain non-coding sequence.
-nogene         the upstream/downstream sequence can only contain non-transcribed sequence.
-maskcoding     all coding sequence is replaced by Ns in the retrieved sequence
-rm		Use the repeat masked version of the genome.
-alltranscripts get sequences for all transcript of genes
-firstintron    with feattype intron, get only first intron sequence
-noncoding      with feattype exon, get only non-coding (part of) exons
-chrom          chromosome name or number (to use with -left and -right)
-left           left limit of sequence to retrieve
-right          right limit of sequence to retrieve
-strand         strand of seauence to retrieve when using -left and -right. Values: 1, -1
-ftfile         feature file
-ftfileformat   feature file format. Supported: ft, gft
End_short_help
  close HELP;
  exit;
}

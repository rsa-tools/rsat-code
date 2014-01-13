#!/usr/bin/perl -w
############################################################
#
# $Id: retrieve-ensembl-seq.pl,v 1.82 2013/09/02 12:36:15 rsat Exp $
#
# Time-stamp
#
############################################################
#use strict;
use DBI();
use Data::Dumper;

BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
}
require "RSA.lib";
require "RSA.seq.lib";
require RSAT::util;
push (@INC, $ENV{bioperl});
push (@INC, $ENV{ensembl});
push (@INC, $ENV{compara});

## EnsEMBL libraries
require Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::DBSQL::SliceAdaptor;
#use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

################################################################
#### main package
package main;
{

  ################################################################
  #### initialise parameters
  local $start_time = &RSAT::util::StartScript();

  local $verbose = 0;
  local $feattype = "mrna";    # other values: gene, intron, exon, cds and utr
  local %supported_feattype = ("gene"=>1,
			       "transcript"=>1,
			       "mrna"=>1,
			       "cds"=>1,
			       "intron"=>1,
			       "exon"=>1,
			       "utr"=>1);
  local $supported_feattypes = join (",", sort(keys %supported_feattype));

  local $type = "upstream";
  local $from = -800;
  local $to = -1;
  local $noorf = 0;
  local $nogene = 0;
  local $rm = 0;
  local $all_transcripts = 0;
  local $first_intron = 0;
  local $utr = "all";    # other values: 5prime and 3prime
  local $non_coding = 0;
  local $all = 0;
  local $query_file;
  local @queries;
  local $left_limit = 1;
  local $right_limit;
  local $strand = 1;
  local $chrom;
  local $ft_file;
  local $ft_file_format = "gft";
  local $mask_coding = 0;
  local $lw = 0;
  local $uniq_seqs = 0;

  local $output_file;
  local $homologs_table;

  local $ortho;
  local $ortho_type;
  local $taxon;

  local $common_name = '';
  local $header = 'scientific';
  local $header_org = '';
  local $label = '';
  local %supported_label = ("query"=>1);
  local $supported_labels = join (",", sort(keys(%supported_label)));

  ## Connection to the EnsEMBL MYSQL database
  local $ensembl_host = $ENV{ensembl_host};
#  local $ensembl_host = 'ensembldb.ensembl.org';  # db at EBI (use outside BIGRE)
#  local $ensembl_host = 'xserve2.bigre.ulb.ac.be';  # Local db (use inside BIGRE)
  local $ensembl_user = "anonymous";
  local $dbname = '';
  local $org = '';
  local $dbversion = '';
  local $port = '';

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
#    &PrintArguments;
    print "\n";
  }

  ################################################################
  ## Get ensembl mysql port from db version
  if ($dbversion && ($dbversion < '48')) {
      $port = '3306';
  } else {
      $port = '5306';
  }

  ################################################################
  ## If option -org is used, connect to ensembldb to get list of 
  ## databases and pick the latest one corresponding to chosen organism
  if ($org) {
      &RSAT::message::TimeWarn (join("\t", "Connecting EnsEMBL to get the dbname for organism ", $org, 
				     "host=".$ensembl_host, 
				     "user=".$ensembl_user )) if ($main::verbose >= 1);
      my $dbh = DBI->connect("DBI:mysql:host=$ensembl_host:port=$port", "$ensembl_user", "", {'RaiseError' => 1});
      my $sth = $dbh->prepare("SHOW DATABASES");
      $sth->execute();
      while (my $ref = $sth->fetchrow_hashref()) {
	  if ($ref->{'Database'} =~ /($org)_core_\d+/) {
	      $dbname = $ref->{'Database'};
	  }
      }
      $sth->finish();
      $dbh->disconnect();
      unless ($dbname) {
	  die "Error: there is no organism named $org in the EnsEMBL database. Use the command supported-organisms-ensembl to obtain a full list of supported organisms.\n";
      }
  } else {  # get organism name from dbname
      $org = $dbname;
      $org =~s/_core_.+//;
  }

  &RSAT::message::Info (join("\t", "dbname = ", $dbname)) if ($main::verbose >= 1);

  ################################################################
  ## Get EnsEMBL db version from db name
  unless ($dbversion) {
      $dbversion = $dbname;
      $dbversion =~ s/($org)_core_//;
      $dbversion =~ s/_.+//;
  }

  &RSAT::message::Info (join("\t", "dbversion", $dbversion)) if ($main::verbose >= 1);

  ################################################################
  ## Open a new connection to EnsEMBL database, but this time we specify the DB name
    &RSAT::message::TimeWarn("Connecting EnsEMBL to retrieve the organism", 
			     "host=".$ensembl_host,
			     "user=".$ensembl_user,
			     "dbname=".$dbname,
			     ) if ($main::verbose >= 1);

  my $registry = "Bio::EnsEMBL::Registry";

  $registry->load_registry_from_db(
				   -host => $ensembl_host,
				   -user => $ensembl_user,
				   -db_version => $dbversion,
                                   -port => $port,
				   -verbose => "0" );

  local $db = Bio::EnsEMBL::Registry->get_DBAdaptor($org, "core");

# local $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $ensembl_host, -user => $ensembl_user, -dbname => $dbname); ## deprecated

    &RSAT::message::Debug("Db", $db) if ($main::verbose >= 3);

  ################################################################
  ### Open output stream
  $fh = &OpenOutputFile($output_file);
#  if ($output_file) {
#    $fh = 'OUT';
#    open $fh, ">".$output_file || die "cannot open file ".$output_file."\n";
#  } else {
#    $fh = *STDOUT;
#  }

  if ($homologs_table) {
      $table_handle = 'TAB';
      open $table_handle, ">".$homologs_table || die "cannot open file ".$homologs_table."\n";
#  } else {
#      $table_handle = *STDOUT;
  }

  #### print verbose
  &Verbose() if ($verbose);

  ################################################################
  ## Get query (one gene, several genes or all genes)

  local $slice_adaptor = $db->get_sliceAdaptor();

  ## Left and right limits
  if ($left_limit && $right_limit && $chrom) {
    local $chromosome = $slice_adaptor->fetch_by_region('chromosome', $chrom);
    ## Get sequence (repeat masked or not)
    $sequence = &GetSequence($left_limit, $right_limit);
    my $size = $right_limit - $left_limit + 1;

    my $rsat_strand;
    if ($strand == 1) {
      $rsat_strand = "D";
    } else {
      $rsat_strand = "R";
    }

    if ($header eq 'none') {
	$header_org = '';
    } else {
	$header_org = ucfirst($org)."-";
    }

    my $fasta_header = ">$header_org$chrom-$left_limit-$right_limit\t$chrom-$left_limit-$right_limit; from 1 to $size; size: $size; location: $chrom $left_limit $right_limit $rsat_strand";

    &PrintSequence($sequence, $fasta_header);

  ## Feature file
  } elsif ($ft_file) {
    open FEAT, $ft_file;
    my $ft_name;
    my $fasta_header;
    while ($line = <FEAT>) {
      chomp($line);
      next if (($line =~/^[#|;]/)||($line eq ""));
      if ($ft_file_format eq "ft") {
	($chrom, $ft_type, $ft_id, $strand, $left_limit, $right_limit,@other_comments) = split (/\t/,$line);
      } elsif ($ft_file_format eq 'gft') {
	($ft_id, $ft_type, $ft_name, $chrom, $left_limit, $right_limit, $strand, @other_comments) = split (/\t/,$line);
      }

      ## Extract only chromosome number if necessary
      $chrom =~ s/chromosome:[\w\.]*?://;
      $chrom =~ s/:.*//;

      ## Tranforms strand in ensembl format
      $strand =~ s/F/1/;
      $strand =~ s/R/-1/;
      $strand =~ s/D/1/;
      $strand =~ s/W/1/;
      $strand =~ s/C/-1/;
      $strand =~ s/>/1/;
      $strand =~ s/</-1/;

      local $chromosome = $slice_adaptor->fetch_by_region('chromosome', $chrom);
      ## Get sequence (repeat masked or not)
      $sequence = &GetSequence($left_limit, $right_limit);
      my $size = $right_limit - $left_limit + 1;

      my $rsat_strand;
      if ($strand == 1) {
	$rsat_strand = "D";
      } elsif ($strand == -1) {
	$rsat_strand = "R";
      }

      if ($header eq 'none') {
	  $header_org = '';
      } else {
	  $header_org = ucfirst($org)."-";
      }

      if ($ft_id) {
	  $fasta_header = ">$header_org$ft_id\t$ft_id; from 1 to $size; size: $size; location: $chrom $left_limit $right_limit $rsat_strand";
      } else {
	  $fasta_header = ">$header_org$chrom-$left_limit-$right_limit\t$chrom-$left_limit-$right_limit; from 1 to $size; size: $size; location: $chrom $left_limit $right_limit $rsat_strand";
      }

      &PrintSequence($sequence, $fasta_header);
    }

  ## All genes
  } elsif ($all) {
      ## from one chromosome
      if ($chrom) {
	  my $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom);
	  my @genes = @{$slice->get_all_Genes()};
	  foreach my $gene (@genes) {
	      &Main($gene, $org);
	  }
      ## From all chromosomes
      } else {
	  my @slices = @{$slice_adaptor->fetch_all("chromosome")};
	  foreach my $slice (@slices) {
	      my @genes = @{$slice->get_all_Genes()};
	      foreach my $gene (@genes) {
		  &Main($gene, $org);
	      }
	  }
      }

  ## List of queries
  } else {
#    my $gene_id = "CG40293"; #gene on D strand
#    my $gene_id = "CG18001"; #gene on R strand
#    @genes = @{$gene_adaptor->fetch_all_by_external_name('BRCA2')};

    my $gene_adaptor = $db->get_GeneAdaptor();

    ## Read queries from input file
    if ($query_file) {
      &RSAT::message::Info("Reading queries from query file", $query_file) if ($main::verbose >= 1);
      my ($in, $input_dir) = &OpenInputFile($query_file);
      while (<$in>) {
	chomp();
	s/\r//g;
	next if ((/^;/) || (/^\#/) || (/^--/)); ## Skip comment and header lines
	if (/(\S+)/) {
	  push @queries, $1;
	}
      }
      close $in;
    }

###      open IN, $query_file;
#       while ($line = <$in>) {
# 	my @genes;
# 	$line =~s/\t.*//;
# 	chomp($line);
# 	if (($line =~ /ENST\d/) || ($line =~ /ENS...T/)) {
# 	  push(@genes, $gene_adaptor->fetch_by_transcript_stable_id($line));
# 	} elsif (($line =~ /ENSP\d/) || ($line =~ /ENS...P/)) {
# 	  push(@genes, $gene_adaptor->fetch_by_translation_stable_id($line));
# 	} elsif (($line =~ /ENSG\d/) || ($line =~ /ENS...G/)) {
# 	  #		my $gene_id = $line;
# 	  push(@genes,$gene_adaptor->fetch_by_stable_id($line));
# 	} else {
# 	  if ($gene_adaptor->fetch_by_stable_id($line)) {
# 	    push(@genes,$gene_adaptor->fetch_by_stable_id($line));
# 	  } else {
# 	    @genes = @{$gene_adaptor->fetch_all_by_external_name($line)};
# 	  }
# 	}
# 	if (@genes) {
# 	  foreach my $gene (@genes) {
# 	    if ($gene) {
# 	      ## get-orthologs if wanted
# 	      if ($ortho) {
# 		&Ortho($gene->stable_id);
# 		## or not
# 	      } else {
# 		&Main($gene, $org);
# 	      }
# 	    } else {
# 	      &RSAT::message::Warning (join("\t", "No sequence for query", $line, "Check validity of your query"));
# 	    }
# 	  }
# 	} else {
# 	  &RSAT::message::Warning (join("\t", "No sequence for query", $line, "Check validity of your query"));
# 	}
#       }
#      close $in;

    ## Treat the list of queries
    my $q=0;
    my $query_nb = scalar(@queries);
    foreach my $id (@queries) {
      $q++;
      &RSAT::message::Info("Treating query", $q."/".$query_nb, $id) if ($main::verbose >= 2);
      my @genes = ();
      if (($id =~ /ENST\d/) || ($id =~ /ENS...T/)) {
	push (@genes, $gene_adaptor->fetch_by_transcript_stable_id($id));
      } elsif (($id =~ /ENSP\d/) || ($id =~ /ENS...P/)) {
	push (@genes, $gene_adaptor->fetch_by_translation_stable_id($id));
      } elsif (($id =~ /ENSG\d/) || ($id =~ /ENS...G/)) {
	push (@genes, $gene_adaptor->fetch_by_stable_id($id));
      } else {
	if ($gene_adaptor->fetch_by_stable_id($id)) {
	  push(@genes,$gene_adaptor->fetch_by_stable_id($id));
#	} elsif ($gene_adaptor->fetch_by_dbID($id)) {
#	  push(@genes,$gene_adaptor->fetch_by_dbID($id));
	} else {
	  @genes = @{$gene_adaptor->fetch_all_by_external_name($id)};
	}
      }

      my $gene_nb = scalar(@genes);
      my $g = 0;
      if ($gene_nb > 0) {
	foreach my $gene (@genes) {
	  $g++;
	  &RSAT::message::Info("Treating gene", $g."/".$gene_nb, $gene) if ($main::verbose >= 3);
	  if ($gene) {
	    ## get-orthologs if wanted
	    if ($ortho) {
	      &Ortho($gene->stable_id);
	      ## or not
	    } else {
	      &Main($gene, $org, $id);
	    }
	  } else {
	    &RSAT::message::Warning (join("\t", "No sequence for query", $id, "Check validity of your query"));
	  }
	}
      } else {
	&RSAT::message::Warning (join("\t", "No sequence for query", $id, "Check validity of your query"));
      }
    }
}



  ################################################################
  ###### Close output stream
  my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
  print $fh $exec_time if ($main::verbose >= 1); ## only report exec time if verbosity is specified
  close $fh if ($output_file);

  close $table_handle if ($homologs_table);

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
      if (&IsNatural($ARGV[$a+1])) {
	$verbose = $ARGV[$a+1];
      } else {
	$verbose = 1;
      }

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
      &RSAT::error::FatalError($feattype, "is not a valid feature type. Supported: ", $supported_feattypes) 
	unless ($supported_feattype{$feattype});

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
    } elsif ($ARGV[$a] eq "-utr") {
      $utr = $ARGV[$a+1];

      ### Non coding
    } elsif ($ARGV[$a] eq "-noncoding") {
      $non_coding = 1;

      ### Line width
    } elsif ($ARGV[$a] eq "-lw") {
      $lw = $ARGV[$a+1];

      ### Unique sequences
    } elsif ($ARGV[$a] eq "-uniqseqs") {
      $uniq_seqs = 1;

      ### Orthologs
  } elsif ($ARGV[$a] eq "-ortho") {
      $ortho  = 1;

      ### Homology type
  } elsif ($ARGV[$a] eq "-ortho_type") {
      $ortho_type  = $ARGV[$a+1];

      ### Taxonomic level
  } elsif ($ARGV[$a] eq "-taxon") {
      $taxon  = $ARGV[$a+1];

      ### Homologs table
  } elsif ($ARGV[$a] eq "-homologs_table") {
      $homologs_table  = $ARGV[$a+1];

      ### Header organism
  }  elsif ($ARGV[$a] eq "-header_org") {
      $header  = $ARGV[$a+1];

      ### Header label
  }  elsif ($ARGV[$a] eq "-label") {
      $label  = $ARGV[$a+1];
      &RSAT::error::FatalError($label, "is not a valid value for the option -label. Supported: ".$supported_labels) 
	  unless $supported_label{$label};
  }
}
}
################################################################
#### verbose message
sub Verbose {
#  print "; retrieve-ensembl-seq.pl ";
#  &PrintArguments();
}
################################################################
### Get sequence(s) relative to a feature
sub Main {
  my ($gene, $main_org, $query) = @_;
  my $gene_id = $gene -> stable_id();

  $db = Bio::EnsEMBL::Registry->get_DBAdaptor($main_org, "core");
  $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($main_org, 'Core', 'Slice');

#  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($org, "core", "Gene");
#  my $gene = $gene_adaptor->fetch_by_stable_id($gene_id);

  my $gene_name = $gene->external_name();
  &RSAT::message::TimeWarn("Getting sequence", $type, $gene_name, $feattype) if ($main::verbose >= 2);
  unless ($gene_name) {
    $gene_name = "";
  }

  unless ($common_name) {
      $common_name = $main_org;
  }

  if ($header eq 'common') {
      $header_org = $common_name."-";
  } elsif ($header eq 'scientific') {
      $header_org = ucfirst($main_org)."-";
  } elsif ($header eq 'none') {
      $header_org = '';
  }

  local $chromosome = $gene->slice;
  local $coord_sys  = $chromosome->coord_system()->name();
  my $chromosome_name = $chromosome->name();

  my $gene_start = $gene->start();
  my $gene_end = $gene->end();
  $strand = $gene->strand();

  my $description = $gene->description();
  unless ($description) {
      $description = "no description";
  }

  my $rsat_strand;
  if ($strand == 1) {
    $rsat_strand = "D";
  } else {
    $rsat_strand = "R";
  }

  if ($main::verbose >= 3) {
    &RSAT::message::Info ("Gene:");
    &RSAT::message::Info (join("\t", "# ID", "Name", "Contig", "Start", "End", "Strand", "Description"));
    &RSAT::message::Info (join("\t", $gene_id, $gene_name, $chromosome_name, $gene_start, $gene_end, $rsat_strand, $description));
  }

  if ($feattype eq "gene") {
    my $sequence;
    my $fasta_header;
    my $size;

    # Output complete gene sequence
    if ($type eq "feature") {
      $sequence = &GetSequence($gene_start, $gene_end);
      $size = $gene_end - $gene_start + 1;
      if ($label eq 'query') {
	$fasta_header = ">$query $header_org$gene_id-$gene_name\t$gene_id; size: $size; location: $chromosome_name $gene_start $gene_end $rsat_strand";
      } else {
	$fasta_header = ">$header_org$gene_id-$gene_name\t$gene_id; size: $size; location: $chromosome_name $gene_start $gene_end $rsat_strand";
      }

    } else { # Output upstream or downstream sequence
      my ($left, $right) = &GetLimits($gene_id, $gene_start, $gene_end);
      # Get sequence (repeat masked or not)
      $sequence = &GetSequence($left, $right);
      $size = $new_to - $new_from + 1;
      if ($label eq 'query') {
	$fasta_header = ">$query $header_org$gene_id-$gene_name\t$gene_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
      } else {
	$fasta_header = ">$header_org$gene_id-$gene_name\t$gene_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
      }
    }

    &PrintSequence ($sequence, $fasta_header);

  } else { # feattype = mrna, cds, introns...
    # Get transcripts
    @transcripts = @{$gene->get_all_Transcripts};

    # Initialize transcripts limits
    my $three_primest_id = $transcripts[0]->display_id();
    my $five_primest_id = $transcripts[0]->display_id();
    my $three_primest_start;
    my $five_primest_end;
    if ($strand == 1) {
      $three_primest_start = $transcripts[0]->start();
      $five_primest_end = $transcripts[0]->end();
    } else {
      $three_primest_start = $transcripts[0]->end();
      $five_primest_end = $transcripts[0]->start();
    }

    %seq_limits = ();

    foreach my $transcript(@transcripts) {
      my $transcript_id = $transcript->display_id();
      my $transcript_start = $transcript->start();
      my $transcript_end = $transcript->end();

      if ($main::verbose >= 3) {
	&RSAT::message::Info ("Transcript:");
	&RSAT::message::Info (join("\t", "# ID", "Start", "End"));
	&RSAT::message::Info (join("\t", $transcript_id, $transcript_start, $transcript_end));
      }

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

      if ($feattype eq 'transcript' && $all_transcripts) {
	# Output complete transcript sequence
	if ($type eq "feature") {
	  $seq_limits{$transcript_id} = [$transcript_start, $transcript_end];

	  unless ($uniq_seqs) {
	    $sequence = &GetSequence($transcript_start, $transcript_end);
	    $size = $transcript_end - $transcript_start + 1;
	    if ($label eq 'query') {
	      $fasta_header = ">$query $header_org$gene_id-$gene_name-$transcript_id\t$gene_id-$transcript_id; size: $size; location: $chromosome_name $transcript_start $transcript_end $rsat_strand";
	    } else {
	      $fasta_header = ">$header_org$gene_id-$gene_name-$transcript_id\t$gene_id-$transcript_id; size: $size; location: $chromosome_name $transcript_start $transcript_end $rsat_strand";
	    }
	    &PrintSequence ($sequence, $fasta_header);
	  }

	} else { # Output upstream or downstream sequence
	  my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $transcript_start, $transcript_end);

	  $seq_limits{$transcript_id} = [$left, $right];

	  unless ($uniq_seqs) {
	    $sequence = &GetSequence($left, $right);
	    my $size = $new_to - $new_from + 1;
	    my $fasta_header;
	    if ($label eq 'query') {
	      $fasta_header = ">$query $header_org$gene_id-$gene_name-$transcript_id\t$gene_id-$transcript_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
	    } else {
	      $fasta_header = ">$header_org$gene_id-$gene_name-$transcript_id\t$gene_id-$transcript_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
	    }
	    &PrintSequence ($sequence, $fasta_header);
	  }
	}
      }

      if ($feattype eq 'mrna' && $all_transcripts) {

	# Output complete mRNA sequence (introns removed)
	if ($type eq "feature") {
	  my $mrna_sequence = '';
	  my $mrna_size = 0;

	  foreach my $exon (@{$transcript->get_all_Exons()}) {
	    my $exon_id = $exon->stable_id();
	    my $exon_start = $exon->start();
	    my $exon_end = $exon->end();

	    if ($main::verbose >= 3) {
	      &RSAT::message::Info ("Exon");
	      &RSAT::message::Info (join("\t", "# ID", "Start", "End"));
	      &RSAT::message::Info (join("\t", $exon_id, $exon_start, $exon_end));
	    }

	    $seq_limits{$exon_id} = [$exon_start,$exon_end];
	    $sequence = &GetSequence($exon_start, $exon_end);
	    my $size = $exon_end - $exon_start + 1;
	    $mrna_size = $mrna_size + $size;
	    if ($strand == 1) {
	      $mrna_sequence = $mrna_sequence.$sequence;
	    } else {
	      $mrna_sequence = $sequence.$mrna_sequence;
	    }
	  }

	  unless ($uniq_seqs) {
	    if ($label eq 'query') {
	      $fasta_header = ">$query $header_org$gene_id-$gene_name-$transcript_id-mRNA\t$gene_id-$transcript_id-mRNA; size: $mrna_size; location: $chromosome_name $transcript_start $transcript_end $rsat_strand";
	    } else {
	      $fasta_header = ">$header_org$gene_id-$gene_name-$transcript_id-mRNA\t$gene_id-$transcript_id-mRNA; size: $mrna_size; location: $chromosome_name $transcript_start $transcript_end $rsat_strand";
	    }
	    &PrintSequence ($mrna_sequence, $fasta_header);
	  }

	} else { # Output upqtream or downstream sequence
	  my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $transcript_start, $transcript_end);

	  $seq_limits{$transcript_id} = [$left, $right];

	  unless ($uniq_seqs) {
	    # Output upstream or downstream sequence
	    $sequence = &GetSequence($left, $right);
	    my $size = $new_to - $new_from + 1;
	    my $fasta_header;
	    if ($label eq 'query') {
	      $fasta_header = ">$query $header_org$gene_id-$gene_name-$transcript_id\t$gene_id-$transcript_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
	    } else {
	      $fasta_header = ">$header_org$gene_id-$gene_name-$transcript_id\t$gene_id-$transcript_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
	    }
	    &PrintSequence ($sequence, $fasta_header);
	  }
	}
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

	foreach my $intron (@{$transcript->get_all_Introns()}) {
	  $i++;
	  my $intron_start = $intron->start();
	  my $intron_end = $intron->end();
	  my $intron_id = $transcript_id."-".$i;

	  if ($main::verbose >= 3) {
	    &RSAT::message::Info ("Intron:");
	    &RSAT::message::Info (join("\t", "# ID", "Start", "End"));
	    &RSAT::message::Info (join("\t", $intron_id, $intron_start, $intron_end));
	  }

	  if ($first_intron) {
	    if ($strand == 1) {
	      if ($intron->start() < $start1) {
		$intron1 = $intron;
		$intron1_id = $transcript_id."-".$i;
		$start1 = $intron->start();
		$end1 = $intron->end();
	      }
	    } else {
	      if ($intron->end() > $end1) {
		$intron1 = $intron;
		$intron1_id = $transcript_id."-".$i;
		$start1 = $intron->start();
		$end1 = $intron->end();
	      }
	    }
	  }

	  unless ($first_intron) {

	    $seq_limits{$intron_id} = [$intron_start, $intron_end];

	    unless ($uniq_seqs) {
	      $sequence = &GetSequence($intron->start(), $intron->end());
	      my $size = ($intron->end() - $intron->start()) + 1;

	      my $fasta_header;
	      if ($label eq 'query') {
		$fasta_header = ">$query $header_org$gene_id-$gene_name-$transcript_id-$i\t$gene_id-$transcript_id-$i; from 1 to $size; size: $size; location: $chromosome_name $intron_start $intron_end $rsat_strand";
	      } else {
		$fasta_header = ">$header_org$gene_id-$gene_name-$transcript_id-$i\t$gene_id-$transcript_id-$i; from 1 to $size; size: $size; location: $chromosome_name $intron_start $intron_end $rsat_strand";
	      }

	      &PrintSequence ($sequence, $fasta_header);
	    }
	  }
	}

	if ($first_intron && $i != 0) {
	  if ($main::verbose >= 3) {
	    &RSAT::message::Info ("First intron:");
	    &RSAT::message::Info (join("\t", "# ID", "Start", "End"));
	    &RSAT::message::Info (join("\t", $intron1_id, $start1, $end1));
	  }

	  $seq_limits{$intron1_id} = [$start1, $end1];

	  unless ($uniq_seqs) {
	    $sequence = &GetSequence($start1, $end1);
	    my $size = ($end1 - $start1) + 1;

	    my $fasta_header;
	    if ($label eq 'query') {
	      $fasta_header = ">$query $header_org$gene_id-$gene_name-$intron1_id\t$gene_id-$intron1_id; from 1 to $size; size: $size; location: $chromosome_name $start1 $end1 $rsat_strand";
	    } else {
	      $fasta_header = ">$header_org$gene_id-$gene_name-$intron1_id\t$gene_id-$intron1_id; from 1 to $size; size: $size; location: $chromosome_name $start1 $end1 $rsat_strand";
	    }
	    &PrintSequence ($sequence, $fasta_header);
	  }
	}

	if ($i == 0) {
	  &RSAT::message::Warning ("Gene $gene_id - transcript $transcript_id has no intron");
	}
      }

      my $coding_region_start = $transcript->coding_region_start();
      my $coding_region_end = $transcript->coding_region_end();

      if (($coding_region_start) && ($coding_region_end)) {
	if ($main::verbose >= 3) {
	  &RSAT::message::Info ("Coding region start: $coding_region_start") ;
	  &RSAT::message::Info ("Coding region end: $coding_region_end");
	}
      } else {
	&RSAT::message::Warning ("Gene $gene_id - transcript $transcript_id has no coding region") unless ($feattype eq "mrna");
      }

      # Exons
      if ($feattype eq "exon") {
	foreach my $exon (@{$transcript->get_all_Exons()}) {
	  my $exon_id = $exon->stable_id();
	  my $exon_start = $exon->start();
	  my $exon_end = $exon->end();

	  if ($main::verbose >= 3) {
	    &RSAT::message::Info ("Exon");
	    &RSAT::message::Info (join("\t", "# ID", "Start", "End"));
	    &RSAT::message::Info (join("\t", $exon_id, $exon_start, $exon_end));
	  }

	  if (($non_coding) && ($coding_region_start) && ($coding_region_end) ) {
	    if ($coding_region_start > $exon_start && $coding_region_start < $exon_end) {

	      $seq_limits{$exon_id} = [$exon_start, $coding_region_start - 1];

	      unless ($uniq_seqs) {
		$sequence = &GetSequence($exon_start, $coding_region_start - 1);
		my $non_coding_exon_right = $coding_region_start - 1;
		my $size = $coding_region_start - $exon_start;

		my $fasta_header;
		if ($label eq 'query') {
		  $fasta_header = ">$query $header_org$gene_id-$gene_name-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $exon_start $non_coding_exon_right $rsat_strand";
		} else {
		  $fasta_header = ">$header_org$gene_id-$gene_name-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $exon_start $non_coding_exon_right $rsat_strand";
		}
		&PrintSequence ($sequence, $fasta_header);
	      }
	    } elsif ($coding_region_end < $exon_end && $coding_region_end > $exon_start) {

	      $seq_limits{$exon_id} = [$coding_region_end + 1, $exon_end];

	      unless ($uniq_seqs) {
		$sequence = &GetSequence($coding_region_end + 1, $exon_end);
		my $non_coding_exon_left = $coding_region_end + 1;
		my $size = $exon->end() - $coding_region_end;

		my $fasta_header;
		if ($label eq 'query') {
		  $fasta_header = ">$query $header_org$gene_id-$gene_name-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $non_coding_exon_left $exon_end $rsat_strand";
		} else {
		  $fasta_header = ">$header_org$gene_id-$gene_name-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $non_coding_exon_left $exon_end $rsat_strand";
		}
		&PrintSequence ($sequence, $fasta_header);
	      }

	    } elsif ($coding_region_start > $exon_end) {

	      $seq_limits{$exon_id} = [$exon_start, $exon_end];

	      unless ($uniq_seqs) {
		$sequence = &GetSequence($exon_start, $exon_end);
		my $size = ($exon_end - $exon_start) + 1;

		my $fasta_header;
		if ($label eq 'query') {
		  $fasta_header = ">$query $header_org$gene_id-$gene_name-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $exon_start $exon_end $rsat_strand";
		} else {
		  $fasta_header = ">$header_org$gene_id-$gene_name-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $exon_start $exon_end $rsat_strand";
		}
		&PrintSequence ($sequence, $fasta_header);
	      }

	    } elsif ($coding_region_end < $exon_start) {

	      $seq_limits{$exon_id} = [$exon_start, $exon_end];

	      unless ($uniq_seqs) {
		$sequence = &GetSequence($exon_start, $exon_end);
		my $size = ($exon_end - $exon_start) + 1;

		my $fasta_header;
		if ($label eq 'query') {
		  $fasta_header = ">$query $header_org$gene_id-$gene_name-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $exon_start $exon_end $rsat_strand";
		} else {
		  $fasta_header = ">$header_org$gene_id-$gene_name-$exon_id-non_coding\t$gene_id-$exon_id-non_coding; from 1 to $size; size: $size; location: $chromosome_name $exon_start $exon_end $rsat_strand";
		}
		&PrintSequence ($sequence, $fasta_header);
	      }
	    }
	  } else {

	    $seq_limits{$exon_id} = [$exon_start, $exon_end];

	    unless ($uniq_seqs) {
	      $sequence = &GetSequence($exon_start, $exon_end);
	      my $size = ($exon_end - $exon_start) + 1;

	      my $fasta_header;
	      if ($label eq 'query') {
		$fasta_header = ">$query $header_org$gene_id-$gene_name-$exon_id\t$gene_id-$exon_id; from 1 to $size; size: $size; location: $chromosome_name $exon_start $exon_end $rsat_strand";
	      } else {
		$fasta_header = ">$header_org$gene_id-$gene_name-$exon_id\t$gene_id-$exon_id; from 1 to $size; size: $size; location: $chromosome_name $exon_start $exon_end $rsat_strand";
	      }
	      &PrintSequence ($sequence, $fasta_header);
	    }
	  }
	}
      }

      # UTR
      if (($feattype eq "utr") && ($coding_region_start) && ($coding_region_end)) {
	my $utr5_start;
	my $utr5_end;
	my $utr3_start;
	my $utr3_end;
	my $utr5_flag = 0;
	my $utr3_flag = 0;
	if ($strand == 1) {
	  $utr5_start = $transcript_start;
	  $utr5_end = $coding_region_start - 1;
	  unless ($transcript_start == $coding_region_start) {
	    $utr5_flag = 1;
	  }
	  $utr3_start = $coding_region_end + 1;
	  $utr3_end = $transcript_end;
	  unless ($transcript_end == $coding_region_end) {
	    $utr3_flag = 1;
	  }
	} else {
	  unless ($transcript_start == $coding_region_start) {
	    $utr3_flag = 1;
	  }
	  $utr3_start = $transcript_start;
	  $utr3_end = $coding_region_start - 1;
	  unless ($transcript_end == $coding_region_end) {
	    $utr5_flag = 1;
	  }
	  $utr5_start = $coding_region_end + 1;
	  $utr5_end = $transcript_end;
	}

	my $utr5_id = $transcript_id."-5prime_UTR";
	my $utr3_id = $transcript_id."-3prime_UTR";

	if (($utr5_flag == 1) && (($utr eq 'all') || ($utr eq '5prime'))) {
	  $seq_limits{$utr5_id} = [$utr5_start, $utr5_end];
	}
	if (($utr3_flag == 1) && (($utr eq 'all') || ($utr eq '3prime'))) {
	  $seq_limits{$utr3_id} = [$utr3_start, $utr3_end];
	}

	unless ($uniq_seqs) {
	  if (($utr5_flag == 1) && (($utr eq 'all') || ($utr eq '5prime'))) {
	    $utr5_sequence = &GetSequence($utr5_start, $utr5_end);
	    my $utr5_size = $utr5_end - $utr5_start + 1;

	    my $utr5_fasta_header;
	    if ($label eq 'query') {
	      $utr5_fasta_header = ">$query $header_org$gene_id-$gene_name-$transcript_id-5prime_UTR\t$gene_id-$transcript_id-5prime_UTR; from 1 to $utr5_size; size: $utr5_size; location: $chromosome_name $utr5_start $utr5_end $rsat_strand";
	    } else {
	      $utr5_fasta_header = ">$header_org$gene_id-$gene_name-$transcript_id-5prime_UTR\t$gene_id-$transcript_id-5prime_UTR; from 1 to $utr5_size; size: $utr5_size; location: $chromosome_name $utr5_start $utr5_end $rsat_strand";
	    }
	    &PrintSequence ($utr5_sequence, $utr5_fasta_header);
	  }
	  if (($utr3_flag == 1) && (($utr eq 'all') || ($utr eq '3prime'))) {
	    $utr3_sequence = &GetSequence($utr3_start, $utr3_end);
	    my $utr3_size = $utr3_end - $utr3_start + 1;

	    my $utr3_fasta_header;
	    if ($label eq 'query') {
	      $utr3_fasta_header = ">$query $header_org$gene_id-$gene_name-$transcript_id-3prime_UTR\t$gene_id-$transcript_id-3prime_UTR; from 1 to $utr3_size; size: $utr3_size; location: $chromosome_name $utr3_start $utr3_end $rsat_strand";
	    } else {
	      $utr3_fasta_header = ">$header_org$gene_id-$gene_name-$transcript_id-3prime_UTR\t$gene_id-$transcript_id-3prime_UTR; from 1 to $utr3_size; size: $utr3_size; location: $chromosome_name $utr3_start $utr3_end $rsat_strand";
	    }
	    &PrintSequence ($utr3_sequence, $utr3_fasta_header);
	  }
	}
      }

      # CDS
      if (($feattype eq 'cds') && ($coding_region_start) && ($coding_region_end)) {
	my $cds_id = $transcript->translation()->stable_id();

	# Output complete CDS sequence
	if ($type eq "feature") {
	  my $cds_sequence = '';
	  my $cds_size = 0;
	  foreach my $exon (@{$transcript->get_all_Exons()}) {
	    my $exon_id = $exon->stable_id();
	    my $exon_start = $exon->start();
	    my $exon_end = $exon->end();

	    if ($main::verbose >= 3) {
	      &RSAT::message::Info ("Exon");
	      &RSAT::message::Info (join("\t", "# ID", "Start", "End"));
	      &RSAT::message::Info (join("\t", $exon_id, $exon_start, $exon_end));
	    }

	    if ($coding_region_start > $exon_start && $coding_region_start <= $exon_end) {
	      if ($coding_region_end > $exon_end) {
		$seq_limits{$exon_id} = [$coding_region_start,$exon_end];
		$sequence = &GetSequence($coding_region_start,$exon_end);
		my $size = $exon_end - $coding_region_start + 1;
		$cds_size = $cds_size + $size;
	      } else {
		$seq_limits{$exon_id} = [$coding_region_start,$coding_region_end];
		$sequence = &GetSequence($coding_region_start,$coding_region_end);
		my $size = $coding_region_end - $coding_region_start + 1;
		$cds_size = $cds_size + $size;
	      }
	      if ($strand == 1) {
		$cds_sequence = $cds_sequence.$sequence;
	      } else {
		$cds_sequence = $sequence.$cds_sequence;
	      }

	    } elsif ($coding_region_end < $exon_end && $coding_region_end >= $exon_start) {
	      if ($coding_region_start < $exon_start) {
		$seq_limits{$exon_id} = [$exon_start,$coding_region_end];
		$sequence = &GetSequence($exon_start, $coding_region_end);
		my $size = $coding_region_end - $exon_start + 1;
		$cds_size = $cds_size + $size;
#		} else {
#		$seq_limits{$exon_id} = [$coding_region_start,$coding_region_end];
#		$sequence = &GetSequence($coding_region_start, $coding_region_end);
	      }
	      if ($strand == 1) {
		$cds_sequence = $cds_sequence.$sequence;
	      } else {
		$cds_sequence = $sequence.$cds_sequence;
	      }

	    } elsif ($coding_region_start <= $exon_start && $coding_region_end >= $exon_end) {
	      $seq_limits{$exon_id} = [$exon_start,$exon_end];
	      $sequence = &GetSequence($exon_start, $exon_end);
	      my $size = $exon_end - $exon_start + 1;
	      $cds_size = $cds_size + $size;
	      if ($strand == 1) {
		$cds_sequence = $cds_sequence.$sequence;
	      } else {
		$cds_sequence = $sequence.$cds_sequence;
	      }
	    }
	  }

	  unless ($uniq_seqs) {
	    if ($label eq 'query') {
	      $fasta_header = ">$query $header_org$gene_id-$gene_name-$transcript_id-$cds_id\t$gene_id-$transcript_id-$cds_id; size: $cds_size; location: $chromosome_name $coding_region_start $coding_region_end $rsat_strand";
	    } else {
	      $fasta_header = ">$header_org$gene_id-$gene_name-$transcript_id-$cds_id\t$gene_id-$transcript_id-$cds_id; size: $cds_size; location: $chromosome_name $coding_region_start $coding_region_end $rsat_strand";
	    }
	    &PrintSequence ($cds_sequence, $fasta_header);
	  }

	} else { # Output upstream or downstream sequence
	  my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $coding_region_start, $coding_region_end);
	  $seq_limits{$cds_id} = [$left, $right];

	  unless ($uniq_seqs) {
	    # Output upstream or downstream sequence
	    $sequence = &GetSequence($left, $right);
	    my $size = $new_to - $new_from + 1;

	    my $fasta_header;
	    if ($label eq 'query') {
	      $fasta_header = ">$query $header_org$gene_id-$gene_name-$transcript_id-$cds_id\t$gene_id-$transcript_id-$cds_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
	    } else {
	      $fasta_header = ">$header_org$gene_id-$gene_name-$transcript_id-$cds_id\t$gene_id-$transcript_id-$cds_id; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
	    }
	    &PrintSequence ($sequence, $fasta_header);
	  }
	}
      }

      if ($feattype eq 'transcript' && !$all_transcripts) {
	# Output complete transcript sequence
	if ($type eq 'feature') {
	  &RSAT::message::Warning ("To get the sequence covering all transcripts of a gene, use -feattype gene");

	} else { # Output upstream or downstream sequence
	  if ($main::verbose >= 3) {
	    &RSAT::message::Info ("Three_primest transcript start: $three_primest_start");
	    &RSAT::message::Info ("Five_primest transcript end: $five_primest_end");
	  }

	  my $ref_transcript;

	  if ($strand == 1) {
	    if ($five_primest_end <= $three_primest_start) {
	      &RSAT::message::Warning ("Gene $gene_id has disjoint alternative transcripts; Retrieve sequences relative to each alternative transcript (-alltranscripts)");
	    } else {
	      my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $three_primest_start, $five_primest_end);
	    }
	  } else {
	    if ($three_primest_start <= $five_primest_end) {
	      &RSAT::message::Warning ("Gene $gene_id has disjoint alternative transcripts; Retrieve sequences relative to each alternative transcript (-alltranscripts)");
	    } else {
	      my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $five_primest_end, $three_primest_start);
	    }
	  }

	  unless (($strand == 1 && $five_primest_end <= $three_primest_start) || ($strand == -1 && $three_primest_start <= $five_primest_end)) {
	    # Output upstream or downstream sequence
	    $sequence = &GetSequence($left, $right);
	    my $size = $new_to - $new_from + 1;
	    if ($type eq "upstream") {
	      $ref_transcript = $three_primest_id;
	    } else {
	      $ref_transcript = $five_primest_id;
	    }

	    my $fasta_header;
	    if ($label eq 'query') {
	      $fasta_header = ">$query $header_org$gene_id-$gene_name-$ref_transcript\t$gene_id-$ref_transcript; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
	    } else {
	      $fasta_header = ">$header_org$gene_id-$gene_name-$ref_transcript\t$gene_id-$ref_transcript; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
	    }
	    &PrintSequence ($sequence, $fasta_header);
	  }
	}
      }

      if ($feattype eq 'mrna' && !$all_transcripts) {
	# Output complete mRNA sequence
	if ($type eq 'feature') {
	  &RSAT::message::Warning ("The option to get the sequence covering all mRNAs of a gene is not implemented yet");

	} else { # Output upstream or downstream sequence (same as for transcript)
	  if ($main::verbose >= 3) {
	    &RSAT::message::Info ("Three_primest transcript start: $three_primest_start");
	    &RSAT::message::Info ("Five_primest transcript end: $five_primest_end");
	  }

	  my $ref_transcript;

	  if ($strand == 1) {
	    if ($five_primest_end <= $three_primest_start) {
	      &RSAT::message::Warning ("Gene $gene_id has disjoint alternative transcripts; Retrieve sequences relative to each alternative transcript (-alltranscripts)");
	    } else {
	      my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $three_primest_start, $five_primest_end);
	    }
	  } else {
	    if ($three_primest_start <= $five_primest_end) {
	      &RSAT::message::Warning ("Gene $gene_id has disjoint alternative transcripts; Retrieve sequences relative to each alternative transcript (-alltranscripts)");
	    } else {
	      my ($left, $right, $new_from, $new_to) = &GetLimits($gene_id, $five_primest_end, $three_primest_start);
	    }
	  }

	  unless (($strand == 1 && $five_primest_end <= $three_primest_start) || ($strand == -1 && $three_primest_start <= $five_primest_end)) {
	    # Output upstream or downstream sequence
	    $sequence = &GetSequence($left, $right);
	    my $size = $new_to - $new_from + 1;
	    if ($type eq "upstream") {
	      $ref_transcript = $three_primest_id;
	    } else {
	      $ref_transcript = $five_primest_id;
	    }

	    my $fasta_header;
	    if ($label eq 'query') {
	      $fasta_header = ">$query $header_org$gene_id-$gene_name-$ref_transcript\t$gene_id-$ref_transcript; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
	    } else {
	      $fasta_header = ">$header_org$gene_id-$gene_name-$ref_transcript\t$gene_id-$ref_transcript; $type from $new_from to $new_to; size: $size; location: $chromosome_name $left $right $rsat_strand";
	    }
	    &PrintSequence ($sequence, $fasta_header);
	  }
	}
      }
    }

    ## Unique sequences
    if ($uniq_seqs && %seq_limits) {
      # Sort retrieved sequences by left limits
      my @ordered_seqs;
      foreach $value (sort {$seq_limits{$a}[0] <=> $seq_limits{$b}[0] } keys %seq_limits) {
	push @ordered_seqs, $value;
      }

      # Initialisation
      $new_seq_index = 1;
      %new_seq_limits = ();
      $new_seq_limits{$new_seq_index} = [$seq_limits{$ordered_seqs[0]}[0], $seq_limits{$ordered_seqs[0]}[1]];

      # Search unique sequences
      foreach $seq (@ordered_seqs) {
	if ($seq_limits{$seq}[0] <= $new_seq_limits{$new_seq_index}[1] + 1) {
	  if ($seq_limits{$seq}[1] > $new_seq_limits{$new_seq_index}[1]) {
	    $new_seq_limits{$new_seq_index}[1] = $seq_limits{$seq}[1];
	  }
	} else {
	  $new_seq_index++;
	  $new_seq_limits{$new_seq_index} = [$seq_limits{$seq}[0], $seq_limits{$seq}[1]];
	}
      }

      # print Dumper(\%new_seq_limits);

      # Get and print unique sequences
      foreach $value (keys(%new_seq_limits)) {
	$sequence = &GetSequence($new_seq_limits{$value}[0], $new_seq_limits{$value}[1]);
	$size = $new_seq_limits{$value}[1] - $new_seq_limits{$value}[0] + 1;
	unless ($feattype eq 'transcript' || $feattype eq 'mrna' || $feattype eq 'cds') {
	  $type = '';
	}

	my $fasta_header;
	if ($label eq 'query') {
	  $fasta_header = ">$query $header_org$gene_id-$gene_name-$value\t$gene_id-$value; $feattype $type unique sequence; size: $size; location: $chromosome_name $new_seq_limits{$value}[0] $new_seq_limits{$value}[1] $rsat_strand";
	} else {
	  $fasta_header = ">$header_org$gene_id-$gene_name-$value\t$gene_id-$value; $feattype $type unique sequence; size: $size; location: $chromosome_name $new_seq_limits{$value}[0] $new_seq_limits{$value}[1] $rsat_strand";
	}
	&PrintSequence ($sequence, $fasta_header);
      }
#      }
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
      $expanded_slice = $slice_adaptor->fetch_by_region($coord_sys, $chromosome->seq_region_name(), $start + $from, $end);
      $closest_neighbour_limit = &GetNeighbours($gene_id, $start, $end, $expanded_slice);

      &RSAT::message::Info ("Closest neighbour limit: $closest_neighbour_limit") if ($main::verbose >= 3);

    }
#    if (($nogene || $noorf) && $closest_neighbour_limit > $start + $from) {
    if (($nogene || $noorf) && $closest_neighbour_limit > $start + $from && $closest_neighbour_limit < $start) {
      $left = $closest_neighbour_limit + 1;
    } else {
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
      $expanded_slice = $slice_adaptor->fetch_by_region($coord_sys, $chromosome->seq_region_name(), $start, $end - $from);
      $closest_neighbour_limit = &GetNeighbours($gene_id, $start, $end, $expanded_slice);
      &RSAT::message::Info ("Closest neighbour limit: $closest_neighbour_limit") if ($main::verbose >= 3);
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
      $expanded_slice = $slice_adaptor->fetch_by_region($coord_sys, $chromosome->seq_region_name(), $start, $end + $to);
      $closest_neighbour_limit = &GetNeighbours($gene_id, $start, $end, $expanded_slice);

      &RSAT::message::Info ("Closest neighbour limit: $closest_neighbour_limit") if ($main::verbose >= 3);

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
      $expanded_slice = $slice_adaptor->fetch_by_region($coord_sys, $chromosome->seq_region_name(), $start - $to, $end);
      $closest_neighbour_limit = &GetNeighbours($gene_id, $start, $end, $expanded_slice);

      &RSAT::message::Info ("Closest neighbour limit: $closest_neighbour_limit") if ($main::verbose >= 3);

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

  &RSAT::message::Info ("Neighbour gene(s):") if ($main::verbose >= 3);
  &RSAT::message::Info (join("\t", "# ID", "Name", "Contig", "Start", "End", "Strand", "Description")) if ($main::verbose >= 3);

  foreach my $neighbour_gene(@neighbour_genes) {
    if ($neighbour_gene->stable_id() eq $gene_id) {    # the query gene itself is in the expanded slice
      next
    }
    my $neighbour_name = $neighbour_gene->external_name();
    unless ($neighbour_name) {
      $neighbour_name = "";
    }

    my $neighbour_description = $neighbour_gene->description();
    unless ($neighbour_description) {
	$neighbour_description = "";
    }

    my $neighbour_strand = $neighbour_gene->seq_region_strand();
    if ($neighbour_strand == 1) {
      $neighbour_strand = "D";
    } else {
      $neighbour_strand = "R";
    }
    my $neighbour_start = $neighbour_gene->seq_region_start();
    my $neighbour_end = $neighbour_gene->seq_region_end();

    &RSAT::message::Info (join("\t", $neighbour_gene->stable_id(), $neighbour_name, $chromosome->name(), $neighbour_start, $neighbour_end, $neighbour_strand, $neighbour_description)) if ($main::verbose >= 3);

    #Find closest neighbour limit
    if ($nogene) {    # neighbour limits are closest gene limits
      if ($strand == 1 && $type eq "upstream") {
	if ($neighbour_gene->seq_region_end() > $closest_neighbour_limit) {
	 if ($neighbour_gene->seq_region_end() < $start + $to) {
	  $closest_neighbour_limit = $neighbour_gene->seq_region_end();
	 } elsif (($neighbour_gene->seq_region_start() < $start + $to) && ($neighbour_gene->seq_region_end() >= $start + $to)) {
	  $closest_neighbour_limit = $start + $to;
	 }
	}
      } elsif ($strand == 1 && $type eq "downstream") {
	if ($neighbour_gene->seq_region_start() < $closest_neighbour_limit) {
	 if ($neighbour_gene->seq_region_start() > $end - $to) {
	  $closest_neighbour_limit = $neighbour_gene->seq_region_start();
	 } elsif (($neighbour_gene->seq_region_start() <= $end - $to) && ($neighbour_gene->seq_region_end() > $end - $to)) {
	  $closest_neighbour_limit = $end - $to;
	 }
	}
      } elsif ($strand == -1 && $type eq "upstream") {
	if ($neighbour_gene->seq_region_start() < $closest_neighbour_limit) {
	 if ($neighbour_gene->seq_region_start() > $end - $to) {
	  $closest_neighbour_limit = $neighbour_gene->seq_region_start();
	 } elsif (($neighbour_gene->seq_region_start() <= $end - $to) && ($neighbour_gene->seq_region_end() > $end - $to)) {
	  $closest_neighbour_limit = $end - $to;
	 }
	}
      } elsif ($strand == -1 && $type eq "downstream") {
	if ($neighbour_gene->seq_region_end() > $closest_neighbour_limit) {
	 if ($neighbour_gene->seq_region_end() < $start + $to) {
	  $closest_neighbour_limit = $neighbour_gene->seq_region_end();
	 } elsif (($neighbour_gene->seq_region_start() < $start + $to) && ($neighbour_gene->seq_region_end() >= $start + $to)) {
	  $closest_neighbour_limit = $start + $to;
	 }
	}
      }
    } elsif ($noorf) {    # neighbour limits are closest CDS limits
      my @transcripts = @{$neighbour_gene->get_all_Transcripts};
      foreach my $transcript(@transcripts) {
       if ($transcript->coding_region_start() && $transcript->coding_region_end()) {
	my $coding_region_start = $transcript->coding_region_start();
	my $coding_region_end = $transcript->coding_region_end();
	if (($strand == 1 && $type eq "upstream") || ($strand == -1 && $type eq "downstream")) {
	  if ($type eq "upstream") {
	  $coding_region_end = $coding_region_end + $start - abs($from) - 1;
	} else {
	  $coding_region_end = $coding_region_end + $start - abs($to) - 1;
	}
	  if ($coding_region_end > $closest_neighbour_limit) {
	   if ($coding_region_end < $start + $to) {
	    $closest_neighbour_limit = $coding_region_end;
	   } elsif (($coding_region_start < $start + $to) && ($coding_region_end >= $start + $to)) {
	    $closest_neighbour_limit = $start + $to;
	   }
	  }
	} elsif (($strand == 1 && $type eq "downstream") || ($strand == -1 && $type eq "upstream")) {
	  $coding_region_start = $coding_region_start + $start - 1;
	  if ($coding_region_start < $closest_neighbour_limit) {
	   if ($coding_region_start > $end - $to) {
	    $closest_neighbour_limit = $coding_region_start;
	   } elsif (($coding_region_start <= $end - $to) && ($coding_region_end > $end - $to)) {
	    $closest_neighbour_limit = $end - $to;
	   }
	  }
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

    my $mini_slice = $slice_adaptor->fetch_by_region($coord_sys, $chromosome->seq_region_name(), $left, $right);

    eval {
	my $masked_slice = $mini_slice->get_repeatmasked_seq();
      };
    if ($@) {
      &RSAT::message::Warning("Masking of repeats not available for this homolog");
    } else {
      my $masked_slice = $mini_slice->get_repeatmasked_seq();

#    my $masked_slice = $mini_slice->get_repeatmasked_seq([],0,{"repeat_class_Dust" => 0});
#    my $masked_slice = $mini_slice->get_repeatmasked_seq(['TRF']);
      $sequence = $masked_slice->seq();
    }
#    my $unmasked_seq = $slice->seq();
#    my $hardmasked_seq = $slice->get_repeatmasked_seq();
#    my $softmasked_seq = $slice->get_repeatmasked_seq(undef, 1);
  }
  if ($mask_coding) {
    my $retrieved_slice = $slice_adaptor->fetch_by_region($coord_sys, $chromosome->seq_region_name(), $left, $right);
    my $transcript_adaptor = $db->get_TranscriptAdaptor();
    my @retrieved_transcripts = @{$transcript_adaptor->fetch_all_by_Slice($retrieved_slice)};
    foreach $retrieved_transcript (@retrieved_transcripts) {
      my @retrieved_exons = @{$retrieved_transcript->get_all_translateable_Exons};
      foreach my $retrieved_exon (@retrieved_exons) {
	&RSAT::message::Info ("Translateable exon start: ".$retrieved_exon->start()."\tTranslateable exon end: ".$retrieved_exon->end()."\tTranslateable exon strand: ".$retrieved_exon->strand()) if ($main::verbose >= 5);
	&RSAT::message::Info ("Translateable exon start: ".$retrieved_exon->seq_region_start()."\tTranslateable exon end: ".$retrieved_exon->seq_region_end()) if ($main::verbose >= 5);
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
	  if ($coding_length > 0) {
	    $masking = "N" x $coding_length;
	    #	  foreach (1..$coding_length) {
	    #	    $masking .= "N";
	    #	  }
##	    &RSAT::message::Debug($retrieved_exon, "seq", $sequence,  "coding_length", $coding_length, "masking", $masking) if ($main::verbose >= 0);
	    substr($sequence, 0, $coding_length) = $masking;
	  }
	}
      }
    }
  }
  if (($strand == -1) && ($sequence)) {
    $sequence = &ReverseComplement($sequence);
  }
  return $sequence;
}


  ################################################################
  ### Get orthologous sequences
  sub Ortho {

# 	my $genome_db_adaptor = $registry->get_adaptor(
#     'Multi', 'compara', 'GenomeDB');

#     ## featch all
#     my $all_genome_dbs = $genome_db_adaptor->fetch_all();
# 	foreach my $this_genome_db (@{$all_genome_dbs}) {
#   		print $this_genome_db->name, "\n";
# 	}

# 	## fetch by slice example
# 	my $gene_adaptor_M  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
# 	#my $gene_id = "OTTMUSG00000018868" ;
# 	my $gene_id_M = "ENSMUSG00000038253";
# 	my $gene_M = $gene_adaptor_M->fetch_by_stable_id($gene_id_M);
# 	my $chromosome_M = $gene_M->slice;

# 	my $gdb = $genome_db_adaptor->fetch_by_Slice($chromosome_M);
# 	print $gdb->name;

# 	die();

#      if (scalar(@queries) > 1) {
#	  &RSAT::message::Warning("Only your first query will be treated");
#      }

      my $ortho_id = shift;

      my $compara_dbname = 'Multi';
#     my $compara_dbname = 'compara'; ## works also...

#     my $ma = Bio::EnsEMBL::Registry->get_adaptor($compara_dbname,'compara','Member');
      my $ma = Bio::EnsEMBL::Registry->get_adaptor($compara_dbname,'compara','GeneMember');

#     Sample Ids to test
#     $ortho_id = 'ENSMUSG00000038253';
#     $ortho_id = 'ENSG00000004059';

      my $gene_member = $ma->fetch_by_source_stable_id('ENSEMBLGENE', $ortho_id);

      # print out some information about the Member
      &RSAT::message::Info("# Chrom_name Chrom_start Chrom_end Description Source_name Taxon_id") if ($main::verbose >= 1);
      &RSAT::message::Info(join (" ", map { $gene_member->$_ } qw(chr_name chr_start chr_end description source_name taxon_id))) if ($main::verbose >= 1);

      my $compara_taxon = $gene_member->taxon;
      &RSAT::message::Info("# Common_name; Genus; Species; Organism; Classification") if ($main::verbose >= 1);
      &RSAT::message::Info(join ("; ", map { $compara_taxon->$_ } qw(common_name genus species binomial classification))) if ($main::verbose >= 1);

      $common_name = $gene_member->taxon->common_name;
      if ($common_name) {
	$common_name =~ s/\s+/_/g;
      }

      my $gene = $gene_member->get_Gene;
      &Main($gene, $org);
#      &Main($ortho_id, $org);

#      my @taxons = split (/ /, $compara_taxon->classification);
#      my @limited_taxons;

#      if ($taxon) {
#	  foreach my $tax (@taxons) {
#	      push(@limited_taxons,$tax);
#	      last if (lc($tax) eq lc($taxon));
#	  }
#      }

      # then you get the homologies where the member is involved
      my $ha = Bio::EnsEMBL::Registry->get_adaptor($compara_dbname,'compara','Homology');
      my $homologies = $ha->fetch_all_by_Member($gene_member);
      # That will return an array reference with all homologies (orthologues, and in some cases paralogues) against other species.
      # Then for each homology, you get all the Members implicated

      if ($homologs_table) {
	  print $table_handle "# Gene_id\tOrganism\tGene_description\tHomology_type\tTaxon_level\tPerc_id\tPerc_pos\tPerc_cov\tQuery_gene\tQuery_organism\n";
      }

      my $taxon_filter_flag = 0;

      ## Get classifications for all EnsEMBL organisms at once
      my %classifications;
      if ($taxon) {
	  ## Get species names - 'new' way
#	  my @dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors();
#	  my @species_;
#	  foreach my $dba (@{$dbas[0]}) {
#	      my @all_species = @{$dba->all_species()};
#	      if ($all_species[0] =~ '_') { # there are not only organism names (multi, Ancestor species)
#		  chomp $all_species[0];
#		  push @species_, $all_species[0];
#	      }
#	  }
	  # Sort because some species come several time
#	  my %hash = map { $_ => 1 } @species_;
#	  my @species = sort keys %hash;

	  ## Get species names - 'old' way
	  my @species;
	  my $dbh = DBI->connect("DBI:mysql:host=ensembldb.ensembl.org:port=5306", "anonymous", "", {'RaiseError' => 0});
	  my $sth = $dbh->prepare("SHOW DATABASES");
	  $sth->execute();
	  my $previous_org = "init";
	  while (my $ref = $sth->fetchrow_hashref()) {
#	      if ($ref->{'Database'} =~ /_core_\d+/) {
	      if ($ref->{'Database'} =~ /_core_$dbversion/) {
		  $ref->{'Database'} =~ s/_core_.+//;
		  if ($ref->{'Database'} ne $previous_org) {
		      push @species, $ref->{'Database'};
		      $previous_org = $ref->{'Database'};
		  }
	      }
	  }
	  $sth->finish();
	  $dbh->disconnect();

	  # Get classifications
	  foreach my $beast (@species) {
	      my $meta_container = Bio::EnsEMBL::Registry->get_adaptor($beast, 'Core', 'MetaContainer');
#	      my $_species = $meta_container->get_Species();  # get_Species is deprecated
	      my $_species = $meta_container->get_classification();
	      my @info = @{$_species};
	      $classifications{$beast} = [@info];
	  }
      }

      foreach my $homology (@{$homologies}) {

      # You will find different kind of description UBRH, MBRH, MBRH, RHS, YoungParalogues see ensembl-compara/docs/docs/schema_doc.html for more details
	&RSAT::message::Info("Homology type:",$homology->description,"Taxon:", $homology->taxonomy_level) if ($main::verbose >= 1);

	# each homology relation have only 2 members, you should find there the initial member used in the first fetching
#	foreach my $member_attribute (@{$homology->get_all_Members}) { # method depreceated from version 68
	foreach my $member (@{$homology->get_all_GeneMembers}) {

	    # for each Member, you get information on the Member specifically and in relation to the homology relation via Attribute object
#	    my ($member, $attribute) = @{$member_attribute};  # Member object has all the attributes itself
	    &RSAT::message::Info("Gene ID:", $member->stable_id, "Taxon ID:", $member->taxon_id) if ($main::verbose >= 1);
#	    &RSAT::message::Info("Perc. id.:", $attribute->perc_id,"Perc. pos.:",$attribute->perc_pos,"Perc. cov.:",$attribute->perc_cov) if ($main::verbose >= 1);
#	    &RSAT::message::Info("Perc. id.:", $member->perc_id,"Perc. pos.:",$member->perc_pos,"Perc. cov.:",$member->perc_cov) if ($main::verbose >= 1);

	    unless ($member->stable_id eq $ortho_id) {

		$bin_name = $member->taxon->binomial;
		$bin_name =~ s/\s+/_/g;

		$common_name = $member->taxon->common_name;
		if ($common_name) {
		    $common_name =~ s/\s+/_/g;
		} else {
		    $common_name = $bin_name;
		}
		&RSAT::message::Info("Common name:", $common_name) if ($main::verbose >= 1);

		my @homolog_classification;
		if ($taxon) {
#		    @homolog_classification = split (/ /, $member->taxon->classification); # this is too slow
		    my $classif_name = lc($bin_name);
		    if ($classif_name eq "canis_lupus_familiaris") {$classif_name = "canis_familiaris"};
		    if ($classif_name eq "gorilla_gorilla_gorilla") {$classif_name = "gorilla_gorilla"};
		    @homolog_classification = @{$classifications{$classif_name}};
		    &RSAT::message::Debug (join(" ", "Homolog classification :", @homolog_classification)) if ($main::verbose >= 3);
		}

		# Prints all homologs to table if asked for
		if ($homologs_table) {
#		    print $table_handle join("\t", $member->stable_id, $bin_name, $member->description, $homology->description, $homology->subtype, $attribute->perc_id, $attribute->perc_pos, $attribute->perc_cov, $ortho_id, $compara_taxon->binomial, "\n");
		    print $table_handle join("\t", $member->stable_id, $bin_name, $member->description, $homology->description, $homology->taxonomy_level, $member->perc_id, $member->perc_pos, $member->perc_cov, $ortho_id, $compara_taxon->binomial, "\n");
		}

		if ($ortho_type) {
		    if ($taxon) {
			if (($homology->description =~ /$ortho_type/) && (lc($taxon) eq lc($homology->taxonomy_level))) {	
			    $taxon_filter_flag = 1;
			    my $gene = $member->get_Gene;
			    &Main($gene, $bin_name);
			} else {
			    foreach my $tax (@homolog_classification) {
				if (($homology->description =~ /$ortho_type/) && (lc($taxon) eq lc($tax))){
				    $taxon_filter_flag = 1;
				    my $gene = $member->get_Gene;
				    &Main($gene, $bin_name);
				}
			    }
			}
		    } else {
			if ($homology->description =~ /$ortho_type/) {
			    my $gene = $member->get_Gene;
			    &Main($gene, $bin_name);
			}
		    }
		} else {
		    if ($taxon) {
			if (lc($taxon) eq lc($homology->subtype)) {
			    $taxon_filter_flag = 1;
			    my $gene = $member->get_Gene;
			    &Main($gene, $bin_name);
			} else {
			    foreach my $tax (@homolog_classification) {
				if (lc($taxon) eq lc($tax)) {
				    $taxon_filter_flag = 1;
				    my $gene = $member->get_Gene;
				    &Main($gene, $bin_name);
				}
			    }
			}
		    } else {
			my $gene = $member->get_Gene;
			&Main($gene, $bin_name);
		    }
		}
	    }
	}
      }
      if ($taxon && ($taxon_filter_flag == 0)) {
	  &RSAT::message::Warning (join("\t", "None of the homologs matches your taxonomic filter", $taxon, "Are you sure of its validity?"));
      }
}

################################################################
#### Print sequence to file
sub PrintSequence {
    my ($sequence, $fasta_header) = @_;

    $sequence = &FoldSequence($sequence, $lw);

    &RSAT::message::Debug("Sequence:") if ($main::verbose >= 3 && !$rm);
    &RSAT::message::Debug("Repeat masked sequence:") if ($main::verbose >= 3 && $rm);
    &RSAT::message::Debug($fasta_header) if ($main::verbose >= 3);
    &RSAT::message::Debug($sequence) if ($main::verbose >= 3);

    # Export sequence to file
    if ($sequence) {
	print $fh "$fasta_header\n";
	print $fh "$sequence\n";
    }
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

	Returns upstream, downstream, intronic, exonic or UTR DNA
	sequences for a list of query genes.

CATEGORY
	genomics
	sequences

REMARK

    This script requires a local installation of the EnsEMBL Perl Core
    and Compara APIs (see RSAT installation guide).

OPTIONS
	-org organism
	        underscore between words (eg 'homo_sapiens')

	        If this option is not used, the option -dbname must be used
	         instead.

	        (type 'supported-organisms | grep EnsEMBL' to obtain the list of supported
	         organisms)

        -ensemblhost
                address of ensembl database server (default is EBI server)

	-dbname	name of EnsEMBL database
		(alternative to organism)

        -dbversion
	        version of ensembl database (example: 47)

	-feattype
		Feature type.
		Supported: $supported_feattypes
		Defaut: $feattype

	-type	sequence type
		Currently supported sequence types
			upstream (default)
			downstream
                        feature

       -utr utr_type

              Type(s) of UTR (untranslated region) to return.

              Supported: all | 5prime | 3prime

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
                Get sequences for all transcripts of genes.
                Use -uniqseqs if you do motif discovery afterwards

        -uniqseqs
                With -alltranscripts, returns only non-redondant sequences

        -firstintron
                With feattype intron, get only first intron sequence

        -noncoding
                With feattype exon, get only non-coding (part of) exons

        -chrom  Chromosome name or number (to use with -left and -right)

        -left   Left limit of sequence to retrieve

        -right  Right limit of sequence to retrieve

        -strand Strand of sequence to retrieve when using -left and -right. Values: 1, -1

        -ftfile Feature file

        -ftfileformat
                Feature file format. Supported: ft, gft

	-ortho  Retrieve homologous sequences from EnsEMBL Compara databases

	-ortho_type Type
                Filter on homology type. (example: ortholog, ortholog_one2one)

        -homologs_table File name
                Prints homology info to a tab delimited file

        -taxon Taxon  Filter on taxonomic level (example: Mammalia)

	-header_org   Type of organism name to use in the fasta header (scientific, common or none).
		      Default is scientific. Common name is only accessible with -ortho.

        -label label_type
	       Information used as sequence label in the fasta header. 

	       Supported label types: 

	       -label query 
	       	      use as sequence label the identifier or name used as query. 

End_help
    close HELP;
    exit;
}

################################################################
#### list of options
sub PrintShortHelp {
  open(HELP, "| less");
  print HELP<<End_short_help;
retrieve-ensembl-seq options
----------------------------
-org		organism
-ensemblhost    address of ensembl database server (default is EBI server)
-dbname         name of ensembl db
-feattype	accepted feature types. Supported: $supported_feattypes
-type		upstream | downstream | feature
-utr utr_type   Type(s) of UTR to return. Supported: all | 5prime | 3prime
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
-uniqseqs       with -alltranscripts, returns only non-redondant sequences
-firstintron    with feattype intron, get only first intron sequence
-noncoding      with feattype exon, get only non-coding (part of) exons
-chrom          chromosome name or number (to use with -left and -right)
-left           left limit of sequence to retrieve
-right          right limit of sequence to retrieve
-strand         strand of seauence to retrieve when using -left and -right. Values: 1, -1
-ftfile         feature file
-ftfileformat   feature file format. Supported: ft, gft
-ortho          retrieve homologous sequences from EnsEMBL Compara databases
-ortho_type     homology type to filter on
-homologs_table file on which to print homology info
-taxon          taxonomic level to filter on
-header_org     type of organism name to use in the fasta header (common, scientific or none)
-label          information used to label sequence in the fasta header (supported: -label query)

End_short_help
  close HELP;
exit;
}

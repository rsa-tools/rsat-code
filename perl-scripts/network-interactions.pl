#!/usr/bin/env perl
# network-interactions version 1

use strict;
use warnings;

BEGIN {
  if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
  }
}
require "RSA.lib";

################################################################
## Main package
package main;
{

  ################################################################
  ## Initialise parameters
  our $start_time = &RSAT::util::StartScript();
  our $program_version = do { my @r = (q$Revision: 1.48 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
  #    $program_version = "0.00";

  our %infile = ();
  our %outfile = ();

  our $verbose = 0;
  our $in = \*STDIN;
  our $out = \*STDOUT;

  # command-line arguments
  our $TFs; # file with TFs
  our $cres; # bed file of regulatory sequences
  our $network; # input network
  our $databases = ""; # matrixes dbs
  our $genome = ""; # working genome version
  our $orgdb = ""; # alternative to databases
  our $report_net = 0; # boolean, 1 if network is provided
  #our $net_type = "not_grn"; # default for network type

  # other variables
  my $enhancer; # filehandle for cres bed file
  my $new_enhancer; # filehandle of filtered bed file
  my @genes_reg = (); # genes regulated
  my %nonRep_genes_reg; # hash of non-repetitive genes regulated
  my $val = 0; # values of has %nonRep_genes_reg
  my $i; # general counter variable

  my %nonRepTFs; # hash of non repetitive TFs
  my $TF; # filehandle from TFs
  my $report_tfs_noRegSeq = ""; # string with TFs not found on bed file
  my $report_tfs_RegSeq = ""; # string with TFs found on bed file

  my (@net_old_allgenes, @net_old_onlyTFs); # networks from input
  my %noRegSeq = (); # genes lacking regulatory sequence
  my $noRegSeq_reporting = ""; # string of genes lacking reg seq
  my %RegSeq = (); # genes with regulatory sequence on input
  my $RegSeq_reporting = ""; # string of genes with reg seq
  my $z = 0; # values of a hash
  my $flag = 0;
  my $ngenes = 0; # number of genes regulated
  my $ntfs = 0; # number of genes regulated

  my @matfiles_path = (); # rsat path to matrixes in use
  my @dbs = ();  # list of databases
  my @fields = (); # general var, list of fields from rows in a file
  my @nonRepTFs_k; # TFs list
  my $fasta; # filehandle of sequences
  my $new_fasta; # filehandle of sequences with headers
  my $cmd = ""; # general variable for &doit
  my $ACs = ""; # list of ACs from matrixes to use
  my @seq_id; # seq_id splitted from matrix-scan results

  my (@net_new_allgenes, @net_new_onlyTFs); # networks
  my $last_gene = ""; # filter
  my $last_tf = ""; # filter
  my @deadend = (); # binary list with length of total TFs, 1 if the TF has no direct interactions
  my $matrixscan_re; # filehandle of matrix-scan results

  my $interaction_gene_row; # filehandle of file with all direct interactions (main output)
  my $intersection; # filehandle  of file with overlap of networks
  my $rejected; # filehandle of file with "rejected" interactions
  my $found; # filehandle of file with found interactions

  my $count = 0;
  my %nonRepTFs_inv; # inverse hash of nonRepTFs

  my $direct; # filehandle of file with GRN direct interactions (main output)
  my $path = ""; # construction of paths in node

  ################################################################
  ## Read argument values
  &ReadArguments();

  ################################################################
  ## Check argument values
  unless ($main::TFs) {
    &RSAT::error::FatalError("File with the network's transcription factors must be provided\n It must be a one-column file with list of TFs");
  }
  unless ($main::cres) {
    &RSAT::error::FatalError("Bed file of enhancer sequences must be provided\n On the 4th column the gene it is referring to must be specified");
  }
  unless ($main::genome) {
    &RSAT::error::FatalError("Genome version must be specified with option -genome (e.g. mm9, hg19)");
  }
  if($report_net){
    unless ($main::network) {
      &RSAT::error::FatalError("Turned on network reporting (-rep_net) but no network (-net) provided");
    }
  }

  # Decide which transfac files to use
  # CODE: add all from organism option

  if($main::databases eq ""){ # default option
    #print("Entered default option of dbs\n");
    $matfiles_path[0] = "rsat/public_html/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.tf";
    #print($matfiles_path[0]."\n");
  } else {
    if($main::databases eq "all"){ # CODE: falta anadir algo que soporte que meta los nombres de dbs como sea
      @dbs = qw(ENCODE epigram Hocomoco_human Hocomoco_mouse Homer hPDI HT_Selex Jolma_2013 Jolma_2015 NCAP-SELEX CAP-SELEX RSAT JASPAR cisBP Yeastract DrosophilaTFs);
    } else { @dbs = split(/ |,|;/, $main::databases) }

    # retrieve paths to databases from db_matrix_files.tab
    open(my $dbs_file, "<", "rsat/public_html/motif_databases/db_matrix_files.tab") or die "Cannot open file $!";
    $i = 0;
    while(my $row = <$dbs_file>){
      chomp($row);
      next if( $row =~ /^;/);
      next if( $row =~ /^#/);

      foreach my $db (@dbs){
        if( $row =~ /^$db/){
          @fields = split(/\t/,$row);
          $matfiles_path[$i++] = "rsat/public_html/motif_databases/".$fields[2];
          #print(" matrix file path number ".$c." : ".$matfiles_path[$c]."\n");
          #$c++;
        }
      }
    }
    close $dbs_file;
  }
  print("# of paths = ".@matfiles_path."\n");

  ################################################################
  ####################################### Save Input Data

  # Rewrite to bedfile without NAs or if no gene is refered to
  # Also, save the genes refered to in bed file of enhancers
  #        I want to preserve the order
  $val = 0; # first value is always 1 in hash
  $i = 1;

  open( $enhancer, "<", $main::cres) or die "Cannot open file $!";
  open( $new_enhancer, ">", "prueba_regulatory_regions.bed") or die "Cannot open file $!";
  print $new_enhancer "# Bed file used for cis-regulatory regions\n";
  #$enhancer = &OpenInputFile($main::cres);
  while(my $enh_row = <$enhancer>){
    chomp($enh_row);
    next unless ($enh_row =~ /\S/); ## Skip empty rows
    next if ($enh_row =~ /^;/); ## Skip comment rows
    next if ($enh_row =~ /^#/); ## Skip header rows

    if($enh_row =~ /^chr\w{1,2}\t\d+\t\d+\t(\w+).*/){
      if($1 ne "NA"){

          # write to new bed file and save order of genes
          print $new_enhancer $enh_row."\n";
          push(@genes_reg,$1);
          # CODE: differentiate genes if they're from TFS or not...?

          # add values in sequential order
          # gene is key, value is a number
          if( $val == 0 ){
            $nonRep_genes_reg{$1} = $val++;
          } elsif(!$nonRep_genes_reg{$1}) {
            $nonRep_genes_reg{$1} = $val++;
          }

      } #else {    # it does filters out only NA
        #print($1."\n");
      #}

    } else {
      print(" You are missing a gene name in bed file at line: ".$i."\n");
      #&RSAT::error::FatalError("Error: not a fourth column of genes in file ".$main::cres."\n");
      #exit(1);
    }
    $i++;
  }

  close $enhancer;
  close $new_enhancer;

  #foreach my $element (sort keys %nonRep_genes_reg){
  #  print("\tkey ".$element."\n");
  #}

  print ("# of genes regulated = ".@genes_reg."\n");
  print ("# of nonRep genes regulated = ".%nonRep_genes_reg."\n");

  #############################################################
  # Save non-repited TFs
  $i = 0;

  open( $TF,"<",$main::TFs) or die "Cannot open file $!";
  while( my $row = <$TF> ){
    next unless ($row =~ /\S/); ## Skip empty rows
    next if ($row =~ /^;/); ## Skip comment rows
    next if ($row =~ /^#/); ## Skip header rows
    chomp($row);

    # distinguish if their sequences were reported
    if ( $nonRep_genes_reg{$row} ){

      # add values in sequential order, beginning in 1
      if( $i == 0 ){
        $nonRepTFs{$row} = $i++;
      } elsif(!$nonRepTFs{$row}) {
        $nonRepTFs{$row} = $i++;
      }

      $report_tfs_RegSeq .= $row.",";

    } else {
      $report_tfs_noRegSeq .= $row.",";
    }

  }
  print("# of nonRepTFs = ".%nonRepTFs."\n");

  # notice user about the lack of regulatory seqences
  print("\nSeems like you haven't reported the regulatory sequences for the following TFs : ".$report_tfs_noRegSeq."\n");
  if( $report_net){
    print(" The program will take into account the rest for further network comparision : ".$report_tfs_RegSeq."\n\n");
  }

  #############################################################
  # Save input network and create input GRN
  $flag = 0;
  $ngenes = %nonRep_genes_reg-1;
  $ntfs = %nonRepTFs -1;

  if($report_net){

    #initialize matrixes to zeros
    # input all-genes net
    #print($val."\n");
    for $i (0..$ngenes){
      for my $j (0..$ngenes){
        $net_old_allgenes[$i][$j] = 0;
      }
    }
    # input only-TFs net
    for $i (0..$ntfs){
      for my $j (0..$ntfs){
        $net_old_onlyTFs[$i][$j] = 0;
      }
    }

    $i = 0;
    # Save input network interactions
    open(my $net, "<", $main::network) or die "Cannot open file $!";
    while (my $net_row = <$net>){
      chomp($net_row);
      next unless ($net_row =~ /\S/); ## Skip empty rows
      next if ($net_row =~ /^;/); ## Skip comment rows
      next if ($net_row =~ /^#/); ## Skip header rows

      @fields = split(/\t/, $net_row);

      # make list of genes in network with regulatory sequences reported and not reported
      if($nonRep_genes_reg{$fields[0]}){
        $RegSeq{$fields[0]} = $z++;
        $flag++;
      } else {
        $noRegSeq{$fields[0]} = $z++;
      }
      if($nonRep_genes_reg{$fields[1]}){
        $RegSeq{$fields[1]} = $z++;
        $flag++;
      } else {
        $noRegSeq{$fields[1]} = $z++;
      }

      # en el reporte de comparacion de redes solo se tomaran en cuenta las interacciones en las que ambos genes tenian la sec reportada en -cre
      # if both genes have a reported regulatory sequence
      if( $flag == 2 ){
        # construct network of all genes
        $net_old_allgenes[ $nonRep_genes_reg{$fields[0]}-1 ][ $nonRep_genes_reg{$fields[1]}-1 ] = 1;
        # if the interaction is between TFs, add this one to the input's only-TFs network
        if( $nonRepTFs{$fields[0]} && $nonRepTFs{$fields[1]} ){
          $net_old_onlyTFs[ $nonRepTFs{$fields[0]}-1 ][ $nonRepTFs{$fields[1]}-1 ] = 1;
        }
      }
      # reinitialize value
      $flag =0;

    }

    # message if sequences of any genes in network were not reported
    if(%noRegSeq){
      foreach my $gene (sort keys %noRegSeq){
        $noRegSeq_reporting .= "".$gene.",";
        #print($gene."\n");
      }
      foreach my $gene (sort keys %RegSeq){
        $RegSeq_reporting .= "".$gene.",";
      }
      # Tell user which regulatory sequences are lacking and which genes are used in the program
      print("\nSeems like you haven't reported the regulatory sequences for the following genes in your input network: ".$noRegSeq_reporting."\n");
      print("The program will use the rest of the genes for further network comparison: ".$RegSeq_reporting."\n\n");
    }

    #for(my $a = 0; $a <$i; $a++){
    #    print($net_old_allgenes[$a][0]."    ".$net_old_allgenes[$a][1]."\n");
    #}


  } # report_net


  ################################################################
  ####################################### RSAT processing

  # Create the file with all the matrixes to use
  #   retrieve-matrix foreach database foreach TF
  @nonRepTFs_k = keys %nonRepTFs;


  foreach my $db (@matfiles_path){
    #print(" Second: ".$db."\n");
    $ACs = &retrieveACsfromTFs($db, @nonRepTFs_k);
    $cmd = "retrieve-matrix -i ".$db." -id ".$ACs." >> prueba_matrixes.tf \n"; # CODE: tengo que quitarle los logs porque se guardan en prueba_matrixes.tf
    &doit($cmd);
  }

  # Retrieve sequences of regulatory elements
  print("\nRetrieving enhancer sequences ...");
  &doit("fetch-sequences -i prueba_regulatory_regions.bed -genome ".$main::genome." -v 1  -header_format galaxy -o prueba_".$main::genome."_cre_seqs.fa");
  print("\tDone!\n");

  # Add regulated genes names to fasta of regulatory sequences
  #     Im using 'external_id' not ensembl_id
  open( $fasta, "<", "prueba_".$main::genome."_cre_seqs.fa") or die "Cannot open file $!";
  open( $new_fasta, ">", "prueba_".$main::genome."_cre_seqs_geneheader.fa") or die "Cannot open file $!";

  $i = 0;
  while(my $fa_row = <$fasta> ){
    chomp($fa_row);

    if($fa_row =~ /^>/){
      print $new_fasta $fa_row."_".$genes_reg[$i++]."\n";# this must be at the identifier line e.g.: >seq1_$genes_reg[$i++]
                                                         #  only so that matrix-scan takes this whole name as seq_id
    } else {
      print $new_fasta $fa_row."\n";
    }
  }
  close $fasta;
  close $new_fasta;

  # Search motifs in all regulatory sequences
  print("\n matrix-scan running ...");
  $cmd = "matrix-scan -quick -v 1 -matrix_format transfac -m prueba_matrixes.tf -bginput -markov 1";
  $cmd.= " -i prueba_".$main::genome."_cre_seqs_geneheader.fa -o prueba_matrix_scan_results.tsv";
  $cmd.= " -id_as_name -lth score 5";
  &doit($cmd);
  print("\tDone!\n");

  ################################################################
  ####################################### Networks comparision

  #if($report_net){
    # initialize new matrixes
    for $i (0..$ngenes){
      for my $j (0..$ngenes){
        $net_new_allgenes[$i][$j] = 0;
      }
    }
    # input only-TFs net
    for $i (0..$ntfs){
      for my $j (0..$ntfs){
        $net_new_onlyTFs[$i][$j] = 0;
      }
      # initialize to zeros
      $deadend[$i] = 0;
    }
  #} # if

  print("Network comparision running ...");
  # Create file of new directed network:
  #     one row per interaction: TF1   TF2

  open( $matrixscan_re, "<", "prueba_matrix_scan_results.tsv") or die "Cannot open file $!";
  # network of genes (complete)
  #open( my $interaction_gene_row, ">", "prueba_direct_interactions_complete.tsv") or die "Cannot open file $!";
  # network of TFs only (GRN)
  #open( my $interaction_TF_row, ">", "prueba_direct_interactions_GRN.tsv") or die "Cannot open file $!";

  $last_gene = "";
  $last_tf = "";

  while(my $results_row = <$matrixscan_re>){
    chomp($results_row);
    next unless ($results_row =~ /\S/); ## Skip empty rows
    next if ($results_row =~ /^;/); ## Skip comment rows
    next if ($results_row =~ /^#/); ## Skip header rows

    @fields = split(/\t/,$results_row);
    @seq_id = split(/_/,$fields[0]); # >seq_id_gene

    # for new interactions
    if ( $last_gene ne $seq_id[$#seq_id] || $last_tf ne $fields[2] ){

      #print $interaction_gene_row "".$fields[2]."\t".$seq_id[$#seq_id]."\n"; #    TFx (regulates) Geney

      # Save new network interactions for comparison
      #if($report_net){
        if($nonRep_genes_reg{$fields[2]} && $nonRep_genes_reg{$seq_id[$#seq_id]}){
            $net_new_allgenes[ $nonRep_genes_reg{$fields[2]}-1 ][ $nonRep_genes_reg{$seq_id[$#seq_id]}-1 ] = 1;
        }
      #}

      # if it is a TF-TF interaction, save it on the GRN
      if ($nonRepTFs{$fields[2]} && $nonRepTFs{$seq_id[$#seq_id]}){
        #print $interaction_TF_row "".$fields[2]."\t".$seq_id[$#seq_id]."\n"; #    TFx (regulates) TFy

        # Save new network interactions for comparison
        # all is based on the regulatory seqs so every gene appears on %nonRep_genes_reg
        if($report_net){
          $net_new_onlyTFs[ $nonRepTFs{$fields[2]}-1 ][ $nonRepTFs{$seq_id[$#seq_id]}-1 ] = 1;
        }

      } # if TF-TF interaction

      $last_gene = $seq_id[$#seq_id];
      $last_tf = $fields[2];
    }

  } #  while
  #close $interaction_gene_row;
  #close $interaction_TF_row;
  close $matrixscan_re;
  print("\t Done\n");


  # Comparisons and write-to-file reporting

  if( $report_net ){

    ## Comparison of complete gene networks
    open( $interaction_gene_row, ">", "prueba_direct_interactions_complete.tsv") or die "Cannot open file $!";
    print $interaction_gene_row "# Complete Network\n#TFx (regulates) Geney\n"; # header
    open( $intersection, ">", "prueba_network_intersection.tsv") or die "Cannot open file $!";
    print $intersection "# Interactions in common from:\n#\t".$main::network." and prueba_direct_interactions_complete.tsv\n#TFx (regulates) Geney\n"; # header
    open( $rejected, ">", "prueba_network_rejected_interactions.tsv") or die "Cannot open file $!";
    print $rejected "# Interactions in:\t".$main::network." not found in prueba_direct_interactions_complete.tsv\n#Genex\tGeney\n"; # header
    open( $found, ">", "prueba_network_found_interactions.tsv") or die "Cannot open file $!";
    print $found "# Newly found interactions in:\tprueba_direct_interactions_complete.tsv not found in ".$main::network."\n#TFx\tGeney\n"; # header

    my %nonRep_genes_reg_inv = reverse %nonRep_genes_reg;

    for $i (0..$ngenes){
      for my $j (0..$ngenes){

        # write interaction to file according to the case: intersection, rejected or newly found
        if( $net_old_allgenes[$i][$j] ){
          if( $net_new_allgenes[$i][$j] ){
            # report direct interaction TFx Geney
            print $interaction_gene_row "".$nonRep_genes_reg_inv{ $i+1 }."\t".$nonRep_genes_reg_inv{ $j+1 }."\n";
            # overlap (only interactions, not lack of them)
            print $intersection "".$nonRep_genes_reg_inv{ $i+1 }."\t".$nonRep_genes_reg_inv{ $j+1 }."\n";
          } else {
            # rejected
            print $rejected "".$nonRep_genes_reg_inv{ $i+1 }."\t".$nonRep_genes_reg_inv{ $j+1 }."\n";
          }
        } elsif ( $net_new_allgenes[$i][$j] ) {
            # report direct interaction TFx Geney
            print $interaction_gene_row "".$nonRep_genes_reg_inv{ $i+1 }."\t".$nonRep_genes_reg_inv{ $j+1 }."\n";
            # newly found
            print $found "".$nonRep_genes_reg_inv{ $i+1 }."\t".$nonRep_genes_reg_inv{ $j+1 }."\n";
        }

      } # for
    } # for

    close $interaction_gene_row;
    close $intersection;
    close $rejected;
    close $found;

    ## Comparison of GRNs ?
  } # report net

  ################################################################
  ## Report all direct interactions (main output)

  open( $direct, ">", "prueba_GRN_direct_interactions.tsv") or die "Cannot open file $!";
  print $direct "#Regulator\tRegulated\n"; # header

  %nonRepTFs_inv = reverse %nonRepTFs;
  $count = 0;

  for my $i (0..$ntfs){
    for my $j (0..$ntfs){

      # write to file direct interactions
      if( $net_new_onlyTFs[$i][$j] ){
        print $direct "".$nonRepTFs_inv{ $i+1 }."\t".$nonRepTFs_inv{ $j+1 }."\n";
      } else {
        # add node to list of deadend nodes
        if($count == $ntfs){
          $deadend[$i] = 1;
        }
        $count++;
      }

    } # for
  } # for

  close $direct;


  ################################################################
  #  CASES OF DIRECT AND INDIRECT INTERACTIONS SIMULTANEOUSLY
  # basic case
  # TF1 TF2
  # TF2 TF3
  # TF1 TF3
  $path = "\n"; # string to concat paths of 3 nodes

  ## Report all paths of three (instead of cases of direct and indirect interactions simoultaneously)
  open( my $indirect, ">", "prueba_GRN_indirect_interactions.tsv") or die "Cannot open file $!";
  print $indirect "# Report of 3-nodes direct interactions\n# cases of no direct interactions are specified.\n#node1 ==> node2 ==> node3\n\n"; #header

  for my $row (0..$ntfs){
    #$path = "";
    # check if the node is a deadend
    if ( !$deadend[$row] ){
      for my $con1 (0..$ntfs){
        # find direct connections
        if( $net_new_onlyTFs[$row][$con1] ){
          # construct path
          $path .= $nonRepTFs_inv{ $row+1 }." ==> ".$nonRepTFs_inv{ $con1+1 };
          # report if its a connection to itself
          if ($row == $con1){
            print $indirect "".$path."\n";
          } elsif ( !$deadend[$con1] ){
            for my $con2 (0..$ntfs){
              # skip it its a connection to itself
              if( $con1 == $con2){ next; }
              if( $net_new_onlyTFs[$con1][$con2] ){
                # report 3-node connection to file
                $path .= " ==> ".$nonRepTFs_inv{ $con2+1 }."\n";
                print $indirect $path;
                # set the 2-node connection for more 3-node connections
                $path = "\n".$nonRepTFs_inv{ $row+1 }." ==> ".$nonRepTFs_inv{ $con1+1 };
              } # direct connection 2 found
            }
          } else {
            # if it was a dead end report 2-node connection
            print $indirect "".$path."\n";
          }
          $path = "\n";

        } # direct connection 1 found
      } # for connection 1
    } else {
      print $indirect "".$nonRepTFs_inv{ $row+1 }." has no targets.\n";
      $path = "\n";
    }

  } # for row

  close $indirect;

  ## Print verbose
  #&Verbose() if ($main::verbose >= 1);

  #&close_and_quit();

} # main

=pod
################################################################
## Close output file and quit
sub close_and_quit {

  ## Report execution time
  my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
  print $main::out $exec_time if ($main::verbose >= 1); ## only report exec time if verbosity is specified

  # Close output file
  if ($outfile{output}) {
    close $main::out;
    &RSAT::message::TimeWarn("Output file", $outfile{output}) if ($main::verbose >= 2);
  }

  ## CLOSE OTHER FILES HERE IF REQUIRED

  exit(0);
}
=cut
################################################################
## Display full help message
sub PrintHelp {
  system "pod2text -c $0";
  exit()
}

## Display short help message
sub PrintOptions {
  &PrintHelp();
}

################################################################
## Read arguments
sub ReadArguments {
  my $arg;
  my @arguments = @ARGV; ## create a copy to shift, because we need ARGV to report command line in &Verbose()
  while (scalar(@arguments) >= 1) {
    $arg = shift (@arguments);
=pod
=head1 OPTIONS
=over 4

=item B<-v #>
Level of verbosity.
=cut
    if ($arg eq "-v") {
      if (&IsNatural($arguments[0])) {
          $main::verbose = shift(@arguments);
      } else {
          $main::verbose = 1;
      }
=item B<-h>
Display full help message.
=cut
    } elsif ($arg eq "-h") {
      &PrintHelp();
=item B<-help>
Same as -h.
=cut
    } elsif ($arg eq "-help") {
      &PrintOptions();
    }# elsif ($arg eq "-i") {
      #$main::infile{input} = shift(@arguments);
    #}
=item B<-tfs TFs_infile>
Mandaroty.
One-column input file with list of TFs in network.
=cut
    elsif ($arg eq "-tfs") {
      $main::TFs = shift(@arguments);
=item B<-cre cre_infile>
Mandaroty.
Bed file with of regulatory sequences per gene in network.
Each sequence must refer to its regulated-gene on the 4th column.
=cut
    } elsif ($arg eq "-cre") {
      $main::cres = shift(@arguments);
=item B<-genome genome_version>
Mandaroty.
Working genome version.
=cut
    } elsif ($arg eq "-genome") {
      $main::genome = shift(@arguments);
=item B<-db databases>
Database(s) to use separated by commas.
=cut
    } elsif ($arg eq "-db") {
      $main::databases = shift(@arguments);
=item B<-net network_infile>
Mandaroty if -report_net is set to 1.
The network must be a two-column file with each row containing a single interaction (direct or indirect).
=cut
    } elsif ($arg eq "-net") {
      $main::network = shift(@arguments);
=item B<-report_net 0|1>
If report_net is set to 1, then network reporting is turned on,
  i.e., differences and overlap between the input network and the new network are outputted
=cut
    } elsif ($arg eq "-report_net") {
      $main::report_net = shift(@arguments);
#=item B<-net_type grn|not_grn>
#Network type specification. Could be GRN (Gene Regulatory Network) containing only TFs,
#of not_grn containing both genes codifying for TFs and genes codifying for not TFs.
#=cut
    } #elsif ($arg eq "-net_type") {
      #$main::net_type = shift(@arguments);
    #}# elsif ($arg eq "-o") {
      #$outfile{output} = shift(@arguments);
    #}
    else {
      &FatalError(join("\t", "Invalid option", $arg));
    }
  } # while

=pod
=back
=cut

} # sub ReadArguments

=pod
################################################################
## Verbose message
sub Verbose {
    printdatabases $out "; template ";
    &PrintArguments($out);
    printf $out "; %-22s\t%s\n", "Program version", $program_version;
    if (%main::infile) {
        print $out "; Input files\n";
        while (my ($key,$value) = each %main::infile) {
            printf $out ";\t%-13s\t%s\n", $key, $value;
        }
    }
    if (%main::outfile) {
        print $out "; Output files\n";
        while (my ($key,$value) = each %main::outfile) {
            printf $out ";\t%-13s\t%s\n", $key, $value;
        }
    }
}
=cut

################################################################
## Obtain list of the accesions from specified TFs to call retrieve-matrix on
#$ACs = &retrieveACsfromTFs($db, @nonRepTFs);

sub retrieveACsfromTFs {
    my ($database, @TFnames) = @_;
    #print(" In function: ".$database."\n");

    # Check in which field the TF name appears
    my $isinAC = 0; # boolean to check if TFname is in AC field or ID
    my @cols; # list of elementes in field AC
    my $retrievingACs = ""; # final "list" with all AC's from the TFnames required

    open( my $db, "<", $database) or die "Cannot open file $!";
    while(my $db_row = <$db>){
        chomp($db_row);
        next unless($db_row =~ /^AC/ || $db_row =~ /^ID/);
        # este paso solo es para no repetir si viene en ambas columnas
        if( $db_row =~ /^AC\s+[^MA{0,1}\d{4,5}(_\d){0,1}(\.\d+){0,1}]/){ # si encuentra algo distinto de MA..
            $isinAC = 1;
            last;
        } elsif( $db_row =~ /^ID\s+[^MA{0,1}\d{4,5}(_\d){0,1}(\.\d+){0,1}]/){
            $isinAC = 0;
            last;
        } # else, no TFname provided
    } #while
    close $db;

    # Generate array of ACs
    open( $db, "<", $database) or die "Cannot open file $!";
    while(my $db_row = <$db>){ # does this continues from the line before or does it starts over again
        chomp($db_row);
        next unless($db_row =~ /^AC/ || $db_row =~ /^ID/);

        if( $db_row =~ /^AC/){
            @cols = split(/\s+/,$db_row); # guarda el AC para cada TF
        }
        foreach my $name (@TFnames){
            if( $db_row =~ /.*$name.*/){
              if( !$isinAC && $db_row =~ /^AC/){ # if TF name is found in AC field
                $retrievingACs .= $cols[1].",";
              } elsif( $isinAC && $db_row =~ /^ID/){ # if TF is found in ID field
                $retrievingACs .= $cols[1].",";
              }
            } # if
        } # foreach
    } # while
    close $db;

    #print("ACs = ".$retrievingACs."\n");

    return $retrievingACs;
} # sub

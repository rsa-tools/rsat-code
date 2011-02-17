#!/usr/bin/perl -w
############################################################
#
# $Id: get-ensembl-genome-no-objects.pl,v 1.23 2011/02/17 05:07:46 rsat Exp $
#
# Time-stamp
#
############################################################
#use strict;
use DBI();

# use Devel::Size qw(size total_size);
# use Scalar::Util qw(weaken);

BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
}
require "RSA.lib";
require RSAT::util;

## parsing libraries
require "RSA.lib";
push @INC, $ENV{RSAT}."/perl-scripts/parsers/" if ($ENV{RSAT});
require "lib/load_classes.pl";
require "lib/parsing_util.pl";

## EnsEMBL libraries
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

################################################################
#### main package
package main;
{
    
    ################################################################
    #### initialise parameters
    local $start_time = &RSAT::util::StartScript();
    
    
    local $slice_type = "chromosome";
    local @chromnames = ();
    local $no_seq = 0;
    local $no_masked = 0;
    local $test = 0;
    local $test_number=20;

    %dir = ();

    ### Working directory
    $dir{main} = `pwd`; #### remember working directory
    chomp($dir{main});

    local %outfile = ();
    $outfile{stats}="get-ensembl-genome_stats.txt";
    $outfile{log}="get-ensembl-genome_log.txt";
    $outfile{err}="get-ensembl-genome_err.txt";
    $outfile{contigs} = "contigs.txt";
    $outfile{contig} = "contig.tab";
    $outfile{mrna} = "mrna.tab";
    $outfile{scrna} = "scrna.tab";
    $outfile{trna} = "trna.tab";
    $outfile{rrna} = "rrna.tab";
    $outfile{misc_rna} = "misc_rna.tab";
    $outfile{repeat_region} = "repeat_region.tab";
    $outfile{cds} = "cds.tab";
    $outfile{gene} = "gene.tab";
    $outfile{mrna_name} = "mrna_names.tab";
    $outfile{scrna_name} = "scrna_names.tab";
    $outfile{trna_name} = "trna_names.tab";
    $outfile{rrna_name} = "rrna_names.tab";
    $outfile{misc_rna_name} = "misc_rna_names.tab";
    $outfile{cds_name} = "cds_names.tab";
    $outfile{gene_name} = "gene_names.tab";
    $outfile{organism} = "organism.tab";
    $outfile{organism_name} = "organism_names.tab";
    $outfile{xreference} = "xreference.tab";
    $outfile{intron} = "intron.tab";
    $outfile{exon} = "exon.tab";
    $outfile{utr} = "utr.tab";
    $outfile{coding_exon} = "coding_exon.tab";
    local $verbose = 0;
    
    ## Connection to the EnsEMBL MYSQL database
    $ensembl_host = $ENV{ensembl_host};
#    $ensembl_host = 'ensembldb.ensembl.org';
#    $ensembl_host = 'xserve2.bigre.ulb.ac.be';
    $ensembl_user = "anonymous";
    $dbname = '';
    $org = '';
    
    #### Options for the exported SQL database
    $host= "localhost";
    $schema="ensembl";
    $user="ensembl";
    $password="ensembl";
    
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


    ### verbose ###
    if ($verbose >= 1) {
	print "; get-ensembl-genome ";
	&PrintArguments();
    }
    
    ################################################################
    ## If option -org is used, connect to ensembldb to get list of 
    ## databases and pick the latest one corresponding to chosen organism
    if ($org) {
	&RSAT::message::TimeWarn (join("\t", "Connecting EnsEMBL to get the dbname for organism ", $org, 
				       "host=".$ensembl_host, 
				       "user=".$ensembl_user )) if ($main::verbose >= 1);
	my $dbh = DBI->connect("DBI:mysql:host=$ensembl_host", "$ensembl_user", "", {'RaiseError' => 1});
	my $sth = $dbh->prepare("SHOW DATABASES");
	$sth->execute();
	while (my $ref = $sth->fetchrow_hashref()) {
	    if ($ref->{'Database'} =~ /($org)_core_\d+/) {
		$dbname = $ref->{'Database'};
	    }
	}
	&RSAT::message::Info (join("\t", "dbname = ", $dbname)) if ($main::verbose >= 1);
	$sth->finish();
	$dbh->disconnect();
    }

    ################################################################
    ### open output streams
    unless ($dir{output}) {
	$dir{output} = join("_", "EnsEMBL", $dbname);
	if (scalar(@chromnames) > 0) {
	    $dir{output} .= "_chrom_";
	    $dir{output} .= join("_", @chromnames);
	}
	if ($test) {
	    $dir{output} .= "_test".$test_number;
	}
    }
    &RSAT::util::CheckOutDir($dir{output});
    chdir($dir{output});
    
    ## log file
    open $log, ">".$outfile{log} || die "cannot open error log file".$outfile{log}."\n";

    ## error file
    open ERR, ">".$outfile{err} || die "cannot open error log file".$outfile{err}."\n";

    ## Contig file
    open CTG, ">", $outfile{contigs} || die "cannot open contig file".$outfile{contigs}."\n"; # file with contig IDs

    ## feature files
    $ORG = &OpenOutputFile($outfile{organism});
    &PrintOrgHeader();
    $ORG_NAME = &OpenOutputFile($outfile{organism_name});
    &PrintFtNameHeader("organism", *ORG_NAME);
    $CONTIG = &OpenOutputFile($outfile{contig});
    &PrintContigHeader;
    unless ($no_rep || $seq_only) {
	$REPEAT = &OpenOutputFile($outfile{repeat_region});
	&PrintFtHeader("repeat_region", *REPEAT);
    }
    unless ($seq_only) {
	$MRNA = &OpenOutputFile($outfile{mrna});
	&PrintFtHeader("mrna", *MRNA);
	$TRNA = &OpenOutputFile($outfile{trna});
	&PrintFtHeader("trna", *TRNA);
	$RRNA = &OpenOutputFile($outfile{rrna});
	&PrintFtHeader("rrna", *RRNA);
	$SCRNA = &OpenOutputFile($outfile{scrna});
	&PrintFtHeader("scrna", *SCRNA);
	$MISC_RNA = &OpenOutputFile($outfile{misc_rna});
	&PrintFtHeader("misc_rna", *MISC_RNA);
	$CDS = &OpenOutputFile($outfile{cds});
	&PrintFtHeader("cds", *CDS);
	$GENE = &OpenOutputFile($outfile{gene});
	&PrintFtHeader("gene", *GENE);
	$INTRON = &OpenOutputFile($outfile{intron});
	&PrintFtHeader("intron", *INTRON);
	$EXON = &OpenOutputFile($outfile{exon});
	&PrintFtHeader("exon", *EXON);
	$UTR = &OpenOutputFile($outfile{utr});
	&PrintFtHeader("utr", *UTR);
	$CODING_EXON = &OpenOutputFile($outfile{coding_exon});
	&PrintFtHeader("coding_exon", *CODING_EXON);

	$MRNA_NAME = &OpenOutputFile($outfile{mrna_name});
	&PrintFtNameHeader("mrna", *MRNA_NAME);
	$TRNA_NAME = &OpenOutputFile($outfile{trna_name});
	&PrintFtNameHeader("trna", *TRNA_NAME);
	$RRNA_NAME = &OpenOutputFile($outfile{rrna_name});
	&PrintFtNameHeader("rrna", *RRNA_NAME);
	$SCRNA_NAME = &OpenOutputFile($outfile{scrna_name});
	&PrintFtNameHeader("scrna", *SCRNA_NAME);
	$MISC_RNA_NAME = &OpenOutputFile($outfile{misc_rna_name});
	&PrintFtNameHeader("misc_rna", *MISC_RNA_NAME);
	$CDS_NAME = &OpenOutputFile($outfile{cds_name});
	&PrintFtNameHeader("cds", *CDS_NAME);
	$GENE_NAME = &OpenOutputFile($outfile{gene_name});
	&PrintFtNameHeader("gene", *GENE_NAME);
    }

    ## xref file
#    open $XREF_TABLE, ">".$outfile{xreference} || die "cannot open error log file".$outfile{xreference}."\n";

    ## cds sequences in fasta format
    unless ($seq_only) {
	@dbsplit = split /_core_/, $dbname;
	$org = $dbsplit[0];
	$outfile{pp} = $org."_aa.fasta";
	open PP, ">", $outfile{pp} || die "cannot open sequence file".$outfile{pp}."\n";
    }

    #### print verbose
    &Verbose() if ($verbose);
    
    ################################################################
    ## Open a new connection to EnsEMBL database, but this time we specify the DB name
    &RSAT::message::TimeWarn("Connecting EnsEMBL to retrieve the organism", 
			     "host=".$ensembl_host,
			     "user=".$ensembl_user,
			     "dbname=".$dbname,
			     ) if ($main::verbose >= 1);
    my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $ensembl_host, -user => $ensembl_user, -dbname => $dbname);

    ## Get Species object
    my $meta_container = $db->get_MetaContainer();
    my $tax_id = $meta_container->get_taxonomy_id();
    my $species = $meta_container->get_Species();
    my $organism_name = $species->binomial();
    my $common_name = $species->common_name();
    my @classification = $species->classification();
    &RSAT::message::Info (join (":", "; Classification ", reverse(@classification)))  if ($main::verbose >= 1);
    &RSAT::message::Info(join ("\t", "; Organism = ", $organism_name, 
			       "common name = ", $common_name, 
			       "NCBI taxid = ", $tax_id))
	if ($main::verbose >= 1);

    &RSAT::message::Debug("Db", $db) if ($main::verbose >= 3);

    print $ORG join ("\t", $tax_id, join (";", reverse(@classification)), $organism_name), "\n";
    print $ORG_NAME join ("\t", $tax_id, $organism_name, "primary"), "\n";
    print $ORG_NAME join ("\t", $tax_id, $common_name, "alternate"), "\n";

    my $slice_adaptor = $db->get_SliceAdaptor();
    &RSAT::message::Debug("Adaptor", $slice_adaptor) if ($main::verbose >= 3);
    

    ## Collect the list of chromosomes
    my @slices;

    &RSAT::message::psWarn("BEFORE GETTING SLICES (ADAPTOR)") if ($main::verbose >= 0);

    if (scalar(@chromnames) > 0) {
	## Get selected chromosomes
	foreach $chromname (@chromnames)  {
	    &RSAT::message::Info(join("\t", "Getting slice", $slice_type, $chromname)) if ($main::verbose >= 1);
	    push @slices, $slice_adaptor->fetch_by_region($slice_type, $chromname);
	}
    } else {
	## get all the chromosomes
	@slices = @{$slice_adaptor->fetch_all($slice_type)};
    }

    my $slices_number = scalar(@slices);
    &RSAT::message::psWarn("AFTER GETTING SLICES (ADAPTOR)") if ($main::verbose >= 0);
    
    my $s=0;
    while (my $slice = shift(@slices)) {

#	&RSAT::message::Info(join("\t", "Slice size", size($slice), "Slice total size", total_size($slice))) if ($main::verbose >= 1);

	$s++;
	my $slice_id = $slice->id();
	my $slice_name = $slice->name();

	################################################################
	## TEMPORARY : a tricky fix for human genome, which contains 109 chromosome slices !
	## We discard slices which do not correspond to full chromosomes.
	## These slices have a different name (they contain _NT_)
	if ($slice_name =~ /_NT_/) {
	    &RSAT::message::Warning(join "\t", "Skipping slice (fragment of chromosome)", 
				    $s."/".$slices_number, 
				    $slice_id, $slice_name) if ($main::verbose >= 1);
	    next;
	}

	&RSAT::message::psWarn("BEFORE GETTING SLICE: ".$s."/".$slices_number) if ($main::verbose >= 0);


	&RSAT::message::TimeWarn("Collecting data for slice", $slice_type, 
				 $s."/".$slices_number, 
				 "name:".$slice->name(), 
				 "seq_region_name:".$slice->seq_region_name(), 
				 ) if ($main::verbose >= 1);

	print $CONTIG join ("\t", $slice_id, $slice->accession_number(), "<NULL>", $slice_type, $slice->length(), $slice->desc()), "\n";

	## Get all repeat regions (features)
	unless (($no_rep) || ($seq_only)) {
	    &RSAT::message::psWarn("Before collecting repeats for slice", $slice_name) if ($main::verbose >= 0);
	    &RSAT::message::TimeWarn("\tGetting repeats for slice", $slice_name) if ($main::verbose >= 1);
	    my $rep = 0;
	    my @ensembl_repeats = @{$slice->get_all_RepeatFeatures()};

	    &RSAT::message::psWarn("After collecting repeats for slice", $slice_name) if ($main::verbose >= 0);
	    my $repeats_number = scalar(@ensembl_repeats);

	    while (my $ensembl_repeat =shift(@ensembl_repeats)) {

#		&RSAT::message::Info(join("\t", "Repeat size", size($ensembl_repeat), "Repeat total size", total_size($ensembl_repeat))) if ($main::verbose >= 1);

		$rep++;
		if (($test) && ($rep > $test_number)) {
		    &RSAT::message::Info(join ("\t","TEST", $test_number, "skipping next repeats for contig", 
					       $slice_name)) if ($main::verbose >= 1);
		    last;
		}
		my $repeat_name = $ensembl_repeat->display_id();
		$repeat_id = $repeat_name.".".$rep;
		my $repeat_start = $ensembl_repeat->hstart();
		my $repeat_end = $ensembl_repeat->hend();
		my $repeat_strand = $ensembl_repeat->hstrand();
		my $repeat_type = "repeat region";
		print $REPEAT join("\t", $repeat_id, $repeat_type, $repeat_name, $slice_id, $repeat_start, $repeat_end, $repeat_strand, $repeat_name), "\n";

		&RSAT::message::Info(join("\t", ";", $slice_type, 
					  $slice->seq_region_name(),  
					  $s."/".$slices_number, 
					  "repeat", $rep."/".$repeats_number, $repeat_name))  if ($main::verbose >= 3);
	    }
	  }

	## Get all Gene objects
	unless ($seq_only) {
	    &RSAT::message::psWarn("Before collecting genes for slice", $slice_name) if ($main::verbose >= 0);
	    &RSAT::message::TimeWarn("\tGetting genes for slice", $slice_name) if ($main::verbose >= 1);
	    my $g = 0;
	    my @ensembl_genes = @{$slice->get_all_Genes()};
	    my $genes_number = scalar(@ensembl_genes);
	    while (my $ensembl_gene = shift(@ensembl_genes)) {

#		&RSAT::message::Info(join("\t", "Gene size", size($ensembl_gene), "Gene total size", total_size($ensembl_gene))) if ($main::verbose >= 1);

		$g++;

		if (($test) && ($g > $test_number)) {

		    &RSAT::message::Info(join ("\t","TEST", $test_number, "skipping next genes for contig", 
					       $slice_name)) if ($main::verbose >= 1);
		    last;
		}

		## Create a new gene object
		my $gene_name = $ensembl_gene->external_name() || $ensembl_gene->stable_id();

		warn join("\t", ";", $slice_type, $slice->seq_region_name(),  $s."/".$slices_number, "gene", 
			  $g."/".$genes_number, $gene_name), "\n" if ($main::verbose >= 3);

		my @feature = &collect_attributes($ensembl_gene);
		my $gene_id = $feature[0];
		$feature[1] = "gene";
		$feature[3] = $slice_id; ## Introduced to remove potential circular reference (test)
		my $gene_description = $feature[7];

		## Exchange start and stop coordinates if on R strand
 #                   if ($feature[6] eq "R") {
 #                       my $xchanger = $feature[4];
 #                       $feature[4] = $feature[5];
 #                       $feature[5] = $xchanger;
 #                   }

		print $GENE join("\t", @feature), "\n";
		print $GENE_NAME join ("\t", $feature[0], $feature[2], "primary"), "\n";
		unless ($feature[0] eq $feature[2]) {
		    print $GENE_NAME join ("\t", $feature[0], $feature[0], "alternate"), "\n";;
		}
#		print_DBEntries($ensembl_gene->get_all_DBLinks());

		## Get all Transcript objects for the current gene
		my $tr = 0;
		my @ensembl_transcript = @{$ensembl_gene->get_all_Transcripts()};

		while (my $trans = shift(@ensembl_transcript)) {

#		    &RSAT::message::Info(join("\t", "Transcript size", size($trans), "Transcript total size", total_size($trans))) if ($main::verbose >= 1);

		    $tr++;
		    warn join("\t", "; Collecting features for transcript", $trans), "\n" if ($main::verbose >= 5);
		    my @feature = &collect_attributes($trans);

		    warn join("\t", "; Attributes collected"), "\n" if ($main::verbose >= 5); ######

		    my $transcript_id = $feature[0];
		    push @feature, $gene_id;
		    $feature[3] = $slice_id;

		    ## Exchange start and stop coordinates if on R strand
#		    if ($feature[6] eq "R") {
#                        my $xchanger = $feature[4];
#                        $feature[4] = $feature[5];
#                        $feature[5] = $xchanger;
#                    }

		    if ($feature[7] eq "<no description>") {
			$feature[7] = $gene_description;
		    }
		    if ($feature[1] eq "protein_coding") {
			$feature[1] = "mRNA";
			print $MRNA join("\t", @feature), "\n";
			print $MRNA_NAME join("\t", $feature[0], $feature[2], "primary"), "\n";
			unless ($feature[0] eq $feature[2]) {
			    print $MRNA_NAME join ("\t", $feature[0], $feature[0], "alternate"), "\n";
			}
			unless ($gene_id eq $feature[0] || $gene_id eq $feature[2]) {
			    print $MRNA_NAME join ("\t", $feature[0], $gene_id, "alternate"), "\n";
			}
			unless ($gene_name eq $feature[0] || $gene_name eq $feature[2] || $gene_name eq $gene_id) {
			    print $MRNA_NAME join ("\t", $feature[0], $gene_name, "alternate"), "\n";
			}
		    } elsif ($feature[1] =~ /tRNA/) {
			print $TRNA join("\t", @feature), "\n";
			print $TRNA_NAME join("\t", $feature[0], $feature[2], "primary"), "\n";
			unless ($feature[0] eq $feature[2]) {
			    print $TRNA_NAME join ("\t", $feature[0], $feature[0], "alternate"), "\n";
			}
			unless ($gene_id eq $feature[0] || $gene_id eq $feature[2]) {
			    print $TRNA_NAME join ("\t", $feature[0], $gene_id, "alternate"), "\n";
			}
			unless ($gene_name eq $feature[0] || $gene_name eq $feature[2] || $gene_name eq $gene_id) {
			    print $TRNA_NAME join ("\t", $feature[0], $gene_name, "alternate"), "\n";
			}
		    } elsif ($feature[1] =~ /rRNA/) {
			print $RRNA join("\t", @feature), "\n";
			print $RRNA_NAME join("\t", $feature[0], $feature[2], "primary"), "\n";
			unless ($feature[0] eq $feature[2]) {
			    print $RRNA_NAME join ("\t", $feature[0], $feature[0], "alternate"), "\n";
			}
			unless ($gene_id eq $feature[0] || $gene_id eq $feature[2]) {
			    print $RRNA_NAME join ("\t", $feature[0], $gene_id, "alternate"), "\n";
			}
			unless ($gene_name eq $feature[0] || $gene_name eq $feature[2] || $gene_name eq $gene_id) {
			    print $RRNA_NAME join ("\t", $feature[0], $gene_name, "alternate"), "\n";
			}
		    } elsif ($feature[1] =~ /scRNA/) {
			print $SCRNA join("\t", @feature), "\n";
			print $SCRNA_NAME join("\t", $feature[0], $feature[2], "primary"), "\n";
			unless ($feature[0] eq $feature[2]) {
			    print $SCRNA_NAME join ("\t", $feature[0], $feature[0], "alternate"), "\n";
			}
			unless ($gene_id eq $feature[0] || $gene_id eq $feature[2]) {
			    print $SCRNA_NAME join ("\t", $feature[0], $gene_id, "alternate"), "\n";
			}
			unless ($gene_name eq $feature[0] || $gene_name eq $feature[2] || $gene_name eq $gene_id) {
			    print $SCRNA_NAME join ("\t", $feature[0], $gene_name, "alternate"), "\n";
			}
		    } elsif ($feature[1] =~ /misc_RNA/) {
			print $MISC_RNA join("\t", @feature), "\n";
			print $MISC_RNA_NAME join("\t", $feature[0], $feature[2], "primary"), "\n";
			unless ($feature[0] eq $feature[2]) {
			    print $MISC_RNA_NAME join ("\t", $feature[0], $feature[0], "alternate"), "\n";
			}
			unless ($gene_id eq $feature[0] || $gene_id eq $feature[2]) {
			    print $MISC_RNA_NAME join ("\t", $feature[0], $gene_id, "alternate"), "\n";
			}
			unless ($gene_name eq $feature[0] || $gene_name eq $feature[2] || $gene_name eq $gene_id) {
			    print $MISC_RNA_NAME join ("\t", $feature[0], $gene_name, "alternate"), "\n";
			}
		    }

		    warn join("\t", "; Collecting Translation"), "\n" if ($main::verbose >= 5); ######

		    ## Get CDS ID and coordinates (relative to chromosome) - there is a strand trick (see API doc)
		    ## Problem: coordinates are strange (really?)
		    my $ensembl_translation = $trans->translation();
		    my $coding_region_start;
		    my $coding_region_end;

		    warn join("\t", "; Translation collected", $ensembl_translation), "\n" if ($main::verbose >= 5); ######

		    if($ensembl_translation) {

		      $coding_region_start = $trans->coding_region_start();
		      $coding_region_end = $trans->coding_region_end();

			## Print the translated sequence for the current CDS
			&PrintNextSequence(PP,"fasta",60,$ensembl_translation->seq(),$ensembl_translation->stable_id());

			## VERIFIER SI CECI EST TOUJOURS UTILE
#			if ($feature[6] eq 'D') {
			    $feature[4] = $coding_region_start;
			    $feature[5] = $coding_region_end;
#			} else {
#			    $feature[4] = $coding_region_end;
#			    $feature[5] = $coding_region_start;
#			}
			$feature[0] = $ensembl_translation->stable_id();
			$feature[1] = "CDS";
			if ($feature[7] eq "<no description>") {
			    $feature[7] = $gene_description;
			}
			$feature[8] = $transcript_id;
			push @feature, $gene_id;
			print $CDS join("\t", @feature), "\n";
			print $CDS_NAME join("\t", $feature[0], $gene_name, "primary"), "\n";
			unless ($gene_id eq $gene_name) {
			    print $CDS_NAME join("\t", $feature[0], $gene_id, "alternate"), "\n";
			}
			unless ($feature[2] eq $gene_name || $feature[2] eq $gene_id) {
			    print $CDS_NAME join("\t", $feature[0], $feature[2], "alternate"), "\n";
			}
			unless ($transcript_id eq $feature[2]) {
			    print $CDS_NAME join("\t", $feature[0], $transcript_id, "alternate"), "\n";
			}
			unless ($feature[0] eq $feature[2]) {
			    print $CDS_NAME join ("\t", $feature[0], $feature[0], "alternate"), "\n";
			}
#		    }

		    ## Getting UTRs
		    my @utrfeature = ();
		    if ($trans->start() < $coding_region_start) {
		      if ($feature[6] eq "D") {
			push @utrfeature, "5'UTR-".$transcript_id;
			push @utrfeature, "5'UTR";
			push @utrfeature, "5'UTR-".$transcript_id;
                      } else {
			push @utrfeature, "3'UTR-".$transcript_id;
			push @utrfeature, "3'UTR";
			push @utrfeature, "3'UTR-".$transcript_id;
		      }
                      push @utrfeature, $slice_id;
		      push @utrfeature, $trans->start;
                      push @utrfeature, $coding_region_start - 1;
		      push @utrfeature, $feature[6];
                      push @utrfeature, "";
                      push @utrfeature, $transcript_id;
                      push @utrfeature, $gene_id;
		      print $UTR join ("\t", @utrfeature), "\n";
		    }
		    @utrfeature = ();
		    if ($trans->end() > $coding_region_end) {
                      if ($feature[6] eq "D") {
                        push @utrfeature, "3'UTR-".$transcript_id;
                        push @utrfeature, "3'UTR";
                        push @utrfeature, "3'UTR-".$transcript_id;
                      } else {
                        push @utrfeature, "5'UTR-".$transcript_id;
                        push @utrfeature, "5'UTR";
                        push @utrfeature, "5'UTR-".$transcript_id;
		      }
                      push @utrfeature, $slice_id;
                      push @utrfeature, $coding_region_end + 1;
		      push @utrfeature, $trans->end;
		      push @utrfeature, $feature[6];
                      push @utrfeature, "";
                      push @utrfeature, $transcript_id;
                      push @utrfeature, $gene_id;
		      print $UTR join ("\t", @utrfeature), "\n";
                    }

		    unless ($no_ons) {
			## Get all Exon objects
			warn join("\t", "; Collecting Exons"), "\n" if ($main::verbose >= 5);
			foreach my $exon (@{$trans->get_all_Exons()}) {
			    my @exonfeature = &get_exonfeature($exon);
			    $exonfeature[3] = $slice_id;
			    push @exonfeature, $transcript_id;

			    ## Coding region (not working for the moment; these methods are under development @ ensembl)
#                           push @exonfeature, $exon->coding_region_start($trans);
#                           push @exonfeature, $exon->coding_region_end($trans);
			    ## instead:
			    push @exonfeature, "";
			    push @exonfeature, "";

			    push @exonfeature, $gene_id;

			    print $EXON join("\t", @exonfeature), "\n";

			    ## Getting coding regions of exons
			    my @codingexonfeature = ();

			    warn join ("\t", "Coding region START=", $coding_region_start), "\n";
			    warn join ("\t", "Coding region END=", $coding_region_end), "\n";
			    warn join ("\t", "EXON START=", $exonfeature[4]), "\n";
			    warn join ("\t", "EXON END=", $exonfeature[5]), "\n";

			    if ($exonfeature[4] < $coding_region_start && $exonfeature[5] > $coding_region_end) {
				push @codingexonfeature, "Coding-".$exonfeature[0];
				push @codingexonfeature, "coding_exon";
				push @codingexonfeature, "Coding-".$exonfeature[0];
				push @codingexonfeature, $slice_id;
				push @codingexonfeature, $coding_region_start;
				push @codingexonfeature, $coding_region_end;
				push @codingexonfeature, $feature[6];
				push @codingexonfeature, "";
				push @codingexonfeature, $transcript_id;
				push @codingexonfeature, $gene_id;
				print $CODING_EXON join ("\t", @codingexonfeature), "\n";
			    } elsif ($exonfeature[4] >= $coding_region_start && $exonfeature[5] <= $coding_region_end) {
				push @codingexonfeature, "Coding-".$exonfeature[0];
                                push @codingexonfeature, "coding_exon";
				push @codingexonfeature, "Coding-".$exonfeature[0];
                                push @codingexonfeature, $slice_id;
				push @codingexonfeature, $exonfeature[4];
				push @codingexonfeature, $exonfeature[5];
                                push @codingexonfeature, $feature[6];
				push @codingexonfeature, "";
                                push @codingexonfeature, $transcript_id;
                                push @codingexonfeature, $gene_id;
				print $CODING_EXON join ("\t", @codingexonfeature), "\n";
                            } elsif ($exonfeature[4] < $coding_region_start && $exonfeature[5] <= $coding_region_end && $exonfeature[5] >= $coding_region_start) {
				push @codingexonfeature, "Coding-".$exonfeature[0];
                                push @codingexonfeature, "coding_exon";
				push @codingexonfeature, "Coding-".$exonfeature[0];
                                push @codingexonfeature, $slice_id;
				push @codingexonfeature, $coding_region_start;
				push @codingexonfeature, $exonfeature[5];
                                push @codingexonfeature, $feature[6];
				push @codingexonfeature, "";
                                push @codingexonfeature, $transcript_id;
                                push @codingexonfeature, $gene_id;
				print $CODING_EXON join ("\t", @codingexonfeature), "\n";
                            } elsif ($exonfeature[4] >= $coding_region_start && $exonfeature[4] <= $coding_region_end && $exonfeature[5] > $coding_region_end) {
				push @codingexonfeature, "Coding-".$exonfeature[0];
                                push @codingexonfeature, "coding_exon";
				push @codingexonfeature, "Coding-".$exonfeature[0];
                                push @codingexonfeature, $slice_id;
				push @codingexonfeature, $exonfeature[4];
				push @codingexonfeature, $coding_region_end;
                                push @codingexonfeature, $feature[6];
				push @codingexonfeature, "";
                                push @codingexonfeature, $transcript_id;
                                push @codingexonfeature, $gene_id;
				print $CODING_EXON join ("\t", @codingexonfeature), "\n";
                            }
			  }

			## Get all Intron objects
			warn join("\t", "; Collecting Introns"), "\n" if ($main::verbose >= 5);
			my $int = 0;
			foreach my $intron (@{$trans->get_all_Introns()}) {
			    $int++;
			    my @intronfeature = &get_intronfeature($intron);
			    $intronfeature[0] = "Intron".$int."-".$transcript_id;
			    $intronfeature[2] = "Intron".$int."-".$transcript_id;
			    $intronfeature[3] = $slice_id;
			    push @intronfeature, $gene_id;

#			    if ($feature[6] eq "R") {
#                                my $xchanger = $intronfeature[4];
#                                $intronfeature[4] = $intronfeature[5];
#                                $intronfeature[5] = $xchanger;
#                            }

			    print $INTRON join("\t", @intronfeature), "\n";
			  }
		      }
		    }
		  }
	      }
	    &RSAT::message::psWarn("After collecting genes for slice", $s, $slice_name) if ($main::verbose >= 0);

	    ## Get encode features
#	    my $enc_regions = $slice->get_all_MiscFeatures('encode');
#	    foreach my $enc_region (@$enc_regions) {
#	      foreach my $attr (@{$enc_region->get_all_Attributes()}) {
#		print $attr->name(), ':', $attr->value(), "\n";
#	      }
#	      print "Analysis: ",$enc_region -> analysis, "\n";
#	      print "Start: ",$enc_region -> start, "\n";
#	      print "End: ",$enc_region -> end, "\n";
#	      print "Strand: ",$enc_region -> strand, "\n";
#	    }

#	    my $misc_features = $slice->get_all_MiscFeatures();
#	    foreach my $misc_feature (@$misc_features) {
#	      foreach my $set (@{$misc_feature->get_all_MiscSets()}) {
#		print $attr->name(), ':', $attr->value(), "\n";
#	      print "Code: ", $set->code(), "\n";
#	      print "Name: ", $set->name(), "\n";
#	      print "Description: ", $set->description(), "\n";
#	      print "\n";
#	      }
#	    }


	  }

	################################################################
	## Export sequence unless otherwise specified
	my $seq_file = $slice_id;
	$seq_file =~ s/\:/_/g;
	$seq_file .= ".raw";
	my $masked_seq_file = $seq_file;
	$masked_seq_file =~ s/\.raw$/_repeat_masked.raw/;
#	if ($slice->is_circular()) {
#	    $slice_shape = 'circular';
#	} else {
#	    $slice_shape = 'linear';
#	}
#	print CTG join ("\t", $seq_file,  $slice_id, $slice_shape), "\n";
	print CTG join ("\t", $seq_file,  $slice_id), "\n";

	unless ($no_seq) {	    
	    ## Export slice sequence (unmasked)
	    &RSAT::message::TimeWarn("Saving sequence for slice", $s."/".$slices_number, 
				     $slice_type, $slice->seq_region_name(), $slice_name) if ($main::verbose >= 1);
	    open SEQ, ">".$seq_file || die "cannot open error log file".$seq_file."\n";
	    print SEQ $slice->seq();
	    close SEQ;
	    ## not sure this is useful, but to try improving garbage collection
	    &RSAT::message::psWarn("After saving sequence for slice", $s, $slice_name) if ($main::verbose >= 0);
	    
	    ## Export slice sequence (hard masked)
	    unless ($no_masked) {
		&RSAT::message::TimeWarn("Getting masked sequence for slice", $s."/".$slices_number, 
					 $slice_type, $slice->seq_region_name(), $slice_name) if ($main::verbose >= 1);
		&RSAT::message::psWarn("Before collecting repeatmasked sequence for slice", $s, $slice_name) if ($main::verbose >= 0);
		my $masked_sequence_slice = $slice->get_repeatmasked_seq();
		open MASKED_SEQ, ">".$masked_seq_file || die "cannot open error log file".$masked_seq_file."\n";
		print MASKED_SEQ $masked_sequence_slice->seq();
		close MASKED_SEQ;
		
		## not sure this is useful, but to try improving garbage collection
 		&RSAT::message::psWarn("After collecting repeatmasked sequence for slice", $s, $slice_name) if ($main::verbose >= 0);
	    }
	}
	&RSAT::message::psWarn("AFTER COLLECTING INFO FOR SLICE: ".$s."/".$slices_number) if ($main::verbose >= 0);
      }

    ## Report execution time and close log file
    my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
    print $log $exec_time if ($main::verbose >= 1); ## only report exec time if verbosity is specified
    print $log &RSAT::util::ReportExecutionTime($start_time) if ($main::verbose >= 1);

    ################################################################
    ###### Close output stream
    close $log if ($outfile{log});
    close ERR if ($outfile{err});
    close CTG if ($outfile{contigs});
    close $ORG if ($outfile{organism});
    close $ORG_NAME if ($outfile{organism_name});
    close $CONTIG if ($outfile{contig});
    unless ($no_rep || $seq_only) {
	close $REPEAT if ($outfile{repeat_region});
    }
    unless ($seq_only) {
	close $MRNA if ($outfile{mrna});
	close $SCRNA if ($outfile{scrna});
	close $TRNA if ($outfile{trna});
	close $RRNA if ($outfile{rrna});
	close $MISC_RNA if ($outfile{misc_rna});
	close $CDS if ($outfile{cds});
	close $GENE if ($outfile{gene});
	close $MRNA_NAME if ($outfile{mrna_name});
	close $SCRNA_NAME if ($outfile{scrna_name});
	close $TRNA_NAME if ($outfile{trna_name});
	close $RRNA_NAME if ($outfile{rrna_name});
	close $MISC_RNA_NAME if ($outfile{misc_rna_name});
	close $CDS_NAME if ($outfile{cds_name});
	close $GENE_NAME if ($outfile{gene_name});
#	close $XREF_TABLE if ($outfile{xreference});
	close PP;
	close $UTR;

	unless ($no_ons) {
	    close $INTRON if ($outfile{intron});
	    close $EXON if ($outfile{exon});
	    close $CODING_EXON if ($outfile{coding_exon});
	}
    }

    ################################################################
    ## Report the output directory
    
    warn "; Results stored in directory\t", $dir{output}, "\n" if ($main::verbose >= 1);
    exit(0);
}
################################################################
################### SUBROUTINE DEFINITION ######################
################################################################


################################################################
#### display full help message 
sub PrintHelp {
  open HELP, "| more";
  print HELP <<End_of_help;
NAME
	get-ensembl-genome.pl

CREATION DATE
        March 2005

AUTHORS
	Olivier Sand (oly\@bigre.ulb.ac.be)
	Jacques van Helden (jvanheld\@bigre.ulb.ac.be)
	
DESCRIPTION

	Retrieve information from EnsEMBL (http://www.ensembl.org) to
	obtain the required data for installing it in RSAT.  

REQUIREMENTS

        This script requires the BioPerl and Bio::EnsEMBL libraries.

CATEGORY
	util

USAGE
        get-ensembl-genome.pl [-i inputfile] [-o outputfile] [-v]

OPTIONS
	-h	display full help message
	-help	display options
	-v	verbose
	-outdir output directory
	-noseq  do not export the sequence (only the features)
	-noons  do not export the introns and exons (this should save memory)
	-norep  do not export the repeats
	-seqonly only exports sequences (no features (God Save the RAM))
	-test #  perform a quick test on # genes (default: $test_number)
	-chrom  import a selected chromosome only
		This option can be used iteratively to specify several
		chromosome names.
			   -chrom 21 -chrom 22 -chrom X

		Multiple chromosomes can also be entered separated by
		commas
			-chrom 21,22,X

   Connection to the EnsEMBL MYSQL server
        -org organism (example: saccharomyces_cerevisiae) ; No caps, underscore separated
	-dbname	EnsEMBL database name (example: saccharomyces_cerevisiae_core_33_1b)

   Options for the automaticaly generated SQL scripts
	-schema database schema (default: $schema)
	-host	database host (default: $host)
	-user	database user (default: $user)
	-pass   database password (default: $password)

SUBSEQUENT STEPS

	This program creates a directory where all the required
	information will be stored (sequences, features, contigs). By
	default, this directory corresponds to the dbname argument,
	and indicates the organism name and its version in EnsEMBL.
	    genus_species_core_version
	(for example EnsEMBL_homo_sapiens_core_33_35f or EnsEMBL_mus_musculus_core_33_34a)

	After this, the genome is still not installed in RSAT. For
	this, you need to use the program install-organism.

	If you want to install the new genome in RSAT, after getting
	the genome from EnsEMBL, you need to

	1) Create a directory in RSAT

    	   mkdir -p $RSAT/public_html/data/genomes/Genus_species/genome
	   
	   where Genus_species is the name of the organism. 

	2) Move the data obtained from ensembl to this directory 

	   mv genus_species_core_version/* \
	       $RSAT/public_html/data/genomes/Genus_species/genome

	   where genus_species_core_version is the output directory
	   generated by this script (by default this is the same
	   string as the dbname argument).

	3) Use install-organism to configure RSAT for this organism

	   install-organism -org Genus_species -task config -up_from XXXX

	   Where XXXX is the default upstream region start relative to the 
	   transcription start site (e.g. -800 for yeast, -2000 for mammals).
	
	4) Use install-organisms to check the start and stop codons.

	   install-organism -org Genus_species -task start_stop
	   
	   more $RSAT/public_html/data/genomes/Genus_species/genome/*_start_codon_frequencies
	   more $RSAT/public_html/data/genomes/Genus_species/genome/*_stop_codon_frequencies

	5) Calculate oligo and dyad frequencies (this takes several hours)

	   install-organism -org Genus_species \
	       -task allup,upstream_freq,oligos,dyads

	6) Clean unnecessary sequence files

	   install-organism -org Genus_species \
	       -task clean

End_of_help
  close HELP;
  exit;
}

################################################################
#### display short help message
sub PrintOptions {
  open HELP, "| more";
  print HELP <<End_short_help;
get-ensembl-genome.pl options
----------------
-h		(must be first argument) display full help message
-help		(must be first argument) display options
-i		input file
-outdir		output directory
-noseq  	do not export the sequence (only the features)
-noons  	do not export introns or exons (this should save memory)
-norep  	do not export repeats
-seqonly	only exports sequences (no features (God Save the RAM))
-chrom  	import a selected chromosome
-test #  	perform a test on # genes (default: $test_number)
-v		verbose
-org		organism (example: saccharomyces_cerevisiae) ; No caps, underscore separated
-dbname		EnsEMBL database name (example: saccharomyces_cerevisiae_core_33_1b)
-schema		database schema (default: $schema)
-host		database host (default: $host)
-user		database user (default: $user)
-pass		database password (default: $password)
End_short_help
  close HELP;
  exit;
}


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
	    &PrintOptions();
	    
	    ### output dir  
	} elsif ($ARGV[$a] eq "-outdir") {
	    $dir{output} = $ARGV[$a+1];

            ### chromosome name
        } elsif ($ARGV[$a] eq "-chrom") {
            push @chromnames,  split(",", $ARGV[$a+1]);
	    
	    ### EnsEMBL database name
	} elsif ($ARGV[$a] eq "-dbname") {
	    $dbname = $ARGV[$a+1];

            ### organism
        } elsif ($ARGV[$a] eq "-org") {
            $org = lc($ARGV[$a+1]); 

	    ################################################################
	    #### SQL database parameters for the export

	    ### schema
	} elsif ($ARGV[$a] eq "-schema") {
	    $schema = $ARGV[$a+1];

	    ### host
	} elsif ($ARGV[$a] eq "-host") {
	    $host = $ARGV[$a+1];

	    ### user
	} elsif ($ARGV[$a] eq "-user") {
	    $user = $ARGV[$a+1];
	    
	    ### password 
	} elsif ($ARGV[$a] =~ /^-pass/) {
	    $password = $ARGV[$a+1];

	    ### do not export the sequence
	} elsif ($ARGV[$a] eq "-noseq") {
	    $no_seq=1;

	    ### do not export repeats
	} elsif ($ARGV[$a] eq "-norep") {
	    $no_rep=1;

	    ### do not export exons and introns
	} elsif ($ARGV[$a] eq "-noons") {
	    $no_ons=1;

	    ### only export the sequence
	} elsif ($ARGV[$a] eq "-seqonly") {
	    $seq_only=1;
	    $no_rep=1;
	    $no_ons=1;

	    ### do not export the masked sequence
	} elsif ($ARGV[$a] eq "-nomask") {
	    $no_masked = 1;

	    ### quick test
	} elsif ($ARGV[$a] eq "-test") {
	    $test=1;
	    if (&IsInteger($ARGV[$a+1])) {
		$test_number = $ARGV[$a+1];
	    }
	}
    }
}
################################################################
#### verbose message
sub Verbose {
    print $log "; get-ensembl-genome.pl ";
    &PrintArguments($log);
    if (%dir) {
	print $log "; Directories\n";
	while (($key,$value) = each %dir) {
	    print $log ";\t$key\t$value\n";
	}
    }
    if (%main::outfile) {
	print $log "; Output files\n";
	while (($key,$value) = each %outfile) {
	    print $log ";\t$key\t$value\n";
	}
    }
}
################################################################
## Convert a feature to a string for export to the feature table
sub collect_attributes {
    my ($ensembl_object) = @_;
    my @feature = ();
   
    ## ID
    my $id = $ensembl_object->stable_id();
    push @feature, $id;

    ## Type
    my $type = $ensembl_object->biotype();
    push @feature, $type;

    ## Name(s)
    my $name = $ensembl_object->external_name();
    unless ($name) {
	$name = $id;
	print ERR join ("\t", "No name for gene", $id, "Using ID instead"), "\n";
    }
    push @feature, $name;

    ## Contig
    my $contig = $ensembl_object->slice->id(); ## =CIRCULAR REFERENCE?
    push @feature, $contig;

    ## Start position
    my $start_pos = $ensembl_object->start();
    push @feature,  $start_pos;

    ## End position
    my $end_pos = $ensembl_object->end();
    push @feature,  $end_pos;

    ## Strand
    my $strand = "D";
    unless ($ensembl_object->strand() == 1) {
	$strand =  "R";
    }
    push @feature, $strand; 

    ## Description
    my $description = $ensembl_object->description();
    unless ($description) {
	$description = "<no description>";
	print ERR join ("\t", "No description for gene", $id), "\n";
    }
    push @feature, $description; 
    
    return @feature;
}
################################################################
## Convert an exon feature to a string for export to the feature table
sub get_exonfeature {
    my ($exon) = @_;
    my @exonfeature = ();
  
    ## ID
    my $id = $exon->stable_id();
    push @exonfeature, $id;

    ## Type
    push @exonfeature, "exon";
    
    ## Exon name
    push @exonfeature, $id;
        
    ## Chromosome name.
    push @exonfeature, $exon->slice->seq_region_name();
#    push @exonfeature, "contig";

    ## Start position
    push @exonfeature, $exon->start();

    ## End position
    push @exonfeature, $exon->end();
    
    ## Strand
    my $strand = "D";
    unless ($exon->strand() == 1) {
        $strand =  "R";
    }
    push @exonfeature, $strand;

    ## Description
    push @exonfeature, "";

    return @exonfeature;
}
################################################################
## Convert an intron feature to a string for export to the feature table
sub get_intronfeature {
    my ($intron) = @_;
    my @intronfeature = ();
 
    ## ID (introns have no ID in EnsEMBL)
    push @intronfeature, "";

    ## Type
    push @intronfeature, "intron";
   
    ## Intron name
    push @intronfeature, "";

    ## Chromosome name.
    push @intronfeature, $intron->slice->seq_region_name(); ## CIRCULAR REFERENCE?
#    push @intronfeature, "contig";

    ## Start position
    push @intronfeature, $intron->start();

    ## End position
    push @intronfeature, $intron->end();
   
    ## Strand
    my $strand = "D";
    unless ($intron->strand() == 1) {
        $strand =  "R";
    }
    push @intronfeature, $strand;

    ## Description
    push @intronfeature, "";
  
    return @intronfeature;
}
################################################################
# Print cross-references  !!!il y a des doublons!!!
#sub print_DBEntries {
#    my $db_entries = shift;
#    foreach my $dbe (@$db_entries) {
#	    print $XREF_TABLE $feature[0],"\t",$dbe->dbname(),"\t",$dbe->display_id(),"\n";
#	    if ($dbe->dbname() eq 'HUGO') {
#    		print $FT_NAME_TABLE $feature[0],"\t",$dbe->display_id(),"\n";
#	    }
#    }
#}
################################################################
## Print header for the feature tables
sub PrintFtHeader {
    local ($feature_type, *filehandle_name) = @_;
    print $filehandle_name "-- dump date	", &AlphaDate(), "\n";
    print $filehandle_name "-- class	EnsEMBL::",$feature_type, "\n";
    print $filehandle_name "-- table	",$feature_type, "\n";
    print $filehandle_name "-- table	main", "\n";
    print $filehandle_name "-- field 1	id", "\n";
    print $filehandle_name "-- field 2	type", "\n";
    print $filehandle_name "-- field 3	name", "\n";
    print $filehandle_name "-- field 4	contig", "\n";
    print $filehandle_name "-- field 5	start_pos", "\n";
    print $filehandle_name "-- field 6	end_pos", "\n";
    print $filehandle_name "-- field 7	strand", "\n";
    print $filehandle_name "-- field 8	description", "\n";
    if ($feature_type eq "gene" || $feature_type eq "repeat_region") {
	print $filehandle_name "-- header", "\n";
	print $filehandle_name "-- id	type	name	contig	start_pos	end_pos	strand	description", "\n";
    } elsif ($feature_type eq "cds" || $feature_type eq "utr" || $feature_type eq "coding_exon") {
	print $filehandle_name "-- field 9	transcript", "\n";
	print $filehandle_name "-- field 10	GeneID", "\n";
	print $filehandle_name "-- header", "\n";
	print $filehandle_name "-- id	type	name	contig	start_pos	end_pos	strand	description	transcript	GeneID", "\n";
    } elsif ($feature_type eq "exon") {
	print $filehandle_name "-- field 9	coding_start", "\n";
	print $filehandle_name "-- field 10	coding_end", "\n";
	print $filehandle_name "-- field 11	GeneID", "\n";
	print $filehandle_name "-- header", "\n";
	print $filehandle_name "-- id	type	name	contig	start_pos	end_pos	strand	description	transcript	coding_start	coding_end	GeneID", "\n";
    } else {
	print $filehandle_name "-- field 9	GeneID", "\n";
	print $filehandle_name "-- header", "\n";
	print $filehandle_name "-- id	type	name	contig	start_pos	end_pos	strand	description	GeneID", "\n";
    }
}

sub PrintOrgHeader {
    print $ORG "-- dump date	", &AlphaDate(), "\n";
    print $ORG "-- class	EnsEMBL::Organism", "\n";
    print $ORG "-- table	organism", "\n";
    print $ORG "-- table	main", "\n";
    print $ORG "-- field 1	id", "\n";
    print $ORG "-- field 2	taxonomy", "\n";
    print $ORG "-- field 3	name", "\n";
    print $ORG "-- header", "\n";
    print $ORG "-- id	taxonomy	name", "\n";
}

sub PrintFtNameHeader {
    local ($feature_type, *filehandle_name) = @_;
    print $filehandle_name "-- dump date	", &AlphaDate(), "\n";
    print $filehandle_name "-- class	EnsEMBL::",$feature_type, "\n";
    print $filehandle_name "-- table	",$feature_type, "_names\n";
    print $filehandle_name "-- table	lateral", "\n";
    print $filehandle_name "-- field 1	id", "\n";
    print $filehandle_name "-- field 2	names", "\n";
    print $filehandle_name "-- field 3	status", "\n";
    print $filehandle_name "-- header", "\n";
    print $filehandle_name "-- id	names	status", "\n";
}

sub PrintContigHeader {
    print $CONTIG "-- dump date	", &AlphaDate(), "\n";
    print $CONTIG "-- class	EnsEMBL::Contig", "\n";
    print $CONTIG "-- table	contig", "\n";
    print $CONTIG "-- table	main", "\n";
    print $CONTIG "-- field 1	id", "\n";
    print $CONTIG "-- field 2	accession", "\n";
    print $CONTIG "-- field 3	version", "\n";
    print $CONTIG "-- field 4	type", "\n";
    print $CONTIG "-- field 5	length", "\n";
    print $CONTIG "-- field 6	description", "\n";
    print $CONTIG "-- header", "\n";
    print $CONTIG "-- id	accession	version	type	length	description", "\n";
}

#!/usr/bin/perl -w
############################################################
#
# $Id: get-ensembl-genome.pl,v 1.12 2005/03/17 15:07:49 jvanheld Exp $
#
# Time-stamp: <2003-07-04 12:48:55 jvanheld>
#
############################################################
#use strict;
use DBI();
BEGIN {
    if ($0 =~ /([^(\/)]+)$/) {
	push (@INC, "$`lib/");
    }
}
require "RSA.lib";
require RSAT::util;

## parsing libraries
require "RSA.lib";
push @INC, $RSA."/perl-scripts/parsers/" if ($RSA);
require "lib/load_classes.pl";
require "lib/parsing_util.pl";

## EnsEMBL libraries
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

## TO DO
## Introns and exons as transcript attributes
## Get CDS and check start/stop codons
## Test with other genomes (Anopheles)
## Check the start and stop codons of Anopheles. Half of them are false.
## Add cross-references to the RSAT objects ($rsat__gene, $rsat_transcript, $rsat_cds, ...)


################################################################
#### Class for EMBL feature
package EMBL::Feature;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "ft_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       names=>"ARRAY",
			       organism=>"SCALAR",
			       type=>"SCALAR",
			       description=>"SCALAR",
			       position=>"SCALAR",
			       contig=>"SCALAR",
			       strand=>"SCALAR",
			       start_pos=>"SCALAR",
			       end_pos=>"SCALAR",
			       xrefs=>"EXPANDED"
			      );
}


################################################################
#### Class for Gene
package EMBL::Gene;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "ft_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       names=>"ARRAY",
			      );
}

################################################################
#### Class for Transcript
package EMBL::Transcript;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "ft_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       names=>"ARRAY",
			      );
}

################################################################
#### Class for EMBL organism
package EMBL::Organism;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "ft_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       names=>"ARRAY",
			       taxonomy=>"SCALAR",
			      );
}

################################################################
#### Class for EMBL contig
package EMBL::Contig;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "ft_";
    @_objects = ();
    %_name_index = ();
    %_id_index = ();
    %_attribute_count = ();
    %_attribute_cardinality = (id=>"SCALAR",
			       accession=>"SCALAR",
			       version=>"SCALAR",
			       type=>"SCALAR",
			       length=>"SCALAR",
			       description=>"SCALAR",
			      );
}


################################################################
#### main package
package main;
{
    
    
    ################################################################
    #### initialise parameters
    my $start_time = &AlphaDate();
    
    
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
    $outfile{feature} = "prev_feature.tab";
    $outfile{xreference} = "prev_xreference.tab";
    $outfile{ft_name} = "prev_feature_name.tab";
    local $verbose = 0;
    
    ## Connection to the EnsEMBL MYSQL database
    $ensembl_host = 'ensembldb.ensembl.org';
    $ensembl_user = "anonymous";
    $dbname = 'homo_sapiens_core_28_35a';
    $org = 'schtroumpf';
    
    #### Options for the exported SQL database
    $host= "localhost";
    $schema="ensembl";
    $user="ensembl";
    $password="ensembl";
    
    ################################################################
    ## Class factories
    
    #### Organism factory
    my $organisms = classes::ClassFactory->new_class(object_type=>"EMBL::Organism",prefix=>"org_");
    $organisms->set_out_fields(qw (id
				   taxonomy
				   name
				   names
				   ));
    
    #### Contig factory
    my $contigs = classes::ClassFactory->new_class(object_type=>"EMBL::Contig",prefix=>"ctg_");
    $contigs->set_out_fields(qw (id
				 accession
				 version
				 form
				 type
				 file
				 length
				 description
				 ));
    ## Gene factory
    my $genes = classes::ClassFactory->new_class(object_type=>"EMBL::Gene",prefix=>"gene_");
    $genes->set_out_fields(qw(id
			      type
			      name
			      contig
			      start_pos
			      end_pos
			      strand
			      description
			      names
			      db_xref
			      introns
			      exons
			      ));

    ## Transcript factory
    my $transcripts = classes::ClassFactory->new_class(object_type=>"EMBL::Transcript",prefix=>"trans_");
    
    ## Feature factory
    my $features = classes::ClassFactory->new_class(object_type=>"EMBL::Feature",prefix=>"ft_");
    $features->set_out_fields(qw( id
			      type
			      name
			      contig
			      start_pos
			      end_pos
			      strand
			      description

			      names
			      db_xref
			      introns
			      exons
                             ));

    #### classes
    @classes = qw (EMBL::Organism EMBL::Contig EMBL::Gene EMBL::Transcript EMBL::Feature);
    

    ################################################################
    ## Read arguments
    &ReadArguments();

 
    ################################################################
    ## Connect to ensembldb to get list of databases and pick the one corresponding to chosen organism
    my $dbh = DBI->connect("DBI:mysql:host=$ensembl_host", "$ensembl_user", "", {'RaiseError' => 1});
    my $sth = $dbh->prepare("SHOW DATABASES");
    $sth->execute();
    while (my $ref = $sth->fetchrow_hashref()) {
        if ($ref->{'Database'} =~ /($org)_core/) {
            $dbname = $ref->{'Database'};
        }
    }
    warn "; dbname = ", $dbname, "\n" if ($main::verbose >= 1);
    $sth->finish();
    $dbh->disconnect();

   
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
    
    ## List of contif sequence files
    open CTG, ">", $outfile{contigs} || die "cannot open contig file".$outfile{contigs}."\n";; # file with contig IDs

    ## feature file
    $FT_TABLE =  &OpenOutputFile($outfile{feature});
    &PrintFtHeader();
    
    ## xref file
    open $XREF_TABLE, ">".$outfile{xreference} || die "cannot open error log file".$outfile{xreference}."\n";
    
    ## feature_name table
    open $FT_NAME_TABLE, ">".$outfile{ft_name} || die "cannot open error log file".$outfile{ft_name}."\n";
    
    
    #### print verbose
    &Verbose() if ($verbose);
    
    ## Connect to EnsEMBL database
    my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $ensembl_host, -user => $ensembl_user, -dbname => $dbname);

    ## Get Species object
    my $meta_container = $db->get_MetaContainer();
    my $tax_id = $meta_container->get_taxonomy_id();
    my $species = $meta_container->get_Species();
    my $organism_name = $species->binomial();
    my $common_name = $species->common_name();
    my @classification = $species->classification();
    warn join (":", "; Classification ", @classification),"\n" if ($main::verbose >= 1);
    warn "; Organism = ", $organism_name, ", common name = ", $common_name, ", NCBI taxid = ", $tax_id, "\n" if ($main::verbose >= 1);

    ## The method species does not return the name of the organism !
    ## TEMPORARY: I use a trick until we find the organism name in EnsEMBL API
    #if (($organism_name eq "") || ($organism_name eq "DEFAULT")) {
    #   $organism_name = ucfirst($dbname);
    #   $organism_name =~ s/_core.*//;
    #}

    my $rsat_organism = $organisms->new_object(name=>$organism_name);
    $rsat_organism->push_attribute("names", $organism_name);
    
    &RSAT::message::Info(join ("\t", "Db", $db)) if ($main::verbose >= 3);
    my $slice_adaptor = $db->get_SliceAdaptor();
    warn join ("\t", "; Adaptor", $slice_adaptor), "\n" if ($main::verbose >= 3);
    

    
    my @slices;
    if (scalar(@chromnames) > 0) {
	## Get selected chromosomes
	foreach $chromname (@chromnames)  {
	    &RSAT::message::Info(join("Getting slice", $slice_type, $chromname)) if ($main::verbose >= 1);
	    push @slices, $slice_adaptor->fetch_by_region($slice_type, $chromname);
	}
    } else {
	## get all the chromosomes
	@slices = @{$slice_adaptor->fetch_all($slice_type)};
    }
    my $s=0;
    foreach my $slice (@slices) {
	$s++;
	my $slice_id = $slice->id();
	my $slice_name = $slice->name();

	## TEMPORARY : a tricky fix for human genome, which contains 109 chromosome slices !
	if ($slice_name =~ /_NT_/) {
	    &RSAT::message::Warning(join "\t", "Skipping slice (fragment of chromosome)", $s."/".scalar(@slices), $slice_id, $slice_name) if ($main::verbose >= 1);
	    next;
	}

	&RSAT::message::TimeWarn("Slice", $slice_type, 
				 $s."/".scalar(@slices), 
				 "name:".$slice->name(), 
				 "seq_region_name:".$slice->seq_region_name(), 
				 $slice->name()) if ($main::verbose >= 2);


	my $rsat_contig = $contigs->new_object();
		
	$rsat_contig->force_attribute("id", $slice_id);
	$rsat_contig->set_attribute("name", $slice_name);
	$rsat_contig->set_attribute("type", $slice_type);
	$rsat_contig->set_attribute("accession", $slice->accession_number());
	$rsat_contig->set_attribute("length", $slice->length());
	$rsat_contig->set_attribute("description", $slice->desc());
	$rsat_contig->set_attribute("chromosome", $slice->seq_region_name());

	## Get all Gene objects
	&RSAT::message::TimeWarn("Getting all genes for slice", $slice_name) if ($main::verbose >= 1);
	my $g = 0;
	my @ensembl_genes = @{$slice->get_all_Genes()};
	foreach my $ensembl_gene (@ensembl_genes) {
	    $g++;

	    if (($test) && ($g > $test_number)) {
		warn "; TEST DONE";
		last;
	    } 

	    my $gene_name = $ensembl_gene->external_name() || $ensembl_gene->stable_id();
	    warn join("\t", ";", $slice_type, $slice->seq_region_name(),  $s."/".scalar(@slices), "gene", $g."/".scalar(@ensembl_genes), $gene_name), "\n" if ($main::verbose >= 3);

	    my $rsat_gene = $genes->new_object(source=>"ensembl", name=>$gene_name);
	    @feature = &collect_attributes($ensembl_gene, $rsat_gene);
	    $rsat_gene->force_attribute("type", "gene");
	    $rsat_gene->force_attribute("organism", $organism_name);
	    
	    print $FT_TABLE join("\t", @feature), "\n";
	    print_DBEntries($ensembl_gene->get_all_DBLinks());

	    ## Get all Transcript objects
	    my $tr = 0;
	    foreach my $trans (@{$ensembl_gene->get_all_Transcripts()}) {
		$tr++;
		warn join("\t", "transcript", $trans), "\n" if ($main::verbose >= 5);
		my $rsat_transcript = $features->new_object(source=>"ensembl");
		my @feature = &collect_attributes($trans, $rsat_transcript);
		$rsat_transcript->force_attribute("type", "mRNA");
		$rsat_transcript->force_attribute("organism", $organism_name);
		my $transcript_name = $rsat_transcript->get_attribute("name");
		$transcript_name .= ".".$tr;
		$rsat_transcript->force_attribute("name", $transcript_name);
		$rsat_transcript->push_attribute("names", $rsat_gene->get_attribute("name"));
		$rsat_transcript->push_attribute("names", $transcript_name);
		$rsat_transcript->set_attribute("gene", $rsat_gene->get_attribute("id"));
		unless($rsat_transcript->get_attribute("description")) {
		    $rsat_transcript->set_attribute("description", $rsat_gene->get_attribute("description"));
		}

		$feature[1] = "transcript";
		print $FT_TABLE join("\t", @feature), "\n";
		$transcriptID = $feature[0];

		## Get CDS ID and coordinates (relative to chromosome) - there is a strand trick (see API doc)
		## Problem: corrdinates are strange
		my $coding_region_start = $trans->coding_region_start();
		my $coding_region_end = $trans->coding_region_end();
		my $ensembl_translation = $trans->translation();
		if($ensembl_translation) {
		    my $rsat_cds = $features->new_object();
		    $rsat_cds->force_attribute("id", $ensembl_translation->stable_id());
		    $rsat_cds->set_attribute("type","CDS");
		    $rsat_cds->set_attribute("name",  $rsat_transcript->get_attribute("name"));
		    $rsat_cds->set_attribute("contig",  $rsat_transcript->get_attribute("contig"));
		    $rsat_cds->set_attribute("start_pos",  $trans->coding_region_start());
		    $rsat_cds->set_attribute("end_pos",  $trans->coding_region_end());
		    $rsat_cds->set_attribute("strand",  $rsat_transcript->get_attribute("strand"));
		    $rsat_cds->set_attribute("description", $rsat_transcript->get_attribute("description"));
		    $rsat_cds->set_attribute("organism", $rsat_transcript->get_attribute("organism"));
		    $rsat_cds->set_attribute("transcript",  $rsat_transcript->get_attribute("id"));
		    $rsat_cds->push_attribute("names", $rsat_gene->get_attribute("name"));
		    $rsat_cds->push_attribute("names", $rsat_transcript->get_attribute("name"));

		    $feature[0] =  $ensembl_translation->stable_id();
		    $feature[1] = "CDS";
		    if ($feature[6] eq 'D') {
			$feature[4] = $coding_region_start;
			$feature[5] = $coding_region_end;
		    } else {
			$feature[4] = $coding_region_end;
			$feature[5] = $coding_region_start;
		    }
		    print $FT_TABLE join ("\t", @feature), "\n"; 
		}

## Tests to make sure a transcripts includes UTRs (are also in first and last exons!)
#        print $feature[0], " : ", $trans->spliced_seq(), "\n";
#        print $feature[0], " : ", $trans->translateable_seq(), "\n";
#        my $fiv_utr = $trans->five_prime_utr();
#        my $thr_utr = $trans->three_prime_utr();
#        print $feature[0], " : ", ($fiv_utr) ? $fiv_utr->seq() : 'No 5 prime UTR', "\n";
#        print $feature[0], " : ", ($thr_utr) ? $thr_utr->seq() : 'No 3 prime UTR', "\n";

		## Get all Exon objects
		foreach my $exon (@{$trans->get_all_Exons()}) {
		    my @exonfeature = &get_exonfeature($exon);
		    print $FT_TABLE join("\t", @exonfeature), "\n";
		}

		## Get all Intron objects
		foreach my $intron (@{$trans->get_all_Introns()}) {
		    my @intronfeature = &get_intronfeature($intron);
		    $intronfeature[0] = "Trnscrpt - ".$transcriptID;
		    print $FT_TABLE join("\t", @intronfeature), "\n";
		}
	    }
	}

	################################################################
	## Export sequence unless otherwise specified
	my $seq_file = $rsat_contig->get_attribute("id");
	$seq_file =~ s/\:/_/g;
	$seq_file .= ".raw";
	my $masked_seq_file = $seq_file;
	$masked_seq_file =~ s/\.raw$/_masked.raw/;
	print CTG join ("\t", $seq_file,  $rsat_contig->get_attribute("id")), "\n";
	unless ($no_seq) {

	    ## Export slice sequence (unmasked)
	    &RSAT::message::TimeWarn("Getting sequence for slice", $s."/".scalar(@slices), 
				     $slice_type, $slice->seq_region_name(), $slice_name) if ($main::verbose >= 1);
	    open SEQ, ">".$seq_file || die "cannot open error log file".$seq_file."\n";
	    print SEQ $slice->seq();
	    close SEQ;

	    ## Export slice sequence (hard masked)
	    unless ($no_masked) {
		&RSAT::message::TimeWarn("Getting masked sequence for slice", $s."/".scalar(@slices), 
					 $slice_type, $slice->seq_region_name(), $slice_name) if ($main::verbose >= 1);
		my $masked_sequence_slice = $slice->get_repeatmasked_seq();
		open MASKED_SEQ, ">".$masked_seq_file || die "cannot open error log file".$masked_seq_file."\n";
		print MASKED_SEQ $masked_sequence_slice->seq();
		close MASKED_SEQ;
	    }
	}
    }

    ################################################################
    ### Save result in tab files
    &RSAT::message::TimeWarn("Exporting collected objects to directory",$dir{output}) if ($main::verbose >= 1);
    chdir $dir{main};
    warn "; Main directory\t", $dir{main}, "\n" if ($main::verbose >= 3);

    foreach $class_factory ($organisms, $contigs,  $features, $genes) {
#    foreach $class_factory ($organisms, $contigs, $features) {
	$class_factory->dump_tables();
	$class_factory->generate_sql(dir=>"$dir{output}/sql_scripts",
				schema=>$schema,
				host=>$host,
				user=>$user,
				password=>$password
				);
    }
    &ExportMakefile(@classes);

    chdir($dir{output});
    &PrintStats($outfile{stats}, @classes);

    


    ################################################################
    ###### finish verbose
    if ($verbose >= 1) {
	my $done_time = &AlphaDate();
	print $log "; Job started $start_time\n";
	print $log "; Job done    $done_time\n";
    }


    ################################################################
    ###### close output stream
    close $log if ($outfile{log});
    close ERR if ($outfile{err});
    close CTG if ($outfile{contigs});
    close $FT_TABLE if ($outfile{feature});
    close $XREF_TABLE if ($outfile{xreference});
    close $FT_NAME_TABLE if ($outfile{ft_name});

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
	Olivier Sand (oly\@scmbb.ulb.ac.be)
	Jacques van Helden (jvanheld\@scmbb.ulb.ac.be)
	
DESCRIPTION

	Retrieve information from EnsEMBL (http://www.ensembl.org) to
	obtain the required data for installing it in RSAT.  

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
	-test #  perform a quick test on # genes (default: $test_number)
	-chrom  import a selected chromosome only
		This option can be used iteratively to specify several
		chromosome names.
			   -chrom 21 -chrom 22 -chrom X

		Multiple chromosomes can also be entered separated by
		commas
			-chrom 21,22,X

   Connection to the EnsEMBL MYSQL server
        -org organism (default: $org) ; No caps, underscore separated
	-dbname	EnsEMBL database name (default: $dbname)

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
	(for example EnsEMBL_homo_sapiens_core_28_35a)

	After this, the genome is still not installed in RSAT. For
	this, you need to use the program install-organisms.

	If you want to install the new genome in RSAT, after getting
	the genome from EnsEMBL, you need to

	1) Create a directory in RSAT

    	   mkdir -p $RSAT/data/genomes/Genus_species/genome
	   
	   where Genus_species is the name of the organism. 

	2) Move the data obtained from ensembl to this directory 

	   mv genus_species_core_version/* \
	       $RSAT/data/genomes/Genus_species/genome

	   where genus_species_core_version is the output directory
	   generated by this script (by default this is the same
	   string as the dbname argument).

	3) Use install-organisms to configure RSAT for this organism

	   install-organism -org Genus_species -task config -up_from XXXX

	   Where XXXX is the default upstream region size (e.g. 800 for
	   yeast, 2000 for mammals).
	
	4) Use install-organisms to check the start ad stop codons.

	   install-organism -org Genus_species -task start_stop
	   
	   more $RSAT/data/genomes/Genus_species/genome/*_start_codon_frequencies
	   more $RSAT/data/genomes/Genus_species/genome/*_stop_codon_frequencies

	5) Calculate oligo and dyad frequencies (this takes several hours)

	   install-organism -org Genus_species \
	       -task allup,clean,upstream_freq,oligos,dyads

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
-chrom  	import a selected chromosome
-test #  	perform a test on # genes (default: $test_number)
-v		verbose
-org		organism (default: $org) ; No caps, underscore separate
-dbname		EnsEMBL database name (default: $dbname)
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
            $org = $ARGV[$a+1]; 

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
    if (defined(%dir)) {
	print $log "; Directories\n";
	while (($key,$value) = each %dir) {
	    print $log ";\t$key\t$value\n";
	}
    }
    if (defined(%outfile)) {
	print $log "; Output files\n";
	while (($key,$value) = each %outfile) {
	    print $log ";\t$key\t$value\n";
	}
    }
}


################################################################
## Convert a feature to a string for export to the feature table
sub collect_attributes {
    my ($ensembl_object, $rsat_object) = @_;
    my @feature = ();
   
    ## ID
    my $id = $ensembl_object->stable_id();
    push @feature, $id;
    $rsat_object->force_attribute("id", $id);

    ## Type
    my $type = $ensembl_object->type();
    push @feature, $type;
    $rsat_object->set_attribute("type", $type);

    ## Gene name
    my $name = $ensembl_object->external_name();
    unless ($name) {
	$name = $id;
	print ERR join ("\t", "No name for gene", $id, "Using ID instead"), "\n";
    }
    push @feature, $name;
    $rsat_object->push_attribute("names", $name);
    $rsat_object->force_attribute("name", $name);

    ## Contig
    my $contig = $ensembl_object->slice->id(); 
    push @feature, $contig;
    $rsat_object->set_attribute("contig", $contig);

    ## Start position
    my $start_pos = $ensembl_object->start();
    push @feature,  $start_pos;
    $rsat_object->set_attribute("start_pos", $start_pos);

    ## End position
    my $end_pos = $ensembl_object->end();
    push @feature,  $end_pos;
    $rsat_object->set_attribute("end_pos", $end_pos);

    ## Strand
    my $strand = "D";
    unless ($ensembl_object->strand() == 1) {
	$strand =  "R";
    }
    push @feature, $strand; 
    $rsat_object->set_attribute("strand", $strand);

    ## Description
    my $description = $ensembl_object->description();
    unless ($description) {
	$description = "<no description>";
	print ERR join ("\t", "No description for gene", $id), "\n";
    }
    push @feature, $description; 
    $rsat_object->set_attribute("description",$description);
    
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
    push @exonfeature, "";
        
    ## Chromosome name.
    push @exonfeature, $exon->slice->seq_region_name();

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
    push @intronfeature, $intron->slice->seq_region_name();

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
sub print_DBEntries {
    my $db_entries = shift;
    foreach my $dbe (@$db_entries) {
	    print $XREF_TABLE $feature[0],"\t",$dbe->dbname(),"\t",$dbe->display_id(),"\n";
	    if ($dbe->dbname() eq 'HUGO') {
    		print $FT_NAME_TABLE $feature[0],"\t",$dbe->display_id(),"\n";
	    }
    }
}


################################################################
## Print header for the feature table
sub PrintFtHeader {
    print $FT_TABLE "-- dump date   	", &AlphaDate(), "\n";
    print $FT_TABLE "-- class       	EnsEMBL feature", "\n";
    print $FT_TABLE "-- table       	feature", "\n";
    print $FT_TABLE "-- table       	main", "\n";
    print $FT_TABLE "-- field 1	id", "\n";
    print $FT_TABLE "-- field 2	type", "\n";
    print $FT_TABLE "-- field 3	name", "\n";
    print $FT_TABLE "-- field 4	contig", "\n";
    print $FT_TABLE "-- field 5	start_pos", "\n";
    print $FT_TABLE "-- field 6	end_pos", "\n";
    print $FT_TABLE "-- field 7	strand", "\n";
    print $FT_TABLE "-- field 8	description", "\n";
    print $FT_TABLE "-- header", "\n";
    print $FT_TABLE "-- id	type	name	contig	start_pos	end_pos	strand	description", "\n";
}

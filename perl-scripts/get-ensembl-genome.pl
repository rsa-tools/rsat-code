#!/usr/bin/perl -w

################################################################
################################################################
################################################################
################            OBSOLETE            ################
################################################################
################################################################
################################################################


############################################################
#
# $Id: get-ensembl-genome.pl,v 1.42 2011/04/10 13:49:57 jvanheld Exp $
#
# Time-stamp: <2003-07-04 12:48:55 jvanheld>
#
############################################################
#use strict;
use DBI();

## use Devel::Peek; ## This requires a specific compilation of Perl

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
			       xrefs=>"EXPANDED",
			       gene=>"SCALAR"
			       );
}

################################################################
#### Class for transcript
package EMBL::Transcript;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "trans_";
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
			       xrefs=>"EXPANDED",
			       gene=>"SCALAR"
			       );
}

################################################################
#### Class for mRNA
package EMBL::mRNA;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "mRNA_";
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
			       xrefs=>"EXPANDED",
			       gene=>"SCALAR"
			       );
}


################################################################
#### Class for scRNA
package EMBL::scRNA;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "scRNA_";
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
			       xrefs=>"EXPANDED",
			       gene=>"SCALAR"
			       );
}


################################################################
#### Class for tRNA
package EMBL::tRNA;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "tRNA_";
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
			       xrefs=>"EXPANDED",
			       gene=>"SCALAR"
			       );
}


################################################################
#### Class for rRNA
package EMBL::rRNA;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "rRNA_";
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
			       xrefs=>"EXPANDED",
			       gene=>"SCALAR"
			       );
}


################################################################
#### Class for miscellaneous RNA
package EMBL::misc_RNA;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "misc_RNA_";
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
			       xrefs=>"EXPANDED",
			       gene=>"SCALAR"
			       );
}


################################################################
#### Class for repeat regions
package EMBL::repeat_region;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "rep_";
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
#### Class for CDS
package EMBL::CDS;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "cds_";
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
			       xrefs=>"EXPANDED",
			       transcript=>"SCALAR"
			       );
}


################################################################
#### Class for sources
package EMBL::Source;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "src_";
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
#### Class for miscellaneous features
package EMBL::misc_feature;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "misc_feature_";
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
    $_prefix = "gn_";
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
#package EMBL::Transcript;
#{
#    @ISA = qw ( classes::DatabaseObject );
#    ### class attributes
#    $_count = 0;
#    $_prefix = "ft_";
#    @_objects = ();
#    %_name_index = ();
#    %_id_index = ();
#    %_attribute_count = ();
#    %_attribute_cardinality = (id=>"SCALAR",
#			       names=>"ARRAY",
#			      );
#}

################################################################
#### Class for EMBL organism
package EMBL::Organism;
{
    @ISA = qw ( classes::DatabaseObject );
    ### class attributes
    $_count = 0;
    $_prefix = "org_";
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
    $_prefix = "ctg_";
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
    $outfile{feature} = "prev_feature.tab";
    $outfile{xreference} = "prev_xreference.tab";
    $outfile{ft_name} = "prev_feature_name.tab";
    local $verbose = 0;
    
    ## Connection to the EnsEMBL MYSQL database
    $ensembl_host = 'ensembldb.ensembl.org';
    $ensembl_user = "anonymous";
    $dbname = '';
    $org = '';
    $port = '5306';
    
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
    $transcripts->set_out_fields(qw( id
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
				     gene
				     ));
    
    ## Feature factory
#    my $features = classes::ClassFactory->new_class(object_type=>"EMBL::Feature",prefix=>"ft_");
#    $features->set_out_fields(qw( id
#			      type
#			      name
#			      contig
#			      start_pos
#			      end_pos
#			      strand
#			      description
#			      names
#			      db_xref
#			      introns
#			      exons
#			      gene
#                             ));

    ## mRNA factory
    my $mRNAs = classes::ClassFactory->new_class(object_type=>"EMBL::mRNA",prefix=>"mRNA_");
    $mRNAs->set_out_fields(qw( id
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
			       gene
			       ));

    ## scRNA factory
    my $scRNAs = classes::ClassFactory->new_class(object_type=>"EMBL::scRNA",prefix=>"scRNA_");
    $scRNAs->set_out_fields(qw( id
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
				gene
				));

    ## tRNA factory
    my $tRNAs = classes::ClassFactory->new_class(object_type=>"EMBL::tRNA",prefix=>"tRNA_");
    $tRNAs->set_out_fields(qw( id
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
			       gene
			       ));

    ## rRNA factory
    my $rRNAs = classes::ClassFactory->new_class(object_type=>"EMBL::rRNA",prefix=>"rRNA_");
    $rRNAs->set_out_fields(qw( id
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
			       gene
			       ));

    ## Repeated regions factory
    my $repeat_regions = classes::ClassFactory->new_class(object_type=>"EMBL::repeat_region",prefix=>"rep_");
    $repeat_regions->set_out_fields(qw( id
					type
					name
					contig
					start_pos
					end_pos
					strand
					));

    ## Miscellaneous RNA factory
    my $misc_RNAs = classes::ClassFactory->new_class(object_type=>"EMBL::misc_RNA",prefix=>"misc_RNA_");
    $misc_RNAs->set_out_fields(qw( id
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
				   gene
				   ));

    ## CDS factory
    my $CDSs = classes::ClassFactory->new_class(object_type=>"EMBL::CDS",prefix=>"cds_");
    $CDSs->set_out_fields(qw( id
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
			      transcript
			      ));

    ## Miscellaneous features factory
    my $misc_features = classes::ClassFactory->new_class(object_type=>"EMBL::misc_feature",prefix=>"misc_feature_");
    $misc_features->set_out_fields(qw( id
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

    ## Sources factory
    my $sources = classes::ClassFactory->new_class(object_type=>"EMBL::Source",prefix=>"src_");
    $sources->set_out_fields(qw( id
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
    @classes = qw (EMBL::Organism
		   EMBL::Contig
		   EMBL::Gene
		   EMBL::mRNA
		   EMBL::scRNA
		   EMBL::tRNA
		   EMBL::rRNA
		   EMBL::misc_RNA
		   EMBL::repeat_region
		   EMBL::CDS
		   EMBL::Feature
		   EMBL::misc_feature
		   EMBL::Source);
    

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
	my $dbh = DBI->connect("DBI:mysql:host=$ensembl_host:port=$port", "$ensembl_user", "", {'RaiseError' => 1});
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
    
    ## List of contig sequence files
    open CTG, ">", $outfile{contigs} || die "cannot open contig file".$outfile{contigs}."\n"; # file with contig IDs

    ## feature file
    $FT_TABLE =  &OpenOutputFile($outfile{feature});
    &PrintFtHeader();

   ## mRNA file
#    $mRNA_TABLE =  &OpenOutputFile($outfile{mRNA});
#    &PrintFtHeader();

   ## scRNA file
#    $scRNA_TABLE =  &OpenOutputFile($outfile{scRNA});
#    &PrintFtHeader();

   ## tRNA file
#    $tRNA_TABLE =  &OpenOutputFile($outfile{tRNA});
#    &PrintFtHeader();

   ## rRNA file
#    $rRNA_TABLE =  &OpenOutputFile($outfile{rRNA});
#    &PrintFtHeader();

   ## misc_RNA file
#    $misc_RNA_TABLE =  &OpenOutputFile($outfile{misc_RNA});
#    &PrintFtHeader();

   ## repeat region file
#    $REP_TABLE =  &OpenOutputFile($outfile{rep});
#    &PrintFtHeader();

   ## Miscellaneous feature file
#    $MISC_FT_TABLE =  &OpenOutputFile($outfile{misc_feature});
#    &PrintFtHeader();

   ## Source file
#    $SRC_TABLE =  &OpenOutputFile($outfile{source});
#    &PrintFtHeader();

    ## xref file
    open $XREF_TABLE, ">".$outfile{xreference} || die "cannot open error log file".$outfile{xreference}."\n";
    
    ## feature_name table
    open $FT_NAME_TABLE, ">".$outfile{ft_name} || die "cannot open error log file".$outfile{ft_name}."\n";
    
    ## cds sequences in fasta format
    @dbsplit = split /_core_/, $dbname;
    $org = $dbsplit[0];
    $outfile{pp} = $org."_aa.fasta";
    open PP, ">", $outfile{pp} || die "cannot open sequence file".$outfile{pp}."\n";
    
    #### print verbose
    &Verbose() if ($verbose);
    
    ################################################################
    ## Open a new connection to EnsEMBL database, but this time we specify the DB name
    &RSAT::message::TimeWarn("Connecting EnsEMBL to retrieve the organism", 
			     "host=".$ensembl_host,
			     "user=".$ensembl_user,
			     "dbname=".$dbname,
			     ) if ($main::verbose >= 1);
    my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $ensembl_host, -user => $ensembl_user, -dbname => $dbname, -port => $port);

    ## Get Species object
    my $meta_container = $db->get_MetaContainer();
    my $tax_id = $meta_container->get_taxonomy_id();
    my $species = $meta_container->get_Species();
    my $organism_name = $species->binomial();
    my $common_name = $species->common_name();
    my @classification = $species->classification();
    &RSAT::message::Info (join (":", "; Classification ", @classification))  if ($main::verbose >= 1);
    &RSAT::message::Info(join ("\t", "; Organism = ", $organism_name, 
			       "common name = ", $common_name, 
			       "NCBI taxid = ", $tax_id))
	if ($main::verbose >= 1);

    ## Instantiate an object for the organism
    my $rsat_organism = $organisms->new_object(name=>$organism_name);
    $rsat_organism->push_attribute("names", $organism_name);
    &RSAT::message::Debug("Db", $db) if ($main::verbose >= 3);

    my $slice_adaptor = $db->get_SliceAdaptor();
    &RSAT::message::Debug("Adaptor", $slice_adaptor) if ($main::verbose >= 3);
    

    ## Collect the list of chromosomes
    my @slices;
#    my $slices_ref;
    if (scalar(@chromnames) > 0) {
	## Get selected chromosomes
	foreach $chromname (@chromnames)  {
	    &RSAT::message::Info(join("\t", "Getting slice", $slice_type, $chromname)) if ($main::verbose >= 1);
	    push @slices, $slice_adaptor->fetch_by_region($slice_type, $chromname);
	}
    } else {
	## get all the chromosomes
	@slices = @{$slice_adaptor->fetch_all($slice_type)};
#	$slices_ref = $slice_adaptor->fetch_all($slice_type);
    }
    
    my $s=0;

#    my $slice;
    while (my $slice = shift(@slices)) {
#    while (scalar(@$slices_ref)) {
#	$slice = shift(@$slices_ref);

#	&RSAT::message::Debug("MEMORY USAGE", &Devel::Peek::mstat()) if ($main::verbose >= 0);

	$s++;
	my $slice_id = $slice->id();
	my $slice_name = $slice->name();

	################################################################
	## TEMPORARY : a tricky fix for human genome, which contains 109 chromosome slices !
	## We discard slices which do not correspond to full chromosomes.
	## These slices have a different name (they contain _NT_)
	if ($slice_name =~ /_NT_/) {
	    &RSAT::message::Warning(join "\t", "Skipping slice (fragment of chromosome)", 
				    $s."/".scalar(@slices), 
				    $slice_id, $slice_name) if ($main::verbose >= 1);
	    next;
	}

	&RSAT::message::psWarn("BEFORE GETTING SLICE: ".$s."/".scalar(@slices)) if ($main::verbose >= 0);


	&RSAT::message::TimeWarn("Collecting data for slice", $slice_type, 
				 $s."/".scalar(@slices), 
				 "name:".$slice->name(), 
				 "seq_region_name:".$slice->seq_region_name(), 
				 ) if ($main::verbose >= 1);


	my $rsat_contig = $contigs->new_object();
	
	$rsat_contig->force_attribute("id", $slice_id);
	$rsat_contig->set_attribute("name", $slice_name);
	$rsat_contig->set_attribute("type", $slice_type);
	$rsat_contig->set_attribute("accession", $slice->accession_number());
	$rsat_contig->set_attribute("length", $slice->length());
	$rsat_contig->set_attribute("description", $slice->desc());
	$rsat_contig->set_attribute("chromosome", $slice->seq_region_name());
#	$rsat_contig->set_attribute("is_circular", $slice->is_circular());

	## Get all repeat regions (features)
	unless (($no_rep) || ($seq_only)) {

	    &RSAT::message::psWarn("Before collecting repeats for slice", $slice_name) if ($main::verbose >= 0);
	    &RSAT::message::TimeWarn("\tGetting repeats for slice", $slice_name) if ($main::verbose >= 1);
	    my $rep = 0;
	    my @ensembl_repeats = @{$slice->get_all_RepeatFeatures()};
#	    my $ensembl_repeats = $slice->get_all_RepeatFeatures();

	    &RSAT::message::psWarn("After collecting repeats for slice", $slice_name) if ($main::verbose >= 0);

	    while (my $ensembl_repeat = shift(@ensembl_repeats)){
		$rep++;
		if (($test) && ($rep > $test_number)) {
		    &RSAT::message::Info(join ("\t","TEST", $test_number, "skipping next repeats for contig", 
						   $rsat_contig->get_attribute("name"))) if ($main::verbose >= 1);
		    last;
		}
		my $repeat_name = $ensembl_repeat->display_id();
		$repeat_id = $repeat_name.".".$rep;

		&RSAT::message::Info(join("\t", ";", $slice_type, 
					  $slice->seq_region_name(),  
					  $s."/".scalar(@slices), 
					  "repeat", $rep."/".scalar(@ensembl_repeats), $repeat_name))  if ($main::verbose >= 3);

		my $rsat_repeat = $repeat_regions->new_object(source=>"ensembl", name=>$repeat_name);

		$rsat_repeat->force_attribute("id", $repeat_id);
		$rsat_repeat->force_attribute("type", "repeat_region");
		$rsat_repeat->set_attribute("contig", $slice_id);
		$rsat_repeat->set_attribute("start_pos", $ensembl_repeat->hstart);
		$rsat_repeat->set_attribute("end_pos", $ensembl_repeat->hend);
		my $ensembl_repeat_strand = $ensembl_repeat->hstrand;
		if ($ensembl_repeat_strand == 1){
		    $rsat_repeat->set_attribute("strand", "D");
		} else {
		    $rsat_repeat->set_attribute("strand", "R");
		}
	    }
	    @ensembl_repeats = undef;
	}

	## Get all Gene objects
	unless ($seq_only) {
	    &RSAT::message::psWarn("Before collecting genes for slice", $slice_name) if ($main::verbose >= 0);
	    &RSAT::message::TimeWarn("\tGetting genes for slice", $slice_name) if ($main::verbose >= 1);
	    my $g = 0;
	    my @ensembl_genes = @{$slice->get_all_Genes()};
	    while (my $ensembl_gene = shift(@ensembl_genes)) {
		$g++;
		
		if (($test) && ($g > $test_number)) {
		    &RSAT::message::Info(join ("\t","TEST", $test_number, "skipping next genes for contig", 
						   $rsat_contig->get_attribute("name"))) if ($main::verbose >= 1);
		    last;
		}
		
		## Create a new gene object
		my $gene_name = $ensembl_gene->external_name() || $ensembl_gene->stable_id();

		warn join("\t", ";", $slice_type, $slice->seq_region_name(),  $s."/".scalar(@slices), "gene", $g."/".scalar(@ensembl_genes), $gene_name), "\n" if ($main::verbose >= 3);
		
		my $rsat_gene = $genes->new_object(source=>"ensembl", name=>$gene_name);
		my @feature = &collect_attributes($ensembl_gene, $rsat_gene);
		$rsat_gene->force_attribute("type", "gene");
		$rsat_gene->force_attribute("organism", $organism_name);
		
		print $FT_TABLE join("\t", @feature), "\n";
#	    print_DBEntries($ensembl_gene->get_all_DBLinks());
		@feature = undef;

		## Get all Transcript objects for the current gene
		my $tr = 0;
		my @ensembl_transcript = @{$ensembl_gene->get_all_Transcripts()};
		while (my $trans = shift(@ensembl_transcript)) {
		    $tr++;
		    warn join("\t", "transcript", $trans), "\n" if ($main::verbose >= 5);
		    my $rsat_transcript = $transcripts->new_object(source=>"ensembl");
		    my @feature = &collect_attributes($trans, $rsat_transcript);
		    $rsat_transcript->force_attribute("organism", $organism_name);
		    my $transcript_name = $rsat_transcript->get_attribute("name");
#		    $transcript_name .= ".".$tr;
		    $rsat_transcript->force_attribute("name", $transcript_name);
		    unless (($rsat_transcript->get_attribute("id") eq $rsat_gene->get_attribute("id")) ||
			    ($rsat_transcript->get_attribute("name") eq $rsat_gene->get_attribute("id"))){
			$rsat_transcript->push_attribute("names", $rsat_gene->get_attribute("id"));
		    }
		    unless (($rsat_transcript->get_attribute("id") eq $rsat_gene->get_attribute("name")) ||
			    ($rsat_transcript->get_attribute("name") eq $rsat_gene->get_attribute("name"))){
			$rsat_transcript->push_attribute("names", $rsat_gene->get_attribute("name"));
		    }
#		    $rsat_transcript->push_attribute("names", $transcript_name);
		    $rsat_transcript->set_attribute("gene", $rsat_gene->get_attribute("id"));
		    if ($rsat_transcript->get_attribute("description") eq "<no description>") {
			$rsat_transcript->force_attribute("description", $rsat_gene->get_attribute("description"));
		    }
		    if ($rsat_transcript->get_attribute("type") eq "protein_coding"){
			my $rsat_mrna = $mRNAs->new_object(source=>"ensembl");
			$rsat_transcript->force_attribute("type", "mRNA");
			$rsat_mrna->force_attribute("id", $rsat_transcript->get_attribute("id"));
			$rsat_mrna->set_attribute("type", $rsat_transcript->get_attribute("type"));
			$rsat_mrna->set_attribute("organism", $organism_name);
			$rsat_mrna->set_attribute("name", $rsat_transcript->get_attribute("name"));
			$rsat_mrna->push_attribute("names", $rsat_transcript->get_attribute("names"));
			$rsat_mrna->set_attribute("gene", $rsat_transcript->get_attribute("gene"));
			$rsat_mrna->set_attribute("description", $rsat_transcript->get_attribute("description"));
			$rsat_mrna->set_attribute("start_pos", $rsat_transcript->get_attribute("start_pos"));
			$rsat_mrna->set_attribute("end_pos", $rsat_transcript->get_attribute("end_pos"));
			$rsat_mrna->set_attribute("strand", $rsat_transcript->get_attribute("strand"));
			$rsat_mrna->set_attribute("contig", $rsat_transcript->get_attribute("contig"));
		    }
		    if ($rsat_transcript->get_attribute("type") eq "tRNA"){
			my $rsat_trna = $tRNAs->new_object(source=>"ensembl");
			$rsat_trna->force_attribute("id", $rsat_transcript->get_attribute("id"));
			$rsat_trna->set_attribute("type", $rsat_transcript->get_attribute("type"));
			$rsat_trna->set_attribute("organism", $organism_name);
			$rsat_trna->set_attribute("name", $rsat_transcript->get_attribute("name"));
			$rsat_trna->push_attribute("names", $rsat_transcript->get_attribute("names"));
			$rsat_trna->set_attribute("gene", $rsat_transcript->get_attribute("gene"));
			$rsat_trna->set_attribute("description", $rsat_transcript->get_attribute("description"));
			$rsat_trna->set_attribute("start_pos", $rsat_transcript->get_attribute("start_pos"));
			$rsat_trna->set_attribute("end_pos", $rsat_transcript->get_attribute("end_pos"));
			$rsat_trna->set_attribute("strand", $rsat_transcript->get_attribute("strand"));
			$rsat_trna->set_attribute("contig", $rsat_transcript->get_attribute("contig"));
		    }
		    if ($rsat_transcript->get_attribute("type") eq "rRNA"){
			my $rsat_rrna = $rRNAs->new_object(source=>"ensembl");
			$rsat_rrna->force_attribute("id", $rsat_transcript->get_attribute("id"));
			$rsat_rrna->set_attribute("type", $rsat_transcript->get_attribute("type"));
			$rsat_rrna->set_attribute("organism", $organism_name);
			$rsat_rrna->set_attribute("name", $rsat_transcript->get_attribute("name"));
			$rsat_rrna->push_attribute("names", $rsat_transcript->get_attribute("names"));
			$rsat_rrna->set_attribute("gene", $rsat_transcript->get_attribute("gene"));
			$rsat_rrna->set_attribute("description", $rsat_transcript->get_attribute("description"));
			$rsat_rrna->set_attribute("start_pos", $rsat_transcript->get_attribute("start_pos"));
			$rsat_rrna->set_attribute("end_pos", $rsat_transcript->get_attribute("end_pos"));
			$rsat_rrna->set_attribute("strand", $rsat_transcript->get_attribute("strand"));
			$rsat_rrna->set_attribute("contig", $rsat_transcript->get_attribute("contig"));
		    }
		    if ($rsat_transcript->get_attribute("type") eq "scRNA"){
			my $rsat_scrna = $scRNAs->new_object(source=>"ensembl");
			$rsat_scrna->force_attribute("id", $rsat_transcript->get_attribute("id"));
			$rsat_scrna->set_attribute("type", $rsat_transcript->get_attribute("type"));
			$rsat_scrna->set_attribute("organism", $organism_name);
			$rsat_scrna->set_attribute("name", $rsat_transcript->get_attribute("name"));
			$rsat_scrna->push_attribute("names", $rsat_transcript->get_attribute("names"));
			$rsat_scrna->set_attribute("gene", $rsat_transcript->get_attribute("gene"));
			$rsat_scrna->set_attribute("description", $rsat_transcript->get_attribute("description"));
			$rsat_scrna->set_attribute("start_pos", $rsat_transcript->get_attribute("start_pos"));
			$rsat_scrna->set_attribute("end_pos", $rsat_transcript->get_attribute("end_pos"));
			$rsat_scrna->set_attribute("strand", $rsat_transcript->get_attribute("strand"));
			$rsat_scrna->set_attribute("contig", $rsat_transcript->get_attribute("contig"));
		    }
		    if ($rsat_transcript->get_attribute("type") eq "misc_RNA"){
			my $rsat_miscrna = $misc_RNAs->new_object(source=>"ensembl");
			$rsat_miscrna->force_attribute("id", $rsat_transcript->get_attribute("id"));
			$rsat_miscrna->set_attribute("type", $rsat_transcript->get_attribute("type"));
			$rsat_miscrna->set_attribute("organism", $organism_name);
			$rsat_miscrna->set_attribute("name", $rsat_transcript->get_attribute("name"));
			$rsat_miscrna->push_attribute("names", $rsat_transcript->get_attribute("names"));
			$rsat_miscrna->set_attribute("gene", $rsat_transcript->get_attribute("gene"));
			$rsat_miscrna->set_attribute("description", $rsat_transcript->get_attribute("description"));
			$rsat_miscrna->set_attribute("start_pos", $rsat_transcript->get_attribute("start_pos"));
			$rsat_miscrna->set_attribute("end_pos", $rsat_transcript->get_attribute("end_pos"));
			$rsat_miscrna->set_attribute("strand", $rsat_transcript->get_attribute("strand"));
			$rsat_miscrna->set_attribute("contig", $rsat_transcript->get_attribute("contig"));
		    }

		    $feature[1] = "transcript";
		    print $FT_TABLE join("\t", @feature), "\n";
		    $transcriptID = $feature[0];
#		    $transcriptID = $rsat_transcript->stable_id();
		    
		    ## Get CDS ID and coordinates (relative to chromosome) - there is a strand trick (see API doc)
		    ## Problem: coordinates are strange
		    my $coding_region_start = $trans->coding_region_start();
		    my $coding_region_end = $trans->coding_region_end();
		    my $ensembl_translation = $trans->translation();
		    if($ensembl_translation) {
			my $rsat_cds = $CDSs->new_object(source=>"ensembl");
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
			unless ($rsat_gene->get_attribute("name") eq $rsat_gene->get_attribute("id")){
			    $rsat_cds->push_attribute("names", $rsat_gene->get_attribute("id"));
			}
			unless (($rsat_gene->get_attribute("id") eq $rsat_transcript->get_attribute("name")) ||
				($rsat_gene->get_attribute("name") eq $rsat_transcript->get_attribute("name"))){
			    $rsat_cds->push_attribute("names", $rsat_transcript->get_attribute("name"));
			}

			## Print the translated sequence for the current CDS
#		    print PP $header, "\n";
			&PrintNextSequence(PP,"fasta",60,$ensembl_translation->seq(),$ensembl_translation->stable_id(),$rsat_gene->get_attribute("description"));

			## VERIFIER SI CECI EST TOUJOURS UTILE
			if ($feature[6] eq 'D') {
			    $feature[4] = $coding_region_start;
			    $feature[5] = $coding_region_end;
			} else {
			    $feature[4] = $coding_region_end;
			    $feature[5] = $coding_region_start;
			}

		    }

		    ## Tests to make sure a transcript includes UTRs (are also in first and last exons!)
#        print $feature[0], " : ", $trans->spliced_seq(), "\n";
#        print $feature[0], " : ", $trans->translateable_seq(), "\n";
#        my $fiv_utr = $trans->five_prime_utr();
#        my $thr_utr = $trans->three_prime_utr();
#        print $feature[0], " : ", ($fiv_utr) ? $fiv_utr->seq() : 'No 5 prime UTR', "\n";
#        print $feature[0], " : ", ($thr_utr) ? $thr_utr->seq() : 'No 3 prime UTR', "\n";


		    unless ($no_ons) {
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
		    @feature = undef;
		}
		@ensembl_transcript = undef;
	    }
	    &RSAT::message::psWarn("After collecting genes for slice", $s, $slice_name) if ($main::verbose >= 0);
	    @ensembl_genes = undef;
	    &RSAT::message::psWarn("After undefining genes for slice", $s, $slice_name) if ($main::verbose >= 0);
	}

	################################################################
	## Export sequence unless otherwise specified
	my $seq_file = $rsat_contig->get_attribute("id");
	$seq_file =~ s/\:/_/g;
	$seq_file .= ".raw";
	my $masked_seq_file = $seq_file;
	$masked_seq_file =~ s/\.raw$/_repeat_masked.raw/;
#	if ($rsat_contig->get_attribute("is_circular")) {
#	    $contig_shape = 'circular';
#	} else {
#	    $contig_shape = 'linear';
#	}
	print CTG join ("\t", $seq_file,  $rsat_contig->get_attribute("id")), "\n";
	unless ($no_seq) {

	    ## Export slice sequence (unmasked)
	    &RSAT::message::TimeWarn("Saving sequence for slice", $s."/".scalar(@slices), 
				     $slice_type, $slice->seq_region_name(), $slice_name) if ($main::verbose >= 1);
	    open SEQ, ">".$seq_file || die "cannot open error log file".$seq_file."\n";
	    print SEQ $slice->seq();
	    close SEQ;
	    ## not sure this is useful, but to try improving garbage collection
	    &RSAT::message::psWarn("After saving sequence for slice", $s, $slice_name) if ($main::verbose >= 0);

	    ## Export slice sequence (hard masked)
	    unless ($no_masked) {
		&RSAT::message::TimeWarn("Getting masked sequence for slice", $s."/".scalar(@slices), 
					 $slice_type, $slice->seq_region_name(), $slice_name) if ($main::verbose >= 1);
		&RSAT::message::psWarn("Before collecting repeatmasked sequence for slice", $s, $slice_name) if ($main::verbose >= 0);
		my $masked_sequence_slice = $slice->get_repeatmasked_seq();
		open MASKED_SEQ, ">".$masked_seq_file || die "cannot open error log file".$masked_seq_file."\n";
		print MASKED_SEQ $masked_sequence_slice->seq();
		close MASKED_SEQ;

		## not sure this is useful, but to try improving garbage collection
		&RSAT::message::psWarn("After collecting repeatmasked sequence for slice", $s, $slice_name) if ($main::verbose >= 0);
#		$masked_sequence_slice->DESTROY();
#		&RSAT::message::psWarn("After destroying repeatmasked sequence for slice", $s, $slice_name) if ($main::verbose >= 0);
	    }
	}
	&RSAT::message::psWarn("AFTER COLLECTING INFO FOR SLICE: ".$s."/".scalar(@slices)) if ($main::verbose >= 0);
	$slice = undef;
	&RSAT::message::psWarn("AFTER UNDEFINING INFO FOR SLICE: ".$s."/".scalar(@slices)) if ($main::verbose >= 0);
    }

    ################################################################
    ### Save result in tab files

    unless ($seq_only) {
	&RSAT::message::TimeWarn("Exporting collected objects to directory",$dir{output}) if ($main::verbose >= 1);
	chdir $dir{main};
	warn "; Main directory\t", $dir{main}, "\n" if ($main::verbose >= 3);
	
	foreach $class_factory ($organisms, 
				$contigs,  
#			    $features, 
				$genes, 
				$repeat_regions, 
				$CDSs, 
				$mRNAs, 
				$tRNAs, 
				$rRNAs, 
				$scRNAs, 
				$misc_RNAs) {
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
    }

    ################################################################
    ## Report execution time and close output stream
    my $exec_time = &RSAT::util::ReportExecutionTime($start_time); ## This has to be exectuted by all scripts
    print $log $exec_time if ($main::verbose >= 1); ## only report exec time if verbosity is specified
    close $log if ($outfile{log});
    close ERR if ($outfile{err});
    close CTG if ($outfile{contigs});
    close $FT_TABLE if ($outfile{feature});
    close $XREF_TABLE if ($outfile{xreference});
    close $FT_NAME_TABLE if ($outfile{ft_name});
    close PP;

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

	6) Clean unnessessary sequence files

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
    my ($ensembl_object, $rsat_object) = @_;
    my @feature = ();
   
    ## ID
    my $id = $ensembl_object->stable_id();
    push @feature, $id;
    $rsat_object->force_attribute("id", $id);

    ## Type
    my $type = $ensembl_object->biotype();
    push @feature, $type;
    $rsat_object->set_attribute("type", $type);

    ## Name(s)
    my $name = $ensembl_object->external_name();
    unless ($name) {
	$name = $id;
	print ERR join ("\t", "No name for gene", $id, "Using ID instead"), "\n";
    }
    push @feature, $name;
    $rsat_object->force_attribute("name", $name);
    $rsat_object->push_attribute("names", $name);
    if ($id ne $name) {
	$rsat_object->push_attribute("names", $id);
    }    

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

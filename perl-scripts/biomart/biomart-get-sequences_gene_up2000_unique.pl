# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
use strict;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

# Select species
my $dataset = "bdistachyon_eg_gene";

#my $confFile = "PATH TO YOUR REGISTRY FILE UNDER biomart-perl/conf/. For Biomart Central Registry navigate to http://www.biomart.org/biomart/martservice?type=registry";
my $confFile = $ENV{biomart_urlconfig};

#
# NB: change action to 'clean' if you wish to start a fresh
# configuration and to 'cached' if you want to skip configuration step
# on subsequent runs from the same registry.
#

my $action='clean';
warn("Initializing Biomart\n");
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

## Generate the query
warn("Generating query\n");
my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');
$query->setDataset($dataset);
$query->addFilter("upstream_flank", ["2000"]);
$query->addAttribute("ensembl_gene_id");
$query->addAttribute("ensembl_transcript_id");
$query->addAttribute("gene_flank");
$query->addAttribute("start_position");
$query->addAttribute("end_position");
$query->addAttribute("external_gene_id");
$query->addAttribute("ensembl_peptide_id");
$query->addAttribute("transcript_start");
$query->addAttribute("transcript_end");
$query->addAttribute("strand");
$query->addAttribute("chromosome_name");
$query->addAttribute("cds_start");
$query->addAttribute("cds_end");
$query->addAttribute("5_utr_start");
$query->addAttribute("5_utr_end");
$query->addAttribute("3_utr_start");
$query->addAttribute("3_utr_end");
$query->addAttribute("description");
$query->addAttribute("uniparc");
$query->addAttribute("uniprot_swissprot_accession");
$query->addAttribute("uniprot_sptrembl");
$query->addAttribute("external_gene_db");

$query->formatter("FASTA");

my $query_runner = BioMart::QueryRunner->new();
############################## GET COUNT ############################
# $query->count(1);
# $query_runner->execute($query);
# print $query_runner->getCount();
#####################################################################


############################## GET RESULTS ##########################
# to obtain unique rows only
$query_runner->uniqueRowsOnly(1);
$query_runner->execute($query);
my $count = $query_runner->getCount();
warn("Query returned ", $count, " results\n");

warn("Printing header\n");
my $header = $query_runner->printHeader();
print($header);

warn("Printing results\n");
my $result = $query_runner->printResults();
print($result);

warn("Printing footer\n");
my $footer = $query_runner->printFooter();
print($footer);

warn("Job done\n");
#####################################################################

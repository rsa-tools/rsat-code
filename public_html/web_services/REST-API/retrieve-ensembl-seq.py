#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/retrieve-ensembl-seq/" ##Returns upstream, downstream, intronic, exonic or UTR DNA sequences for a list of query genes.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "org" : "Saccharomyces_cerevisiae_GCF_000146045.2_R64", ##String. Query underscore between words (eg. ‘homo_sapiens’)
        "ensemblhost" : None, ##String. Address of ensembl database server (default is EBI server)
        "dbname" : None, ##String. Name of EnsEMBL database (alternative to organism).
        "dbversion" : None, ##Integer. Version of ensembl database (eg. 47
        "feattype" : None, ##String. Feature type. Supported - cds,exon,gene,intron,mrna,transcript,utr
        "type" : None, ##SString. Sequence type. Currently supported sequence types - upstream, downstream, feature.
        "utr" : None, ##String. Type(s) of UTR (untranslated region) to return. Supported - all, 5prime, 3prime
        "q" : None, ##String query. The query should be an EnsEMBL gene identifier (eg ‘ENSG00000177799’). Multiple queries can be entered, separated by ‘,’
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.
        "all" : None, ##Boolean. Return all genomic upstream regions. 
        "from" : None, ##Integer. Limits of the region to extract, relative to orf start. Default is organism dependant (example: Saccharomyces cerevisiae = -800).
        "to" : None, ##Superior limit of the region to retrieve. Default is '-1'.
        "noorf" : None, ##Boolean. The upstream/downstream sequence can only contain non-coding sequence.
        "nogene" : None, ##Boolean. the upstream/downstream sequence can only contain non-transcribed sequence.
        "maskcoding" : None, ##Boolean. all coding sequence is replaced by N in the retrieved sequence
        "labelsep" : None, ##Separator between the label fields. Default: | (pipe character).
        "rm" : None, ##Boolean. Use the repeat masked version of the genome. Attention - repeated regions are annotated for some genomes only.
        "alltranscripts" : None, ##Boolean. Get sequences for all transcripts of genes. Use -uniqseqs if you do motif discovery afterwards
        "uniqseqs" : None, ##Boolean. With -alltranscripts, returns only non-redondant sequences
        "firstintron" : None, ##Boolean. With feattype intron, get only first intron sequence
        "noncoding" : None, ##Boolean. With feattype exon, get only non-coding (part of) exons.
        "chrom" : None, ##String. Chromosome name or number (to use with -left and -right)
        "left" : None, ##Integer. Left limit of sequence to retrieve
        "right" : None, ##Integer. Right limit of sequence to retrieve
        "strang" : None, ##Integer. Strand of sequence to retrieve when using -left and -right. Supported: 1, -1
        "ftfile_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "ftfile_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping. Default value : text
        "ftfileformat" : None, ##String. Feature file format. Supported - ft, gft
        "ortho" : None, ##Boolean. Retrieve homologous sequences from EnsEMBL Compara databases
        "ortho_type" : None, ##String. Filter on homology type. (eg. ortholog, ortholog_one2one)
        "homologs_table" : None, ##String. File_name. Prints homology info to a tab delimited file
        "taxon" : None, ##String. Filter on taxonomic level (eg. Mammalia)
        "header_org" : None, ##String. Type of organism name to use in the fasta header (scientific, common or none). Default is scientific. Common name is only accessible with -ortho
        "label" : None, ##String. Information used as sequence label in the fasta header. Supported -label query
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/retrieve-seq/" ##Returns upstream, downstream or coding DNA sequences for list of query genes.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "org" : "Saccharomyces_cerevisiae_GCF_000146045.2_R64", ##String. Query organism
        "seq_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server
        "seq_string_type" : "text", ##Type of information provided by the input string (URL, piping, text)
        "informat" : None, ##Input sequence format. Supported: IG (Intelligenetics), WC (wconsensus), raw, FastA
        "prefix" : None, ##Prefix for sequence identifier
        "feattype" : None, ##Feature type. Supported - gene,mRNA,tRNA,rRNA,scRNA,misc_RNA,CDS,start_codon,stop_codon,exon
        "type" : None, ##Sequence type. Currently supported sequence types - upstream, downstream, orf, random
        "n" : None, ##Integer. Number of sequecnes (only with -type random).
        "q" : "DAL5", ##String query. The query should be an orf identifier (eg ‘metR’). The query is case-insensitive.
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.
        "ids_only" : None, ##Boolean. Use this option if the queries are provided as a list of IDs. 
        "all" : None, ##Boolean. Return all genomic upstream regions. 
        "oft" : None, ##Boolean. Output features file. Available values : true, false
        "from" : None, ##Integer. Limits of the region to extract, relative to orf start. Default is organism dependant (example: Saccharomyces cerevisiae = -800).
        "to" : None, ##Superior limit of the region to retrieve. Default is '-1'.
        "format" : None, ##Allows to select different output formats. Available values : IG, wconsensus, multi, fasta
        "lw" : None, ##Line width, A newline character will be inserted in the sequence every lw bases (0 for whole sequence on one line).
        "label" : None, ##Field(s) to be used in the sequence label. Multiple fields can be specified, separated by commas. Supported: id, name, organism_name, sequence_type, current_from, current_to, ctg, orf_strand, reg_left, reg_right. Default: name. 
        "labelsep" : None, ##Separator between the label fields. Default: | (pipe character).
        "noorf" : None, ##Boolean. Prevent overlap with neighbout genes. Available values : true, false
        "rm" : None, ##Boolean. Use the repeat masked version of the genome. Available values : true, false
        "nocom" : None, ##Prevents warning when a gene cannot be identified. By default, the comment indicates the ORF and upstream sequence coordinates.
        "imp_pos" : None, ##Admit imprecise positions. 
        "nowarn" : None, ##Prevents warning when a gene cannot be identified. Available values : true, false
        "randsels" : None, ##Integer. Select a random set of N genes in the genome annotations.
        "if_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "lf_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping. Default value : text
        "features_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "features_string_type" : "text" ##Type of information provided by the input string. Available values : text, url, piping. Default value : text
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
##print(r.text)
# print (repr(r.json))
 
name_of_file = input("What is the name of the file: ")
completeName = name_of_file + ".html"
f = open(completeName, "w+")
f.write(r.text)
f.close()
folder = '/Users/VanVourdalak/Desktop/RSAT/REST-API/APIsResults'
savedfile = os.path.join(folder, completeName)
(savedfile)

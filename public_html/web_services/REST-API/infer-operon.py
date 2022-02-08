#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/infer-operon/" ##Given a list of input genes, infer the operon to which each of these genes belong.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.        
        "format" : "id", ##String. Input format. Supported formats = id, varBed, bed  
        "org" : "Saccharomyces_cerevisiae_GCF_000146045.2_R64", ##String. Query organism
        "all" : None, ##Boolean. Infer operons for all the genes of the query organism.
        "q" : None, ##String. Query gene. This option can be used iteratively on the same command line to specify several query genes.
        "dist" : None, ##Integer. Distance threshold.
        "sep" : None, ##String. Specify the separator for multi-value fields (e.g. genes) in the output table.
        "min_gene_nb" : None ##Integer. Specify a threshold on the number of genes in the operon.
        "return" : None, ##String. List of fields to return. Separated by ',’. Supported - leader,trailer,operon,query,name,upstr_dist,q_info,up_info,down_info
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
print(r.text)
# print (repr(r.json))
 

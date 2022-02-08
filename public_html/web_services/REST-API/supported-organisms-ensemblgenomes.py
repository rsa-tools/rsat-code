#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/supported-organisms-ensemblgenomes/" ##Get the list of organisms supported at the Ensembl Genomes database, using the Perl API.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "q" : "", ##String. Query string. Multiple query strings can be specified, separated by ‘,’
        "query_type" : None, ##String. Can be - branch, name, species_taxid
        "query_file_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "query_file_string_type" : None, ##Type of information provided by the input string. Available values : text, url, piping.
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
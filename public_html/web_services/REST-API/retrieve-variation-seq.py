#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/retrieve-variation-seq/" ##Given a set of IDs for polymorphic variations, retrieve the corresponding variants

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##String describing the matrix collection. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows). Example = http://rsat.sb-roscoff.fr/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.tf
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "org" : None, ##String. Organism name.
        "format" : None, ##String. Variation format. Format of the input file
        "mml" : None, ##Integer. Length of the longest Matrix
        "col" : None, ##String. Column containing the variation IDs with the input format ‘id’
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

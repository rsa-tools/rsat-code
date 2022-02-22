#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/retrieve-matrix/" ##Retrieve a subset of matrices from a collection of position-specific scoring matrices.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat.sb-roscoff.fr/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.tf", ##String describing the matrix collection. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows). Example = http://rsat.sb-roscoff.fr/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.tf
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "v" : 1, ##Integer. Verbosity.
        "id" : "MA0265.1", ##Query IDs or ACs. The search is case-insensitive. Several queries can be provided, separated by commas without spaces, e.g. FOXA1,pou2f2.
        "id_file_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "id_file_string_type" : "text" ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content.
    }
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()


print(r.text)
# print (repr(r.json))

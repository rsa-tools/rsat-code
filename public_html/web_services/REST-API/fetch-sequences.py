#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/fetch-sequences/" ##Retrieve genome sequences for a set of coordinates specified in a bed file.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "genome" : "mm9", ##String. Genome version (e g mm9, hg19)
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "header_format" : "galaxy", ##Format for sequence headers. Supported header formats. UCSC, galaxy
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))


#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/sequence-lengths/" ##Return the length of each sequence for a user-speficied sequence file. Optionally, return the sum of lengths.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows)
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "unit" : None, ##String. Units for sequence lengths. Supported - bp,kb,mb,gb
        "in_format" : None, ##Input format. The input file can contain either sequences or genomic coordinates (-in_format bed).
        "sum" : None, ##Boolean. Only return sum of sequene lengths
    }
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()


print(r.text)
# print (repr(r.json))

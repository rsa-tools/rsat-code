#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/variation-info/" ##Taking as input a set of either variation IDs (rs numbers) or genomic regions (bed format), retrieve information about matching variants. The collected information is reported in varBed format, which is a specific type of bed format for variations (see convert-variations for a detailed description of varBed and for interconversions between variation formats).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.
        "format" : "id", ##String. Input format. Supported formats = id, varBed, bed
        "org" : "Saccharomyces_cerevisiae_GCF_000146045.2_R64", ##String. Query organism
        "type" : None, ##String. Specify one or several accepted types of variation. Supported types = SNV, insertion, deletion, substitution.
        "col" : None, ##Integer. Number of the column containing the variation IDs with the input format “id”
        "release" : None, ##Integer. The Ensembl release number of database (e.g. 84)
        "skip" : None, ##Integer. Skip the N first variations of the input file. This option is useful to run quick tests, or to split an analysis in separate tasks.
        "last" : None ##Integer. Stop after the N first variations of the list. This option is useful to run0 quick tests, or to split an analysis in separate tasks.
    }
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()


print(r.text)
# print (repr(r.json))

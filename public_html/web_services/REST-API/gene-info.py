#Begin by importing the Requests module
import requests, sys ## sys module provides information about constants, functions and methods of the Python interpreter.
import json
import os.path

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/gene-info/" ##Returns the information about genes (CDS, mRNA, …) specified either by their identifier, name, or by any supported synonym.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "org" : "Saccharomyces_cerevisiae_GCF_000146045.2_R64", ##String. Query organism
        "q" : "ARG3", ##String query. The query should be an orf identifier (eg ‘metR’). The query is case-insensitive.
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.
        "full" : None, ##Boolean. Full match only (no substring matching)
        "noquery" : None, ##Boolean. Do not print the query at the beginning of each line
        "descr" : None, ##Boolean. Match queries against the description (in addition to gene ID and names)
        "feattype" : None ##Boolean. Feature type. Supported - gene,mRNA,tRNA,rRNA,scRNA,misc_RNA,CDS,start_codon,stop_codon,exon
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()

 
#print(r.text)  
#print(repr(r.json)) 

name_of_file = input("What is the name of the file: ")
completeName = name_of_file + ".html"
f = open(completeName, "w+")
f.write(repr(r.json))
f.close()


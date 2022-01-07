#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/add-gene-info/" ##Takes as input a tab-delimited file with one ore more columns containing gene IDs, and adds columns with information about the corresponding genes.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "col" : None, ##Integer. Column containing gene IDs.
        "org" : None, ##String. Organism
        "info" : None, ##String. Information type (supported - id,ctg,strand,left,right,name,descr,names,upstr_neighb_name,upstr_neighb_id,upstr_limit,upstr_size,downstr_neighb_name,downstr_neighb_id,downstr_limit,downstr_size,right_neighb_name,right_neighb_id,right_limit,right_size,left_neighb_name,left_neighb_id,left_limit,left_size)
        "before" : None, ##Boolean. Add the information before the input line (by default, the info is added at the end of each input line).
        "null" : None, ##Boolean. String to display for undefined values (default ).
        "feattype" : None, ##String. Feature type. Supported - gene,mRNA,tRNA,rRNA,scRNA,misc_RNA,CDS,start_codon,stop_codon,exon
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

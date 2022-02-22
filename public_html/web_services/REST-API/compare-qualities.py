#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/compare-qualities/" ##Compare the empirical distributioins of weight score obtained with different matrices in a given sequence set.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "quality_list_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "quality_list_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.        
        "cluster" : None, ##Integer. Cluster matrix-quality results
        "img_format" : None, ##String. Image format for the plot comparin CV_FPR_50% / Matrix_sites_FPR_50%. To display the supported formats, type the following command XYgraph -h.
    }
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 
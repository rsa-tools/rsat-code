#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/variation-scan/" ##Scan variation sequences with position-specific scroring matrices (PSSM)

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.        
        "org" : "", ##Organism name for background model. Use only if bg is not specified
        "m_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "m_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.        
        "top_matrices" : None, ##Integer. Only work with a given number of top matrices. This option is useful for quick tests and debugging.
        "top_variations" : None, ##Integer. Only work with a given number of top matrices. This option is useful for quick tests and debugging.
        "m_format" : None, ##String. Matrix file format. Supported = transfac or tab.
        "markov_order" : None, ##Integer. Alternative to the bg option. Automatically choose a background file on the server based on the species name and assembly.
        "bg_string" : None, ##String. Input string specifying the query.
        "bg_string_type" : None, ##String. Available values : text, url, piping
        "uth_pval" : None, ##String. Upper threshold for p value
        "lth_pval_ratio" : None, ##Integer. Lower threshold for p value ratio
        "lth_score" : None, ##String. Weight of predicted sites
        "lth_w_diff" : None, ##String. Weight difference between variants
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

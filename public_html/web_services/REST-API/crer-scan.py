#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/crer-scan/" ##IScanning of predicted sites on sequence. And Detection of putative cis-regulatory enriched regions (CRERs).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "in_format" : None, ##String. Available values : ft, bed
        "autoparam" : None, ##Boolean. Extract some input parameters from the commented rows (starting with ‘;’) of the input file.
        "s" : None, ##Boolean. Sort the list of sites.
        "return_limits" : None, ##Boolean. Return every limits of sequences
        "return_limits_filtered" : None, ##Boolean. Return the limits filtered of the sequence
        "uth_site_pval" : None, ##String. Maximal p-value of sites to be considered.
        "number_of matrix" : None, ##Integer. number of matrix used for the discovery of transcription factor binding sites.
        "lth_score" : None, ##Number. Minimal site score to be considered
        "uth_score" : None, ##Number. Maximal site score to be considered
        "lth_crer_size" : None, ##Number. Minimal size of the enriched region (in bp).
        "uth_crer_size" : None, ##String. Maximal size of the enriched region (in bp).
        "lth_crer_sites" : None, ##Number. Minimal number of sites covered by the enriched region
        "uth_crer_sites" : None, ##Number. Maximal number of sites covered by the enriched region
        "lth_crer_sites_distance" : None, ##Number. Distance between successive sites to be considered.
        "uth_crer_sites_distance" : None, ##Number. Distance between successive sites to be considered.
        "uth_crer_pval" : None, ##String. Maximal binomial p-value
        "uth_crer_eval" : None, ##String. Maximal e-value
        "lth_crer_sig" : None, ##Number. Minimal binomial significance
        "uth_overlap" : None, ##Integer. Maximal overlap to define two distinct sites
        "nopval" : None, ##Boolean. Compute crer without p value
        "pre_table" : None, ##Boolean. Compute a table where is all possible p_value
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/compare-features/" ##Compare two or more sets of features. This program takes as input several feature files (two or more), and calculates the intersection, union and difference between features. It also computes contingency tables and comparison statistics.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "filelist_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "filelist_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "ref_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "ref_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "iformat" : None, ##String. Input feature format. bed, dnapat, ft, galaxy_seq, gft, gff, gff3bed, swembl, ucsc_seq
        "self" : None, ##Boolean. Also perform comparison between features in the same file (self-comparison). This can be useful to detect redundancy between annotated features.
        "return" : "stats,inter", ##String. Specify the output type(s). Supported output types. stats,inter,diff.
        "lth_inter_len" : None, ##Number. Minimal overlap (bp). Lower threshold on Length (in residues) of the intersection between two features.
        "lth_inter_cov" : None, ##Number. Intersection coverage (0-1). Lower threshold on Coverage of the intersection between two features.
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

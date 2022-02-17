#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/pattern-assembly/" ##Assemble a set of oligonucleotides or dyads into groups of overlapping patterns (assemblies).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##nput string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows). .pat files
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "2str" : None, ##Boolean. Oligonucleotide occurrences found on both stands are summed. Available values : true, false
        "1str" : None, ##Boolean. Inactivates the summation of occurrences on both strands. Available values : true, false
        "sc" : 2, ##Integer. Score column
        "cc" : None, ##Integer. Define a column containing cluster names or numbers.
        "maxfl" : None, ##Integer. Maximum flanking segment size (default 1).
        "subst" : None, ##Integer. Maximum allowed substitutions (default 0).
        "maxpat" : None, ##Integer. Maximum number of allowed patterns (default 0).
        "toppat" : 100, ##Integer. Maximum number of patterns to assemble.
        "match" : None, ##Integer. Minimum number of matching residues to include a pattern in an assembly (default 0).
        "weight" : None, ##Number. Minimum matching weight to include a pattern in an assembly (default 0)
        "max_asmb_nb" : None, ##Integer. Maximal number of assemblies (default 5)
        "max_asmb_per_cluster" : None, ##Integer. Maximal number of assemblies per cluster (default 2).
        "max_asmb_size" : None, ##Integer. Maximal assembly size, i.e. the number of patterns per alignment group (default 50)
        "max_asmb_width" : None, ##Integer. Maximal width for an assembly (default 0)
        "single_sep" : None, ##Boolean. Report the isolated words (i.e. words that do not match any other words) separately.
    }
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()


print(r.text)
# print (repr(r.json))

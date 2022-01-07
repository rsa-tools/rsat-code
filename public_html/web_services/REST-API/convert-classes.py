#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/convert-classes/" ##Interconversions between different formats of class/cluster files.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "names_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "names_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "from" : None, ##String. Input format. Supported - tab, mcl, profiles, mcode, rnsc
        "to" : None, ##String. Output format. Supported - tab, mcl,profiles
        "mcol" : None, ##Integer. Member column. Column containing the member names in the tab format (default 1).
        "ccol" : None, ##Integer. Class column. Column containing the class names in the tab format (default 1).
        "scol" : None, ##String. Score column. Column containing the scores in tab format. If not specified, scores are not defined.
        "null" : None, ##Boolean. Null string used as score in the profile output for the undefined class memberships (default 0).
        "ing" : None, ##String. Value to display as replacement for the infinite values (obtained e.g. from log(0)).
        "all_scores" : None, ##Boolean. Assign each node to all the clusters with the ERMG format, with the posterior probability as score.
        "min_scores" : None ##Integer. Minimal score value for member to class assignation.
} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/count-words/" ##Calculates oligomer frequencies from a set of sequences.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.        
        "l" : None, ##Integer. Set oligomer length to l (monad size when using dyads)
        "2str" : None, ##Boolean. Add reverse complement.
        "1str" : None, ##Boolean. Do not add reverse complement.
        "sp" : None, ##String. Spacing between the two parts of the dyads from N to M. Ex.0-20
        "noov" : None, ##Boolean. Do not allow overlapping occurrences.
        "grouprc" : None, ##Boolean. Group reverse complement with the direct sequence.
        "nogrouprc" : None, ##Boolean. Do not group reverse complement with the direct sequence.
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 
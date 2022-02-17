#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/info-gibbs/" ##Gibbs sampling algorithm for motifs discovery. Searches for highly conserved motifs in a set of DNA sequences.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "w" : 20, ##Integer. Set the motif width to w
        "maxspacing" : None, ##Integer. Set maximal spacing between motif monad to maxspacing (only for dyadic motif).
        "minspacing" : None, ##Integer. Set minimal spacing between motif monad to minspacing (only for dyadic motif).
        "strand" : None, ##String. Search in foward strand + or in both strands ±
        "n" : 1000, ##Integer. Maximum number of Gibbs sampling iterations.
        "sites" : None, ##Integer. Number of motif occurrences that are expected to be found (incompatible with -e)
        "e" : None, ##Number. Mean number of motif occurrences (sites) expected per sequence that are expected to be found (incompatible with --sites)
        "zoops" : None, ##Boolean. Try to find 0 or 1 site per sequence
        "m" : None, ##Integer. Number of motifs to extract (one by default)
        "b" : None, ##Integer. Use b predefined INCLUSive background model
        "d" : None, ##Integer. Set minimal distance between 2 motif occurrences to d
        "t" : None, ##Number. Set the temperature (should be in range [0.6 1.4])
        "r" : 5, ##Integer. Try to run the Gibbs sampling seach r times
        "collect" : None, ##Integer. Try to collect the N best sites using their weight scores
        "seedmatrix" : None ##String. Start sampling form sites collected by scanning the sequences with matrix seedmatrix
    }
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()


print(r.text)
# print (repr(r.json))

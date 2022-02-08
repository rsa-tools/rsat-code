#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/matrix-distrib/" ##Computes the theoretical distribution of score probabilities of a given PSSM. Score probabilities can be computed according to bernoulli as well as markov-chain background models

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "m_string" : "", ##String. Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "m_string_type" : "text", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "top" : None, ##Integer. Top_matrix_nb. Restrict the analysis to the N top matrices of the input file.
        "mlist_string" : "http://rsat-tagc.univ-mrs.fr/rsat/motif_databases/cisBP2/cisBP_Saccharomyces_cerevisiae_2019.tf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "mlist_string_type" : "url", ##String. Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_format" : None, ##String. Matrix format. Default is tab.
        "pseudo" : 1, ##Integer. Pseudo-count for the matrix (default 1).
        "org" : "Saccharomyces_cerevisiae", ##String. Organism for background model file
        "markov_order" : None ##Integer. Markov order for the background model file
        "bgfile_string" : None, ##String. Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bgfile_string_type" : "text", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "bg_format" : None, ##String. Background_format. Supported formats - all the input formats supported by convert-background-model.
        "bg_pseudo" : 0.01, ##Number. Pseudo frequency for the background models. Value must be a real between 0 and 1 (default 0)
        "decimals" : None, ##Integer. Number of decimals to print or the transition probabilities.
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
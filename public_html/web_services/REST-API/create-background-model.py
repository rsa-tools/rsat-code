#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/create-background-model/" ##Create a background model by computing oligonucleotide frequencies in a user-specified sequence file (he “background” file).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat//demo_files/ChIP-seq_peaks/Oct4_peaks_top1000.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "out_format" : "transitions", ##String. Output format. Supported transitions, tab, oligo-analysis, oligos, meme, motifsampler, ms, inclusive, patser, tables.
        "markov" : 2 ##Integer. Markov order.
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
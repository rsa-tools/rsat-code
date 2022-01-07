#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/convert-features/" ##Interconversions between various formats of feature description.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "from" : "dnapat", ##String. Input format. Supported - galaxy_seq,ucsc2_seq,bed,gff,gff3,swembl,getfasta_seq,dnapat,gtf,ft,bed3col,gft
        "to" : "ft", ##String. output format. Supported - fasta,bed3col,gft,ft,dnapat,great,gff3,bed,gff.
        "add_chrm" : None, ##Boolean. only for BED output - add “chr” in front of the chromosome name and change MT into chrM
        "remove_chrm" : None, ##Boolean. Remove the “chr” in front of the chromosone and change chrM into MT
        "yeast_to_roman" : None, ##Boolean. Convert arabic to roman numbers to denote chromosome numbers according to S.cerevisiae specifications. Only works for chromosomes I to XVI.
        "featname" : None, ##Boolean. Set a name for all features of the file. This option can be convenient for conversions to bed files.
        "summits" : None, ##Boolean. only valid for SWEMBL input - replace start and end coordinates by peak summit position
        "extend" : None, ##Integer. Extend peak coordinates on both sides (start and end).
        "extend_start" : None, ##Integer. Extend start coordinate leftwise
        "extend_end" : None, ##Integer. Extend end coordinate rightwise.
        "coord" : None, ##Boolean. bedfile with absolute coordinate of the sequence relative to which the features were defined (e.g. features from promoter-wise to genome-wise coordinates).
        "origin" : None, ##String. Origin of coordinates relative to sequence fragment. This option is only valid when combined with the option coord. Supported: start, end, center
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

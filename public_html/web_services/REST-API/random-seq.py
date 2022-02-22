#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/random-seq/" ##Generate random DNA or protein sequences according to various probabilistic models (independently distributed residues or Markov models).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "l" : 1000, ##Integer. Sequence length
        "n" : 10, ##Integer. Number of sequences
        "lw" : 50, ##Integer. Line width. A newline character will be inserted in the sequence every lw base. Default is 70. A value of 0 will prevent newline insertion.
        "type" : None, ##String. Type of sequence(s) to generate. Protein, DNA, or other
        "seed" : None ##Integer. Number of seed for the random generator
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
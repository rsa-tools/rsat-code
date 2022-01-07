#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/IUPAC-to-regular/" ##Convert a pattern described with the IUPAC code for ambiguous nucleotides into an equivalent regular expression. This expression can be used to search for complex patterns with general string search program like grep or gais.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "q" : None, ##Query String. 
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

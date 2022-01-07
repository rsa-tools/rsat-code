#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/classfreq/" ##Computes frequency distribution of numerical values found in a column of a tab-delimited text file. Class intervals can be specified as a fixed value (-ci) or automatically derived from the data.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.        
        "max" : None, ##Integer. Numbers strictly greater than max are not taken into account.
        "min" : None, ##Integer. Numbers strictly smaller than min are not taken into account. min also determines the lower limit of the first class.
        "ci" : None, ##Integer. Class interval. If not specified, takes the value (max - min)/20 so that 21 classes are specified.
        "col" : None, ##Integer. Column to which apply the program. This option can be used iteratively.
        "from" : None, ##Integer. Inferior limit for the classes to display values lower than this limit are however taken into account in the calculation of statistics (mean, variance, …) and of class frequencies (In contrast with the -min option).
        "to" : None, ##Integer. Superior limit for the classes to display values higher than this limit are however taken into account in the calculation of statistics (mean, variance, …) and of class frequencies (In contrast with the -max option).
        "thr" : None, ##Integer. Threshold. Only display classes with absolute frequency higher than or equal to the threshold. This option is useful for sparse data, where many classes do not contain any observation (-thr 1).
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

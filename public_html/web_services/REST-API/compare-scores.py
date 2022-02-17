#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/compare-scores/" ##Compare the score associated to keys in different input files (basically, this amounts to join different tables on the basis of a unique identifier).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "i_2_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_2_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "files" : None, ##String. List of files specified on the command line.
        "filelilst_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "filelist_string_type" : "text", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "format" : None, ##String. Export format (default profiles) Supported formats - classes,profiles
        "sc" : None, ##Integer. Score column.
        "ic" : None, ##Integer. Identifier column (default 1)
        "header" : None, ##Boolean. Use the first line of each input file as column headers.
        "lc" : None, ##Boolean. By default, the comparison is case-insensitive, but the ID case is maintained in the output. This can however b modified with the options '-lc’ (IDs converted to lowercases) and '-uc’ (IDs converted to uppercases).
        "uc" : None, ##Boolean. See -lc.
        "null" : None, ##Boolean. Null string (default ) displayed when one file contains no value for a given key)
        "numeric" : None, ##Boolean. Sort IDs numerically rather than alphabetically
        "decreasing" : None, ##Boolean. Sort IDs numerically in a decreasing order
        "basename" : None, ##Boolean. Remove path (directory) from file names in the header
        "suppress" : None, ##String. Suppress a given substring from file names in the header. Ex. -suppress ‘.tab’ -suppress ‘oligos_’
        "subst" : None, ##String. Substitute a given substring from file names in the header by a specified substring. Ex. -subst ‘oligo_’ ‘ol’
} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
#print (repr(r.json))
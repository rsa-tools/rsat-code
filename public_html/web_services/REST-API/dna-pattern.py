#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/dna-pattern/" ##Searches all occurrences of a pattern within DNA sequences. The pattern can be entered as a simple nucleotide sequence, but can also include degenerate nucleotide codes, or regular expressions.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.        
        "mask" : None, ##String. Mask lower or uppercases, respecively, i.e. replace selected case by N characters. Supported: upper, lower
        "format" : None, ##String. Input sequence format. The accepted formats are fasta, IG, raw, multi, filelist
        "pl_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "pl_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.        
        "subst" : None, ##Integer. Allow subst substitutions.
        "noIUPAC" : None, ##Boolean. The pattern is considered as a standard regular expression.
        "sc" : None, ##Integer. Score column
        "noid" : None, ##Boolean. Do not search pattern identifier in the second column of pattern file.
        "noov" : None, ##Boolean. Do not count overlapping matches for self-overlapping patterns.
        "2str" : None, ##Boolean. Search matches on both strands (direct and reverse complement)
        "1str" : None, ##Boolean. Search matches only on the direct strand.
        "R" : None, ##Boolean. Search matches only on the reverse complement strand
        "id" : None, ##Boolean. Pattern identifier (one word).
        "return" : None, ##String. List of fields to return. Multiple fields can be entered separated by commas, or by using iteratively the option.
        "match_format" : None, ##String. Format for returning matches (supported - fasta, table)
        "th" : None, ##Integer. Threshold. Return match count only for sequences with greater than or equal number of matches
        "merge" : None, ##Boolean. Merge mutually overlapping matches.
        "N" : None, ##Integer. Return matching sequences with N flanking nucleotides.
        "NL" : None, ##Integer. Return matching sequences with NL left flanking nucleotides.
        "NR" : None, ##Integer. Return matching sequences with NR right flanking nucleotides.
        "origin" : None, ##Integer. Define origin as the origin for the calculation of positions.
        "window" : None, ##Integer. Sliding window size.
        "top" : None, ##Boolean. (With sliding window only). only return the top score obtained with the sliding window for each sequence.
        "sort" : None, ##Boolean. (With -top only). sort sequences according to their top score
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 
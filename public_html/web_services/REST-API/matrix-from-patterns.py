#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/matrix-from-patterns/" ##Build PSSMs from a sequence set using as seeds a set of patterns (oligos, dyads) or an assembly.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "seq_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "seq_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "format" : None, ##String. Sequence format
        "pl_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "pl_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "toppat" : None, ##Integer. Max number of patterns to assemble. This argument is passed to pattern-assembly.
        "max_asmb_nb" : None, ##Integer. This parameter is passed to pattern-assembly, to indicate the maximal number of assemblies to return.
        "top_seq" : None, ##Integer. Max number of sequences to scan for building final matrices.
        "sc" : None, ##Integer. Column containing the pattern scores in the input pattern file. This argument is passed to pattern-assembly.
        "cc" : None, ##Integer. Column indicating the pattern clusters in the input pattern file (default 1). This argument is passed to pattern-assembly.
        "max_asmb_per_cluster" : None, ##Integer. This parameter is passed to pattern-assembly, to indicate the maximal number of assemblies to return per cluster.
        "subst" : None, ##Integer. Maximum number of allowed substitution for pattern assembly.
        "maxfl" : None, ##Integer. Maximum number of flanking residues for pattern assembly.
        "match" : None, ##Integer. Minimum number of matching residues for pattern assembly.
        "weight" : None, ##Integer. Minimum matching weight for pattern assembly.
        "max_asmb_size" : None, ##Integer. Maximum assembly size (number of patterns per assembly).
        "max_asmb_width" : None, ##Integer. Maximum assembly width.
        "asmb_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "asmb_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "1str" : None, ##Boolean. Use a single strand to build the motifs
        "2str" : None, ##Boolean. Use a both strands to build the motifs
        "prefix" : None, ##String. Matrix_prefix
        "cluster" : None, ##String. Run matrix-clustering to filter out redundant matrices, on the significance matrices (-cluster sig), on the count matrices (-cluster counts), on both (i<-cluster both>) or none (-cluster none).
        "sites" : None, ##String. Export the sites used to build the count matrix.
        "collect_method" : None, ##String. Method for converting sig matrices into count matrices. Supported - matrix-scan-quick (Default), info-gibbs (slow), matrix-scan (slow, obsolete)
        "gibbs_msps" : None, ##Integer. Mean number of sites per sequences passed to info-gibbs for converting significance matrices into count matrices
        "gibbs_iter" : None, ##Integer. Number of iterations for info-gibbs.
        "flanks" : None, ##Integer. Number of flanking residues to be added on each side of the significance matrix in order to extend the motif size
        "min_weight" : None, ##Integer. Minimal weight
        "gibbs_final" : None, ##Boolean. Run the final cycle with info-gibbs to collect the best sites.
        "logo" : None, ##Boolean. Export the sequence logos representing the count matrix.
        "links" : None, ##Boolean. Return HTML links in the convert-matrix result, to send the matrices to external tools (TOMTOM) for comparigon with motif collections.
        "scan_param" : None, ##String. The next argument is passed to matrix-scan (this will raise an error if these arguments are not supported). Ex. -scan_param '-uth Pval 1e-3 -uth rank 40
} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

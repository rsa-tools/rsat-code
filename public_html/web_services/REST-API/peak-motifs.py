#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"

#The endpoint indicates which RSAT resource you are interested in
ext = "/peak-motifs/" ##Workflow combining various algorithms to discover motifs from set of peak sequences, e.g. genomic regions obtained from ChIP-seq or related experiments (STARR-seq, ChIP-chip, ChIP-PET).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "ctrl_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "ctrl_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "max_seq_len" : None, ##Intger. msl. Maximal sequence length.
        "ref_motifs_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "ref_motifs_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "motifs_db" : None, ##String. db_name db_format db_file. File containing a database of transcription factor biding motifs (e.g. JASPAR, TRANSFAC, RegulonDB, etc.). Ex. -motif_db TRANSFAC transfac transfac_download_dr/cgi-bin/data/matrix.dat
        "title" : None, ##String. graph_title. Title displayed on top of the graphs
        "str" : "-2str", ##String. Single-strand (-1str) or double-strand (-2str) analysis
        "img_format" : None, ##String. Image format
        "task" : None, ##String. Default: purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan
        "disco" : None, ##String. Specify the software tool(s) that will be used for motif discovery.
        "nmotifs" : None, ##Integer. max_motif_number
        "minol" : None, ##Integer. Minimal lengths of oligonucleotide for word-counting approaches
        "maxol" : None, ##Integer. Maximal length of oligonucleotide for word-counting approaches
        "markov" : "auto", ##String. Order of the Markov model used to estimate expected oligonucleotide frequencies
        "min_markov" : None, ##Integer. Minimal value can be specified for the Markov order
        "max_markov" : None, ##Integer. Maximal value can be specified for the Markov order
        "noov" : None, ##Boolean. Treatment of self-overlapping words for motif discovery
        "r_plot" : None ##Boolean. Use R rather than the Perl GD library to generate plots.
    } 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

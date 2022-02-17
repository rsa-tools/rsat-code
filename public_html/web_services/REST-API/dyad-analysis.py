#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/dyad-analysis/" ##Detects overrepresented dyads (spaced pairs) in a set of DNA sequences. A dyad is defined here as a pair of oligonucleotides of the same size separated by a fixed number of bases.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "mask" : None, ##String. Mask lower or uppercases, respecively, i.e. replace selected case by N characters. Supported: upper, lower, none, non-dna
        "format" : "fasta", ##String. Input sequence format. Various standards are supported - raw, multi, ig, fasta, wconsensus, ncbi, tab.
        "l" : 3, ##Integer. Oligonucleotide size (default 3). This is the size of a single element (a half dyad).
        "spacing" : None, ##String. Default 0-20. The spacing is the number of bases between the end of the first element and the start of the second one.
        "type" : "any", ##String. dyad_type. In order to fasten execution, the program can be asked to restrict its analysis to symmetric dyads. Supported:  dr, ir, any, rep
        "accept_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "accept_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "groupsp" : None, ##Boolean. Group dyads made of the same words (monads) but with different spacings.
        "2str" : None, ##Boolean. Count on both strands
        "1str" : None, ##Boolean. Single strand count
        "prot" : None, ##Boolenan. Input sequence is proteic. In this case, the analysis concerns pairs of oligopeptides instead of oligonucleotides
        "expfreq_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "expfreq_string_type" : None, ##Type of ping: result file from other tool; text: input content
        "bg" : "upstream-noorf", ##Background model. Type of sequences used as background model for estimating expected dyad frequencies. Supported: upstream, upstream-noorf, intergenic, monads, input
        "org" : "Saccharomyces_cerevisiae", ##String. Query organism, to which the query genes belong.
        "taxon" : None, ##String. Reference taxon, in which orthologous genes have to be collected.
        "lth_occ" : 1, ##Number. Lower threshold on some parameters.
        "lth_occ_P" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_E" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_sig" : 0, ##Number. Lower threshold on some parameters.
        "lth_observed_freq" : None, ##Number. Lower threshold on some parameters.
        "lth_exp_freq" : None, ##Number. Lower threshold on some parameters.
        "lth_zscore" : None, ##Number. Lower threshold on some parameters.
        "lth_mseq" : None, ##Number. Lower threshold on some parameters.
        "lth_ms_P" : None, ##Number. Lower threshold on some parameters.
        "lth_ms_E" : None, ##Number. Lower threshold on some parameters.
        "lth_ms_sig" : None, ##Number. Lower threshold on some parameters.
        "lth_ratio" : None, ##Number. Lower threshold on some parameters.
        "lth_rank" : None, ##Number. Lower threshold on some parameters.
        "uth_occ" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_P" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_E" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_sig" : None, ##Number. Upper threshold on some parameters.
        "uth_observed_freq" : None, ##Number. Upper threshold on some parameters.
        "uth_exp_freq" : None, ##Number. Upper threshold on some parameters.
        "uth_zscore" : None, ##Number. Upper threshold on some parameters.
        "uth_mseq" : None, ##Number. Upper threshold on some parameters.
        "uth_ms_P" : None, ##Number. Upper threshold on some parameters.
        "uth_ms_E" : None, ##Number. Upper threshold on some parameters.
        "uth_ms_sig" : None, ##Number. Upper threshold on some parameters.
        "uth_ratio" : None, ##Number. Upper threshold on some parameters.
        "uth_rank" : 50, ##Number. Upper threshold on some parameters.
        "sort" : None, ##Boolean. Sort results by decreasing order of significance.
        "return" : "occ,proba,rank", ##String. output_fields. Output fields may contain one or several of the following words - freq, occ, proba, zscore, ratio, rank
        "under" : None, ##Boolean. Detect under-represented instead of over-represented dyads (left-tail significance test).
        "two_tails" : None, ##Boolean. Detect under-represented and over-represented dyads (two-tail significance test).
        "dry" : None, ##Boolenan. Dry run - print the commands but do not execute them.
        "zeroocc" : None, ##Boolean. Report also dyads with zero occurrences (provided they fit the other thresholds).
        "quick" : None, ##Boolean. Quick count mode -delegate the counting of word occurrences to count-words, a program written in C by Matthieu Defrance.
        "noov" : None, ##Boolean. Do not allow overlapping matches of the same word.
        "seqtype" : None ##Boolean. Input sequence type. Available values : dna, prot, other
}
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()


print(r.text)
# print (repr(r.json))

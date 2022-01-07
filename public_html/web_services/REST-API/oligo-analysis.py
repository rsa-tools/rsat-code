#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/oligo-analysis/" ##Calculates oligomer frequencies in a set of sequences, and detects overrepresented oligomers.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "mask" : None, ##String. Mask lower or uppercases, respecively, i.e. replace selected case by N characters. Supported: upper, lower
        "format" : "fasta", ##String. Input sequence format. Various standards are supported - fasta, wconsensus, IG, filelist, raw
        "l" : None, ##Integer. Oligonucleotide size.
        "expfreq_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "expfreq_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "calibN_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "calibN_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "calib1_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "calib1_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "bg" : None, ##String. Type of sequences used as background model for estimating expected oligonucleotide frequencies. Supported: upstream, upstreamL, upstream-noorf, intergenic, input
        "org" : None, ##String. Organism name for background model. Only when bgfile is not specified
        "taxon" : None, ##String. Organism or taxon that used as reference for the estimation of a background model based on a genome subset (option -bg). Either -org or -taxon is required with the option -bg.
        "markov" : None, ##Integer. Markov chain - the frequency expected for each word is calculated on basis of subword frequencies observed in the input set.
        "lexicon" : None, ##String. Expected word frequencies are calculated on the basis of subword frequencies, in a similar (but not identical) way to the “dictionary” approach developed by Harmen Bussemaker.
        "pseudo" : None, ##Number. Pseudo-frequency for the background model, where pseudo must be a real value between 0 and 1.
        "noov" : None, ##Boolean. No overlapping.
        "2str" : None, ##Boolean. Oligonucleotide occurrences found on both stands are summed.
        "1str" : None, ##Boolean. Inactivates the summation of occurrences on both strands.
        "seqtype" : None, ##String. Input sequence type. Supported: dna, prot, other
        "return" : None, ##String. List of statistics to return. Supported - occ,mseq,freq,proba,ratio,zscore,like,pos,rank
        "pal" : None, ##Boolean. Only return reverse palindroms
        "table" : None, ##Boolean. Return a table where rows represents input sequences, and columns the counts of occurrences for each possible oligo
        "distrib" : None, ##Boolean. Return occurrence distributions (one row per pattern)
        "grouprc" : None, ##Boolean. Group reverse complement with the direct sequence in the output file.
        "nogrouprc" : None, ##Boolean. Inactivates grouping of reverse complement pairs.
        "oneN" : None, ##Boolean. Group oligonucleotides by neighborhood, where one neighborhood is defined as a set of oligo differing by one mismatch at a common position.
        "onedeg" : None, ##Boolean. Sucessively insert one ambiguous nucleotide code at each position of each pattern
        "sort" : None, ##Boolean. Sort oligomers according to overrepresentation.
        "under" : None, ##Boolean. Detect under-represented instead of over-represented words (left-tail significance test, see below for details).
        "two_tails" : None, ##Boolean. Detect under-represented and over-represented words (two-tails significance test, see below for details).
        "zeroocc" : None, ##Boolean. Report also patterns with zero occurrences (provided they fit the other thresholds).
        "quick" : None, ##Boolean. Quick count mode - delegate the counting of word occurrences to count-words, a program written in C by Matthieu Defrance.
        "quick_if_possible" : None, ##Boolean. Evaluate if the quick mode is compatible with the selected output parameters, and, if so, run in this mode.
        "accept_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "accept_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "lth_occ" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_P" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_E" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_sig" : None, ##Number. Lower threshold on some parameters.
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
        "uth_rank" : None, ##Number. Upper threshold on some parameters. 
} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

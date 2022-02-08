#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/position-analysis/" ##Calculates the positional distribution of oligonucleotides in a set of sequences, and detects those which significantly discard from a homogeneous repartition.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "seqtype" : None, ##String. Input sequence type. Supported: dna, any.
        "last" : None, ##Integer. Stop after N sequences (for quick testing)
        "skip" : None, ##Integer. Skip the first N sequences.
        "first" : None, ##Integer. First sequence to analyze.
        "seqnb" : None, ##Integer. Maximal number of sequences to analyze.
        "mask" : None, ##String. Mask lower or uppercases, respecively, i.e. replace selected case by N characters. Supported: upper, lower.
        "format" : None, ##String. Input file format. Must be followed by one of the following options - fasta, wconsensus, IG, filelist, raw
        "l" : None, ##Integer. Oligomer length.
        "ci" : None, ##Integer. Alphabet. window interval (default 20 bases).
        "origin" : None, ##String. Reference for calculating positions. Supported: start, center, end
        "offset" : None, ##Integer. Add an offset to site positions.
        "noov" : None, ##Boolean. No overlapping.
        "2str" : None, ##Boolean. Oligonucleotide occurrences found on both stands are summed.
        "1str" : None, ##Boolean. Inactivates the summation of occurrences on both strands.
        "grouprc" : None, ##Boolean. Group reverse complement with the direct sequence in the output file.
        "nogrouprc" : None, ##Boolean. Inactivates grouping of reverse complement pairs.
        "sort" : None, ##Boolean. Sort oligomers according to overrepresentation.
        "return" : None, ##String. List of statistics to return. Supported - html,distrib,exp_occ,chi,rank,graphs,clusters
        "task" : None, ##String. Supported tasks - pos, clusters, matrices, graphs, index, all
        "markov" : None, ##Integer. Order for the Markov model use to compute position-specific expected word frequencies.
        "max_graphs" : None, ##Integer. Maximal number of graphs to export
        "pl_string" : None, ##String. Input string specifying the query.
        "pl_string_type" : None, ##String. Available values : text, url, piping
        "sc" : None, ##String. Score column. (only valid whith the option -pl)
        "minpos" : None, ##Integer. Minimal position to take into account for the chi-square calculation This value must be a multiple of the window interval.
        "maxpos" : None, ##Integer. Maximal position to take into account for the chi-square calculation This value must be a multiple of the window interval.
        "max_asmb_per_cluster" : None, ##Integer. 
        "nocheck" : None, ##Boolean. Do not check the applicability condition on the chi-square.
        "nofilter" : None, ##Boolean. Do not discard oligos which do not fit the condition of applicability
        "header" : None, ##String. Information to display in column headers of the distributions. Available values : mid, midfloor, min, max, interval
        "top_seq_for_matrices" : None, ##Integer. Select the top N sequences for building position-specific scoring matrices (PSSM).
        "img_format" : None, ##String. Image format (this parameter is passed to XYgraph).
        "title" : None, ##String. Title for the index table and position profile plots.
        "clust_method" : None, ##String. Agglomeration rule for the hierarchical clustering. Supported - complete, average, single, ward
        "clust_nb" : None, ##Integer. Number of clusters (default 8).
        "clust_suffix" : None, ##String. Suffix to append to the cluster file and the directory contianing cluster graphics. Default 'clusters’
        "lth_chi" : None, ##Number. Lower threshold on chi2
        "lth_sig" : None, ##Number. Lower threshold on significance
        "lth_occ" : None, ##Integer. Lower threshold on occurrences
        "uth_rank" : None, ##Integer. Upper threshold on rank
} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/matrix-scan/" ##Scan sequences with one or several position-specific scoring matrices (PSSM) to identify instances of the corresponding motifs (putative sites). This program supports a variety of background models (Bernoulli, Markov chains of any order).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "seq_format" : None, ##String. Sequence format.
        "m_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "m_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_format" : None, ##String. Matrix suffix. This argument is mandatory.
        "first_matrix" : None, ##Integer. Start scanning with the Nth matrix (kip the N-1 first matrices of the matrix file)
        "last_matrix" : None, ##Integer. Only scan with the top N matrices per matrix file.
        "matrix_ac" : None, ##String. Select one or more matrices specified by their ID. Ex. -matrix_ac MA0049.1,MA0221.1
        "matrix_name" : None, ##String. Select one or more matrices specified by their ID. Ex. -matrix_name eve,hb.
        "matrix_id" : None, ##String. Select one or more matrices specified by their ID. Ex. -matrix_id M00010,M00271
        "first_seq" : None, ##Integer. Start scanning at the Nth sequence.
        "last_seq" : None, ##Integer. Only scan with the top N sequences
        "consensus_name" : None, ##Boolenan. Use the motif (degenerate) consensus as matrix name.
        "id_as_name" : None, ##Boolean. Use the motif identifier as matrix name.
        "ac_as_name" : None, ##Boolean. Use the motif accession number as matrix name.
        "mlist_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "mlist_string_type" : None, ##Type of information provided by the input string. Available values : text, url, piping
        "mask" : None, ##String. Mask specific types of characters (lowercases, uppercases, non-dna), i.e. replace them by N characters. Supported. upper, lower, non-dna
        "n" : None, ##String. Treatment of N characters. These characters are often used in DNA sequences to represent undefined or masked nucleotides. Supported: skip, score.
        "pseudo" : None, ##Number. Pseudo-count for the matrix (default 1).
        "equi_pseudo" : None, ##Number. If this option is called, the pseudo-weight is distributed in an equiprobable way between residues.
        "org" : None, ##String. Organism name for background model. Only when bgfile is not specified
        "markov_order" : None, ##Integer. Markov order for background model. Only when bgfile is not specified.
        "bgfile_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bgfile_string_type" : None, ##Type of ping: result file from other tool; text: input content
        "bg_format" : None, ##String. Format of background model file. For supported formats see convert-background-model -h
        "bginput" : None, ##Boolean. Calculate background model from the input sequence set.
        "bg_pseudo" : None, ##Number. Pseudo frequency for the background model. Value must be a real between 0 and 1
        "markov" : None, ##Integer. Order of the markov chain for the background model.
        "window" : None, ## Integer. Size of the sliding window for the background model calculation.
        "origin" : None, ##String. Specify the origin for the calculation of positions. Supported: start, end, center, chrom
        "seq_source" : None, ## String. Sequence source for genomic coordinates. Supported. galaxy, getfasta, ucsc
        "offset" : None, ##Integer. Add a given number to site positions (change the reference point).
        "2str" : None, ##Boolean. Scan both strands for DNA sequences
        "1str" : None, ##Boolean. Single-strand search for DNA sequences
        "return" : None, ##String. lists of fields to return. Supported fields - sites, p_score, pval, seq_scores, rank, normw, proba_BM, limits,weight_limits, distrib, occ_proba, bg_model,bg_residues, matrix, freq_matrix, weight_matrix,crer
        "mth_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "mth_string_type" : None, ##Type of information provided by the input string. Available values : text, url, piping
        "crer_ids" : None, ##Boolean. Assign one separate feature ID per CRER.
        "sort_distrib" : None, ##Boolean. Sort score distributions by decreasing values of significance.
        "bg_distrib_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bg_distrib_string_type" : None, ##Type of information provided by the input string. Available values : text, url, piping
        "base" : None, ##Integer. Base for the logarithms (Default exp(1))
        "recursive" : None, ##Boolean. Run matrix-scan separately for each sequence.
        "batch" : None, ##Integer. Dispatch matrix-scan jobs on a cluster. Number of sequences to be analyzed by job (= on each node of the cluster)
        "lth_score" : None, ##Number. Lower threshold on some parameters.
        "lth_pval" : None, ##Number. Lower threshold on some parameters.
        "lth_sig" : None, ##Number. Lower threshold on some parameters.
        "lth_normw" : 0, ##Number. Lower threshold on some parameters.
        "lth_proba_M" : None, ##Number. Lower threshold on some parameters.
        "lth_proba_B" : None, ##Number. Lower threshold on some parameters.
        "lth_rank" : None, ##Number. Lower threshold on some parameters.
        "lth_rank_pm" : None, ##Number. Lower threshold on some parameters.
        "lth_crer_sig" : None, ##Number. Lower threshold on some parameters.
        "lth_crer_pval" : None, ##Number. Lower threshold on some parameters.
        "lth_crer_sites" : None, ##Number. Lower threshold on some parameters.
        "lth_crer_size" : None, ##Number. Lower threshold on some parameters.
        "lth_crer_site_distance" : None, ##Number. Lower threshold on some parameters.
        "lth_occ" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_sum" : None, ##Number. Lower threshold on some parameters.
        "lth_inv_cum" : None, ##Number. Lower threshold on some parameters.
        "lth_exp_occ" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_pva" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_eval" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_sig" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_sig_rank" : None, ##Number. Lower threshold on some parameters.
        "uth_score" : None, ##Number. Upper threshold on some parameters.
        "uth_pval" : None, ##Number. Upper threshold on some parameters.
        "uth_sig" : None, ##Number. Upper threshold on some parameters.
        "uth_normw" : 0, ##Number. Upper threshold on some parameters.
        "uth_proba_M" : None, ##Number. Upper threshold on some parameters.
        "uth_proba_B" : None, ##Number. Upper threshold on some parameters.
        "uth_rank" : None, ##Number. Upper threshold on some parameters.
        "uth_rank_pm" : None, ##Number. Upper threshold on some parameters.
        "uth_crer_sig" : None, ##Number. Upper threshold on some parameters.
        "uth_crer_pval" : None, ##Number. Upper threshold on some parameters.
        "uth_crer_sites" : None, ##Number. Upper threshold on some parameters.
        "uth_crer_size" : None, ##Number. Upper threshold on some parameters.
        "uth_crer_site_distance" : None, ##Number. Upper threshold on some parameters.
        "uth_occ" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_sum" : None, ##Number. Upper threshold on some parameters.
        "uth_inv_cum" : None, ##Number. Upper threshold on some parameters.
        "uth_exp_occ" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_pva" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_eval" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_sig" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_sig_rank" : None ##Number. Upper threshold on some parameters.
}
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()


print(r.text)
# print (repr(r.json))

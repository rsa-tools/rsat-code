#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/compare-classes/" ##Compare two class/cluster files (the query file and the reference file) and report the intersection between each pair of classes/clusters + some statistics about this intersection (hypergeometric P-value, E-value, …).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "r_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "r_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "q_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "q_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "sc" : None, ##Integer. Score column. Specify a column of the input file containing a score associated to each member.
        "mames_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "mames_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "qmames_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "qmames_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "max_lines_q" : None, ##Integer. Max number of lines to read for the query file.
        "max_lines_r" : None, ##Integer. Max number of lines to read for the reference file.
        "max_lines" : None, ##Integer. Max number of lines to read for both query and reference files
        "cat" : None, ##String. Compare the query file to pre-defined catalogs (e.g. GO, MIPS functional classes, …).
        "org" : None, ##String. Organism (for pre-defined catalogs)
        "return" : None, ##String. Return fields. Supported - Q_only,R_only,common,dotprod,entropy,freq,jac_sim,members,occ,proba,rank,sor_sim
        "pop" : None, ##Integer. Population size
        "sort" : None, ##String. Sort on the basis of the specified key.
        "rep" : None, ##String. Replacement. Sampling was performed with replacement, i.e. a given element can appear several times in the same class.
        "sym" : None, ##Boolean. Symmetric comparison.
        "distinct" : None, ##Boolean. Prevent to compare each class with itself (when the reference and query files contain the same classes).
        "triangle" : None, ##Boolean. Do not perform the reciprocal comparisons - if reference A has already been compared to query B, then reference B does not need to be compared to query A.
        "matrix" : None, ##String. Return a pairwise matrix, where each row corresponds to a reference class, each column to a query class, and each cell contains a comparison between the two classes
        "margins" : None, ##Boolean. 
        "null" : None, ##Boolean. Null string (default NA) displayed for undefined values.
        "base" : None, ##Number. Logarithm base (Default - 2.71828182845905) used for entropy-based metrics
        "dot" : None, ##String. Export a graph with the associations in a dot file.
        "gml_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "gml_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "multi_cor" : None, ##String. Factor used for the multi-testing correction. Supported - nq, nr, nc
        "lth_DPbits" : None, ##Number. Lower threshold on some parameters.
        "lth_E(QR)" : None, ##Number. Lower threshold on some parameters.
        "lth_E_val" : None, ##Number. Lower threshold on some parameters.
        "lth_F(!Q!R)" : None, ##Number. Lower threshold on some parameters.
        "lth_F(Q!R)" : None, ##Number. Lower threshold on some parameters.
        "lth_F(Q)" : None, ##Number. Lower threshold on some parameters.
        "lth_F(QR)" : None, ##Number. Lower threshold on some parameters.
        "lth_F(Q!R)" : None, ##Number. Lower threshold on some parameters.
        "lth_F(R)" : None, ##Number. Lower threshold on some parameters.
        "lth_H(Q)" : None, ##Number. Lower threshold on some parameters.
        "lth_H(Q,R)" : None, ##Number. Lower threshold on some parameters.
        "lth_H(Q|R)" : None, ##Number. Lower threshold on some parameters.
        "lth_I(Q,R)" : None, ##Number. Lower threshold on some parameters.
        "lth_IC" : None, ##Number. Lower threshold on some parameters.
        "lth_P(QR)" : None, ##Number. Lower threshold on some parameters.
        "lth_P(Q|R)" : None, ##Number. Lower threshold on some parameters.
        "lth_P(R|Q)" : None, ##Number. Lower threshold on some parameters.
        "lth_P_val" : None, ##Number. Lower threshold on some parameters.
        "lth_Q" : None, ##Number. Lower threshold on some parameters.
        "lth_QR" : None, ##Number. Lower threshold on some parameters.
        "lth_Q_only" : None, ##Number. Lower threshold on some parameters.
        "lth_QvR" : None, ##Number. Lower threshold on some parameters.
        "lth_R" : None, ##Number. Lower threshold on some parameters.
        "lth_R_only" : None, ##Number. Lower threshold on some parameters.
        "lth_U(Q|R)" : None, ##Number. Lower threshold on some parameters.
        "lth_U(R|Q)" : None, ##Number. Lower threshold on some parameters.
        "lth_common" : None, ##Number. Lower threshold on some parameters.
        "lth_dH(Q,R)" : None, ##Number. Lower threshold on some parameters.
        "lth_dotprod" : None, ##Number. Lower threshold on some parameters.
        "lth_jac_sim" : None, ##Number. Lower threshold on some parameters.
        "lth_rDPbits" : None, ##Number. Lower threshold on some parameters.
        "lth_rank" : None, ##Number. Lower threshold on some parameters.
        "lth_sig" : None, ##Number. Lower threshold on some parameters.
        "lth_sor_sim" : None, ##Number. Lower threshold on some parameters.
        "lth_sqrt_dp" : None, ##Number. Lower threshold on some parameters.
        "uth_DPbits" : None, ##Number. Upper threshold on some parameters.
        "uth_E(QR)" : None, ##Number. Upper threshold on some parameters.
        "uth_E_val" : None, ##Number. Upper threshold on some parameters.
        "uth_F(!Q!R)" : None, ##Number. Upper threshold on some parameters.
        "uth_F(Q!R)" : None, ##Number. Upper threshold on some parameters.
        "uth_F(Q)" : None, ##Number. Upper threshold on some parameters.
        "uth_F(QR)" : None, ##Number. Upper threshold on some parameters.
        "uth_F(Q!R)" : None, ##Number. Upper threshold on some parameters.
        "uth_F(R)" : None, ##Number. Upper threshold on some parameters.
        "uth_H(Q)" : None, ##Number. Upper threshold on some parameters.
        "uth_H(Q,R)" : None, ##Number. Upper threshold on some parameters.
        "uth_H(Q|R)" : None, ##Number. Upper threshold on some parameters.
        "uth_I(Q,R)" : None, ##Number. Upper threshold on some parameters.
        "uth_IC" : None, ##Number. Upper threshold on some parameters.
        "uth_P(QR)" : None, ##Number. Upper threshold on some parameters.
        "uth_P(Q|R)" : None, ##Number. Upper threshold on some parameters.
        "uth_P(R|Q)" : None, ##Number. Upper threshold on some parameters.
        "uth_P_val" : None, ##Number. Upper threshold on some parameters.
        "uth_Q" : None, ##Number. Upper threshold on some parameters.
        "uth_QR" : None, ##Number. Upper threshold on some parameters.
        "uth_Q_only" : None, ##Number. Upper threshold on some parameters.
        "uth_QvR" : None, ##Number. Upper threshold on some parameters.
        "uth_R" : None, ##Number. Upper threshold on some parameters.
        "uth_R_only" : None, ##Number. Upper threshold on some parameters.
        "uth_U(Q|R)" : None, ##Number. Upper threshold on some parameters.
        "uth_U(R|Q)" : None, ##Number. Upper threshold on some parameters.
        "uth_common" : None, ##Number. Upper threshold on some parameters.
        "uth_dH(Q,R)" : None, ##Number. Upper threshold on some parameters.
        "uth_dotprod" : None, ##Number. Upper threshold on some parameters.
        "uth_jac_sim" : None, ##Number. Upper threshold on some parameters.
        "uth_rDPbits" : None, ##Number. Upper threshold on some parameters.
        "uth_rank" : None, ##Number. Upper threshold on some parameters.
        "uth_sig" : None, ##Number. Upper threshold on some parameters.
        "uth_sor_sim" : None, ##Number. Upper threshold on some parameters.
        "uth_sqrt_dp" : None, ##Number. Upper threshold on some parameters.
} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

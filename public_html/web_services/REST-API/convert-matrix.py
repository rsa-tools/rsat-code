#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/convert-matrix/" ##Performs inter-conversions between various formats of position-specific scoring matrices (PSSM).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "residue_type" : None, ##String. Supported: dna, cytomod
        "mlist_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "mlist_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "mlist_name_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "mlist_name_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_id_file_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_id_file_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "split" : None, ##Boolean. Split a single multi-matrices input file in a set of separate files. The output file names start with the prefix specificed by the option -o, followed by a suffix indicating the order of the matrix in the input file (m1, m2, …).
        "bg_pseude" : 0.01, ##Number. Pseudo frequency for the background models. Value must be a real between 0 and 1.
        "from" : "tab", ##String. Available values : alignace, assembly, cb, cis-bp, clustal, cluster-buster, consensus, encode, feature, footprintdb, gibbs, homer, info-gibbs, infogibbs, jaspar, meme, meme_block, motifsampler, mscan, sequences, stamp, stamp-transfac, tab, tf, transfac, uniprobe, yeastract
        "to" : "tab", ##String. Output matrix format. Supported - cb,cluster-buster,consensus,infogibbs,jaspar,param_table,patser,stamp,tab,tf,tomtom,transfac
        "return" : "counts,counts,consensus,parameters,logo", ##String. Return type. Supported - consensus,counts,frequencies,header,info,information,links,logo,logo_matrix,logo_table,margins,parameters,profile,sites,wdistrib,weights
        "sort" : None, ##String. Sort by key. Available values : desc, asc, alpha
        "top" : None, ##Integer. Maximal number of matrices to return.
        "skip" : None, ##Integer. Skip the first N matrices.
        "pseudo" : 1, ##Integer. Pseudo-weight used for the calculation of the weight matrix
        "equi_pseudo" : None, ##Boolean. If this option is called, the pseudo-weight is distributed in an equiprobable way between residues.
        "multiply" : 1, ##Integer. Multiply all the values of the input matrices by the number
        "rescale" : None, ##Integer. Scale the matrix to a fixed value for the sums per columns.
        "insert_col_left" : None, ##Integer. Insert columns on the left flank of the count matrix. The inserted columns are filled with zeros.
        "insert_col_right" : None, ##Integer. Insert columns on the right flank of the count matrix. The inserted columns are filled with zeros.
        "insert_col" : None, ##Integer. Insert columns on the both flanks of the count matrix. The inserted columns are filled with zeros.
        "trim_col_left" : None, ##Integer. Remove columns on the left flank of the count matrix.
        "trim_col_right" : None, ##Integer. Remove columns on the right flank of the count matrix.
        "trim_col" : None, ##Integer. Remove columns on the both flanks of the count matrix.
        "base" : None, ##Integer. Base for the logarithms used in the scores involving a log-likelihood (weight and information content).
        "decimals" : 1, ##Integer. Number of decimals to print for real matrices (frequencies, weights, information) or to compute score distributions.
        "prefix" : None, ##String. Prefix to be added before identifier(s) and name(s) of the input matrix/matrices.
        "attr" : None, ##String. key value. Force an attribute of the matrix (matrices) to have a given value.
        "perm" : 0, ##Integer. Number of permuted matrices to return.
        "max_profile" : None, ##Integer. Maximal width of the profile histogram (units equal number of characters).
        "rc" : None, ##Boolean. Convert the matrix to its reverse complement.
        "logo_program" : "weblogo", ##String. External program used to generate logo drawings. Supported: weblogo, seqlogo.
        "logo_format" : "png", ##String. Format for logo image file. Supported - eps,jpeg,logodata,pdf,png,png_print,svg
        "logo_file" : None ##String. Specifies the name of the logo file.
}
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()


print(r.text)
#print (repr(r.json))

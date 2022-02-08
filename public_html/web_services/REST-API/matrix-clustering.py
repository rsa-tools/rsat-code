#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/matrix-clustering/" ##Apply hierarchical clustering to identify clusters of similar motifs.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "matrix_1_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_1_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_title_1" : None, ##String. The matrix_title will be concatenated to each motif ID in order to create unique motif IDs.
        "matrix_format_1" : None, ##String. Since the program takes several matrices as input, it only accepts matrices in formats supporting several matrices per file (transfac, tf, tab, cluster-buster, cb, infogibbs, meme, stamp, uniprobe).
        "matrix_2_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_2_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_title_2" : None, ##String. The matrix_title will be concatenated to each motif ID in order to create unique motif IDs.
        "matrix_format_2" : None, ##String. Since the program takes several matrices as input, it only accepts matrices in formats supporting several matrices per file (transfac, tf, tab, cluster-buster, cb, infogibbs, meme, stamp, uniprobe).
        "matrix_3_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_3_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_title_3" : None, ##String. The matrix_title will be concatenated to each motif ID in order to create unique motif IDs.
        "matrix_format_3" : None, ##String. Since the program takes several matrices as input, it only accepts matrices in formats supporting several matrices per file (transfac, tf, tab, cluster-buster, cb, infogibbs, meme, stamp, uniprobe).
        "matrix_file_table_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_file_table_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "max_matrices" : None, ##Integer. This option specifies how many matrices can be clustered in the same analysis. If there are more matrices than the specified number, the program restrics the analyses to the first X matrices, and issues a warning.
        "title" : None, ##String. Title displayed on top of the report page
        "ID_link_color_table" : None, ##String. his option allows to add a link to a any website specified by the user and can be used to visualize complete databases (e.g. Jaspar), thus each motif in the logo tree will point to its respective link in the Jaspar website
        "label_in_tree" : None, ##String. Option to select the labels displayed in the logo tree.
        "task" : None, ##String. Specify one or several tasks to be run. Supported tasks - all, comparison, clustering, report
        "hclust_method" : None, ##String. Option to select the agglomeration rule for hierarchical clustering. Available values : complete, average, single
        "metric_build_tree" : None, ##String. Select the metric which will be used to cluster the motifs.based in one metric of to measure motif similarity. Available values : cor, Ncor, dEucl, NdEucl, logocor, logoDP, Nlogocor, Icor, NIcor, SSD, rank_mean, mean_zscore
        "lth_w" : None, ##Number. Lower threshold.
        "lth_cor" : None, ##Number. Lower threshold.
        "lth_Ncor" : None, ##Number. Lower threshold.
        "uth_w" : None, ##Number. Upper threshold.
        "Uth_cor" : None, ##Number. Upper threshold.
        "uth_Ncor" : None, ##Number. Upper threshold.
        "calc" : None, ##String. Specify the operator used to merge matrices (argument passed to merge-matrices). Available values : mean, sum
        "quick" : None, ##Boolean. With this option the motif comparison is done with the program compare-matrices-quick (implemented in C) rather than the program compare-matrices (implemented in Perl)
        "top_matrices" : None, ##Integer. Only analyze the first X motifs of the input file. This options is convenient for quick testing before starting the full analysis.
        "skip_matrices" : None, ##Integer. Skip the first X motifs of the input file. This options is convenient for testing the program on a subset of the motifs before starting the full analysis.
        "return" : "json,heatmap", ##String. List of fields to return. Supported fields - heatmap,json,newick,root_matrices.
        "o" : None, ##String. Prefix for the output files.

} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 
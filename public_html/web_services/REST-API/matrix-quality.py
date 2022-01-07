#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/matrix-quality/" ##Evaluate the quality of a Position-Specific Scoring Matrix (PSSM), by comparing score distributions obtained with this matrix in various sequence sets.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "m_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "m_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "ms_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "ms_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "top" : None, ##Integer. Maximal number of matrices to analyze.
        "matrix_format" : None, ##String. Matrix_format
        "seq_type" : None, ##String. The type of the sequence (which will appear in the leend of the plots). This option is mandatory.
        "seq_file_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "seq_file_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "seq_type_2" : None, ##String. The type of the sequence (which will appear in the leend of the plots). This option is mandatory.
        "seq_file_2_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "seq_file_2_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "seq_format" : None, ##String. Sequence format.
        "scanopt" : None, ##String. seq_type "option1 option2 …". Sequence set-specific options for matrix-scan. These options are added at the end of the matrix-scan command for scanning the specified sequence set.
        "no_cv" : None, ##Boolean. Do not apply the leave-one-out (LOO) test on the matrix site sequences.
        "kfold" : None, ##Integer. k-fold cross-validation. Divide the matrix sites in k chunks for cross-validation
        "noperm" : None, ##Boolean. Skip the matrix permutation step.
        "noscan" : None, ##Boolean. Skip the matrix-scan step.
        "nocompa" : None, ##Boolean. Skip the step of comparisons between distributions.
        "nograph" : None, ##Boolean. Skip the step of drawing comparison graphs.
        "noicon" : None, ##Boolean. Do not generate the small graphs (icons) used for the galleries in the indexes.
        "export_hits" : None, ##Boolean. Return matrix-scan scores in addition to the distribution of scores
        "perm_sep" : None, ##Boolean. Calculate the distributions for each permuted matrix separately.
        "perm" : None, ##Integer. Number of permutations for a specific set (default 0).
        "perm_2" : None, ##Integer. Number of permutations for a specific set (default 0). Optional.
        "pseudo" : None, ##Integer. Pseudo-counts. The pseudo-count reflects the possibility that residues that were not (yet) observed in the model might however be valid for future observations
        "org" : None, ##String. Organism name for background model. Only if bgfile is not specified.
        "markov_order" : None, ##Integer. Markov order for background model. Only if bgfile is not specified.
        "bgfile_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bgile_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "bg_format" : None, ##String. Format for the background model file.
        "bg_pseudo" : None, ##Number. From 0 to 1
        "decimal" : None, ##Integer. Number of decimals for computing weight scores (default 2).
        "graph_option" : None, ##String. Specify options that will be passed to the program XYgraph for generating the distributions and the ROC curves.
        "roc_ref" : None, ##String. Reference distribution for the ROC curve.
        "roc_option" : None, ##String. Specify options that will be passed to the program XYgraph for generating the ROC curves (ot the distribution curves)
        "distrib_option" : None, ##String. Specify options that will be passed to the program XYgraph for generating the distribution curves (not the ROC curves).
        "img_format" : None, ##String. Image format for the plots (ROC curve, score profiles, …).
        "r_plot" : None, ##Boolean. Generate plots using R instead of the Perl GD module.
        "logo_format" : None, ##String. Image format for the sequence logos.
        "plot" : None, ##String. Additions plots will be drawn to compare - a) The enrichment of scores in a set of sequences for different matrices b) The enrichment of scores in different sequence sets for one matrix. Supported - nwd, occ_proba
        "plot_2" : None, ##String. Additions plots will be drawn to compare - a) The enrichment of scores in a set of sequences for different matrices b) The enrichment of scores in different sequence sets for one matrix. Supported - nwd, occ_proba
        "archive" : None, ##Boolean. Available values : true, false
        "html_title" : None, ##String. Get a title for the html page.
        "task" : None, ##String. Specify one or several tasks to be run. Supported - scan,theor,loo,theor_cv,permute,compare,graphs,synthesis,plot
} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

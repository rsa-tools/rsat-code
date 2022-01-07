#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/footprint-scan/" ##Scan promoters of orthologous genes with one or several position-specific scoring matrices (PSSM) in order to detect enriched motifs, and thereby predict phylogenetically conserved target genes.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "q" : "lexA", ##String. Query gene. Can be multiple genes, separated by ‘,’
        "genes_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server
        "genes_string_type" : "text", ##Type of information provided by the input string (URL, piping, text)
        "m_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "m_string_type" : "text", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_format" : None, ##Matrix suffix. This argument is mandatory. 
        "org_list_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "org_list_string_type" : None, ##Type of information provided by the input string (URL, piping, text)
        "org" : "Escherichia_coli_GCF_000005845.2_ASM584v2", ##String. Query organism, to which the query genes belong.
        "taxon" : "Gammaproteobacteria", ##String. Reference taxon, in which orthologous genes have to be collected.
        "orthologs_list_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "orthologs_list_string_type" : None, ##Type of information provided by the input string (URL, piping, text)
        "tf" : None, ##String. Transcription_factor.
        "pseudo" : None, ##Integer. Pseudo-count for the matrix (default 1). See matrix-scan for details.
        "bgfile_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bgfile_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_table_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_table_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "bg_format" : None, ##String. Format of background model file. For supported formats see convert-background-model -h
        "bginput" : None, ##Boolean. Calculate background model from the input sequence set.
        "bg_pseudo" : None, ##Number. Pseudo frequency for the background model. Value must be a real between 0 and 1.
        "markov" : None, ##Integer. Order of the markov chain for the background model.
        "window" : None, ##Integer. Size of the sliding window for the background model calculation.
        "tree" : None, ##String. Phylogenetic tree use for the Bayesian Branch Length Scoring task.
        "bbls_draw" : None, #String. Fomat for the output image. Supported. png, pdf, none.
        "all_genes" : None, ##Boolenan. Automatically analyze all the genes of a query genome, and store each result in a separate folder (the folder name is defined automatically)
        "unique_species" : None, ##Boolean. Retain at most one organism per species. This enables to filter out the numerous strains sequences for some species of particular interest. (e.g. Escherichia coli, Bacillus subtilis, …).
        "unique_genus" : None, ##Boolean. Retain at most one organism per genus. Same filter as for -unique_species, but at the level of the genus.
        "bg_model" : None, ##String. Background model. Allow the user to choose among alternative background model (see Janky & van Helden, 2008). Supported: monads, taxfreq, org_list, file.
        "pval" : None, ##String. Set the threshold on site p-value to report only the evaluated over-representations of binding sites whenever the individual sites crossed it. The default is set to 1e-4.
        "occ_th" : None, ##Number. Threshold set on the occurrence significance (over-representation) for scores that have p-value equal or smaller thant the one given as threshold in the option -pval.
        "plot_format" : None, ##String. Format for the occurrence plots (occurrence frequencies, occurrence sinificance). Supported - all formats supported by the program XYgraph
        "skip_m" : None, ##Integer. Skip the first N matrices in the matrix_table (useful for quick testing and for resuming interrupted tasks when using a matrix_table or when several matrices are entered with the option -m ).
        "last_m" : None, ##Number. Stop after having treated the first N matrices in the matrix table (useful for quick testing when using a matrix_table or when several matrices are entered with the option -m ).
        "dist_thr" : None, ##Integer. Specify here the intergenic distance threshold in base pairs.
        "task" : None, ##String. Specify a subset of tasks to be executed. Supported - all, operons, query_seq, orthologs, ortho_seq, purge, gene_index, index, orthologs_tf, occ_sig, occ_sig_graph, scan, map.
        "no_purge" : None, ##Boolean. This option can only be used combined with the -org_list option, this gives the posibility to analyse a given set of sequences managing sequence redundancy using a list of “no redundant” organisms.
        "use_tree_org" : None, ##Boolean. Only uses organisms in the phylogenetic tree for orthologs search.
        "batch" : None, ##Boolean. Generate one command per query gene, and post it on the queue of a PC cluster.
        "dry" : None, ##Boolenan. Dry run - print the commands but do not execute them.
        "nodie" : None, ##Boolean. Do not die in case a sub-program returns an error.
        "sep_genes" : None, ##Boolean. Search footprints for each query gene separately.
        "info_lines" : None, ##Boolean. Draw reference lines on the significance profile plots, to highlight some particular values.
        "infer_operons" : None, ##Boolean. Infer operons in order to retrieve the promoters of the predicted operon leader genes rather than those located immediately upstream of the orthologs.
        "batch_matrix" : None, ##Boolean. Generate one footprint-scan command per matrix and post it on the queue of a PC cluster.
        "occ_sig_opt" : None, ##String. Additional options passed to matrix-scan for the test of over-representation of matrix hits. Ex. -occ_sig_opt '-uth rank 1’
        "occ_sig_graph_opt" : None, ##Boolean. Additional options passed to XYgraph for drawing the occurrence significance graph.
        "scan_opt" : None, ##String. Additional options passed to matrix-scan for site detection and feature-map drawing. Ex. -scan_opt '-uth pval 0.001’
        "map_opt" : None, ##String. Additional options passed to feature-map for feature-map drawing. Ex. -map_opt '-mapthick 12’
        "filter_bgfile_string" : None, ##String. Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "filter_bgfile_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "filter_pval" : None, ##Number. Set the threshold to filter out TF-interactions that are not present on the query organism.
        "rand" : None, ##Boolenan. When the option -rand is activated, the program replaces each ortholog by a gene selected at random in the genome where this ortholg was found.
        "crer" : None, ##Boolean. Return Cis-Regulatory elements Enriched-Regions (CRER).
        "diamond" : None, ##Boolean. Use ranks_dmnd.tab from diamond blast computed in genome-blast.
        "synthesis" : None, ##Boolean. This option generated synthetic tables (in tab-delimited text and html) for all the results.
        "map_format" : None ##String. Format for the feature map. Supported: jpg, gif, png, ps.
} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
# print (repr(r.json))
 

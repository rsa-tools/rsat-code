#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/footprint-discovery/" ##Detect phylogenetic footprints by applying dyad-analysis in promoters of a set of orthologous genes.

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "q" : "lexA", ##String. Query gene. Can be multiple genes, separated by ‘,’
        "genes_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server
        "genes_string_type" : "text", ##Type of information provided by the input string (URL, piping, text)
        "org_list_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "org_list_string_type" : None, ##Type of information provided by the input string (URL, piping, text)
        "orthologs_list_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "orthologs_list_string_type" : None, ##Type of information provided by the input string (URL, piping, text)
        "org" : "Escherichia_coli_GCF_000005845.2_ASM584v2", ##String. Query organism, to which the query genes belong.
        "taxon" : "Gammaproteobacteria", ##String. Reference taxon, in which orthologous genes have to be collected.
        "max_dyad_degree" : None, ##Integer. Maximal dyad degree for network inference. Default 20.
        "max_genes" : None, ##Integer. Maximal number of genes to analyze.
        "skip" : None, ##Integer. Skip the first N genes (useful for quick testing and for resuming interrupted tasks). 
        "last" : None, ##Integer. Stop after having treated the first N genes (useful for quick testing).
        "use_tree_org" : None, ##Boolean. Only uses organisms in the phylogenetic tree for orthologs search.
        "all_genes" : None, ##Boolenan. Automatically analyze all the genes of a query genome, and store each result in a separate folder (the folder name is defined automatically)
        "unique_species" : None, ##Boolean. Retain at most one organism per species. This enables to filter out the numerous strains sequences for some species of particular interest. (e.g. Escherichia coli, Bacillus subtilis, …).
        "unique_genus" : None, ##Boolean. Retain at most one organism per genus. Same filter as for -unique_species, but at the level of the genus.
        "bg_model" : None, ##String. Background model. Allow the user to choose among alternative background model (see Janky & van Helden, 2008). Supported: monads, taxfreq, org_list, file.
        "bgfile_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bgfile_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "lth_occ" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_P" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_E" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_sig" : None, ##Number. Lower threshold on some parameters.
        "lth_observed_freq" : None, ##Number. Lower threshold on some parameters.
        "lth_zscore" : None, ##Number. Lower threshold on some parameters.
        "lth_ratio" : None, ##Number. Lower threshold on some parameters.
        "lth_rank" : None, ##Number. Lower threshold on some parameters.
        "uth_occ" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_P" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_E" : None, ##Number. Upper threshold on some parameters.
        "uth_occ_sig" : None, ##Number. Upper threshold on some parameters.
        "uth_observed_freq" : None, ##Number. Upper threshold on some parameters.
        "uth_zscore" : None, ##Number. Upper threshold on some parameters.
        "uth_ratio" : None, ##Number. Upper threshold on some parameters.
        "uth_rank" : None, ##Number. Upper threshold on some parameters. 
        "dist_thr" : None, ##Integer. Specify here the intergenic distance threshold in base pairs.
        "task" : None, ##String. Specify a subset of tasks to be executed. Supported - all, operons, query_seq, orthologs, ortho_seq, purge, gene_index, index, filter_dyads, dyads, map, network
        "return" : "occ,proba,rank", ##String. Return fields for dyad-analysis. This argument is passed to dyad-analysis for the discovery of dyads in promoters of orthologous genes. output fields may contain one or several of the following words - freq, occ, proba, zscore, ratio, rank
        "no_purge" : None, ##Boolean. This option can only be used combined with the -org_list option, this gives the posibility to analyse a given set of sequences managing sequence redundancy using a list of “no redundant” organisms.
        "batch" : None, ##Boolean. Generate one command per query gene, and post it on the queue of a PC cluster.
        "dry" : None, ##Boolenan. Dry run - print the commands but do not execute them.
        "nodie" : None, ##Boolean. Do not die in case a sub-program returns an error.
        "sep_genes" : None, ##Boolean. Search footprints for each query gene separately.
        "infer_operons" : None, ##Boolean. Infer operons in order to retrieve the promoters of the predicted operon leader genes rather than those located immediately upstream of the orthologs.
        "filter" : None, ##Boolean. Only accept dyads found in the promoter of the query gene, in the query organism. (option selected by default)
        "no_filter" : None, ##Boolean. Accept all dyads, even if they are not found in the promoter of the query gene, in the query organism. (will cancel -filter option if selected).
        "rand" : None, ##Boolenan. When the option -rand is activated, the program replaces each ortholog by a gene selected at random in the genome where this ortholg was found.
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
 

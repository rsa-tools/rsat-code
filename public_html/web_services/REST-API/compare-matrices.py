#Begin by importing the Requests module
import requests, sys

#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/compare-matrices/" ##Compare two collections of position-specific scoring matrices (PSSM), and return various similarity statistics + matrix alignments (pairwise, one-to-n).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "file1_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "file1_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "file2_string" : "http://rsat-tagc.univ-mrs.fr/rsat/motif_databases/JASPAR/Jaspar_2020/nonredundant/JASPAR2020_CORE_vertebrates_non-redundant_pfms.tf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "file2_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "file_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "file_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "format1" : None, ##String. Matrix_format1. Specify the matrix format for the first input file only (requires format2)
        "format2" : None, ##String. Matrix_format2. Specify the matrix format for the second input file only (requires format1).
        "top1" : None, ##Integer. Only analyze the first X motifs of the first file. This options is convenient for quick testing before starting the full analysis.
        "top2" : None, ##Integer. Only analyze the first X motifs of the second file. This options is convenient for quick testing before starting the full analysis.
        "skip1" : None, ##Integer. Skip the first X motifs of the first input matrix file. This option can be combined with the option -top1 in order to restrict the analysis to a given subset of a large library.
        "skip2" : None, ##Integer. Skip the first X motifs of the second input matrix file. This options can be combined with the option -top2 in order to restrict the analysis to a given subset of a large library.
        "mode" : None, ##String. matches or scores or profiles or scan
        "distinct" : None, ##Boolenan. Skip comparison between a matrix and itself
        "stand" : "DR", ##String. Perform matrix comparisons in direct (D) reverse complementary ® or both orientations (DR, default option). When the R or DR options are activated, all matrices of the second matrix file are converted to the reverse complementary matrix.
        "return" : "cor,Ncor,logoDP,NsEucl,NSW,match_rank,matrix_id,matrix_name,width,strand,offset,consensus,alignments_1ton", ##String. return_fields. List of fields to return (only valid for the formats “profiles” and “matches”). Supported return fields - offset, cor, Ncor, Ncor1, Ncor2, NcorS, cov, SSD, NSW, SW, dEucl, NdEucl, NsEucl, dKL, matrix_number, matrix_id, matrix_name, matrix_label, matrix_ac, width, strand, offset, pos, consensus, offset_rank, match_rank, graph, alignments_pairwise, alignments_lton, alignments, logos, matrix_desc, all
        "labels" : None, ##String. Attributes to inclute in the matrix labels of the description table.
        "lth_w" : None, ##Number. Lower threshold on some parameters. Width = number of aligned columns
        "lth_cor" : None, ##Number. Lower threshold on some parameters. Pearson correlation (computed on residue occurrences in aligned columns)
        "lth_Ncor" : None, ##Number. Lower threshold on some parameters. Relative width-normalized Pearson correlation
        "lth_logoDP" : None, ##Number. Lower threshold on some parameters. Dot product of sequence logos
        "lth_logocor" : None, ##Number. Lower threshold on some parameters. Correlation computed on sequence logos
        "lth_Nlogocor" : None, ##Number. Lower threshold on some parameters. Relative width-normalized logocor
        "lth_Icor" : None, ##Number. Lower threshold on some parameters. Pearson correlation computed on Information content
        "lth_NIcor" : None, ##Number. Lower threshold on some parameters. Relative width-normalized Icor
        "lth_cov" : None, ##Number. Lower threshold on some parameters. Covariance between residues in aligned columns
        "lth_dEucl" : None, ##Number. Lower threshold on some parameters. Euclidian distance between residue occurrences in aligned columns
        "lth_NdEucl" : None, ##Number. Lower threshold on some parameters. Relative width-normalized dEucl
        "lth_NsEucl" : None, ##Number. Lower threshold on some parameters. Similarity derived from Relative width-normalized Euclidian distanc
        "lth_SSD" : None, ##Number. Lower threshold on some parameters. Sum of square deviations
        "lth_SW" : None, ##Number. Lower threshold on some parameters. Sandelin-Wasserman
        "lth_NSW" : None, ##Number. Lower threshold on some parameters. Relative width-normalized Sandelin-Wasserman
        "lth_match_rank" : None, ##Number. Lower threshold on some parameters. Rank of current match among all sorted matches
        "lth_offset" : None, ##Number. Lower threshold on some parameters. Offset between first and second matrices
        "uth_w" : None, ##Number. Upper threshold on some parameters. Width = number of aligned columns
        "uth_cor" : None, ##Number. Upper threshold on some parameters. Pearson correlation (computed on residue occurrences in aligned columns)
        "uth_Ncor" : None, ##Number. Upper threshold on some parameters. Relative width-normalized Pearson correlation
        "uth_logoDP" : None, ##Number. Upper threshold on some parameters. Dot product of sequence logos
        "uth_logocor" : None, ##Number. Upper threshold on some parameters. Correlation computed on sequence logos
        "uth_Nlogocor" : None, ##Number. Upper threshold on some parameters. Relative width-normalized logocor
        "uth_Icor" : None, ##Number. Upper threshold on some parameters. Pearson correlation computed on Information content
        "uth_NIcor" : None, ##Number. Upper threshold on some parameters.  Relative width-normalized Icor
        "uth_cov" : None, ##Number. Upper threshold on some parameters. Covariance between residues in aligned columns
        "uth_dEucl" : None, ##Number. Upper threshold on some parameters. Euclidian distance between residue occurrences in aligned columns
        "uth_NdEucl" : None, ##Number. Upper threshold on some parameters. Relative width-normalized dEucl
        "uth_NsEucl" : None, ##Number. Upper threshold on some parameters. Similarity derived from Relative width-normalized Euclidian distanc
        "uth_SSD" : None, ##Number. Upper threshold on some parameters. Sum of square deviations
        "uth_SW" : None, ##Number. Upper threshold on some parameters. Sandelin-Wasserman
        "uth_NSW" : None, ##Number. Upper threshold on some parameters. Relative width-normalized Sandelin-Wasserman
        "uth_match_rank" : None, ##Number. Upper threshold on some parameters. Rank of current match among all sorted matches
        "uth_offset" : None ##Number. Upper threshold on some parameters. Offset between first and second matrices
}
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()


##print(r.text)
# print (repr(r.json))
name_of_file = input("What is the name of the file: ")
completeName = name_of_file + ".html"
f = open(completeName, "w+")
f.write(r.text)
f.close()

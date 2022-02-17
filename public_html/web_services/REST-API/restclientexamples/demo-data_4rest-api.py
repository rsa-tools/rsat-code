## fungi server
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
## prokaryotes server
server = "http://embnet.ccg.unam.mx/rest.wsgi"
## metazoa server
server = "http://rsat.sb-roscoff.fr/rest.wsgi"

##classfreq

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://embnet.ccg.unam.mx/rsat/demo_files/allup500_Saccharomyces_cerevisiae_some_pattern_counts.tab", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "max" : None, ##Integer. Numbers strictly greater than max are not taken into account.
        "min" : None, ##Integer. Numbers strictly smaller than min are not taken into account. min also determines the lower limit of the first class.
        "ci" : 1, ##Integer. Class interval. If not specified, takes the value (max - min)/20 so that 21 classes are specified.
        "col" : 4, ##Integer. Column to which apply the program. This option can be used iteratively.
        "from" : None, ##Integer. Inferior limit for the classes to display values lower than this limit are however taken into account in the calculation of statistics (mean, variance, …) and of class frequencies (In contrast with the -min option).
        "to" : None, ##Integer. Superior limit for the classes to display values higher than this limit are however taken into account in the calculation of statistics (mean, variance, …) and of class frequencies (In contrast with the -max option).
        "thr" : None, ##Integer. Threshold. Only display classes with absolute frequency higher than or equal to the threshold. This option is useful for sparse data, where many classes do not contain any observation (-thr 1).
    }

## compare-features

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/Blow_2010_GSM348064_forebrain_p300_peaks.bed", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "filelist_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "filelist_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "ref_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/Blow_2010_GSM559653_midbrain_p300_peaks.bed", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "ref_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "iformat" : "bed", ##String. Input feature format. bed, dnapat, ft, galaxy_seq, gft, gff, gff3bed, swembl, ucsc_seq
        "self" : None, ##Boolean. Also perform comparison between features in the same file (self-comparison). This can be useful to detect redundancy between annotated features.
        "return" : "stats,inter", ##String. Specify the output type(s). Supported output types. stats,inter,diff.
        "lth_inter_len" : 1, ##Number. Minimal overlap (bp). Lower threshold on Length (in residues) of the intersection between two features.
        "lth_inter_cov" : 0.25, ##Number. Intersection coverage (0-1). Lower threshold on Coverage of the intersection between two features.
    }

## compare-matrices

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "file1_string" : "http://rsat-tagc.univ-mrs.fr/rsat/motif_databases/JASPAR/Jaspar_2020/nonredundant/JASPAR2020_CORE_non-redundant_pfms.tf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
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
        "lth_w" : 5, ##Number. Lower threshold on some parameters. Width = number of aligned columns
        "lth_cor" : 0.7, ##Number. Lower threshold on some parameters. Pearson correlation (computed on residue occurrences in aligned columns)
        "lth_Ncor" : 0.4, ##Number. Lower threshold on some parameters. Relative width-normalized Pearson correlation
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

## convert-features

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/convert-features_demo.pat", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "from" : "dnapat", ##String. Input format. Supported - galaxy_seq,ucsc2_seq,bed,gff,gff3,swembl,getfasta_seq,dnapat,gtf,ft,bed3col,gft
        "to" : "ft", ##String. output format. Supported - fasta,bed3col,gft,ft,dnapat,great,gff3,bed,gff.
        "add_chrm" : None, ##Boolean. only for BED output - add “chr” in front of the chromosome name and change MT into chrM
        "remove_chrm" : None, ##Boolean. Remove the “chr” in front of the chromosone and change chrM into MT
        "yeast_to_roman" : None, ##Boolean. Convert arabic to roman numbers to denote chromosome numbers according to S.cerevisiae specifications. Only works for chromosomes I to XVI.
        "featname" : None, ##Boolean. Set a name for all features of the file. This option can be convenient for conversions to bed files.
        "summits" : None, ##Boolean. only valid for SWEMBL input - replace start and end coordinates by peak summit position
        "extend" : None, ##Integer. Extend peak coordinates on both sides (start and end).
        "extend_start" : None, ##Integer. Extend start coordinate leftwise
        "extend_end" : None, ##Integer. Extend end coordinate rightwise.
        "coord" : None, ##Boolean. bedfile with absolute coordinate of the sequence relative to which the features were defined (e.g. features from promoter-wise to genome-wise coordinates).
        "origin" : None, ##String. Origin of coordinates relative to sequence fragment. This option is only valid when combined with the option coord. Supported: start, end, center
    }

## convert-matrix

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/Dmelanogaster_segmentation_12matrices.tf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "residue_type" : "dna", ##String. Supported: dna, cytomod
        "mlist_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "mlist_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "mlist_name_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "mlist_name_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_id_file_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_id_file_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "split" : None, ##Boolean. Split a single multi-matrices input file in a set of separate files. The output file names start with the prefix specificed by the option -o, followed by a suffix indicating the order of the matrix in the input file (m1, m2, …).
        "bg_pseude" : 0.01, ##Number. Pseudo frequency for the background models. Value must be a real between 0 and 1.
        "from" : "tf", ##String. Available values : alignace, assembly, cb, cis-bp, clustal, cluster-buster, consensus, encode, feature, footprintdb, gibbs, homer, info-gibbs, infogibbs, jaspar, meme, meme_block, motifsampler, mscan, sequences, stamp, stamp-transfac, tab, tf, transfac, uniprobe, yeastract
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

## convert-seq

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/dna-pattern_demo_seq.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "mask" : None, ##String. Masked characters are replaced by by N characters, or by a dot (option -dotmask). Supported. upper, lower, non-dna
        "noempty" : None, ##Boolean. Remove empty sequences from the set (same as -skip_short 1).
        "mask_short" : None, ##Integer. Mask (replace by N characters) sequences shorter than the specified length.
        "skip_short" : 30, ##Integer. Skip sequences shorter than the specified length. Same functionality as -mask_short, except that short sequences are not returned at all in the output.
        "skip_long" : None, ##Integer.  Skip sequences longer than the specified length.
        "last" : None, ##Integer. Stop after the Nth sequence.
        "top" : None, ##Integer. Same as -last N
        "first" : None, ##Integer. Start at the Nth sequence (skip the N-1 first sequences).
        "skip" : None, ##Integer. Skip the N first sequences (start at sequence N+1).
        "from" : "fasta", ##String. Input format. Supported - embl, fasta, filelist, ft, gcg, genbank, ig, maf,multi, ncbi, raw,tab,wc, wconsensus
        "to" : "wconsensus", ##String. Output format. Supported - fasta, fastq, filelist, ft, ig, multi, raw, tab, wc, wconsensus
        "id_col" : None, ##Integer. Column containing sequence identifiers in tab format.
        "seq_col" : None, ##Integer. Column containing sequence sequences in tab format.
        "comment_col" : None, ##Integer. Column containing sequence comments (sequence description) in tab format.
        "lw" : 60, ##Integer. Line width. A carriage return is inserted every lw characters within the output sequence.
        "addrc" : None, ##Boolean. Adds the reverse complement of each input sequence to the output file.
        "lc" : None, ##Boolean. Lowercase. the sequence is printed in lowercase.
        "uc" : None, ##Boolean. Uppercase. the sequence is printed in uppercase.
        "dna" : None, ##Boolean. Convert any non-acgt character into “n” characters.
        "dotmask" : None, ##Boolean. Convert masked characters into dots.
        "id" : None, ##String. Identifier. sequence identifier (useful for converting a raw sequence from the STDIN)
        "prefix" : None, ##String. Sequence prefix (useful for converting from a multi sequence)
        "nocheckid" : None, ##Boolean. Prevent to check sequence IDs for conversion to file list.
}

## convert-variations

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs.gvf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "from" : "gvf", ##String. Format of the input file. Supported format vcf, gvf, varBed
        "to" : "varBed", ##String. Format of the output file. Supported format vcf, gvf, varBed
        "phased" : True, ##Boolean. Retrieve phasing information of the vcf genotypes in order to write each entry phased.
    }

## create-background-model

data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat//demo_files/ChIP-seq_peaks/Oct4_peaks_top1000.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "out_format" : "transitions", ##String. Output format. Supported transitions, tab, oligo-analysis, oligos, meme, motifsampler, ms, inclusive, patser, tables.
        "markov" : 2 ##Integer. Markov order.
    }

## count-words

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/Dmelanogaster_eve_up5000.fasta", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "l" : 6, ##Integer. Set oligomer length to l (monad size when using dyads)
        "2str" : None, ##Boolean. Add reverse complement.
        "1str" : None, ##Boolean. Do not add reverse complement.
        "sp" : None, ##String. Spacing between the two parts of the dyads from N to M. Ex.0-20
        "noov" : None, ##Boolean. Do not allow overlapping occurrences.
        "grouprc" : None, ##Boolean. Group reverse complement with the direct sequence.
        "nogrouprc" : None, ##Boolean. Do not group reverse complement with the direct sequence.
    }

## crer-scan

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/Drosophila_melanogaster_eve_segmentation_sites_pval0.001_nocomments.ft", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "in_format" : "ft", ##String. Available values : ft, bed
        "autoparam" : None, ##Boolean. Extract some input parameters from the commented rows (starting with ‘;’) of the input file.
        "s" : None, ##Boolean. Sort the list of sites.
        "return_limits" : False, ##Boolean. Return every limits of sequences
        "return_limits_filtered" : None, ##Boolean. Return the limits filtered of the sequence
        "uth_site_pval" : 0.001, ##String. Maximal p-value of sites to be considered.
        "number_of matrix" : None, ##Integer. number of matrix used for the discovery of transcription factor binding sites.
        "lth_score" : None, ##Number. Minimal site score to be considered
        "uth_score" : 500, ##Number. Maximal site score to be considered
        "lth_crer_size" : 100, ##Number. Minimal size of the enriched region (in bp).
        "uth_crer_size" : None, ##String. Maximal size of the enriched region (in bp).
        "lth_crer_sites" : 2, ##Number. Minimal number of sites covered by the enriched region
        "uth_crer_sites" : None, ##Number. Maximal number of sites covered by the enriched region
        "lth_crer_sites_distance" : 1, ##Number. Distance between successive sites to be considered.
        "uth_crer_sites_distance" : 100, ##Number. Distance between successive sites to be considered.
        "uth_crer_pval" : None, ##String. Maximal binomial p-value
        "uth_crer_eval" : None, ##String. Maximal e-value
        "lth_crer_sig" : 0.1, ##Number. Minimal binomial significance
        "uth_overlap" : 1, ##Integer. Maximal overlap to define two distinct sites
        "nopval" : None, ##Boolean. Compute crer without p value
        "pre_table" : None, ##Boolean. Compute a table where is all possible p_value
    }

## dna-pattern

data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/dna-pattern_demo_seq.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "mask" : None, ##String. Mask lower or uppercases, respecively, i.e. replace selected case by N characters. Supported: upper, lower
        "format" : "fasta", ##String. Input sequence format. The accepted formats are fasta, IG, raw, multi, filelist
        "pl_string" : "CACGTG CACGTT", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "pl_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.
        "subst" : None, ##Integer. Allow subst substitutions.
        "noIUPAC" : None, ##Boolean. The pattern is considered as a standard regular expression.
        "sc" : None, ##Integer. Score column
        "noid" : None, ##Boolean. Do not search pattern identifier in the second column of pattern file.
        "noov" : False, ##Boolean. Do not count overlapping matches for self-overlapping patterns.
        "2str" : True, ##Boolean. Search matches on both strands (direct and reverse complement)
        "1str" : None, ##Boolean. Search matches only on the direct strand.
        "R" : None, ##Boolean. Search matches only on the reverse complement strand
        "id" : None, ##Boolean. Pattern identifier (one word).
        "return" : None, ##String. List of fields to return. Multiple fields can be entered separated by commas, or by using iteratively the option.
        "match_format" : None, ##String. Format for returning matches (supported - fasta, table)
        "th" : None, ##Integer. Threshold. Return match count only for sequences with greater than or equal number of matches
        "merge" : None, ##Boolean. Merge mutually overlapping matches.
        "N" : 4, ##Integer. Return matching sequences with N flanking nucleotides.
        "NL" : None, ##Integer. Return matching sequences with NL left flanking nucleotides.
        "NR" : None, ##Integer. Return matching sequences with NR right flanking nucleotides.
        "origin" : None, ##Integer. Define origin as the origin for the calculation of positions.
        "window" : None, ##Integer. Sliding window size.
        "top" : None, ##Boolean. (With sliding window only). only return the top score obtained with the sliding window for each sequence.
        "sort" : None, ##Boolean. (With -top only). sort sequences according to their top score
    }

## dyad-analysis

data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/dyad_analysis_demo_seq.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "mask" : None, ##String. Mask lower or uppercases, respecively, i.e. replace selected case by N characters. Supported: upper, lower, none, non-dna
        "format" : "fasta", ##String. Input sequence format. Various standards are supported - raw, multi, ig, fasta, wconsensus, ncbi, tab.
        "l" : 3, ##Integer. Oligonucleotide size (default 3). This is the size of a single element (a half dyad).
        "spacing" : None, ##String. Default 0-20. The spacing is the number of bases between the end of the first element and the start of the second one.
        "type" : "any", ##String. dyad_type. In order to fasten execution, the program can be asked to restrict its analysis to symmetric dyads. Supported:  dr, ir, any, rep
        "accept_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "accept_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "groupsp" : None, ##Boolean. Group dyads made of the same words (monads) but with different spacings.
        "2str" : None, ##Boolean. Count on both strands
        "1str" : None, ##Boolean. Single strand count
        "prot" : None, ##Boolenan. Input sequence is proteic. In this case, the analysis concerns pairs of oligopeptides instead of oligonucleotides
        "expfreq_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "expfreq_string_type" : None, ##Type of ping: result file from other tool; text: input content
        "bg" : "upstream-noorf", ##Background model. Type of sequences used as background model for estimating expected dyad frequencies. Supported: upstream, upstream-noorf, intergenic, monads, input
        "org" : "Saccharomyces_cerevisiae", ##String. Query organism, to which the query genes belong.
        "taxon" : None, ##String. Reference taxon, in which orthologous genes have to be collected.
        "lth_occ" : 1, ##Number. Lower threshold on some parameters.
        "lth_occ_P" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_E" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_sig" : 0, ##Number. Lower threshold on some parameters.
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
        "uth_rank" : 50, ##Number. Upper threshold on some parameters.
        "sort" : None, ##Boolean. Sort results by decreasing order of significance.
        "return" : "occ,proba,rank", ##String. output_fields. Output fields may contain one or several of the following words - freq, occ, proba, zscore, ratio, rank
        "under" : None, ##Boolean. Detect under-represented instead of over-represented dyads (left-tail significance test).
        "two_tails" : None, ##Boolean. Detect under-represented and over-represented dyads (two-tail significance test).
        "dry" : None, ##Boolenan. Dry run - print the commands but do not execute them.
        "zeroocc" : None, ##Boolean. Report also dyads with zero occurrences (provided they fit the other thresholds).
        "quick" : None, ##Boolean. Quick count mode -delegate the counting of word occurrences to count-words, a program written in C by Matthieu Defrance.
        "noov" : None, ##Boolean. Do not allow overlapping matches of the same word.
        "seqtype" : None ##Boolean. Input sequence type. Available values : dna, prot, other
}

## feature-map

data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/feature-map_demo_data.ft", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "format" : "png", ##String. Output image format. Supported:  jpg, gif, png, ps
        "from" : -800, ##Integer. Lower limit of the positions represented on the graph
        "to" : 0, ##Integer. Upper limit of the positions represented on the graph
        "title" : "Motifs discovered in upstream sequences of 19 MET genes", ##String. Generic Title for the feature map.
        "label" : "id", ##String. Keylist define the info to display for each feature. Valid keys are id, strand, descr, pos
        "boxlabels" : None, ##String. Writes the label within the feature box (by default, the label is written outside of this box).
        "symbol" : False, ##Boolean. Associates a graphical symbol (i.e. rectangle, circle, buterfly, …) to each feature.
        "dot" : False, ##Boolean. A color dot is associated to each feature.
        "mlen" : 500, ##Integer. Map length (in pixels). Default is 600.
        "mapthick" : 25, ##Integer. Thickness refers to either width (for vertical maps) or height (horizintal maps).
        "mspacing" : 2, ##Integer. The size of the border between maps (in pixel).
        "border" : None, ##Integer. Image border (default=10 pixels)
        "origin" : 0, ##Integer. All coordinates are recalculated relative to origin.
        "legend" : True, ##Boolean. Draws a legend on the graph, showing the symbol associated to each distinct feature.
        "scalebar" : True, ##Boolean. Draw a scale bar on the left of the graph.
        "scalestep" : 50, ##Integer. Step between annotations of the scale bar.
        "no_name" : True, ##Boolean. Do not print sequence names besides each sequence.
        "scorethick" : True, ##Boolean. Each feature is displayed with a thickness proportional to its score. Only positive scores are represented.
        "maxscore" : None, ##Integer. (Only valid when -scorethick is active). Maximal allowed score value. Higher score values are clipped for the drawing.
        "minscore" : None, ##Integer. (Only valid when -scorethick is active). minimal allowed score value. Features with smaller score are not displayed.
        "maxfthick" : None, ##Integer. Max feature thickness.
        "minfthick" : None, ##Integer. Min feature thickness.
        "htmap" : None, ##Boolean. HTML map.
        "mono" : False, ##Boolean. Monochrome palette (for printing on black/white printer).
        "aacolors" : None, ##Boolean. Amino acid specific colors.
        "colors_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "colors_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "export_colors" : None, ##String. Export the feature-specific colors in a separate file as a color file.
        "bgcolor" : "220,220,220", ##String. Background color in R,G,B format.
        "horizontal" : True, ##Boolean. Horizontal map (default).
        "vertical" : None, ##Boolean. Vertical map (default is horizontal).
        "select" : None, ##String. Id_list. Only display the features whose ID is in id_list. -select ‘gataag’,’gattag’
        "seq_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "seq_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "seqformat" : None, ##String. Format of the reference sequence file.
    }

## fetch-sequences

data =  {
        "genome" : "mm9", ##String. Genome version (e g mm9, hg19)
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "header_format" : "galaxy", ##Format for sequence headers. Supported header formats. UCSC, galaxy
    }

## footprint-discovery

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
        "bg_model" : "taxfreq", ##String. Background model. Allow the user to choose among alternative background model (see Janky & van Helden, 2008). Supported: monads, taxfreq, org_list, file.
        "bgfile_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bgfile_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "lth_occ" : 1, ##Number. Lower threshold on some parameters.
        "lth_occ_P" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_E" : None, ##Number. Lower threshold on some parameters.
        "lth_occ_sig" : 0, ##Number. Lower threshold on some parameters.
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
        "uth_rank" : 50, ##Number. Upper threshold on some parameters.
        "dist_thr" : None, ##Integer. Specify here the intergenic distance threshold in base pairs.
        "task" : None, ##String. Specify a subset of tasks to be executed. Supported - all, operons, query_seq, orthologs, ortho_seq, purge, gene_index, index, filter_dyads, dyads, map, network
        "return" : "occ,proba,rank", ##String. Return fields for dyad-analysis. This argument is passed to dyad-analysis for the discovery of dyads in promoters of orthologous genes. output fields may contain one or several of the following words - freq, occ, proba, zscore, ratio, rank
        "no_purge" : None, ##Boolean. This option can only be used combined with the -org_list option, this gives the posibility to analyse a given set of sequences managing sequence redundancy using a list of “no redundant” organisms.
        "batch" : None, ##Boolean. Generate one command per query gene, and post it on the queue of a PC cluster.
        "dry" : None, ##Boolenan. Dry run - print the commands but do not execute them.
        "nodie" : None, ##Boolean. Do not die in case a sub-program returns an error.
        "sep_genes" : True, ##Boolean. Search footprints for each query gene separately.
        "infer_operons" : None, ##Boolean. Infer operons in order to retrieve the promoters of the predicted operon leader genes rather than those located immediately upstream of the orthologs.
        "filter" : None, ##Boolean. Only accept dyads found in the promoter of the query gene, in the query organism. (option selected by default)
        "no_filter" : None, ##Boolean. Accept all dyads, even if they are not found in the promoter of the query gene, in the query organism. (will cancel -filter option if selected).
        "rand" : None, ##Boolenan. When the option -rand is activated, the program replaces each ortholog by a gene selected at random in the genome where this ortholg was found.
        "diamond" : None, ##Boolean. Use ranks_dmnd.tab from diamond blast computed in genome-blast.
        "synthesis" : None, ##Boolean. This option generated synthetic tables (in tab-delimited text and html) for all the results.
        "map_format" : None ##String. Format for the feature map. Supported: jpg, gif, png, ps.
}

## footprint-scan

data =  {
        "q" : "lexA,recA", ##String. Query gene. Can be multiple genes, separated by ‘,’
        "genes_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server
        "genes_string_type" : "text", ##Type of information provided by the input string (URL, piping, text)
        "m_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/LexA.2nt_upstream-noorf-ovlp-2str.20.tf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "m_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_format" : None, ##Matrix suffix. This argument is mandatory.
        "org_list_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "org_list_string_type" : None, ##Type of information provided by the input string (URL, piping, text)
        "org" : "Escherichia_coli_GCF_000005845.2_ASM584v2", ##String. Query organism, to which the query genes belong.
        "taxon" : "Gammaproteobacteria", ##String. Reference taxon, in which orthologous genes have to be collected.
        "orthologs_list_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "orthologs_list_string_type" : None, ##Type of information provided by the input string (URL, piping, text)
        "tf" : None, ##String. Transcription_factor.
        "pseudo" : 1, ##Integer. Pseudo-count for the matrix (default 1). See matrix-scan for details.
        "bgfile_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bgfile_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_table_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_table_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "bg_format" : None, ##String. Format of background model file. For supported formats see convert-background-model -h
        "bginput" : None, ##Boolean. Calculate background model from the input sequence set.
        "bg_pseudo" : None, ##Number. Pseudo frequency for the background model. Value must be a real between 0 and 1.
        "markov" : 0, ##Integer. Order of the markov chain for the background model.
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

## gene-info

data =  {
        "org" : "Saccharomyces_cerevisiae", ##String. Query organism
        "q" : None, ##String query. The query should be an orf identifier (eg ‘metR’). The query is case-insensitive.
        "i_string" : "ARG3,PHO,YBL05[\d]W", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.
        "full" : None, ##Boolean. Full match only (no substring matching)
        "noquery" : None, ##Boolean. Do not print the query at the beginning of each line
        "descr" : None, ##Boolean. Match queries against the description (in addition to gene ID and names)
        "feattype" : "gene" ##Boolean. Feature type. Supported - gene,mRNA,tRNA,rRNA,scRNA,misc_RNA,CDS,start_codon,stop_codon,exon
    }

## infer-operon

data =  {
        "i_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "text", ##Type of information provided by the input string. Available values : text, url, piping.
        "format" : "id", ##String. Input format. Supported formats = id, varBed, bed
        "org" : "Escherichia_coli_GCF_000005845.2_ASM584v2", ##String. Query organism
        "all" : False, ##Boolean. Infer operons for all the genes of the query organism.
        "q" : "bioD,trpB,trpE,hisB,metB", ##String. Query gene. This option can be used iteratively on the same command line to specify several query genes.
        "dist" :55, ##Integer. Distance threshold.
        "sep" : ",", ##String. Specify the separator for multi-value fields (e.g. genes) in the output table.
        "min_gene_nb" : 2, ##Integer. Specify a threshold on the number of genes in the operon.
        "return" : "query,name,leader,operon,upstr_dist", ##String. List of fields to return. Separated by ',’. Supported - leader,trailer,operon,query,name,upstr_dist,q_info,up_info,down_info
    }

## info-gibbs

data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/info-gibbs_demo_seq.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "w" : 20, ##Integer. Set the motif width to w
        "maxspacing" : None, ##Integer. Set maximal spacing between motif monad to maxspacing (only for dyadic motif).
        "minspacing" : None, ##Integer. Set minimal spacing between motif monad to minspacing (only for dyadic motif).
        "strand" : None, ##String. Search in foward strand + or in both strands ±
        "n" : 1000, ##Integer. Maximum number of Gibbs sampling iterations.
        "sites" : None, ##Integer. Number of motif occurrences that are expected to be found (incompatible with -e)
        "e" : 1, ##Number. Mean number of motif occurrences (sites) expected per sequence that are expected to be found (incompatible with --sites)
        "zoops" : None, ##Boolean. Try to find 0 or 1 site per sequence
        "m" : 1, ##Integer. Number of motifs to extract (one by default)
        "b" : None, ##Integer. Use b predefined INCLUSive background model
        "d" : None, ##Integer. Set minimal distance between 2 motif occurrences to d
        "t" : None, ##Number. Set the temperature (should be in range [0.6 1.4])
        "r" : 5, ##Integer. Try to run the Gibbs sampling seach r times
        "collect" : None, ##Integer. Try to collect the N best sites using their weight scores
        "seedmatrix" : None ##String. Start sampling form sites collected by scanning the sequences with matrix seedmatrix
    }

## IUPAC-to-regular

data =  {
        "q" : "ACGGWN[ACGT]YCGT", ##Query String.
    }

## matrix-clustering

data =  {
        "matrix_1_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/peak-motifs_Oct4_matrices.tf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_1_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_title_1" : "Oct4_peak_motifs", ##String. The matrix_title will be concatenated to each motif ID in order to create unique motif IDs.
        "matrix_format_1" : "transfac", ##String. Since the program takes several matrices as input, it only accepts matrices in formats supporting several matrices per file (transfac, tf, tab, cluster-buster, cb, infogibbs, meme, stamp, uniprobe).
        "matrix_2_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_2_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_title_2" : None, ##String. The matrix_title will be concatenated to each motif ID in order to create unique motif IDs.
        "matrix_format_2" : None, ##String. Since the program takes several matrices as input, it only accepts matrices in formats supporting several matrices per file (transfac, tf, tab, cluster-buster, cb, infogibbs, meme, stamp, uniprobe).
        "matrix_3_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_3_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_title_3" : None, ##String. The matrix_title will be concatenated to each motif ID in order to create unique motif IDs.
        "matrix_format_3" : None, ##String. Since the program takes several matrices as input, it only accepts matrices in formats supporting several matrices per file (transfac, tf, tab, cluster-buster, cb, infogibbs, meme, stamp, uniprobe).
        "matrix_file_table_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "matrix_file_table_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "max_matrices" : None, ##Integer. This option specifies how many matrices can be clustered in the same analysis. If there are more matrices than the specified number, the program restrics the analyses to the first X matrices, and issues a warning.
        "title" : None, ##String. Title displayed on top of the report page
        "ID_link_color_table" : None, ##String. his option allows to add a link to a any website specified by the user and can be used to visualize complete databases (e.g. Jaspar), thus each motif in the logo tree will point to its respective link in the Jaspar website
        "label_in_tree" : None, ##String. Option to select the labels displayed in the logo tree.
        "task" : None, ##String. Specify one or several tasks to be run. Supported tasks - all, comparison, clustering, report
        "hclust_method" : "average", ##String. Option to select the agglomeration rule for hierarchical clustering. Available values : complete, average, single
        "metric_build_tree" : "Ncor", ##String. Select the metric which will be used to cluster the motifs.based in one metric of to measure motif similarity. Available values : cor, Ncor, dEucl, NdEucl, logocor, logoDP, Nlogocor, Icor, NIcor, SSD, rank_mean, mean_zscore
        "lth_w" : None, ##Number. Lower threshold.
        "lth_cor" : None, ##Number. Lower threshold.
        "lth_Ncor" : None, ##Number. Lower threshold.
        "uth_w" : None, ##Number. Upper threshold.
        "Uth_cor" : None, ##Number. Upper threshold.
        "uth_Ncor" : None, ##Number. Upper threshold.
        "calc" : "sum", ##String. Specify the operator used to merge matrices (argument passed to merge-matrices). Available values : mean, sum
        "quick" : None, ##Boolean. With this option the motif comparison is done with the program compare-matrices-quick (implemented in C) rather than the program compare-matrices (implemented in Perl)
        "top_matrices" : None, ##Integer. Only analyze the first X motifs of the input file. This options is convenient for quick testing before starting the full analysis.
        "skip_matrices" : None, ##Integer. Skip the first X motifs of the input file. This options is convenient for testing the program on a subset of the motifs before starting the full analysis.
        "return" : "json,heatmap", ##String. List of fields to return. Supported fields - heatmap,json,newick,root_matrices.
        "o" : None, ##String. Prefix for the output files.

}


## matrix-distrib

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "m_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/peak-motifs_Oct4_matrices.tf", ##String. Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "m_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "top" : None, ##Integer. Top_matrix_nb. Restrict the analysis to the N top matrices of the input file.
        "mlist_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "mlist_string_type" : None, ##String. Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_format" : "transfac", ##String. Matrix format. Default is tab.
        "pseudo" : 1, ##Integer. Pseudo-count for the matrix (default 1).
        "org" : "Saccharomyces_cerevisiae", ##String. Organism for background model file
        "markov_order" : 0, ##Integer. Markov order for the background model file
        "bgfile_string" : None, ##String. Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bgfile_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "bg_format" : None, ##String. Background_format. Supported formats - all the input formats supported by convert-background-model.
        "bg_pseudo" : None, ##Number. Pseudo frequency for the background models. Value must be a real between 0 and 1 (default 0)
        "decimals" : 1, ##Integer. Number of decimals to print or the transition probabilities.
    }

## matrix-quality

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "m_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/matrix_quality_demo_matrix.tf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "m_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "ms_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "ms_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "top" : None, ##Integer. Maximal number of matrices to analyze.
        "matrix_format" : "transfac", ##String. Matrix_format
        "seq_type" : "fasta", ##String. The type of the sequence (which will appear in the leend of the plots). This option is mandatory.
        "seq_file_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/matrix_quality_demo_seq1.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "seq_file_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "seq_type_2" : "fasta", ##String. The type of the sequence (which will appear in the leend of the plots). This option is mandatory.
        "seq_file_2_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/matrix_quality_demo_seq2.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "seq_file_2_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
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
        "perm" : 1, ##Integer. Number of permutations for a specific set (default 0).
        "perm_2" : 1, ##Integer. Number of permutations for a specific set (default 0). Optional.
        "pseudo" : None, ##Integer. Pseudo-counts. The pseudo-count reflects the possibility that residues that were not (yet) observed in the model might however be valid for future observations
        "org" : "Escherichia_coli_GCF_000005845.2_ASM584v2", ##String. Organism name for background model. Only if bgfile is not specified.
        "markov_order" : 1, ##Integer. Markov order for background model. Only if bgfile is not specified.
        "bgfile_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bgile_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "bg_format" : None, ##String. Format for the background model file.
        "bg_pseudo" : 0.01, ##Number. From 0 to 1
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


## matrix-scan

data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/Dmelanogaster_eve_up5000.fasta", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "seq_format" : "fasta", ##String. Sequence format.
        "m_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/Dmelanogaster_segmentation_12matrices.tf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "m_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "matrix_format" : "transfac", ##String. Matrix suffix. This argument is mandatory.
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
        "n" : "score", ##String. Treatment of N characters. These characters are often used in DNA sequences to represent undefined or masked nucleotides. Supported: skip, score.
        "pseudo" : 1, ##Number. Pseudo-count for the matrix (default 1).
        "equi_pseudo" : None, ##Number. If this option is called, the pseudo-weight is distributed in an equiprobable way between residues.
        "org" : None, ##String. Organism name for background model. Only when bgfile is not specified
        "markov_order" : 1, ##Integer. Markov order for background model. Only when bgfile is not specified.
        "bgfile_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bgfile_string_type" : None, ##Type of ping: result file from other tool; text: input content
        "bg_format" : None, ##String. Format of background model file. For supported formats see convert-background-model -h
        "bginput" : True, ##Boolean. Calculate background model from the input sequence set.
        "bg_pseudo" : 0.01, ##Number. Pseudo frequency for the background model. Value must be a real between 0 and 1
        "markov" : 1, ##Integer. Order of the markov chain for the background model.
        "window" : None, ## Integer. Size of the sliding window for the background model calculation.
        "origin" : "end", ##String. Specify the origin for the calculation of positions. Supported: start, end, center, chrom
        "seq_source" : None, ## String. Sequence source for genomic coordinates. Supported. galaxy, getfasta, ucsc
        "offset" : 0, ##Integer. Add a given number to site positions (change the reference point).
        "2str" : True, ##Boolean. Scan both strands for DNA sequences
        "1str" : None, ##Boolean. Single-strand search for DNA sequences
        "return" : "sites,pval,limits", ##String. lists of fields to return. Supported fields - sites, p_score, pval, seq_scores, rank, normw, proba_BM, limits,weight_limits, distrib, occ_proba, bg_model,bg_residues, matrix, freq_matrix, weight_matrix,crer
        "mth_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "mth_string_type" : None, ##Type of information provided by the input string. Available values : text, url, piping
        "crer_ids" : None, ##Boolean. Assign one separate feature ID per CRER.
        "sort_distrib" : None, ##Boolean. Sort score distributions by decreasing values of significance.
        "bg_distrib_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "bg_distrib_string_type" : None, ##Type of information provided by the input string. Available values : text, url, piping
        "base" : None, ##Integer. Base for the logarithms (Default exp(1))
        "recursive" : None, ##Boolean. Run matrix-scan separately for each sequence.
        "batch" : None, ##Integer. Dispatch matrix-scan jobs on a cluster. Number of sequences to be analyzed by job (= on each node of the cluster)
        "lth_score" : 1, ##Number. Lower threshold on some parameters.
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
        "uth_pval" : 0.0001, ##Number. Upper threshold on some parameters.
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


## oligo-analysis

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/oligo_analysis_demo_seq.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "mask" : None, ##String. Mask lower or uppercases, respecively, i.e. replace selected case by N characters. Supported: upper, lower
        "format" : "fasta", ##String. Input sequence format. Various standards are supported - fasta, wconsensus, IG, filelist, raw
        "l" : 6, ##Integer. Oligonucleotide size.
        "expfreq_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "expfreq_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "calibN_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "calibN_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "calib1_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "calib1_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "bg" : "upstream-noorf", ##String. Type of sequences used as background model for estimating expected oligonucleotide frequencies. Supported: upstream, upstreamL, upstream-noorf, intergenic, input
        "org" : "Saccharomyces_cerevisiae", ##String. Organism name for background model. Only when bgfile is not specified
        "taxon" : None, ##String. Organism or taxon that used as reference for the estimation of a background model based on a genome subset (option -bg). Either -org or -taxon is required with the option -bg.
        "markov" : 2, ##Integer. Markov chain - the frequency expected for each word is calculated on basis of subword frequencies observed in the input set.
        "lexicon" : None, ##String. Expected word frequencies are calculated on the basis of subword frequencies, in a similar (but not identical) way to the “dictionary” approach developed by Harmen Bussemaker.
        "pseudo" : 0.01, ##Number. Pseudo-frequency for the background model, where pseudo must be a real value between 0 and 1.
        "noov" : True, ##Boolean. No overlapping.
        "2str" : True, ##Boolean. Oligonucleotide occurrences found on both stands are summed.
        "1str" : None, ##Boolean. Inactivates the summation of occurrences on both strands.
        "seqtype" : "dna", ##String. Input sequence type. Supported: dna, prot, other
        "return" : "occ,proba,rank", ##String. List of statistics to return. Supported - occ,mseq,freq,proba,ratio,zscore,like,pos,rank
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
        "quick" : True, ##Boolean. Quick count mode - delegate the counting of word occurrences to count-words, a program written in C by Matthieu Defrance.
        "quick_if_possible" : True, ##Boolean. Evaluate if the quick mode is compatible with the selected output parameters, and, if so, run in this mode.
        "accept_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "accept_string_type" : None, ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
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

## pattern-assembly

data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/convert-features_demo.pat", ##nput string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "2str" : True, ##Boolean. Oligonucleotide occurrences found on both stands are summed. Available values : true, false
        "1str" : None, ##Boolean. Inactivates the summation of occurrences on both strands. Available values : true, false
        "sc" : 8, ##Integer. Score column
        "cc" : None, ##Integer. Define a column containing cluster names or numbers.
        "maxfl" : 1, ##Integer. Maximum flanking segment size (default 1).
        "subst" : 1, ##Integer. Maximum allowed substitutions (default 0).
        "maxpat" : 100, ##Integer. Maximum number of allowed patterns (default 0).
        "toppat" : 100, ##Integer. Maximum number of patterns to assemble.
        "match" : None, ##Integer. Minimum number of matching residues to include a pattern in an assembly (default 0).
        "weight" : None, ##Number. Minimum matching weight to include a pattern in an assembly (default 0)
        "max_asmb_nb" : None, ##Integer. Maximal number of assemblies (default 5)
        "max_asmb_per_cluster" : None, ##Integer. Maximal number of assemblies per cluster (default 2).
        "max_asmb_size" : None, ##Integer. Maximal assembly size, i.e. the number of patterns per alignment group (default 50)
        "max_asmb_width" : None, ##Integer. Maximal width for an assembly (default 0)
        "single_sep" : None, ##Boolean. Report the isolated words (i.e. words that do not match any other words) separately.
    }

## peak-motifs

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat.france-bioinformatique.fr/rsat//demo_files/ChIP-seq_peaks/Oct4_peaks_top1000.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "ctrl_string" : "", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "ctrl_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "max_seq_len" : 6, ##Intger. msl. Maximal sequence length.
        "ref_motifs_string" : "http://rsat-tagc.univ-mrs.fr/rsat/motif_databases/JASPAR/Jaspar_2020/nonredundant/JASPAR2020_CORE_vertebrates_non-redundant_pfms.tf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "ref_motifs_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "motifs_db" : None, ##String. db_name db_format db_file. File containing a database of transcription factor biding motifs (e.g. JASPAR, TRANSFAC, RegulonDB, etc.). Ex. -motif_db TRANSFAC transfac transfac_download_dr/cgi-bin/data/matrix.dat
        "title" : None, ##String. graph_title. Title displayed on top of the graphs
        "str" : "-2str", ##String. Single-strand (-1str) or double-strand (-2str) analysis
        "img_format" : None, ##String. Image format
        "task" : "purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan", ##String. Default: purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,motifs_vs_db,scan
        "disco" : "oligos,positions", ##String. Specify the software tool(s) that will be used for motif discovery.
        "nmotifs" : 5, ##Integer. max_motif_number
        "minol" : None, ##Integer. Minimal lengths of oligonucleotide for word-counting approaches
        "maxol" : None, ##Integer. Maximal length of oligonucleotide for word-counting approaches
        "markov" : "auto", ##String. Order of the Markov model used to estimate expected oligonucleotide frequencies
        "min_markov" : None, ##Integer. Minimal value can be specified for the Markov order
        "max_markov" : None, ##Integer. Maximal value can be specified for the Markov order
        "noov" : None, ##Boolean. Treatment of self-overlapping words for motif discovery
        "r_plot" : None ##Boolean. Use R rather than the Perl GD library to generate plots.
    }


## position-analysis

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/dyad_analysis_demo_seq.fa", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "seqtype" : "dna", ##String. Input sequence type. Supported: dna, any.
        "last" : None, ##Integer. Stop after N sequences (for quick testing)
        "skip" : None, ##Integer. Skip the first N sequences.
        "first" : None, ##Integer. First sequence to analyze.
        "seqnb" : None, ##Integer. Maximal number of sequences to analyze.
        "mask" : None, ##String. Mask lower or uppercases, respecively, i.e. replace selected case by N characters. Supported: upper, lower.
        "format" : None, ##String. Input file format. Must be followed by one of the following options - fasta, wconsensus, IG, filelist, raw
        "l" : 6, ##Integer. Oligomer length.
        "ci" : None, ##Integer. Alphabet. window interval (default 20 bases).
        "origin" : None, ##String. Reference for calculating positions. Supported: start, center, end
        "offset" : None, ##Integer. Add an offset to site positions.
        "noov" : True, ##Boolean. No overlapping.
        "2str" : True, ##Boolean. Oligonucleotide occurrences found on both stands are summed.
        "1str" : None, ##Boolean. Inactivates the summation of occurrences on both strands.
        "grouprc" : True, ##Boolean. Group reverse complement with the direct sequence in the output file.
        "nogrouprc" : None, ##Boolean. Inactivates grouping of reverse complement pairs.
        "sort" : None, ##Boolean. Sort oligomers according to overrepresentation.
        "return" : "distrib,chi,rank,graphs,clusters", ##String. List of statistics to return. Supported - html,distrib,exp_occ,chi,rank,graphs,clusters
        "task" : "pos,matrices", ##String. Supported tasks - pos, clusters, matrices, graphs, index, all
        "markov" : 1, ##Integer. Order for the Markov model use to compute position-specific expected word frequencies.
        "max_graphs" : None, ##Integer. Maximal number of graphs to export
        "pl_string" : None, ##String. Input string specifying the query.
        "pl_string_type" : None, ##String. Available values : text, url, piping
        "sc" : None, ##String. Score column. (only valid whith the option -pl)
        "minpos" : None, ##Integer. Minimal position to take into account for the chi-square calculation This value must be a multiple of the window interval.
        "maxpos" : None, ##Integer. Maximal position to take into account for the chi-square calculation This value must be a multiple of the window interval.
        "max_asmb_per_cluster" : None, ##Integer.
        "nocheck" : False, ##Boolean. Do not check the applicability condition on the chi-square.
        "nofilter" : None, ##Boolean. Do not discard oligos which do not fit the condition of applicability
        "header" : None, ##String. Information to display in column headers of the distributions. Available values : mid, midfloor, min, max, interval
        "top_seq_for_matrices" : None, ##Integer. Select the top N sequences for building position-specific scoring matrices (PSSM).
        "img_format" : None, ##String. Image format (this parameter is passed to XYgraph).
        "title" : None, ##String. Title for the index table and position profile plots.
        "clust_method" : "average", ##String. Agglomeration rule for the hierarchical clustering. Supported - complete, average, single, ward
        "clust_nb" : None, ##Integer. Number of clusters (default 8).
        "clust_suffix" : None, ##String. Suffix to append to the cluster file and the directory contianing cluster graphics. Default 'clusters’
        "lth_chi" : 1, ##Number. Lower threshold on chi2
        "lth_sig" : None, ##Number. Lower threshold on significance
        "lth_occ" : 5, ##Integer. Lower threshold on occurrences
        "uth_rank" : None, ##Integer. Upper threshold on rank
}


## random-seq

data =  {
        "l" : 1000, ##Integer. Sequence length
        "n" : 10, ##Integer. Number of sequences
        "lw" : 50, ##Integer. Line width. A newline character will be inserted in the sequence every lw base. Default is 70. A value of 0 will prevent newline insertion.
        "type" : "DNA", ##String. Type of sequence(s) to generate. Protein, DNA, or other
        "seed" : 777 ##Integer. Number of seed for the random generator
    }

## retrieve-matrix

data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/motif_databases/JASPAR/Jaspar_2020/nonredundant/JASPAR2020_CORE_vertebrates_non-redundant_pfms.tf", ##String describing the matrix collection. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows). Example = http://rsat.sb-roscoff.fr/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.tf
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "v" : 1, ##Integer. Verbosity.
        "id" : "MA0265.1", ##Query IDs or ACs. The search is case-insensitive. Several queries can be provided, separated by commas without spaces, e.g. FOXA1,pou2f2.
        "id_file_string" : None, ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "id_file_string_type" : "text" ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content.
    }

## supported-organisms

data =  {
        "format" : "tab", ##Output format. supported tab, tree, html_tree, newick
        "return" : "ID,name,taxid,source,last_update,nb,seq_format,up_from,up_to,taxonomy,data,genome,genome_assembly,genome_version,download_date,variant_available,variant_source,path_to_variant_files", ##String. Output fields. Supported ID,name,taxid,source,last_update,nb,seq_format,up_from,up_to,taxonomy,data,genome,genome_assembly,genome_version,download_date,variant_available,variant_source,variant_date,path_to_variant_files
        "taxon" : None, ##String. Selected_taxon. Only returns organisms belonging to a selected taxon.
        "group" : "Fungi", ##String. Selected_group. Only returns organisms belonging to a selected taxonomy group.
        "source" : None, ##String. Selected_source. Only return organisms from a user-selected source(s).
        "depth" : None, ##Integer. Depth for exploring the taxonomic tree.
        "unique_species" : True, ##Boolean. Select at most one organism per species. This option aims at avoiding to be submerged by hundreds of strains sequenced for the same species.
        "unique_genus" : True, ##Boolean. Select at most one organism per genus.
        "source" : None ##String. Return the list of organisms supported on a remote RSAT server, via the Web services interfaces.
    }

## variation-info

data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs_IDs.txt", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "format" : "id", ##String. Input format. Supported formats = id, varBed, bed
        "org" : "Homo_sapiens_GRCh37", ##String. Query organism
        "type" : None, ##String. Specify one or several accepted types of variation. Supported types = SNV, insertion, deletion, substitution.
        "col" : None, ##Integer. Number of the column containing the variation IDs with the input format “id”
        "release" : None, ##Integer. The Ensembl release number of database (e.g. 84)
        "skip" : None, ##Integer. Skip the N first variations of the input file. This option is useful to run quick tests, or to split an analysis in separate tasks.
        "last" : None, ##Integer. Stop after the N first variations of the list. This option is useful to run0 quick tests, or to split an analysis in separate tasks.
    }

## variation-scan

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs.varseq", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "org" : "", ##Organism name for background model. Use only if bg is not specified
        "m_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/variation_demo_set_MWeirauch_cell_2014_15SNPs_TFs.tf", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "m_string_type" : "url", ##Type of information provided by the input string. Available values : text, url, piping.
        "top_matrices" : None, ##Integer. Only work with a given number of top matrices. This option is useful for quick tests and debugging.
        "top_variations" : None, ##Integer. Only work with a given number of top matrices. This option is useful for quick tests and debugging.
        "m_format" : "transfac", ##String. Matrix file format. Supported = transfac or tab.
        "markov_order" : 2, ##Integer. Alternative to the bg option. Automatically choose a background file on the server based on the species name and assembly.
        "bg_string" : "http://rsat.sb-roscoff.fr//demo_files/all_human_ENCODE_DNAse_mk1_bg.ol", ##String. Input string specifying the query.
        "bg_string_type" : "url", ##String. Available values : text, url, piping
        "uth_pval" : 0.001, ##String. Upper threshold for p value
        "lth_pval_ratio" : 10, ##Integer. Lower threshold for p value ratio
        "lth_score" : 1, ##String. Weight of predicted sites
        "lth_w_diff" : 1, ##String. Weight difference between variants
    }

## XYgraph

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : "http://rsat-tagc.univ-mrs.fr/rsat/demo_files/XYgraph_demo_data.ft", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server
        "i_string_type" : "url", ##Type of information provided by the input string (URL, piping, text)
        "format" : "png", ##String. Output image format. Supported: gif, eps, pdf, png, jpg
        "title1" : "XYgraph dmeo", ##String. first graph title. The title string should be embedded in double quotes if is contains spaces or special chars.
        "title2" : "score comparisons", ##String. Second graph title
        "xleg1" : "positional bias", ##String. First X legend.
        "xleg2" : "(chi-square value)", ##String. Second X legend.
        "yleg1" : "z-score", ##String. First Y legend.
        "yleg2" : "(Markov models, different orders)", ##String. Second Y legend.
        "xmax" : None, ##Integer. Maximal value represented on X axis.
        "ymax" : None, ##Integer. Maximal value represented on Y axis.
        "xmin" : 0, ##Integer. Minimal value represented on X axis.
        "ymin" : None, ##Integer. Minimal value represented on Y axis.
        "same_limits" : None, ##Boolean. Use same limits (max, min) for the X and Y axes.
        "min" : None, ##Integer. Minimal value represented on both X and Y axis (combinates the effects of -xmin and -ymin).
        "max" : None, ##Integer. Maximal value represented on both X and Y axis (combinates the effects of -xmin and -ymin).
        "xgstep1" : None, ##Integer. 1st step value for the grid across X axis.
        "ygstep1" : None, ##Integer. 1st step value for the grid across Y axis.
        "gstep1" : None, ##Integer. 1st step value for the grid across both X and Y axis.
        "xgstep2" : None, ##Integer. 2nd step value for the grid across X axis.
        "ygstep2" : None, ##Integer. 2nd step value for the grid across Y axis
        "gstep2" : None, ##Integer. 2nd step value for the grid across both X and Y axis.
        "xsize" : 500, ##Integer. Size of the X axis (in pixels). Default is 400.
        "ysize" : 400, ##Integer. Size of the Y axis (in pixels). Default is 400.
        "size" : None, ##Integer. Size of the X and Y axis (in pixels).
        "pointsize" : None, ##Integer. Point size (in pixels).
        "columns" : True, ##Boolean. Data fields are in columns (default)
        "rows" : True, ##Boolean. Data fields are in row
        "xcol" : "6", ##String. Column containing data for the X axis.
        "ycol" : "2,3,4,5", ##String. Column containing data for the Y axis.
        "xlog" : None, ##Integer. X data are displayed on a logarithmic scale
        "ylog" : None, ##Integer. Y data are displayed on a logarithmic scale
        "log" : None, ##Integer. Same as combining -xlog N -ylog M
        "lines" : None, ##Boolean. Points are jointed by lines.
        "force_lines" : None, ##Boolean. Points are jointed by lines.
        "line" : None, ##Integer. Same as lines, but for the Nth column only
        "header" : True, ##Boolean. First line of the data file contains a column header
        "legend" : True, ##Boolean. Use the content of the first line from input file as legend for Y data.
        "histo" : None, ##Boolean. Histogram. The X data should in this case contain the middle position of ach class, and Y data the frequencies
        "fhisto" : None, ##Boolean. filled histogram.
        "hbox" : None, ##String. Highlight box. the coordinates of the highlighted box in pixels (left, top, right, bottom respectively).
        "tbox" : None, ##String. Threshold box. units of X and Y data (low_x, high_x, low_y, high_y respectively).
        "bg" : None, ##String. Background color.
        "mono" : None, ##Boolean. Monochrome. All dots are drawn in black, and a specific symbol is associated to each.
        "htmap" : True, ##Boolean. An HTML document is automatically generated
        "lc" : None, ##Boolean. Label column
        "colors_string" : None, ##String. Input string specifying the query.
        "colors_string_type" : None, ##Type of information provided by the input string. Available values : text, url, piping
        "export_colors" : None, ##String. Color_file. Only working with eps and pdf formats
        "gp" : None, ##String. gnuplot additional commands. Ex. -gp ‘set size ratio 0.5’
        "hline" : None, ##String. Draw an horizontal line on the indicated position. Ex. -hline red 0 or -hline red 0,10,20,30
        "null" : None, ##Boolean. Indicate the NULL or NA character in the table to omit them.
        "vline" : None ##String. Draw a vertical line on the indicated position. Ex. -vline green 0 or -vline green 1,2,4,8,16
}

---
title: "[RSAT](RSAT_home.cgi) - matrix-clustering manual"
output:
  html_document:
    toc: yes
    toc_depth: 3
  pdf_document:
    toc: yes
    toc_depth: 3
css: course.css
---

### NAME NAME

matrix-clustering

## VERSION

$program\_version

## DESCRIPTION

Taking as input a set of position-specific scoring matrices, identify
clusters of similar matrices and build consensus motifs by merging the
matrices that belong to the same cluster.

## DEPENDENCIES

The clustering step relies on _MCL_, the graph-based clustering
algorithm developed by Stijn Van Dongen. MCL must be installed and its
path indicated in the RSAT configuration file
($RSAT/RSAT\_config.props). The installation of MCL can be done with a
RSAT makefile:

    cd $RSAT
    make -f makefiles/install_software.mk install_mcl

Various R packages are required in _matrix-clustering_ to convert the 
hierarchical tree into different output formats and to manipulate the
dendrogram which is exported.

    RJSONIO : http://cran.r-project.org/web/packages/RJSONIO/index.html
    ctc : http://www.bioconductor.org/packages/release/bioc/html/ctc.html
    dendextend : http://cran.r-project.org/web/packages/dendextend/index.html

For visualize the logo forest it is required the JavaScript _D3_ 
(Data Driven Documents) library, the user can select an option to connect
 directly with the server to load the functions of this library (see option _-d3\_base_). 

    D3 : http://d3js.org/

## AUTHORS

### Implementation

- Jacques.van-Helden@univ-amu.fr
- Jaime Castro <jcastro@lcg.unam.mx>

### Conception

- Jacques van Helden

    The following collaborator contributed to the definition of
    requirements for this program.

- Carl Herrmann
- Denis Thieffry
- Morgane Thomas-Chollier

## CATEGORY

util

## USAGE

matrix-clustering \[-i inputfile\] \[-o outputfile\] \[-v \] \[...\]

## OUTPUT FORMAT

## SEE ALSO

- _compare-matrices_

    The program _compare-matrices_ is used by _cluster-matrices_ to
    measure pairwise similarities and define the best alignment (offset,
    strand) between each pair of matrices.

## WISH LIST

## OPTIONS

- **-v #**

    Level of verbosity (detail in the warning messages during execution)

- **-h**

    Display full help message

- **-help**

    Same as -h

- **-i input matrix file**

    The input file contains a set of position-specific scoring
    matrices.

- **-matrix\_format matrix\_format**

    Specify the input matrix format.

    **Supported matrix formats**

    Since the program takes several matrices as input, it only accepts
    matrices in formats supporting several matrices per file (transfac,
    tf, tab, clusterbuster, cb, infogibbs, meme, stamp, uniprobe).

    For a description of these formats, see the help of _convert-matrix_.

- **-title graph\_title**

    Title displayed on top of the report page.

- **-display\_title**

    If it is selected. The title is displayed in the trees and in the result table.
    This is ideal when the user wants to compare motifs from different sources (files).

- **-root\_matrices\_only**

    When this option is selected. matrix-clustering returns a file with the 
    motifs at the root of each cluster. This save time and memory consumption because
    the branch-motifs, heatmaps, and trees are not exported.

- **-o output\_prefix**

    Prefix for the output files.

    Mandatory option: since the program _cluster-matrices_ returns a
    list of output files (pairwise matrix comparisons, matrix clusters).

- **-heatmap**

    Display consensus of merged matrices on the internal branches of the
    tree.

- **-export format**

    Specify format for the output tree.

    The hierarchical tree in JSON format is always exported, since it is required to display
    the logo tree with the d3 library. Additional formats are proposed in
    option to enable visualization with classical phylogeny analysis
    tools.

    **Supported trees formats**

    (JSON, newick)

    - _JSON_ (default)

        File format used for D3 library to visualize the logo forest in HTML.

    - _newick_ (optional)

        Widely used textual format to describe phylogenetic trees.

- **-task tasks**

    Specify one or several tasks to be run. If this option is not
    specified, all the tasks are run.

    Note that some tasks depend on other ones. This option should thus be
    used with caution, by experimented users only.

    Supported tasks: (all, comparison, clustering)

    - **all**

        Execute all the parts of the program (default)

    - **clustering**

        Skip the matrix comparison step and only executes the clustering step.

        Assumes the users already have the description table and comparison table 
        exported from the program _compare-matrices_.

        This option is ideal to saving time once all comparison beteen the input motifs had been done. 

- **-label**

    Option to select the matrix label fields displayed in the html tree

    **Supported labels**

        (name, consensus, id)
        

- **-quick**

    With this option the motif comparison step is done with  the program _compare-matrices-quick_ 
    (implemented in C) rather than the classic version compare-matrices (implemented in Perl).
    The quick version runs x100 times faster, but has not all implemented options as in the Perl version.

    We suggest use this option for a big set of input motifs > 300 motifs. 

    **NOTE:** By the moment the only threshold used in quick version is Ncor. 

- **-clone\_input**

    If this option is selected, the input motif database is exported
    in the results folder.

    NOTE: take into account the input file size. 

- **-max\_matrix**

    This option specify how many matrices can be clustered in the same analysis. If there are more matrices than
    the specified number, the program reports an error.

    This parameter can be useful when the user analyse a big dataset of matrices.

- **-hclust\_method**

    Option to select the agglomeration rule for hierarchical clustering.

    Supported agglomeration rules:

    - _complete_ (default)

        Compute inter-cluster distances based on the two most distant nodes.

    - _average_

        Compute inter-cluster distances as the average distance between nodes
        belonging to the relative clusters.

    - _single_

        Compute inter-cluster distances based on the closest nodes.

- **-top X**

    Only analyze the first X motifs of the input file. This options is
    convenient for quick testing before starting the full analysis.

- **-skip X**

    Skip the first X motifs of the input file. This options is convenient
    for testing the program on a subset of the motifs before starting the
    full analysis.

- **-consensus\_labels**

    Option to select the labels displayed in the consensus
    alignment picture

    Default: consensus, id, strand

    **Supported labels**

        (consensus, id, strand, number)

- **-lth param lower\_threshold**
- **-uth param upper\_threshold**

    Threshold on some parameter (-lth: lower, -uth: upper threshold).

    Threshold parameters are passed to compare-classes. 

    In addition, if a threshold is defined in the (unique) metrics used as
    clustering score (option _-score_), this threshold will be used to
    decide whether motifs should be aligned or not. If two motifs have a
    similarity score lower (or distance score higher) than the selected
    threshold, their aligment will be skipped. The status of each motif
     (Aligned or Non-aligned) is reported in the file
    prefix\_matrix\_alignment\_table.tab

    Suggested thresholds:

        cor >= 0.7

        Ncor >= 0.4

- **-score metric**

    Select the metric which will be used to cluster the motifs.

    Supported metrics : cor, Ncor

    Default: Ncor 

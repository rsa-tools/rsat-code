/*************************************************************************/
/* Authors     : Jerome Callut (jerome.callut@uclouvain.be)              */
/*               Pierre Dupont (pierre.dupont@uclouvain.be)              */
/* Institution : Departement d'Ingenierie Informatique                   */
/*               Universite catholique de Louvain                        */
/* Ref.        : P. Dupont, J. Callut, G. Dooms, J.-N. Monette, and      */
/*               Y. Deville. Relevant subgraph extraction from random    */
/*               walks in a graph. Technical Report RR 2006-07, INGI,    */
/*               UCL, 2006.                                              */
/* Version     : 1.00                                                    */
/* Modified    : 29/08/2008                                              */
/* License     : GPL 3.0                                                 */ 
/*               GNU GENERAL PUBLIC LICENSE                              */
/*               Version 3, 29 June 2007                                 */
/*               For details, see the file gpl-3.0.txt in this bundle    */
/*************************************************************************/

1. DIRECTORY CONTENT
2. INSTALLATION
3. COMMAND LINE SYNTAX
4. DATA FORMAT
5. MATLAB FUNCTION
6. REFERENCES


1. DIRECTORY CONTENT
--------------------

kwalks/src      : the directory containing the C source files
kwalks/bin      : the directory containing the binary file and scripts
kwalks/matlab   : the directory containing the matlab sources for unbounded walks
kwalks/examples : the directory containing some example graphs

2. INSTALLATION
---------------

This piece of software is developped to run under UNIX systems (Linux,OS X, ...).
To install the lkwalk tool, go in the directory src/ and type "make".
This will generate the binary file bin/lkwalk. Then, you should add the bin/ 
directory to your $PATH variable.

3. COMMAND LINE SYNTAX
----------------------

Usage     : lkwalk -g graphFile -l walkLen -k relNodes -o outFile [options]

Function  : Compute the limited K-walk Expected Passage Times
Arguments : -g graphFile  the file containing the (sparse) adjacency matrix
            -l walkLen    the length of the limited k-walk
            -k relNodes   the list of nodes of interest
                          groups are separated by # and nodes are separated by :
            -o outFile    the output file containing the expected times
Options   : -u            consider walks up to length walkLen
            -p relProbas  initial probas for nodes of interest
                          groups are separated by # and nodes are separated by :
            -d            generate a dot file (outFile.dot)
            -v verbose    verbosity level (default 1)
            -h            display this help

Version   : 1.00
Author    : J. Callut and P. Dupont (2007)

Example   : lkwalk -g examples/graph10_bin_sp -l 5 -k 1:2#3:5 -o out

Remark    : The node indexing starts at 1

4. DATA FORMAT
--------------

4.1 Input

The input *connected* graph must be provided as sparse adjacency matrix in a text
file. The first line indicates the total number of nodes in the graph. Each next
line corresponds to a graph node. The first entry indicates the node index l and
the next entries are pairs c:w meaning that the edge l->c has a weight w.
Lines must be provided in increasing order with respect to the node index l.
In each line, the node indices c must appear in increasing order.
For instance, here follows a file for an unweighted graph having 10 nodes.

10
1 2:1.0 3:1.0 6:1.0 8:1.0 9:1.0 10:1.0
2 1:1.0 3:1.0 10:1.0
3 2:1.0 3:1.0 7:1.0 8:1.0 9:1.0 10:1.0
4 1:1.0 2:1.0 3:1.0 5:1.0 6:1.0
5 2:1.0 4:1.0
6 1:1.0 3:1.0 4:1.0 5:1.0 6:1.0 8:1.0 10:1.0
7 1:1.0 2:1.0 4:1.0 5:1.0 10:1.0
8 1:1.0 2:1.0 5:1.0 6:1.0 9:1.0 10:1.0
9 1:1.0 4:1.0 6:1.0 7:1.0
10 2:1.0 3:1.0 4:1.0 5:1.0 6:1.0 8:1.0

4.2 Outputs

The program outputs 4 files : outFile.N outFile.E [outFile.dif] [outFile.dot]
1. outFile.N contains the Expected Node Passage Times as a sparse vector of floats.
2. outFile.E contains the Expected Edge Passage Times as a sparse matrix of floats.
3. outFile.dif contains a sparse symmetric matrix with the absolute Expected Edge Passage
   Times difference *IF* the input graph is undirected otherwise no file is provided.
4. outFile.dot contains a graphviz dot file for vizualization. If the graph is
   undirected the edges weights are the absolute Expected Edge Passage Times difference,
   otherwise the edges weights are the normalized Expected Edge Passage Times.

Additional files with the .norm extension are simply normalized by the expected
walk length.
   
5. MATLAB FUNCTION
------------------

In addition of the C module computing limited K-walks, a matlab function
for computing standard K-walks is also provided. The syntax is the following

graph=readMat_sparse(fname);
[N,E]=kwalk(graph,knodes,kprobas);

fname   : the file name containing the sparse adjacency matrix as in 4.1
knodes  : an array containing the relevant nodes
kprobas : an array containing initial probas for the nodes of interest (optional)
N       : an array containing the Expected Node Passage Times
E       : a sparse matrix containing the Expected Edge Passage Times

example : [N,E]=kwalk('../examples/graph10_bin_sp',[1 2]);

Remark : Specifying the nodes of interrest by group is not yet supported.

This module can also be called from the command line using the scripts/kwalk
Python script. Running this script requires to have the Python 2.4, or later,
interpreter installed. The matlab/ directory has to be in your matlab path and
you have to configure the paths by adapting the following lines in scripts/kwalk: 

	MATLAB_CMD  = 'matlab -nojvm -nosplash -nodesktop'
	BASE_PATH   = '~/UCL/kwalks/matlab/'

Usage: kwalk -g graphFile -k relNodes -o outFile [-p relProbas]

(see section 3 for details about the arguments and section 4 about the inputs/outputs)


6. REFERENCES
-------------

[1] P. Dupont, J. Callut, G. Dooms, J.-N. Monette, and Y. Deville. Relevant 
subgraph extraction from random walks in a graph. Technical Report RR 
2006-07, INGI, UCL, 2006. 




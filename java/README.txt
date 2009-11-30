%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% README
%
% Graphtools command line - version 1.1
%
% Author: Karoline Faust
%
% Last modification: 26/03/2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


1. Installation
===============

For the installation of the Graphtools, see the RSAT install guide in the section "Adding RSAT to your path".

For the installation of kWalks and REA, see the RSAT install guide in the section "Installing third-party programs".

2. Use
=======

2.1. List available tools
---------------------------

You can list available tools by typing:
java graphtools.util.ListTools

2.2 Getting help
----------------

All tools provide a -h option to display help.

2.3. Simplify names using aliases
---------------------------------

The command line tool names may be simplified by setting aliases.
Example:
alias Pathfinder="java graphtools.algorithms.Pathfinder"

allows to type:
Pathfinder -h
instead of:
java graphtools.algorithms.Pathfinder -h

2.4. Increasing memory
----------------------

For large graphs, you may need to increase the memory allocated to the java virtual machine.

You can do so by specifying the -Xmx option.

Example:
java -Xmx800m graphtools.algorithms.Pathfinder -h

2.5. Update data files
----------------------

The java web applications depend on four
data files. All these files can be updated using the graphtools.

1) Kegg_organisms_list.txt
Usage:
This file is displayed in the KEGG network provider web interface.
Location:
The updated file should replace the file located at $RSAT/data/KEGG.
Update:
A list of organisms supported by the current KEGG PATHWAY version
can be obtained using the -O option of the MetabolicGraphProvider,
which queries KEGG via its API.

2) rpairs.tab
Usage:
This file is used by the MetabolicGraphProvider and its client, the
KEGG network provider, to generate RPAIR networks.
Location:
The updated file should replace the file located at $RSAT/data/KEGG.
Update:
The recent version of the KEGG LIGAND reaction file can be downloaded
and the list of rpairs can be generated from it using the -p option of the KeggLigandDataManager.

3) metabolicdb_dump_day_month_year.backup
Usage:
This file contains the metabolic database, which is used by the MetabolicGraphProvider and
the metabolic Pathfinder web application.
Location:
The file can be downloaded from the data section of the official NeAT website:
http://rsat.bigre.ulb.ac.be/neat/
It should be loaded into a local postgres database. See the NeAT web server install guide for details.
Update:
Current KEGG LIGAND data can be parsed into a postgres metabolic database using
the KeggLigandDataManager options -A, -l and -z.

4) preloadedNetworks.tgz
Usage:
This file contains the networks that are used by the metabolic Pathfinder web application.
Location:
The updated networks should replace the networks located at $RSAT/data/Stored_networks.
Update:
The preloaded networks for the metabolic pathfinder can be generated with
the KeggLigandDataManager options -M, -S and -T.

2.6. Tests
----------

************* Path finding *************

You can try the following to test your installation:
java graphtools.algorithms.Pathfinder -g $RSAT/public_html/demo_files/yeast_string_database_graph_converted_weights.tab -f flat -s RAS2 -t TEC1 -r 2

After a few seconds (depending on your computing power), the following should be listed:

#start node     end node        path index      rank    distance        steps   path
RAS2    TEC1    1       1       6.25    6       RAS2    CDC42   STE20   FUS3    STE12   TEC1
RAS2    TEC1    2       2       7.5     7       RAS2    CDC42   STE20   STE5    FUS3    STE12   TEC1
RAS2    TEC1    3       2       7.5     7       RAS2    CDC42   STE20   STE11   FUS3    STE12   TEC1
RAS2    TEC1    4       2       7.5     7       RAS2    CDC42   SHO1    STE20   FUS3    STE12   TEC1
RAS2    TEC1    5       2       7.5     7       RAS2    CDC42   STE20   FUS3    DIG2    STE12   TEC1
RAS2    TEC1    6       2       7.5     7       RAS2    CDC42   STE20   FUS3    DIG1    STE12   TEC1
RAS2    TEC1    7       2       7.5     7       RAS2    CDC42   BEM1    STE20   FUS3    STE12   TEC1
RAS2    TEC1    8       2       7.5     7       RAS2    CDC42   STE20   STE7    FUS3    STE12   TEC1
RAS2    TEC1    9       9       8.75    8       RAS2    CDC42   STE20   STE5    STE11   FUS3    STE12   TEC1

If it takes too long, specify a temp directory (option -p). This will greatly speed up path finding the next time you run it on the same graph.

************* Pathway inference *************

Type the following:
java graphtools.algorithms.Pathwayinference -g $RSAT/public_html/demo_files/Ecoli_metabolic_network.tab -f flat -b -s 'R02412>/R02412<#R00986>/R00986<#R01715>/R01715<' -y con -o aromatic_amino_acid_biosynthesis.txt

After a few seconds, the following should be listed in aromatic_amino_acid_biosynthesis.txt (the order of nodes and edges may vary):

; export date: Fri Oct 31 15:02:45 CET 2008
; graph id: path_data_0
; comments: ; Date=Fri Oct 31 15:02:45 CET 2008
; INPUT
; ===============================
; Seed group name(s)=seeds
; Seeds:
; ; members	groups	weights
; R02412>	R02412>G	NaN
; R02412<	R02412>G	NaN
; R00986>	R00986>G	NaN
; R00986<	R00986>G	NaN
; R01715>	R01715>G	NaN
; R01715<	R01715>G	NaN
; R01715>G	seeds	NaN
; R02412>G	seeds	NaN
; R00986>G	seeds	NaN

; Graph=Ecoli_metabolic_network.tab
; Directed=false
; Metabolic=true
; CONFIGURATION
; ===============================
; Algorithm=hybrid
; Weight policy=con
; Inflation factor=1.0
; Exclusion attribute=ExclusionAttribute
; Seed preprocessing=false
; Pathway postprocessing=false
; Hybrid algorithm: subgraph percentage=0.0050
; Iteration number=1
; KWalks algorithm: pruning=false
; PROCESSING
; ===============================
; Runtime in ms=5143
; All seeds connected by preprocessing=false
; Edges introduced by preprocessing=[]
; Edges introduced by postprocessing=[]

;NODES	rgb_color	color	ObjectType	ExclusionAttribute
R01715>	#0000FF	blue	Reaction	R01715
C01269	#FFD700	gold	Compound	C01269
R01715<	#0000FF	blue	Reaction	R01715
C03175	#FFD700	gold	Compound	C03175
R01714>	#FFD700	gold	Reaction	R01714
R03460>	#FFD700	gold	Reaction	R03460
R01714<	#FFD700	gold	Reaction	R01714
R03460<	#FFD700	gold	Reaction	R03460
R00986>	#0000FF	blue	Reaction	R00986
R02412>	#0000FF	blue	Reaction	R02412
R00986<	#0000FF	blue	Reaction	R00986
R02412<	#0000FF	blue	Reaction	R02412
C00251	#FFD700	gold	Compound	C00251
;ARCS	rgb_color	color
R01714>	C01269	#FFD700	gold
R01714<	C01269	#FFD700	gold
R02412>	C03175	#FFD700	gold
R02412<	C03175	#FFD700	gold
R01715>	C00251	#FFD700	gold
R01715<	C00251	#FFD700	gold
C00251	R01714>	#FFD700	gold
C00251	R01714<	#FFD700	gold
C03175	R03460>	#FFD700	gold
R00986<	C00251	#FFD700	gold
C03175	R03460<	#FFD700	gold
C00251	R00986>	#FFD700	gold
R03460>	C01269	#FFD700	gold
R03460<	C01269	#FFD700	gold

2.7. Tool-specific remarks
--------------------------

2.7.1. MetabolicGraphProvider

If MetabolicGraphProvider is run in KGML modus, it downloads KGML files into the $RSAT/public_html/data/KEGG directory.
In the same directory, it also creates a folder to store organism-specific metabolic graphs for later use.

The rpairs.tab file required for some options can be found in $RSAT/java/misc
In the same directory, the lists of supported KEGG organisms (for KGML modus) are stored.

The options -a, -r, -R and the DB modus of MetabolicGraphProvider are only available
if the metabolic database has been installed.
Check the NeAT web server install guide for installation of postgres and the metabolic database.

2.7.2. Pathwayinference

The program "compSteiner" relying on Klein and Ravi's algorithm and implemented by Nadja Betzler is not freely distributed,
thus algorithms 'compSteiner' and 'steinerhybrid' are not available. All other algorithms are available with this distribution.

2.7.3. GraphAnnotator

The GraphAnnotator tool is only functional if the metabolic database has been installed.

3. Disclaimer
=============

REA has been developed by Jimenez and Marzal [1]. The Java wrapper around REA has been inspired by a Python wrapper by Jean-Noel Monette and Pierre Schaus (INGI - UCL [2]).
Jean-Noel Monette and Pierre Schaus also improved the REA code and made this modification available.
The multiple-to-multiple end extension of path finding relies on a graph transformation suggested by Olivier Hubaut (former aMAZE team member) and also mentioned in [3].
The backtracking algorithm has been developed by Fabian Couche and Didier Croes [4,5].
The graph library used is the work of the former aMAZE team [6].
KWalks has been developped mainly by Pierre Dupont and Jerome Callut [7].

MetabolicGraphProvider makes use of data stored in KEGG [8].
For the use of KEGG data, see: http://www.genome.jp/kegg/legal.html

The postgres backup file metabolicdb_dump_day_month_year.backup contains data taken
from KEGG and MetaCyc [9].
For the use of MetaCyc data, see: http://biocyc.org/download-flatfiles.shtml

For licenses of individual jar files see the licenses directory.

4. Literature/Links
===================

[1] V. Jimenez and A. Marzal (1999) Computing the K Shortest Paths: a New Algorithm and an Experimental Comparison, 3rd Workshop on Algorithm Engineering
[2] http://www.uclouvain.be/ingi.html
[3] Duin, C.W., Volgenant, A., and Vo§, S. (2004). "Solving group Steiner problems as Steiner problems." European Journal of Operational Research 154, 323-329.
[4] D. Croes, F. Couche, S. Wodak and J. van Helden (2005) Metabolic PathFinding: inferring relevant pathways in biochemical networks. Nucleic Acids Research 33: W326-W330.
[5] D. Croes, F. Couche, S. Wodak and J.van Helden (2006) Inferring Meaningful Pathways in Weighted Metabolic Networks.
[6] http://www.scmbb.ulb.ac.be/amaze/
[7] P. Dupont, J. Callut, G. Dooms, J.-N. Monette and Y. Deville (2006) Relevant subgraph extraction from random walks in a graph . Research Report UCL/FSA/INGI RR 2006-07.
[8] Kanehisa, M., Araki, M., Goto, S., Hattori, M., Hirakawa, M., Itoh, M., Katayama, T., Kawashima, S., Okuda, S., Tokimatsu, T., and Yamanishi, Y.; KEGG for linking genomes to life and the environment. Nucleic Acids Res. 36, D480-D484 (2008).
[9] Caspi, R., Foerster, H., Fulcher, C. A., Kaipa, P., Krummenacker, M., Latendresse, M., Paley, S., Rhee, S. Y., Shearer, A. G., Tissier, C., Walk, T. C., Zhang, P. & Karp, P. D. (2008). The MetaCyc Database of metabolic pathways and enzymes and the BioCyc collection of Pathway/Genome Databases. Nucleic Acids Res36, D623-D631.
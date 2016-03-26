This folder contains some sample files to generate flowcharts with graphviz. 

The information is contained in dot-formatted files (graphviz
description of graphs). The layout and figures are generated with graphviz. 

The makefile provides the basic rules to generate figures in various
formats (pdf, png, eps, svg).


The simplest way to generate the graphs is to run the command:

    make all


Files
=====

test.dot

	Illustration of the simplest way to describe a directed graph
	in dot format.


template_flowchart.dot

	A more elaborate example with specific box shapes. 

rna-seq_anlaysis.dot

	A conceptual flowchart of the generic operations that can be
	generated with graphviz.

rsat_toolmap.dot

	A very fragmentary tool map of some tools from the software
	suite Regulatory Sequence Analysis Tools (RSAT,
	http://rsat.eu/). This illustrates the use of clickable links
	in the pdf file.

regulatory_variations.dot

	An example of specific tool map for the analysis of regulatory
	variants, with links to the Web forms of the resources used.



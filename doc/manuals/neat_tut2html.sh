#!/bin/bash
#
# Task: Compile Neat tutorial html pages from latex files using latex2html
# Author: Karoline Faust
# Date: 17/01/2008

OUTDIR=../../public_html/tutorials/neat_tutorial
STYLEFILE=../../public_html/main_grat.css

# compile tex files and create library
latex neat_tutorial.tex
bibtex neat_tutorial.aux
latex neat_tutorial.tex
latex neat_tutorial.tex
# convert tex to bib
latex2html neat_tutorial.tex -split 3 -dir $OUTDIR -long_titles 3 -mkdir
# make sure we have correct style
cp $STYLEFILE $OUTDIR/neat_tutorial.css
# remove temp files
rm neat_tutorial.aux
rm neat_tutorial.bbl
rm neat_tutorial.dvi
rm neat_tutorial.blg
rm neat_tutorial.idx
rm neat_tutorial.log
rm neat_tutorial.toc
rm rsat_latex_commands.aux

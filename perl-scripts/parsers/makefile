############################################################
#
# $Id: makefile,v 1.7 2002/07/02 16:32:40 jvanheld Exp $
#
# Time-stamp: <2002-07-02 18:31:06 jvanheld>
#
############################################################

RSYNC_OPT=-e ssh
RSYNC = rsync -ruptvl ${RSYNC_OPT}

PARSER		= .
RUBENS_PARSER	= jvanheld@rubens.ulb.ac.be:/rubens/dsk4/jvanheld/amaze/amaze_programs/amaze_parser

PARSED_DATA	= /win/amaze/amaze_data/parsed_data/
RUBENS_PARSED_DATA = jvanheld@rubens.ulb.ac.be:/rubens/dsk4/jvanheld/amaze/parsed_data

OPT	= -v 1

usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^(\S+):/){ print "\t$$1\n" }' makefile

all:
	make kegg
	make genes
	make swissprot
	make glaxo
	make gaurab
	make rsync_to_rubens

swissprot:
	perl parse_swissprot.pl -data swissprot -data trembl ${OPT}

kegg:
	perl parse_kegg.pl ${OPT}

genes:
	perl parse_genes.pl -org S.cerevisiae -org E.coli -org H.sapiens -suffix "" ${OPT}

glaxo:
	perl parse_glaxo.pl ${OPT}

gaurab:
	perl parse_gaurab.pl ${OPT}

skeletons:
	perl parse_pathway_skeleton.pl ${OPT}

rsync_to_rubens:
	${RSYNC} ${PARSER} ${RUBENS_PARSER}

#DELIVERY_DIR=glaxo_parsed
#DELIVERY_DIR=gaurab_parsed
DELIVERY_DIR=kegg_ligand
DELIVERY_DATE=`date +%Y%m%d`
deliver:
	${RSYNC} --delete ${PARSED_DATA}/${DELIVERY_DIR}/${DELIVERY_DATE} ${RUBENS_PARSED_DATA}/${DELIVERY_DIR}

#### TEMPORARY ####
aliases:
	alias stats='zgrep field Process.tab.gz; for col in 6 7 8 9 10 11 12; do echo $${col} `gunzip -c Process.tab.gz | grep -v "^--" | cut -f $${col} | classfreq -ci 1 -v | grep sum`; done'

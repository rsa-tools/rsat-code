############################################################
#
# $Id: makefile,v 1.6 2002/07/02 13:27:10 jvanheld Exp $
#
# Time-stamp: <2002-07-02 15:24:44 jvanheld>
#
############################################################

RSYNC_OPT=-e ssh
RSYNC = rsync -ruptvl ${RSYNC_OPT}

PARSER		= .
PARSED_DATA	= /win/amaze/amaze_data/parsed_data/

SPICA		= 193.62.198.80
MOZART_DIR	= jvanheld@mozart.ebi.ac.uk:aMAZE_team/jvanheld/
SPICA_DIR	= jvanheld@${SPICA}:~/

RUBENS_DELIVERY	= jvanheld@rubens.ulb.ac.be:/rubens/dsk4/jvanheld/amaze/parsed_data
RUBENS_PARSER	= jvanheld@rubens.ulb.ac.be:/rubens/dsk4/jvanheld/amaze/amaze_programs/amaze_parser

OPT	= -v 1

usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^(\S+:)/){ print "\t"; print }' makefile

all:
	make kegg
	make genes
	make swissprot
	make glaxo
	make gaurab
	make rsync_to_mozart
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

rsync_to_mozart:
	${RSYNC} --delete ${PARSER} ${MOZART_DIR}/programs/parsing/amaze_parsers

rsync_to_rubens:
	${RSYNC} ${PARSER} ${RUBENS_PARSER}

#DELIVERY_DIR=glaxo_parsed
#DELIVERY_DIR=gaurab_parsed
DELIVERY_DIR=kegg_parsed
DELIVERY_DATE=`date +%Y%m%d`
deliver:
	${RSYNC} --delete ${PARSED_DATA}/${DELIVERY_DIR}/${DELIVERY_DATE} ${RUBENS_DELIVERY}/${DELIVERY_DIR}

#### TEMPORARY ####
aliases:
	alias stats='zgrep field Process.tab.gz; for col in 6 7 8 9 10 11 12; do echo $${col} `gunzip -c Process.tab.gz | grep -v "^--" | cut -f $${col} | classfreq -ci 1 -v | grep sum`; done'

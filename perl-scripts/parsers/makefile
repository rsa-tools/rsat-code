############################################################
#
# $Id: makefile,v 1.4 2000/12/12 02:52:01 jvanheld Exp $
#
# Time-stamp: <2000-12-12 03:48:51 jvanheld>
#
############################################################

PARSER		= .
PARSED_DATA	= /win/amaze/amaze_data/parsed_data/

SPICA		= 193.62.198.80
BACH_DIR	= jvanheld@bach.ebi.ac.uk:aMAZE_team/jvanheld/
SPICA_DIR	= jvanheld@${SPICA}:~/

RUBENS_DELIVERY	= jvanheld@rubens.ulb.ac.be:/rubens/dsk3/genomics/delivery/internal/
RUBENS_PARSER	= jvanheld@rubens.ulb.ac.be:/rubens/dsk3/genomics/amaze/aMAZE_team/jvanheld/programs/parsing/amaze_parsers
RUBENS_PARSED_DATA	= jvanheld@rubens.ulb.ac.be:/rubens/dsk3/genomics/amaze/aMAZE_team/jvanheld/parsed_data

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
	make rsync_to_bach
	make rsync_to_rubens

swissprot:
	perl parse_swissprot.pl -w 1 -data swissprot -data trembl

kegg:
	perl parse_kegg.pl -w 1

genes:
	perl parse_genes.pl -w 1 -org S.cerevisiae -org E.coli -org H.sapiens -suffix "" 

glaxo:
	perl parse_glaxo.pl -w 1

gaurab:
	perl parse_gaurab.pl -w 1

rsync_to_bach:
	rsync -ruptvl --delete -e ssh ${PARSER} ${BACH_DIR}/programs/parsing/amaze_parsers
	rsync -ruptvl --delete -e ssh ${PARSED_DATA} ${BACH_DIR}

rsync_to_rubens:
	rsync -avl -e ssh ${PARSER} ${RUBENS_PARSER}
	rsync -avl -e ssh ${PARSED_DATA} ${RUBENS_PARSED_DATA}

deliver:
	rsync -ruptvln -e ssh ${PARSED_DATA}/kegg_parsed/20001212 ${RUBENS_DELIVERY}/kegg_parsed
	rsync -ruptvln -e ssh ${PARSED_DATA}/swissprot_parsed/20001212 ${RUBENS_DELIVERY}/swissprot_parsed


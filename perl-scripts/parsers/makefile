############################################################
#
# $Id: makefile,v 1.1 2000/11/21 08:29:42 jvanheld Exp $
#
# Time-stamp: <2000-11-21 09:29:18 jvanheld>
#
############################################################

PARSER		= .
PARSED_DATA	= /win/amaze/amaze_data/parsed_data

SPICA		= 193.62.198.80
BACH_DIR	= jvanheld@bach.ebi.ac.uk:aMAZE_team/jvanheld/
SPICA_DIR	= jvanheld@${SPICA}:~/
BXL_DIR		= jvanheld@rubens.ulb.ac.be:/rubens/dsk3/genomics/amaze/aMAZE_team/jvanheld/

usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^(\S+:)/){ print "\t"; print }' makefile

all:
	make ligand
	make glaxo
	make genes
	make gaurab
	make rsync_to_bach
	make rsync_to_rubens

ligand:
	perl parse_ligand.pl -w 1

genes:
	perl parse_genes.pl -w 1

glaxo:
	perl parse_glaxo.pl -w 1

gaurab:
	perl parse_gaurab.pl -w 1

rsync_to_bach:
	rsync -avl -e ssh ${PARSER} ${BACH_DIR}/programs/parsing/amaze_parsers
	rsync -avl --delete -e ssh ${PARSED_DATA} ${BACH_DIR}

rsync_to_spica:
	rsync -avl --delete -e ssh ${PARSED_DATA} ${SPICA_DIR}

rsync_to_rubens:
	rsync -avl -e ssh ${PARSER} ${BXL_DIR}/programs/parsing/amaze_parsers
	rsync -avl -e ssh ${PARSED_DATA} ${BXL_DIR}


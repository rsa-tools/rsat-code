############################################################
#
# $Id: makefile,v 1.10 2011/11/03 03:27:54 jvanheld Exp $
#
# Time-stamp: <2002-12-23 15:49:11 jvanheld>
#
############################################################

################################################################
#
# Protocol for parsing all the data required to load amaze
#

#### programs
RSYNC_OPT=-e ssh
RSYNC = rsync -ruptvl ${RSYNC_OPT}
PERL=perl

################################################################
#### parameters
################################################################


SERVER=paulus.ulb.ac.be

# directory containing the parser scripts
#PARSER=/rubens/dsk2/jvanheld/amaze_programs/amaze_parsers
PARSER=/win/amaze/amaze_programs/amaze_parsers

# Root of the directory where the parsed data will be saved.
PARSED_ROOT=${HOME}/parsed_for_amaze
PARSED_DATA	= /win/amaze/amaze_data/parsed_data/
SCMB_PARSED_DATA = jvanheld@${SERVER}:/rubens/dsk4/jvanheld/amaze/parsed_data

PARSING_DATE = `date +%Y%m%d`

# When the makefile is started with the option TEST=-test, only a
# small ssubset of the data file is parsed.
TEST=

# level of verbosity
VERBOSE=1

# general options for all the parsers
OPT=
OPTIONS=-clean -parseddir ${PARSED_ROOT} -v ${VERBOSE} ${TEST} ${OPT}

################################################################
#### target definitions
################################################################

usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^(\S+):/){ print "\t$$1\n" }' makefile

all: parse_all load_all

parse_all: parse_kegg parse_genes parse_swissprot parse_regulation parse_pathway_skeletons

load_all: load_kegg load_genes load_swissprot load_regulation load_pathway_skeletons


################################################################
#
# Polypeptides are parsed from SWALL, i.e. the non-redundant
# compilation of Swiss-prot and TREMBL. For this, the two databases
# have to be specified separately with the option -data
SWISSPROT_ORGANISMS=-org yeast -org ecoli -org human
parse_swissprot:
	perl ${PARSER}/parse_swissprot.pl	\
		-data swissprot -data trembl	\
		${SWISSPROT_ORGANISMS}			\
		${OPTIONS}

SWISSPROT_PARSED=${PARSED_ROOT}/swissprot/${PARSING_DATE}
load_swissprot:
	(cd ${SWISSPROT_PARSED}/sql_scripts ;			\
	make -f polypeptide.mk drop; make -f polypeptide.mk all	\
	)

################################################################
#
# Compounds, reactions, EC numbers, generic and organism-specific
# pathways from KEGG
parse_kegg:
	perl ${PARSER}/parse_kegg.pl ${OPTIONS}

KEGG_PARSED=${PARSED_ROOT}/kegg_ligand/${PARSING_DATE}
load_kegg:
	(cd ${KEGG_PARSED}/sql_scripts ;				\
	make -f compound.mk drop; make -f compound.mk all;		\
	make -f reaction.mk drop; make -f reaction.mk all;		\
	make -f reactant.mk drop; make -f reactant.mk all;		\
	make -f ecset.mk drop; make -f ecset.mk all;			\
	make -f genericpathway.mk drop; make -f genericpathway.mk all;	\
	make -f pathway.mk drop; make -f pathway.mk all			\
	)

################################################################
#
# Genes from KEGG. 
#
KEGG_ORGANISMS= -org S.cerevisiae -org E.coli -org H.sapiens
parse_genes:
	perl ${PARSER}/parse_genes.pl ${KEGG_ORGANISMS} ${OPTIONS}

KEGG_GENES_PARSED=${PARSED_ROOT}/kegg_genes/${PARSING_DATE}
load_genes:
	(cd ${KEGG_GENES_PARSED}/sql_scripts ;		\
	make -f gene.mk drop; make -f gene.mk all	\
	)

################################################################
#
# Metabolic regulation annotated by aMAZE team
#
parse_regulation:
	perl ${PARSER}/parse_regulation.pl ${OPTIONS}

REGULATION_PARSED=${PARSED_ROOT}/regulation/${PARSING_DATE}
load_regulation:
	(cd ${REGULATION_PARSED}/sql_scripts ;				\
	make -f transcriptregul.mk drop; make -f transcriptregul.mk all	\
	)

################################################################
#
# Pathway skeletons annotated by aMAZE team
#
parse_pathway_skeletons:
	perl ${PARSER}/parse_pathway_skeletons.pl -seeds -mirror -d ${OPTIONS}

PATHWAY_SKELETONS_PARSED=${PARSED_ROOT}/pathway_skeletons/${PARSING_DATE}
load_pathway_skeletons:
	(cd ${PATHWAY_SKELETONS_PARSED}/sql_scripts ;		\
	make -f process.mk drop; make -f process.mk all ;	\
	make -f processleaf.mk drop; make -f processleaf.mk all	\
	)

################################################################
#
# BRENDA compounds cleaned by Glaxo
#
parse_glaxo:
	perl ${PARSER}/parse_glaxo.pl ${OPTIONS}

################################################################
#
# Synchronize the parsers between brol and BIGRE
#
BIGRE_PROGRAMS=jvanheld@${SERVER}:amaze_programs/
BIGRE_PARSER=${BIGRE_PROGRAMS}/amaze_parser
rsync_to_bigre:
	${RSYNC} --delete --exclude SWISS --exclude swiss ${PARSER} ${BIGRE_PROGRAMS}

rsync_from_bigre:
	rsync -uptvl ${BIGRE_PARSER}/'*' . 

################################################################
#
# Deliver the result of the parsing
#
#DELIVERY_DIR=regulation
#DELIVERY_DIR=kegg_ligand
#DELIVERY_DIR=kegg_genes
DELIVERY_DIR=swissprot
#DELIVERY_DIR=genbank
DELIVERY_DATE=`date +%Y%m%d`
deliver:
	${RSYNC} --delete							\
		${PARSED_DATA}/${DELIVERY_DIR}/${DELIVERY_DATE}${SUFFIX}	\
		${BIGRE_PARSED_DATA}/${DELIVERY_DIR}

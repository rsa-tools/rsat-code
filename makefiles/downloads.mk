############################################################
#
# $Id: downloads.mk,v 1.2 2003/11/06 21:49:57 jvanheld Exp $
#
# Time-stamp: <2003-10-09 14:02:21 jvanheld>
#
############################################################


################################################################
#### proograms
MAKEFILE=${RSAT}/makefiles/downloads.mk
MAKE = make -sk -f ${MAKEFILE}

DATE = `date +%Y%m%d_%H%M%S`
LOGFILE=-o logs/wget_${DATE}_log.txt
WGET=wget -np -rNL ${LOGFILE}
# #WGET = wget -rNL -o logs/wget_${DATE}_log.txt
RSYNC = rsync -ruptvl -e ssh

### tags
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}

all:
	${MAKE} kegg
	${MAKE} expasy

################################################################
#
# Complete genomes at EBI
#
EBI_GENOMES=ftp://ftp.ebi.ac.uk/pub/databases/genomes/
ebi_genomes:
	@echo importing genomes from EBI
	${WGET} --accept=embl.Z,con,fasta.Z ${EBI_GENOMES}


################################################################
#
# Genbank genome repository
#
GENBANK_DIRS =					\
	genomes					\
	genbank/genomes				\
	refseq 


GENBANK_GENOMES=ftp://ftp.ncbi.nih.gov
GB_DIR=genomes/Saccharomyces_cerevisiae
one_genbank_dir:
	@mkdir -p logs
	@echo "${DATE}	updating dir	$${GB_DIR}" >> wget_updates.txt
	${WGET}							\
		--exclude-directories 'Bacteria.OLD'		\
		--exclude-directories ARCHIVE			\
		--accept=gbk,README,gbff,gaa,faa		\
		--accept=gbk.gz,README.gz,gbff.gz,gaa.gz,faa.gz	\
		"${GENBANK_GENOMES}/${GB_DIR}" 
	@echo "${DATE}	updated dir	$${GB_DIR}" >> wget_updates.txt

genbank:
	@echo "${DATE}	starting to update	${GENBANK_DIRS}" >> wget_updates.txt; 
	@for dir in ${GENBANK_DIRS} ; do		\
		make one_genbank_dir GB_DIR=$${dir} ;	\
	done
	@${WGET} ${GENBANK_GENOMES}/acc
	@${WGET} ${GENBANK_GENOMES}/README
	@${WGET} ${GENBANK_GENOMES}/Bacteria/README
	@${WGET} ${GENBANK_GENOMES}/Bacteria/accessions
	@echo "${DATE}	finished to update	${GENBANK_DIRS}" >> wget_updates.txt; 




################################################################
#
# Download KEGG databases
#
# does not take the full genome sequences

kegg: kegg_ligand kegg_genomes kegg_pathways


KEGG_FTP=ftp://ftp.genome.ad.jp/pub/kegg/

KEGG_GENOMES=${KEGG_FTP}/genomes/
kegg_genomes:
	@mkdir -p logs
	${WGET} -X sequences.old,sequences,sequences.weekly.last.tar.Z,genes.weekly.last.tar.Z ${KEGG_GENOMES}


KEGG_LIGAND=${KEGG_FTP}/ligand/
kegg_ligand:
	@mkdir -p logs
	${WGET} ${KEGG_LIGAND}

KEGG_PATHWAYS=${KEGG_FTP}/pathways/
kegg_pathways:
	@mkdir -p logs
	${WGET} ${KEGG_PATHWAYS}


################################################################
# Databases at Expasy :
# - swiss-prot
# - trembl
# - ENZYME
################################################################
EXPASY=ftp://ftp.expasy.org/databases
EXPASY_DIRS=						\
	README						\
	enzyme						\
	sp_tr_nrdb					\
	swiss-prot/release_compressed		
expasy:
	@mkdir -p logs
	for dir in ${EXPASY_DIRS} ; do		\
		${WGET} ${EXPASY}/$${dir} ;	\
	done

################################################################
#
# Gene ontology
#
################################################################


GO=http://www.godatabase.org/dev/database/archive/latest/
go:
	${WGET} ${GO}

################################################################
#
# Bacilus subtilis transcription factor database

DBTBS=http://elmo.ims.u-tokyo.ac.jp/dbtbs/tfac/
dbtbs:
	${WGET} ${DBTBS}


################################################################
#
# Saccharomyces cerevisiae database
#
SGD=ftp://genome-ftp.stanford.edu/pub/yeast/data_download/
sgd:
	${WGET} --exclude-directories obsolete_files ${SGD}
#	${WGET} ${SGD}/genome_seq/
#	${WGET} ${SGD}/tables/ORF_Locations/ORF_table.txt.gz
#	${WGET} ${SGD}/tables/ORF_Descriptions


################################################################
#
# Alternative yeast genomes
#
yeasts: yeast_mti

################################################################
# Several Saccharomyces genomes from the MIT 
#

#MIT_YEAST=http://www-genome.wi.mit.edu/personal/manoli/yeasts/
MIT_YEAST_DIRS=											\
	http://www-genome.wi.mit.edu/ftp/pub/annotation/fungi/comp_yeasts/			\
	http://www-genome.wi.mit.edu/annotation/fungi/comp_yeasts/downloads.html	
yeast_mit:
	@for dir in ${MIT_YEAST_DIRS} ; do	\
		echo "getting $${dir}";		\
		${WGET} $${dir} ;		\
	done


################################################################
#
# Yeast genome from the MIPS
#
################################################################
MIPS_YEAST = ftp://ftpmips.gsf.de/yeast/
DIRS =						\
	CYGD					\
	catalogues				\
	chri					\
	chrii					\
	chriii					\
	chriv					\
	chrix					\
	chrv					\
	chrvi					\
	chrvii					\
	chrviii					\
	chrx					\
	chrxi					\
	chrxii					\
	chrxiii					\
	chrxiv					\
	chrxv					\
	chrxvi					\
	eurofan					\
	mito

yeast_from_mips:
	for d in ${DIRS}; do					\
		${MAKE} DIR=$${d} yeast_from_mips_one_dir ;	\
	done

DIR=CYGD
yeast_from_mips_one_dir:
	${WGET} ${MIPS_YEAST}/${DIR}

yeast_from_mips_all:
	${WGET} ${MIPS_YEAST}

################################################################
#
# Caenorhabditis genome from euGenes
#
################################################################

WORM_EUGENE=ftp://iubio.bio.indiana.edu/eugenes/worm/features/
worm_from_eugene:
	${WGET} ${WORM_EUGENE}

WORM_WORMBASE =						\
	ftp://ftp.wormbase.org/pub/wormbase/DNA_DUMPS/	\
	ftp://ftp.wormbase.org/pub/wormbase/GENE_DUMPS/
worm_from_wormbase:
	${WGET} ${WORM_WORMBASE}

WORM_SANGER = ftp://ftp.sanger.ac.uk/pub/C.elegans_sequences/ 
worm_from_sanger:
	${WGET} ${WORM_SANGER}

################################################################
#
# Candida albicans genome from Stanford
#
################################################################
CANDIDA = ftp://cycle.stanford.edu/pub/projects/candida/
candida:
	${WGET} ${CANDIDA}

################################################################
#
# Schizosaccharomyces pombe
#
################################################################
POMBE_CONTIGS = ftp://ftp.sanger.ac.uk/pub/yeast/Pombe/CONTIGS
POMBE_GENOME=ftp://ftp.sanger.ac.uk/pub/yeast/sequences/pombe/pombe.dbs.Z
POMBE_PROTEINS=ftp://ftp.sanger.ac.uk/pub/yeast/sequences/pombe/
POMBE_PAULUS = jvanheld@paulus:/permeke/dsk2/genomics/S_pombe/ftp.sanger.ac.uk .
pombe:
	${WGET} ${POMBE_CONTIGS}
	${WGET} ${POMBE_GENOME}
	${WGET} ${POMBE_PROTEINS}
#	${RSYNC} ${POMBE_PAULUS}

################################################################
#
# Ralstonia eutropha
#
################################################################
RALSTONIA = http://genome.ornl.gov/microbial/reut/04dec00/
ralstonia:
	${WGET} ${RALSTONIA}


scpd:
	${WGET} http://cgsigma.cshl.org/jian/


################################################################
#
# Homo sapiens from ENSEMBL

HUMAN_ENSEMBL=ftp://ftp.ensembl.org/pub/current_human/data
human_from_ensembl:
	${WGET} ${HUMAN_ENSEMBL}/golden_path/
	${WGET} ${HUMAN_ENSEMBL}/flatfiles/genbank

ENSEMBL_DIRS=						\
	ftp://ftp.ensembl.org/pub/current_mouse/data/	\
	ftp://ftp.ensembl.org/pub/current_human/data/
ensembl:
	@for dir in ${ENSEMBL_DIRS} ; do				\
		echo "downloading directory $${dir} from ENSEMBL";	\
		${WGET} $${dir} ;					\
	done

################################################################
#
# Homo sapiens from UCSC

HUMAN_UCSC=ftp://genome.ucsc.edu/goldenPath/14nov2002/
#HUMAN_UCSC=http://genome.cse.ucsc.edu/goldenPath/14nov2002/
HUMAN_UCSC_FILES=				\
	bigZips/chromFa.zip			\
	bigZips/contigFa.zip			\
	bigZips/chromFaMasked.zip		\
	bigZips/contigFaMasked.zip		\
	bigZips/upstream5000.zip		\
	database/genename.sql			\
	database/refGene.txt.gz			\
	database/refGene.sql			\
	database/refLink.sql			\
	database/refLink.txt.gz			\
	database/geneName.sql			\
	database/geneName.txt.gz		\
	database/geneid.sql			\
	database/geneid.txt.gz			\
	database/genieAlt.txt			\
	database/genieAltPep.txt 

human_from_ucsc:
	echo "downloading ${HUMAN_UCSC}"
	${WGET} ${HUMAN_UCSC}

#human_from_ucsc:
#	@for f in ${HUMAN_UCSC_FILES} ; do	\
#		echo "downloading ${HUMAN_UCSC}/$${f}" ;	\
#		${WGET} ${HUMAN_UCSC}/$${f} ;	\
#	done


################################################################
#
# Rice 

RICE_GENOME=ftp://ftp.dna.affrc.go.jp/pub/RiceGAAS/current/*.tar.gz
rice:
	${WGET} ${RICE_GENOME}



################################################################
#
# List new bacteria from a wget file
#

UPDATE=200305
new_bacteria:
	cat logs/wget_${UPDATE}*_log.txt		\
		| grep saved			\
		| grep gbk			\
		| grep Bacteria			\
		| perl -pe 's|Bacteria/|\t|'	\
		| cut -f 2			\
		| perl -pe 's|/|\t|'		\
		| cut -f 1			\
		| sort -u

################################################################
#
# Microarray data
#

CELL_CYCLE=http://genome-www.stanford.edu/cellcycle/data/rawdata/individual.html
spellman_cell_cycle:
	${WGET} --accept=.html,.out		\
		${CELL_CYCLE}


################################################################
#
# NCBI taxonomy
#

TAXONOMY=ftp://ftp.ncbi.nih.gov/pub/taxonomy/
taxonomy:
	${WGET} ${TAXONOMY}
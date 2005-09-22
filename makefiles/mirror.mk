############################################################
#
# $Id: mirror.mk,v 1.24 2005/09/22 06:29:40 rsat Exp $
#
# Time-stamp: <2003-10-01 12:05:45 jvanheld>
#
############################################################

RSA=${HOME}/rsa-tools
RSA_SERVER=rsat.scmbb.ulb.ac.be
RSA_SERVER_DIR=rsa-tools
RSA_SERVER_LOGIN=rsat

DATE = `date +%Y%m%d_%H%M%S`


#################################################################
# programs
OPT=
MAKEFILE=${RSAT}/makefiles/mirror.mk
MAKE=make -sk -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT} --exclude '*~'
SSH=-e ssh
RSYNC = rsync ${RSYNC_OPT} ${SSH}

################################################################
#
# Server
RSAT_SERVER = ${RSA_SERVER_LOGIN}@rsat.scmbb.ulb.ac.be:/home/rsat/rsa-tools

SERVER=${RSAT_SERVER}

### tags
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^(\S+):/){ print "\t$$1\n" }' ${MAKEFILE}

################################################################
#### from local machine to servers
################################################################
DIR=perl-scripts
DIRS=perl-scripts public_html doc

rsync_to_server: script_to_server pub_to_server

scripts_to_server:
	echo "Synchronizing perl-scripts to server ${SERVER}"
	${RSYNC} --exclude perllib --exclude perl-scripts/lib/arch --exclude qd.pl perl-scripts ${SERVER}/

pub_to_server:
	echo "Synchronizing public_html to server ${SERVER}"
	${RSYNC} --exclude logs --exclude tmp --exclude data public_html ${SERVER}/  

data_to_server:
	echo "Synchronizing data to server ${SERVER}"
	${RSYNC} ${EXCLUDED_FILES}  public_html/data/* ${SERVER}/public_html/data/

genomes_to_server:
	echo "Synchronizing genomes to server ${SERVER}"
	${RSYNC} ${EXCLUDED_FILES} public_html/data/genomes ${SERVER}/public_html/data/
	${RSYNC}  public_html/data/supported*.pl ${SERVER}/public_html/data/

doc_to_server:
	${MAKE} dir_to_server DIR=doc


dir_to_server:
	echo "Synchronizing dir ${DIR} to server ${SERVER}" 
	${RSYNC} ${DIR} ${SERVER}/ 


################################################################
#### From server to local machine
################################################################
all_from_server: scripts_from_server pub_from_server data_from_server logs_from_server

DIR=doc
TARGET_DIR=${RSA}/
RSYNC_FROM_SERVER_CMD=${RSYNC} ${RSA_SERVER_LOGIN}@${RSA_SERVER}:${RSA_SERVER_DIR}/${DIR} ${TARGET_DIR}
dir_from_server:
	@echo ${RSYNC_FROM_SERVER_CMD}
	${RSYNC_FROM_SERVER_CMD}

doc_from_server:
	${MAKE} dir_from_server DIR=doc

logs_from_server:
	${MAKE} dir_from_server DIR='public_html/logs' RSYNC_OPT='-ruptvl ${OPT}' TARGET_DIR=${RSA}/public_html/

scripts_from_server:
	${MAKE} dir_from_server DIR=perl-scripts

pub_from_server:
	${RSYNC}								\
		--exclude data							\
		--exclude logs							\
		--exclude tmp							\
		${RSA_SERVER_LOGIN}@${RSA_SERVER}:${RSA_SERVER_DIR}/public_html ${RSA}/

EXCLUDED_GENOMES=				\
		--exclude Danio_rerio		\
		--exclude Mus_musculus*		\
		--exclude Gallus_gallus		\
		--exclude Canis_familiaris	\
		--exclude Pan_troglodytes	\
		--exclude Rattus_norvegicus*	\
		--exclude Homo_sapiens*

EXCLUDED_FILES=					\
		--exclude 'blast_db'		\
		--exclude 'blast_hits'		\
		--exclude '*.wc'		\
		--exclude '*.wc.gz'		
#		--exclude '*.fasta'		\
#		--exclude '*.fasta.gz'

EXCLUDED_DIRS=					\
		--exclude embl_genomes		\
		--exclude previous_version	\
		--exclude tmp			\
		--exclude upstream_calibrations	\
		--exclude comparative_genomics

EXCLUDED=${EXCLUDED_GENOMES} ${EXCLUDED_DIRS} ${EXCLUDED_FILES}
data_from_server:
	${RSYNC} ${EXCLUDED}							\
		${RSA_SERVER_LOGIN}@${RSA_SERVER}:${RSA_SERVER_DIR}/public_html/data/*	\
		${RSA}/public_html/data/

################################################################
#### Server in finland
################################################################
ORGS=Saccharomyces_cerevisiae Escherichia_coli_K12 Bacillus_subtilis
medicel:
	${RSYNC} config/medicel.config ${MEDICEL}/config/
	${RSYNC} doc/*.pdf ${MEDICEL}/doc/
	rsync ${SSH} -ruptvL distrib/* ${MEDICEL}/perl-scripts
	for org in ${ORGS}; do				\
		${RSYNC} data/$${org} ${MEDICEL}/data/;	\
	done

################################################################
#### from servers to brol
################################################################

from_rsat:
	${MAKE} dir_from_rsat DIR=perl-scripts
	${MAKE} dir_from_rsat DIR=public_html
	${MAKE} data_from_rsat

dir_from_rsat:
	${RSYNC} ${EXCLUDED} ${RSAT_SERVER}/${DIR} .

data_from_rsat:
	${RSYNC} --exclude Homo_sapiens* --exclude Mus_musculus --exclude Oryza_sativa ${RSAT_SERVER}/public_html/data public_html/ 

SRC=perl-scripts
COMPIL=compil/
PROGRAMS=					\
	random-seq				\
	oligo-analysis				\
	retrieve-seq				\
	dyad-analysis				\
	dna-pattern				\
	orf-info
LIB_TO_COMPILE=					

LIB_TO_COPY=					\
	RSA.stat.lib				\
	RSA.classes				\
	RSA.seq.lib				\
	RSA.cgi.lib				\
	RSA.lib 	
compile:
	@mkdir -p ${COMPIL}/lib
	@mkdir -p ${COMPIL}/bin
	@(cd  ${COMPIL}/bin; ln -fs ../lib)
#	@cp -f config/default.config ${COMPIL}/RSA.config

#### temporary : compiled libraries seem to make problems for the programs
#	@for lb in ${LIB_TO_COMPILE}; do						\
#		echo "compiling library $${lb}";					\
#		cp -f ${SRC}/lib/$${lb} ${COMPIL}/lib/$${lb}.pl ;			\
#		(cd ${COMPIL}/lib; pwd; perlcc -o $${lb} $${lb}.pl && rm -f $${lb}.pl);	\
#	done

	@for lb in ${LIB_TO_COPY}; do				\
		echo "copying library $${lb}";			\
		cp -f ${SRC}/lib/$${lb} ${COMPIL}/lib/ ;	\
	done

	@for pgm in ${PROGRAMS}; do								\
		echo "compiling program $${pgm}";						\
		cp -f ${SRC}/$${pgm} ${COMPIL}/bin/$${pgm}.pl ;					\
		(cd ${COMPIL}/bin; pwd; perlcc -o $${pgm} $${pgm}.pl && rm -f $${pgm}.pl);	\
	done

# more ${GENBANK_DIR}/wget-log_20020108.txt | grep saved| grep .gbk | perl -pe s/.*\`\(.*\)\'.*/\$1/ | perl -pe 's|ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/||'| perl -pe 's|/|\t|' | cut -f 1 | sort -u | xargs
BACTERIA = `ls -1 ${GENBANK_DIR}/Bacteria | grep _ | xargs `
BACT=Mycoplasma_genitalium

install_all_bacteria:
	@for bact in ${BACTERIA}; do			\
		${MAKE} install_one_bacteria BACT=$${bact} ;	\
	done

install_one_bacteria:
	echo
	echo "Installing bacteria ${BACT}"
	${MAKE} install_organism ORGANISM=${BACT}			\
		ORGANISM_DIR=${GENBANK_DIR}/Bacteria/${BACT}

EUKARYOTES=					\
	Encephalitozoon_cuniculi		\
	Plasmodium_falciparum			\
	Saccharomyces_cerevisiae		\
	Schizosaccharomyces_pombe		\
	Caenorhabditis_elegans			\
	Arabidopsis_thaliana			\
	Drosophila_melanogaster
install_all_eukaryotes:
	@for orf in ${EUKARYOTES}; do			\
		${MAKE} install_organism ORGANISM=$${org};	\
	done

link_eukaryotes:
	mkdir -p ${GENBANK_DIR}/Drosophila_melanogaster
	cd ${GENBANK_DIR}/Drosophila_melanogaster ; \
	ln -fs ../../genbank/genomes/D_melanogaster/Scaffolds/LARGE/*.gbk .

	mkdir -p ${GENBANK_DIR}/Caenorhabditis_elegans
	cd ${GENBANK_DIR}/Caenorhabditis_elegans ;	\
	ln -fs ../C_elegans/*/*.gbk .

	mkdir -p ${GENBANK_DIR}/Schizosaccharomyces_pombe
	cd ${GENBANK_DIR}/Schizosaccharomyces_pombe ;	\
	ln -fs ../S_pombe/*/*.gbk .

	mkdir -p ${GENBANK_DIR}/Arabidopsis_thaliana
	cd ${GENBANK_DIR}/Arabidopsis_thaliana ;	\
	ln -fs ../A_thaliana/*/*.gbk .

	mkdir -p ${GENBANK_DIR}/Saccharomyces_cerevisiae
	cd ${GENBANK_DIR}/Saccharomyces_cerevisiae ; \
	ln -fs ../../genbank/genomes/S_cerevisiae/*/*.gbk .

	mkdir -p ${GENBANK_DIR}/Plasmodium_falciparum
	cd ${GENBANK_DIR}/Plasmodium_falciparum ; \
	ln -fs ../../genbank/genomes/P_falciparum/*/*.gbk .


################################################################
#### Install an organism in rsa-tools

################################################################
#### install an organism in rsa-tools after parsing
ORGANISM=Mycoplasma_genitalium
#INSTALL_STEPS= parse,config,start_stop,ncf,allup,oligos,dyads,clean
INSTALL_STEPS=all
install_organism:
	@echo "install log	${INSTALL_LOG}"
	echo "Parsing organism ${ORGANISM}" 
	install-organism -v 1								\
		-org ${ORGANISM}							\
		-step ${INSTALL_STEPS}

################################################################
#### parse a genome from the genbank genome release
#ORGANISM=Plasmodium_falciparum
#ORGANISM=Homo_sapiens
GENBANK_DIR=/lin/genomics/genbank/ftp.ncbi.nih.gov/genomes
ORGANISM_DIR=${GENBANK_DIR}/${ORGANISM}
PARSE_COMMAND=	parse-genbank.pl -v 1 -i ${ORGANISM_DIR}
parse_organism:
	${PARSE_COMMAND}


POMBE_DIR=/win/databases/downloads/ftp.sanger.ac.uk/pub/yeast/Pombe/CONTIGS/
install_pombe:
	echo "Parsing organism Schizosaccharomyces pombe" ;
#	parse-embl.pl -i ${POMBE_DIR} -org 'Schizosaccharomyces pombe' -v 1
	install-organism -v 1											\
		-org Schizosaccharomyces_pombe									\
		-features ${RSA}/data/Schizosaccharomyces_pombe/genome/Gene_Schizosaccharomyces_pombe.tab	\
		-genome ${RSA}/data/genome/Contigs_Schizosaccharomyces_pombe.txt				\
		-format filelist										\
		-source genbank											\
		-step config -step start_stop -step ncf -step oligos -step dyads;

install_gd:
	(cd ${RSA}/lib-sources/GD-1.38;			\
	perl Makefile.PL INSTALLDIRS=site		\
		INSTALLSITELIB=${RSA}/extlib		\
		INSTALLSITEARCH=${RSA}/extlib/arch ;	\
	make ;						\
	make install ;					\
	)


################################################################
#### clean temporary directory
CLEAN_DATE=3
clean_tmp:
	@echo "Before cleaning	" `du -sk public_html/tmp`
	find ${RSA}/public_html/tmp/ -mtime +${CLEAN_DATE} -type f -exec rm -f {} \;	
	@echo "After cleaning	" `du -sk public_html/tmp`
	@echo "Cleaned temporary directory" | mail -s 'cleaning tmp' jvanheld@scmbb.ulb.ac.be


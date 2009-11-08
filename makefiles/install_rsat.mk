############################################################
#
# $Id: install_rsat.mk,v 1.8 2009/11/08 23:58:06 jvanheld Exp $
#
# Time-stamp: <2003-05-23 09:36:00 jvanheld>
#
############################################################

GENBANK_DIR=${RSAT}/downloads/ftp.ncbi.nih.gov/genbank/genomes
NCBI_DIR=${RSAT}/downloads/ftp.ncbi.nih.gov/genomes
V=1
DATE = `date +%Y%m%d_%H%M%S`


#################################################################
# programs

WGET = wget -np -rNL 
MAKEFILE=${RSAT}/makefiles/install_rsat.mk
MAKE=nice -n 19 make -s -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT}
SSH=-e 'ssh -x'
RSYNC = rsync ${RSYNC_OPT} ${SSH}

################################################################
#
# Servers

CIFN = jvanheld@embnet.cifn.unam.mx:rsa-tools/
PAULUS = jvanheld@paulus.ulb.ac.be:rsa-tools/
UPPSALA = jvanheld@bioinformatics.bmc.uu.se:rsa-tools
LIV = jvanheld@liv.bmc.uu.se:rsa-tools
RSAT_SERVER = jvanheld@164.15.61.35:${RSAT}
MEDICEL = root@grimsel.co.helsinki.fi:/work/programs/rsa-tools
SERVERS = ${CIFN} ${UPPSALA} ${LIV}

#SERVER=${RSAT_SERVER}
SERVER=${UPPSALA}

### tags
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^(\S+):/){ print "\t$$1\n"}' ${MAKEFILE}

################################################################
#### from brol to servers
################################################################
DIR=perl-scripts
DIRS=perl-scripts public_html doc
rsync_servers:
	for server in ${SERVERS} ; do					\
		${MAKE} rsync_server SERVER=$${server} DIR=$${dir} ;	\
	done

rsync_server:
	@for dir in ${DIRS}; do							\
		${MAKE} rsync_dir  SERVER=${SERVER} DIR=$${dir} ;	\
	done

rsync_dir:
	echo "Synchronizing dir ${DIR} to server ${SERVER}"
	${RSYNC}  ${OPT} --exclude data --exclude tmp --exclude logs --exclude perl-scripts/lib/arch --exclude qd.pl ${DIR} ${SERVER}/

RSYNC_DATA_CMD=${RSYNC} --exclude 'Mus_musculus*' --exclude 'Homo_sapiens*' public_html/data ${SERVER}/public_html/ 
rsync_data:
	@for server in ${SERVERS} ; do					\
		${MAKE} rsync_data_one_server SERVER=$${server} ;	\
	done

rsync_data_one_server:
	@echo "Synchronizing data to server ${SERVER}" 
	@echo ${RSYNC_DATA_CMD} ;			
	${RSYNC_DATA_CMD};				


ORGS=Saccharomyces_cerevisiae Escherichia_coli_K12 Bacillus_subtilis
medicel:
	${RSYNC} config/medicel.config ${MEDICEL}/config/
	${RSYNC} doc/*.pdf ${MEDICEL}/doc/
	rsync ${SSH} -ruptvL distrib/* ${MEDICEL}/perl-scripts
	for org in ${ORGS}; do				\
		${RSYNC} data/$${org} ${MEDICEL}/data/;	\
	done

rsync_archives:
	@for server in ${SERVERS} ; do				\
		${RSYNC} archives/* $${server}/archives/ ;	\
	done

################################################################
#### from servers to brol
################################################################
rsync_logs:
	@for server in ${SERVERS} ; do					\
		echo "${RSYNC} $${server}/logs/log-file_* logs/" ;	\
		${RSYNC} $${server}/logs/log-file_* logs/ ;		\
	done

rsync_config:
	@for server in ${SERVERS} ; do \
		${RSYNC} $${server}/config/*.config config/ ;\
	done

FOLDERS=data pdf_files
from_ucmb:
	@for folder in ${FOLDERS}; do \
		${RSYNC} jvanheld@${PAULUS}:rsa-tools/$${folder} . ; \
	done

FOLDERS=data 
from_cifn:
	@for folder in ${FOLDERS}; do \
		${RSYNC} jvanheld@${CIFN}:rsa-tools/$${folder}/* ./$${folder} ; \
	done

SRC=perl-scripts
COMPIL=compil/
PROGRAMS=	\
	oligo-analysis	\
	retrieve-seq	\
	dyad-analysis	\
	dna-pattern	\
	orf-info
LIBRARIES=\
	RSA.stat.lib 	\
#	RSA.classes	\
#	RSA.seq.lib	\
#	RSA.cgi.lib	\
#	RSA.lib 	
compile:
	@mkdir -p ${COMPIL}/lib
	@mkdir -p ${COMPIL}/bin
	@(cd  ${COMPIL}/bin; ln -fs ../lib)
	@cp -f config/default.config ${COMPIL}/RSA.config

	@for lb in ${LIBRARIES}; do \
		echo "compiling library $${lb}"; \
		cp -f ${SRC}/lib/$${lb} ${COMPIL}/lib/$${lb}.pl ; \
		(cd ${COMPIL}/lib; pwd; perlcc $${lb}.pl && rm -f $${lb}.pl); \
	done

	@for pgm in ${PROGRAMS}; do \
		echo "compiling program $${pgm}"; \
		cp -f ${SRC}/$${pgm} ${COMPIL}/bin/$${pgm}.pl ; \
		(cd ${COMPIL}/bin; pwd; perlcc $${pgm}.pl && rm -f $${pgm}.pl); \
	dgone

BACTERIA = `ls -1 ${NCBI_DIR}/Bacteria | grep _ | sort -u | grep -v bacteria | xargs `
BACT=Mycoplasma_genitalium

list_bacteria:
	@echo "Bacteria to install"
	@echo ${BACTERIA}

TASK=install_one_bacteria
iterate_all_bacteria:
	@for bact in ${BACTERIA}; do				\
		${MAKE} ${TASK} BACT=$${bact} ;	\
	done

install_all_bacteria:
	${MAKE} iterate_all_bacteria TASK=install_one_bacteria

install_one_bacteria:
	@echo
	@echo "${DATE}	Installing bacteria ${BACT}"
	@${MAKE} install_organism ORGANISM=${BACT}		\
		ORGANISM_DIR=${NCBI_DIR}/Bacteria/${BACT}



#ORGANISM=Plasmodium_falciparum
ORGANISM=Homo_sapiens
ORGANISM_DIR=${NCBI_DIR}/${ORGANISM}
INSTALL_TASK=allup,clean,config,dyads,ncf,intergenic_freq,oligos,parse,start_stop,upstream_freq
install_organism:
	@echo "install log	${INSTALL_LOG}"
	echo "Parsing organism ${ORGANISM}" 
	install-organism -v ${V}								\
		-org ${ORGANISM}							\
		-task  ${INSTALL_TASK}

POMBE_DIR=/win/databases/downloads/ftp.sanger.ac.uk/pub/yeast/Pombe/CONTIGS/
install_pombe:
	echo "Parsing organism Schizosaccharomyces pombe" ;
#	parse-embl.pl -i ${POMBE_DIR} -org 'Schizosaccharomyces pombe' -v ${V}
	install-organism -v ${V}											\
		-org Schizosaccharomyces_pombe									\
		-features ${RSAT}/data/Schizosaccharomyces_pombe/genome/Gene_Schizosaccharomyces_pombe.tab	\
		-genome ${RSAT}/data/genome/Contigs_Schizosaccharomyces_pombe.txt				\
		-format filelist										\
		-source genbank											\
		-step config -step start_stop -step ncf -step oligos -step dyads;

################################################################
#### parse a genome from the genbank genome release
ORGANISM=Plasmodium_faciparum
#ORGANISM=Homo_sapiens
ORGANISM_DIR=${NCBI_DIR}/${ORGANISM}
PARSE_COMMAND=parse-genbank.pl -v ${V} -i ${ORGANISM_DIR} ${OPT}
parse_organism:
	${PARSE_COMMAND}

parse_one_bacteria:
	${MAKE} parse_organism ORGANISM=${BACT} NCBI_DIR=${NCBI_DIR}/Bacteria ${OPT}

parse_all_bacteria:
	${MAKE} iterate_all_bacteria TASK=parse_one_bacteria

################################################################
#### installation of the GD graphical library
#### (obsolete)
GD_DISTRIB_DIR=stein.cshl.org/WWW/software/GD/
GD_DISTRIB= http://${GD_DISTRIB_DIR}/GD.pm.tar.gz
get_gd:
	mkdir -p lib-sources
	(cd lib-sources; \
	${WGET} ${GD_DISTRIB})


GD_VERSION=2.06
uncompress_gd:
	(cd lib-sources/${GD_DISTRIB_DIR};		\
	gunzip -c GD.pm.tar.gz | tar -xpf - ;		\

install_gd:
	(cd lib-sources/${GD_DISTRIB_DIR}/GD-${GD_VERSION} ;	\
	perl Makefile.PL INSTALLDIRS=site			\
		INSTALLSITELIB=${RSAT}/extlib			\
		INSTALLSITEARCH=${RSAT}/extlib/arch ;		\
	make ;							\
	make install ;						\
	)


APP_DIR=${RSAT}/applications
PROGRAM=consensus
PROGRAM_DIR=${APP_DIR}/${PROGRAM}
PROGRAM_ARCHIVE=`ls -1t ${RSAT}/app_sources/${PROGRAM}* | head -1`
uncompress_program:
	@echo installing ${PROGRAM_ARCHIVE} in dir ${PROGRAM_DIR}
	@mkdir -p ${PROGRAM_DIR}
	@(cd ${PROGRAM_DIR} ;				\
	gunzip -c ${PROGRAM_ARCHIVE} | tar -xf - )

INSTALLED_PROGRAM=`ls -1t ${APP_DIR}/${PROGRAM}/${PROGRAM}*`
install_program:
	(cd ${PROGRAM_DIR}; make ${INSTALL_OPT})
	(cd bin; ln -fs ${INSTALLED_PROGRAM} ./${PROGRAM})


install_consensus:
	${MAKE} uncompress_program PROGRAM=consensus
	${MAKE} install_program PROGRAM=consensus INSTALL_OPT='CPPFLAGS=""'

GIBBS_DIR=${APP_DIR}/gibbs/gibbs9_95
install_gibbs:
	${MAKE} uncompress_program PROGRAM=gibbs
	(cd ${GIBBS_DIR}; ./compile; cd ${GIBBS_DIR}/code; make clean)
	(cd bin; ln -fs ${GIBBS_DIR}/gibbs ./gibbs)

install_ext_apps:
	${MAKE} install_consensus
	${MAKE} install_patser
	${MAKE} install_gibbs

################################################################
## Install the BioPerl library
## For this example, we install Bioperl and EnsEMBL libraries 
## in $RSAT/lib, but you can install it in some other place
### (password is 'cvs')
bioperl:
	@mkdir -p ${RSAT}/lib
	@echo "Password is 'cvs'"
	@cvs -d :pserver:cvs@code.open-bio.org:/home/repository/bioperl login
	(cd ${RSAT}/lib;  cvs -d :pserver:cvs@code.open-bio.org:/home/repository/bioperl checkout bioperl-live)


################################################################
## Install the EnsEMBL Perl API
ENSEMBL_BRANCH=56
ensembl_api:	
	@echo  "Password is 'CVSUSER'"
	@cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl login
	@(cd ${RSAT}/lib; \
		cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl \
		checkout -r branch-ensembl-${ENSEMBL_BRANCH} ensembl ; \
		cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl \
		checkout -r branch-ensembl-${ENSEMBL_BRANCH} ensembl-compara)
	@echo "Don't forget to adapt the following lines in the file ${RSAT}/RSAT_config.props"
	@echo "ensembl=${RSAT}/lib/ensembl/modules"
	@echo "compara=${RSAT}/lib/ensembl-compara/modules"
	@echo "bioperl=${RSAT}/lib/bioperl-live"


################################################################
## Get and install patser (matrix-based pattern matching)
#PATSER_TAR=patser-v3e.1.tar.gz
PATSER_VERSION=patser-v3b.5
PATSER_TAR=${PATSER_VERSION}.tar.gz
PATSER_URL=ftp://www.genetics.wustl.edu/pub/stormo/Consensus/
PATSER_DIR=ext/patser/${PATSER_VERSION}
PATSER_APP=`cd ${PATSER_DIR} ; ls -1tr patser-v* | grep -v .tar | tail -1 | xargs`
get_patser:
	@mkdir -p ${PATSER_DIR}
	@echo "Getting patser using ${WGET}"
	(cd ${PATSER_DIR}; ${WGET} -nv ${PATSER_URL}/${PATSER_TAR}; tar -xpzf ${PATSER_TAR})
	@echo "patser dir	${PATSER_DIR}"

install_patser:
	@echo "Installing patser"
#	(cd ${PATSER_DIR}; rm *.o; make)
	rsync -ruptvl ${PATSER_DIR}/${PATSER_APP} ${RSAT}/bin/
	(cd ${RSAT}/bin; ln -fs ${PATSER_APP} patser)
	@echo "ls -ltr ${RSAT}/bin/patser*"
#	${MAKE} uncompress_program PROGRAM_DIR=${PATSER_DIR} PROGRAM=patser
#	${MAKE} install_program PROGRAM=patser

################################################################
## Get and install the program seqlogo
SEQLOGO_URL=http://weblogo.berkeley.edu/release
SEQLOGO_TAR=weblogo.2.8.2.tar.gz
SEQLOGO_DIR=${RSAT}/ext/seqlogo
get_seqlogo:
	@mkdir -p ${SEQLOGO_DIR}
	@echo "Getting seqlogo using ${WGET}"
	(cd ${SEQLOGO_DIR}; ${WGET} -nv ${SEQLOGO_URL}/${SEQLOGO_TAR}; tar -xpzf ${SEQLOGO_TAR})
	@echo "seqlogo dir	${SEQLOGO_DIR}"

install_seqlogo:
	@echo "Installing seqlogo"
	@rsync -ruptl ${SEQLOGO_DIR}/weblogo/seqlogo ${RSAT}/bin/
	@rsync -ruptl ${SEQLOGO_DIR}/weblogo/template.* ${RSAT}/bin/
	@rsync -ruptl ${SEQLOGO_DIR}/weblogo/logo.pm ${RSAT}/bin/

################################################################
## Get and install the program gnuplot
GNUPLOT_URL=http://sourceforge.net/projects/gnuplot/files/gnuplot/4.2.6/
GNUPLOT_VER=gnuplot-4.2.5
GNUPLOT_TAR=${GNUPLOT_VER}.tar.gz
GNUPLOT_DIR=${RSAT}/ext/gnuplot
get_gnuplot:
	@mkdir -p ${GNUPLOT_DIR}
	@echo "Getting gnuplot using ${WGET}"
	(cd ${GNUPLOT_DIR}; ${WGET} -nv ${GNUPLOT_URL}/${GNUPLOT_TAR}; tar -xpzf ${GNUPLOT_TAR})
	@echo "gnuplot dir	${GNUPLOT_DIR}"

install_gnuplot:
	@echo "Installing gnuplot"
	(cd ${GNUPLOT_DIR}/${GNUPLOT_VER}; ./configure && make)

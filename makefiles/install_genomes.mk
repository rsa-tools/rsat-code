############################################################
#
# $Id: install_genomes.mk,v 1.48 2013/04/08 22:35:33 rsat Exp $
#
# Time-stamp: <2003-10-10 22:49:55 jvanheld>
#
############################################################

include ${RSAT}/makefiles/util.mk

DATE = `date +%Y%m%d_%H%M%S`

################################################################
#### Directories
GENBANK_DIR=${RSAT}/downloads/ftp.ncbi.nih.gov/genbank/genomes
NCBI_DIR=${RSAT}/downloads/ftp.ncbi.nih.gov/genomes

#################################################################
#### Programs

WGET = wget -np -rNL 
MAKEFILE=${RSAT}/makefiles/install_genomes.mk
#MAKE=nice -n 19 make -f ${MAKEFILE}
RSYNC_OPT = -ruptvl ${OPT}
SSH=-e 'ssh -x'
RSYNC = rsync ${RSYNC_OPT} ${SSH}

V=1 

################################################################
### Targets

### Install one organism
ORG=Saccharomyces_cerevisiae
MASKING=''
ORG_DIR=${NCBI_DIR}/${ORG}
PARSE_TASK=config,parse
CALIBRATE_TASK=allup,seq_len_distrib,dyads,oligos,start_stop,upstream_freq,genome_segments,protein_freq
INSTALL_TASK=${PARSE_TASK},${CALIBRATE_TASK}
INSTALL_CMD=install-organism -v ${V}		\
		-genbank ${NCBI_DIR}		\
		-org ${ORG} ${MASKING}		\
		-task ${INSTALL_TASK}		\
		${OPT}

parse_one_organism:
	@echo "Parsing organism ${ORG}" 
	@${MAKE}  one_install_command WHEN=now INSTALL_TASK=${PARSE_TASK}

calibrate_one_organism:
	@echo "Calibrating organism ${ORG}" 
	@${MAKE}  one_install_command INSTALL_TASK=${CALIBRATE_TASK}

install_one_organism:
	@echo "Installing organism ${ORG}" 
	@${MAKE} one_install_command 

## Install a set of selected organuisms
SELECTED_ORGANISMS=Mycoplasma_genitalium_G37_uid57707 Escherichia_coli_K_12_substr__MG1655_uid57779
install_selected_organisms:
	@for org in ${SELECTED_ORGANISMS}; do				\
		${MAKE} install_one_organism ORG=$${org};	\
	done

one_install_command:
	${MAKE} my_command MY_COMMAND="${INSTALL_CMD}" JOB_PREFIX=${ORG}_install

### Prokaryote with a small genome for quick testing
PRO=Mycoplasma_genitalium_G37_uid57707

### All the prokaryote in NCBI genome directory
PROKARYOTES = `ls -1 ${NCBI_DIR}/Bacteria | grep _ | sort -u | xargs `
list_prokaryotes:
	@echo "Prokaryotes"
	@echo "${PROKARYOTES}"


### Install all prokaryotes genomes on RSAT
install_all_prokaryotes:
	@for pro in ${PROKARYOTES}; do				\
		${MAKE} install_one_prokaryote PRO=$${pro};	\
	done

### Install all prokaryotes genomes on RSAT
NEW_PRO_FILE=${RSAT}/new_bact.txt
NEW_PRO=`cat ${NEW_PRO_FILE}`
list_new_prokaryotes:
	@echo
	@echo "New prokaryotes from file	${NEW_PRO_FILE}"
	@echo "${NEW_PRO}"


install_new_prokaryotes:
	@for pro in ${NEW_PRO}; do				\
		${MAKE} install_one_prokaryote PRO=$${pro};	\
	done


### Install a single prokaryote genome
install_one_prokaryote:
	@echo
	@echo "${DATE}	Installing prokaryote ${PRO}"
	@${MAKE} install_one_organism ORG=${PRO}		\
		NCBI_DIR=${NCBI_DIR}/Bacteria

################################################################
## Install phages from NCBI
PHAGE=Enterobacteria_phage_lambda
PHAGE_DIR=/bio/databases/Phages4rsat
install_one_phage:
	@echo
	@echo "${DATE}	Installing prokaryote ${PHAGE}"
	@${MAKE} install_one_organism ORG=${PHAGE}		\
		NCBI_DIR=${PHAGE_DIR}

### All the phage in NCBI genome directory
PHAGES = `ls -1 ${PHAGE_DIR} | grep _ | sort -u | xargs `
list_phages:
	@echo "Phages"
	@echo "${PHAGES}"

### Install all phages genomes on RSAT
install_all_phages:
	@for pro in ${PHAGES}; do				\
		${MAKE} install_one_phage PHAGE=$${pro};	\
	done

### All the fungi in NCBI genome directory
NCBI_FUNGI = `ls -1 ${NCBI_DIR}/Fungi | grep _ | grep -v Candida_albicans | sort -u | xargs `
OTHER_FUNGI=					\
	Aspergillus_nidulans \
	Aspergillus_oryzae \
	Aspergillus_terreus \
	Candida_dubliniensis \
	Candida_guilliermondii \
	Candida_lusitaniae \
	Candida_tropicalis \
	Chaetomium_globosum \
	Coccidioides_immitis \
	Fusarium_graminearum \
	Fusarium_oxysporum \
	Fusarium_verticillioides \
	Histoplasma_capsulatum \
	Kluyveromyces_waltii \
	Lodderomyces_elongisporus \
	Magnaporthe_grisea \
	Phanerochaete_chrysosporium \
	Puccinia_graminis \
	Rhizopus_oryzae \
	Saccharomyces_bayanus \
	Saccharomyces_castellii \
	Saccharomyces_kluyveri \
	Saccharomyces_kudriavzevii \
	Saccharomyces_mikatae \
	Saccharomyces_paradoxus \
	Schizosaccharomyces_japonicus_yFS275 \
	Sclerotinia_sclerotiorum \
	Stagonospora_nodorum \
	Trichoderma_reesei \
	Uncinocarpus_reesii

### List the Fungi in the NCBI directory
FUNGI= ${NCBI_FUNGI} ${OTHER_FUNGI}
list_fungi:
	@echo "Fungi"
	@echo "${FUNGI}"

list_ncbi_fungi:
	@echo "NCBI Fungi"
	@echo "${NCBI_FUNGI}"

### Pare fungal genomes from NCBI
parse_ncbi_fungi:
	@for fungus in ${NCBI_FUNGI}; do						\
		${MAKE} install_one_fungus FUNGUS=$${fungus} INSTALL_TASK=parse;	\
	done

config_ncbi_fungi:
	@for fungus in ${NCBI_FUNGI}; do						\
		${MAKE} install_one_fungus FUNGUS=$${fungus} INSTALL_TASK=config;	\
	done

### Install all fungi genomes on RSAT
install_all_fungi:
	@for fungus in ${FUNGI}; do				\
		${MAKE} install_one_fungus FUNGUS=$${fungus} ;	\
	done

### Install a single fungus genome
FUNGUS=Saccharomyces_cerevisiae_uid128
install_one_fungus:
	@echo
	@echo "${DATE}	Installing fungus ${FUNGUS}"
	@${MAKE} install_one_organism ORG=${FUNGUS}		\
		OPT='-prefid CDS locus_tag' \
		NCBI_DIR=${NCBI_DIR}/Fungi 

### List the Protozoa in the NCBI directory
NCBI_PROTOZOA = `ls -1 ${NCBI_DIR}/Protozoa | grep _ | sort -u | xargs `
PROTOZOA=${NCBI_PROTOZOA}
list_protozoa:
	@echo "Protozoa"
	@echo "${PROTOZOA}"


### Install all protozoa genomes on RSAT
install_all_protozoa:
	@for proto in ${PROTOZOA}; do				\
		${MAKE} install_one_proto PROTO=$${proto} ;	\
	done

#PROTO=Plasmodium_falciparum
PROTO=Dictyostelium_discoideum
### Install a single proto genome
install_one_proto:
	@echo
	@echo "${DATE}	Installing proto ${PROTO}"
	@${MAKE} install_one_organism ORG=${PROTO}		\
		NCBI_DIR=${NCBI_DIR}/Protozoa 

### Parse one organism
parse_organism:
	@echo "Parsing organism ${ORG}"
	${MAKE} install_one_organism INSTALL_TASK=parse


ENSEMBL_DIR=${RSAT}/downloads/ftp.ensembl.org/pub/current_worm/data/flatfiles/genbank
parse_organism_ensembl:
	parse-genbank.pl -v 1 -source ensembl -ext dat -i ${ENSEMBL_DIR} -org ${ORG}_ENSEMBL

install_organism_ensembl:
	${MAKE} install_one_organism  ORG=${ORG}_ENSEMBL OPT='-source ensembl -ensembl ${ENSEMBL_DIR}'

################################################################
#### Install all eukaryote genomes
#### Genomes are selected manually because NCBI directories are
#### a bit messy for eukaryotes.
EUKARYOTES=					\
	Saccharomyces_cerevisiae		\
	Plasmodium_falciparum			\
	Schizosaccharomyces_pombe		\
	Encephalitozoon_cuniculi		\
	Apis_mellifera				\
	Drosophila_melanogaster			\
	Arabidopsis_thaliana			\
	Caenorhabditis_elegans

### List selected eukaroytes
list_eukaryotes:
	@echo "Eukaryote to install	${EUKARYOTES}"
	@echo "Reference Eukaryote to install	${REF_EUKARYOTES}"

install_all_eukaryotes:
	@for org in ${EUKARYOTES} ; do \
		${MAKE} install_one_organism ORG=$${org} INSTALL_TASK=${INSTALL_TASK},clean; \
	done
	@for org in ${REF_EUKARYOTES} ; do \
		${MAKE} install_one_ref_eukaryote ORG=$${org} INSTALL_TASK=${INSTALL_TASK},clean; \
	done
	@${MAKE} install_human
	@${MAKE} install_rat
	@${MAKE} install_zebrafish
	@${MAKE} install_mouse

################################################################
## Install organisms found in the root of the NCBI genome distribution. 
FULL_ORG=${ORG}
LINK_DIR=${RSAT}/genome_installations
LINK_DIR_ORG=${LINK_DIR}/${FULL_ORG}
LINK_PATTERN=CHR_*/*.gbk.gz
install_one_eukaryote:
	@echo 'Installing eukaryote ${ORG}	${FULL_ORG}'
	@mkdir -p ${LINK_DIR_ORG}
	@rm -rf ${LINK_DIR_ORG}/*
	@(cd ${LINK_DIR_ORG}; ln -s ${ORG_DIR}/${LINK_PATTERN} .)
	@echo ${LINK_DIR_ORG}
	${MAKE} install_one_organism NCBI_DIR=${LINK_DIR} ORG=${FULL_ORG}

################################################################
## use only the _ref_ files (e.g. for vertebrate genomes)
REF_EUKARYOTES= \
	Gallus_gallus				\
	Canis_familiaris			\
	Pan_troglodytes
install_one_ref_eukaryote:
	${MAKE} install_one_eukaryote LINK_PATTERN=CHR_*/*_ref_*.gbk.gz


################################################################
## Specific treatment for some organisms, because the folder name ($ORG)
## differs from the organism name in the NCBI distribution ($FULL_ORG)
## (e.g. ORG=H_sapiens FULL_ORG=Homo_sapiens)
install_mouse:
	${MAKE} install_one_ref_eukaryote ORG=M_musculus FULL_ORG=Mus_musculus

install_zebrafish:
	${MAKE} install_one_ref_eukaryote ORG=D_rerio FULL_ORG=Danio_rerio

install_rat:
	${MAKE} install_one_ref_eukaryote ORG=R_norvegicus FULL_ORG=Rattus_norvegicus

install_human:
	${MAKE} install_one_ref_eukaryote ORG=H_sapiens FULL_ORG=Homo_sapiens



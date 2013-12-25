################################################################
## Import and convert motif databases

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/motif_databases.mk


################################################################
## Importe and convert a collection from uniprobe
#UNIPROBE_PWMS=http://thebrain.bwh.harvard.edu/uniprobe/downloads/All/all_pwms.zip
uniprobe_one_collection: uniprobe_get_pwms uniprobe_convert_pwms uniprobe_db_descr

UNIPROBE_DIR=public_html/data/motif_databases/uniprobe
UNIPROBE_DOWNLOAD_DIR=${UNIPROBE_DIR}/downloads
UNIPROBE_REF=Cell09
UNIPROBE_DATE=2009
UNIPROBE_AUTHOR=Groves
UNIPROBE_COLLECTION=worm
UNIPROBE_ARCHIVE=worm_pwm_all
UNIPROBE_URL=http://thebrain.bwh.harvard.edu/uniprobe/downloads/${UNIPROBE_REF}/${UNIPROBE_ARCHIVE}.zip
uniprobe_get_pwms:
	@mkdir -p ${UNIPROBE_DOWNLOAD_DIR}
	(cd  ${UNIPROBE_DOWNLOAD_DIR}; wget -NL ${UNIPROBE_URL}; unzip ${UNIPROBE_ARCHIVE}.zip);
	@echo "Collection ${UNIPROBE_COLLECTION} downloaded"
	@ls -l ${UNIPROBE_DOWNLOAD_DIR}/${UNIPROBE_ARCHIVE}.*

UNIPROBE_CONVERTED=${UNIPROBE_DIR}/uniprobe_${UNIPROBE_COLLECTION}_${UNIPROBE_REF}.tf
uniprobe_convert_pwms:
	@echo
	@echo "Converting Uniprobe collection ${UNIPROBE_COLLECTION}"
	convert-matrix -v 0 -from uniprobe -to transfac \
		-i ${UNIPROBE_DOWNLOAD_DIR}/${UNIPROBE_ARCHIVE}.txt \
		-multiply 100 \
		-o ${UNIPROBE_CONVERTED}
	@ls -l ${UNIPROBE_CONVERTED}

uniprobe_db_descr:
	@echo
	@echo "Database description line"
	@echo "Uniprobe_${UNIPROBE_COLLECTION}	uniprobe	${UNIPROBE_CONVERTED}	Uniprobe ${UNIPROBE_COLLECTION} from ${UNIPROBE_AUTHOR} (${UNIPROBE_REF})	2009	${UNIPROBE_URL}"

uniprobe_worm:
	@${MAKE} uniprobe_one_collection \
		UNIPROBE_REF=Cell09 \
		UNIPROBE_AUTHOR=Groves \
		UNIPROBE_COLLECTION=worm \
		UNIPROBE_DATE=2009 \
		UNIPROBE_ARCHIVE=worm_pwm_all

uniprobe_human:
	@${MAKE} uniprobe_get_pwms \
		UNIPROBE_REF=NAR10 \
		UNIPROBE_AUTHOR=Alibes \
		UNIPROBE_COLLECTION=human \
		UNIPROBE_DATE=2010 \
		UNIPROBE_ARCHIVE=NAR10_pwm

################################################################
## Import JASPAR database
JASPAR_COLLECTIONS=all fungi insects nematodes plants urochordates vertebrates
JASPAR_COLLECTION=fungi
JASPAR_MATRICES_URL=http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/redundant
JASPAR_MARICES_ORI=pfm_${JASPAR_COLLECTION}.txt
JASPAR_COLLECTION_URL=${JASPAR_MATRICES_URL}/${JASPAR_MARICES_ORI}

## Download and convert jaspar matrices
jaspar_matrices:
	@echo
	@echo "Importing JASPAR matrices"
	@echo
	@echo "#DB_NAME	FORMAT	FILE	DESCR   VERSION	URL" > new_jaspar_db_matrix_files.tab
	@for col in ${JASPAR_COLLECTIONS}; do \
		${MAKE} jaspar_matrices_one_collection JASPAR_COLLECTION=$${col}; \
	done
	@echo "JASPAR matrices have been parsed"
	@echo
	@echo "You can now check the content of DB matrix file"
	@echo "	new_jaspar_db_matrix_files.tab"
	@echo "and paste its content in the RSAT DB matrix file"
	@echo "	public_html/data/motif_databases/db_matrix_files.tab"

JASPAR_MATRICES_DIR=${RSAT}/public_html/data/motif_databases/JASPAR
JASPAR_RELEASE=2013-11
JASPAR_COLLECTION_TF=${JASPAR_MATRICES_DIR}/jaspar_core_${JASPAR_COLLECTION}_${JASPAR_RELEASE}.tf
jaspar_matrices_one_collection:
	@echo
	@echo "Downloading and converting JASPAR collection	${JASPAR_COLLECTION}"
	(cd ${JASPAR_MATRICES_DIR}; wget --no-verbose --timestamping --no-directories ${JASPAR_COLLECTION_URL}; )
	@echo "	${JASPAR_MATRICES_DIR}/${JASPAR_MARICES_ORI}"
	convert-matrix -v 0 -from jaspar -to transfac \
		-i ${JASPAR_MATRICES_DIR}/${JASPAR_MARICES_ORI} \
		-pseudo 1 -bg_pseudo 0.01 \
		-return counts,consensus,parameters \
		-o ${JASPAR_COLLECTION_TF}
	@echo "	${JASPAR_COLLECTION_TF}"
	@echo "jaspar_core_${JASPAR_COLLECTION}	tf	JASPAR/jaspar_core_${JASPAR_COLLECTION}_${JASPAR_RELEASE}.tf	 JASPAR core ${JASPAR_COLLECTION}	${JASPAR_RELEASE}	${JASPAR_COLLECTION_URL}" >> new_jaspar_db_matrix_files.tab

JASPAR_SITES_URL=http://jaspar.genereg.net/html/DOWNLOAD/sites/
JASPAR_SITES_DIR=
jaspar_sites_download:
	(cd ${JASPAR_SITES_DIR}; wget --no-verbose -rNL --no-directories ${JASPAR_SITES_URL})

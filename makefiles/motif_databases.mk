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

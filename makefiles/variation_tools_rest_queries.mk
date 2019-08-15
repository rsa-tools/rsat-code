## Load site-specific options for this RSAT instance
include ${RSAT}/makefiles/util.mk

################################################################
## Variables
V=1
MAKEFILE=${RSAT}/makefiles/variation_tools_rest_queries.mk
MAKE=make -s -f ${MAKEFILE}
DATE=`date +%Y-%M-%d_%H:%M:%S`
DAY=`date +%Y%m%d`
TIME=`date +%Y%m%d_%H%M%S`
SSH_OPT = -e ssh 
RSYNC_OPT= -ruptvlz  ${SSH_OPT} 
RSYNC = rsync  ${RSYNC_OPT}
WGET=wget --passive-ftp -np -rNL

REST_SERVER=http://rsat-tagc.univ-mrs.fr/rsat/
REST_ROOT=${REST_SERVER}/rest.wsgi/
#VAR_IDS=rs554219 rs1542725 rs2981578 rs35054928 rs7948643
VAR_IDS=rs554219 rs1542725 SNP195 rs1859961 rs7895676 rs6983267 rs35054928 rs57095329 rs12740374 rs11986220 rs2981578 rs75915166 rs7948643 rs8072254 rs4784227
#VAR_ID_STRING=rs554219%2Crs1542725%2Crs2981578%2Crs35054928%2Crs7948643
VAR_ID_STRING=rs554219%2Crs1542725%2CSNP195%2Crs1859961%2Crs7895676%2Crs6983267%2Crs35054928%2Crs57095329%2Crs12740374%2Crs11986220%2Crs2981578%2Crs75915166%2Crs7948643%2Crs8072254%2Crs4784227
RESULT_DIR=test/vartools_rest
VAR_PREFIX=${RESULT_DIR}/some_snps
VAR_ID_FILE=${VAR_PREFIX}_ids.txt
SPECIES=Homo_sapiens
ASSEMBLY=GRCh38
VAR_INFO_GET=${VAR_PREFIX}_varinfo_GET.varBed
VAR_INFO_POST=${VAR_PREFIX}_varinfo_POST.varBed
VAR_SEQ_POST=${VAR_PREFIX}_varseq_POST.varSeq
list_param:
	@echo "Parameters"
	@echo "	SPECIES		${SPECIES}"
	@echo "	ASSEMBLY	${ASSEMBLY}"
	@echo "	REST_SERVER	${REST_SERVER}"
	@echo "	REST_ROOT	${REST_ROOT}"
	@echo "	VAR_ID_STRING	${VAR_ID_STRING}"
	@echo "	RESULT_DIR	${RESULT_DIR}"
	@echo "	VAR_ID_FILE	${VAR_ID_FILE}"
	@echo "	VAR_INFO_GET	${VAR_INFO_GET}"
	@echo "	VAR_INFO_POST	${VAR_INFO_POST}"
	@echo "	VAR_SEQ_POST	${VAR_SEQ_POST}"
	@echo "	SOME_MOTIFS	${SOME_MOTIFS}"

## Run all the queries below
all: dir varinfo_get write_snp_query_file varinfo_post

dir:
	@mkdir -p ${RESULT_DIR}
	@echo "	RESULT_DIR	${RESULT_DIR}"


################################################################
## Select relevant motifs for this study case
SELECTED_TFS=CEBPA,CEBPB,ELF1,ERG,ETS1,ETV4,FOXA1,FOXA2,GABPA,GATA2,GFI1B,POU2F2,RUNX2,RUNX3,ZNF384
MATRIX_COLLECTION=public_html/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.tf
select_matrices:
	retrieve-matrix -v ${V} -i ${MATRIX_COLLECTION} -id ${SELECTED_TFS}


################################################################
## Send a GET request to RSAT variation-info
GET_STRING=${REST_ROOT}variation-info/${SPECIES}/${ASSEMBLY}?i_string=${VAR_ID_STRING}&i_string_type=text&format=id&col=1&content-type=text%2Fplain
varinfo_get:
	curl -X GET  '${GET_STRING}' -H "accept: text/plain" > ${VAR_INFO_GET}
	@echo "	VAR_INFO_GET	${VAR_INFO_GET}"

################################################################
## Write a small query file for the POST query
write_snp_query_file:
	@rm -f ${VAR_ID_FILE}
	@for snp in ${VAR_IDS}; do \
		echo $${snp} >> ${VAR_ID_FILE} ; \
	done
	@echo "	VAR_ID_FILE	${VAR_ID_FILE}"

################################################################
## Sent the REST query as POST
varinfo_post: write_snp_query_file
	curl -X POST "${REST_ROOT}/variation-info/${SPECIES}/${ASSEMBLY}" \
		-H "accept: text/plain" \
		-H "Content-Type: multipart/form-data" \
		-F "i=@${VAR_ID_FILE};type=text/plain" \
		-F "i_string_type=text" \
		-F "format=id" \
		-F "col=1" > ${VAR_INFO_POST}
	@echo "	VAR_INFO_POST	${VAR_INFO_POST}"

## NOTE: one of the SNPs is a deletion, but the variation-info option -type does not report it
#		-F "type=SNV" \
# I should add a specific report for variants filtered out because they dont belong to the selected type

################################################################
## Retrieve variation sequences via a POST query
varseq_post:
	@echo "	VAR_INFO_POST	${VAR_INFO_POST}"
	curl -X POST "${REST_ROOT}/retrieve-variation-seq/${SPECIES}/${ASSEMBLY}" \
		-H "accept: text/plain" \
		-H "Content-Type: multipart/form-data" \
		-F "i=@some_snps_varinfo_GET.varBed;type=text/plain" \
		-F "i_string_type=text" \
		-F "format=varBed" \
		-F "mml=30" \
		-F "col=1" > ${VAR_SEQ_POST}
	@echo "	VAR_INFO_POST	${VAR_SEQ_POST}"

################################################################
## Scan variations via a POST query
varscan_post:
	@echo "	VAR_SEQ_POST	${VAR_SEQ_POST}"
	curl -X POST "${REST_ROOT}/retrieve-variation-seq/${SPECIES}/${ASSEMBLY}" \
		-H "accept: text/plain" \
		-H "Content-Type: multipart/form-data" \
		-F "i=@some_snps_varinfo_GET.varBed;type=text/plain" \
		-F "i_string_type=text" \
		-F "format=varBed" \
		-F "mml=30" \
		-F "col=1" > ${VAR_SCAN_POST}
	@echo "	VAR_INFO_POST	${VAR_SCAN_POST}"

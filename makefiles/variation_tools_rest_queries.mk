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
VARINFO_GET=${VAR_PREFIX}_varinfo_GET.varBed
VARINFO_POST=${VAR_PREFIX}_varinfo_POST.varBed
VARSEQ_POST=${VAR_PREFIX}_varseq_POST.varSeq
VARSCAN_POST=${VAR_PREFIX}_varscan_POST.tsv
SELECTED_MATRICES=${RESULT_DIR}/selected_matrices.tf
SELECTED_MATRICES_GET=${RESULT_DIR}/selected_matrices_GET.tf
list_param:
	@echo "Parameters"
	@echo "	SPECIES			${SPECIES}"
	@echo "	ASSEMBLY		${ASSEMBLY}"
	@echo "	REST_SERVER		${REST_SERVER}"
	@echo "	REST_ROOT		${REST_ROOT}"
	@echo "	VAR_ID_STRING		${VAR_ID_STRING}"
	@echo "	RESULT_DIR		${RESULT_DIR}"
	@echo "	VAR_ID_FILE		${VAR_ID_FILE}"
	@echo "	VARINFO_GET		${VARINFO_GET}"
	@echo "	VARINFO_POST		${VARINFO_POST}"
	@echo "	VARSEQ_POST		${VARSEQ_POST}"
	@echo "	SELECTED_MATRICES	${SELECTED_MATRICES}"
	@echo "	SELECTED_MATRICES_GET	${SELECTED_MATRICES_GET}"

## Run all the queries below
all: dir select_matrices varinfo_get write_snp_query_file varinfo_post varscan_post

################################################################
## Create result directory
dir:
	@mkdir -p ${RESULT_DIR}
	@echo "	RESULT_DIR	${RESULT_DIR}"


################################################################
## Select relevant motifs for this study case
SELECTED_TFS=CEBPA,CEBPB,ELF1,ERG,ETS1,ETV4,FOXA1,FOXA2,GABPA,GATA2,GFI1B,POU2F2,RUNX2,RUNX3,ZNF384
SELECTED_TF_STRING=CEBPA%2CCEBPB%2CELF1%2CERG%2CETS1%2CETV4%2CFOXA1%2CFOXA2%2CGABPA%2CGATA2%2CGFI1B%2CPOU2F2%2CRUNX2%2CRUNX3%2CZNF384
MATRIX_COLLECTION=public_html/motif_databases/JASPAR/Jaspar_2018/nonredundant/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.tf
MATRIX_FILE=
select_matrices:
	retrieve-matrix -v 0 -i ${MATRIX_COLLECTION} -id ${SELECTED_TFS} -o ${SELECTED_MATRICES}
	@echo "	SELECTED_MATRICES	${SELECTED_MATRICES}"

## Select relevant motifs via a GET request to REST web services
MOTIFDB_URL=http%3A%2F%2Frsat.sb-roscoff.fr%2Fmotif_databases%2FJASPAR%2FJaspar_2018%2Fnonredundant%2FJASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.tf
select_matrices_get:
	curl -X GET "${REST_ROOT}/retrieve-matrix/?i_string=${MOTIFDB_URL}&i_string_type=url&v=1&id=${SELECTED_TF_STRING}&content-type=text%2Fplain" -H "accept: application/json" > ${SELECTED_MATRICES_GET}
	@echo "	SELECTED_MATRICES_GET	${SELECTED_MATRICES_GET}"


################################################################
## Send a GET request to RSAT variation-info
VARINFO_GET_STRING=${REST_ROOT}variation-info/${SPECIES}/${ASSEMBLY}?i_string=${VAR_ID_STRING}&i_string_type=text&format=id&col=1&content-type=text%2Fplain
varinfo_get:
	curl -X GET  '${VARINFO_GET_STRING}' -H "accept: text/plain" > ${VARINFO_GET}
	@echo "	VARINFO_GET	${VARINFO_GET}"

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
		-F "col=1" > ${VARINFO_POST}
	@echo "	VARINFO_POST	${VARINFO_POST}"

## NOTE: one of the SNPs is a deletion, but the variation-info option -type does not report it
#		-F "type=SNV" \
# I should add a specific report for variants filtered out because they dont belong to the selected type

################################################################
## Retrieve variation sequences via a POST query
varseq_post:
	@echo "	VARINFO_POST	${VARINFO_POST}"
	curl -X POST "${REST_ROOT}/retrieve-variation-seq/${SPECIES}/${ASSEMBLY}" \
		-H "accept: text/plain" \
		-H "Content-Type: multipart/form-data" \
		-F "i=@${VARINFO_POST};type=text/plain" \
		-F "i_string_type=text" \
		-F "format=varBed" \
		-F "mml=30" \
		-F "col=1" > ${VARSEQ_POST}
	@echo "	VARINFO_POST	${VARSEQ_POST}"

################################################################
## Get a local copy of the background model
BG_FILE=2nt_upstream-noorf_${SPECIES}_${ASSEMBLY}-ovlp-1str.freq.gz
BG_URL=http://rsat.sb-roscoff.fr/data/genomes/Homo_sapiens_GRCh38/oligo-frequencies/${BG_FILE}
BG_FILE=${RESULT_DIR}/

################################################################
## Scan variations via a POST query
VARSCAN_PVAL=0.001
VARSCAN_PVAL_RATIO=10
VARSCAN_TOP_MATRICES=0OB
VARSCAN_TOP_VARIATIONS=0
VARSCAN_MARKOV=2
varscan_post:
	@echo "	VARSEQ_POST		${VARSEQ_POST}"
	@echo "	VARSCAN_MARKOV		${VARSCAN_MARKOV}"
	@echo "	VARSCAN_TOP_MATRICES	${VARSCAN_TOP_MATRICES}"
	@echo "	VARSCAN_TOP_VARIATIONS	${VARSCAN_TOP_VARIATIONS}"
	curl -X POST "${REST_ROOT}/variation-scan/${SPECIES}/${ASSEMBLY}" \
		-H "accept:text/plain" \
		-H "Content-Type: multipart/form-data" \
		-F "v=1" \
		-F "i=@${VARSEQ_POST};type=text/plain" \
		-F "m=@${SELECTED_MATRICES};type=text/plain" \
		-F "m_string_type=url" \
		-F "m_format=transfac" \
		-F "markov_order=${VARSCAN_MARKOV}" \
		-F "top_matrices=${VARSCAN_TOP_MATRICES}" \
		-F "top_variation=${VARSCAN_TOP_VARIATIONS}" \
		-F "uth_pval=${VARSCAN_PVAL}" \
		-F "lth_pval_ratio=${VARSCAN_PVAL_RATIO}" \
		-F "lth_score=1" \
		-F "lth_w_diff=1" > ${VARSCAN_POST}
	@echo "	VARINFO_POST	${VARSCAN_POST}"

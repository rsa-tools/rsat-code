################################################################
## Compute taxon-wide background models

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/taxfreq.mk
V=1

################################################################
## The folowing targets were used to test dyad frequencies directly
## counted on the collection of all promoters for the taxon of
## interest. This takes several hours, so it is not the recommended
## practice, the main goal was to validate the models obtained with
## the program taxon-frequencies

################################################################
## Concatenate all the ortholog sequences in a single file.
## This will be used for estimating bg model + random sequence
## selections.
ALL_ORTHO_SEQ_DIR=data/from_rekins/${TAXON}/orthologs/all_ortho_seq
ALL_ORTHO_SEQ=${ALL_ORTHO_SEQ_DIR}/${TAXON}_all_ortho.fasta
all_ortho_seq:
	@mkdir -p ${ALL_ORTHO_SEQ_DIR}
	cat ${ORTHO_SEQ_DIR}/*.fasta > ${ALL_ORTHO_SEQ}
	@echo ${ALL_ORTHO_SEQ}

################################################################
## Retrieve all upstream sequences from gammaproteobacteria
TAXON=Gammaproteobacteria
TAXON_SPECIES=`supported-organisms -format full | grep ${TAXON} | cut -f 1 | xargs `
ALLUP_DIR=data/sequences/all_promoters
ALLUP_TAXON=${ALLUP_DIR}/${TAXON}_allup${NOORF}.fasta
NOORF=-noorf
taxon_up:
	@mkdir -p ${ALLUP_DIR}
	@echo ${TAXON_SPECIES}
	@rm -f ${ALL_UP_TAXON}
	@for org in ${TAXON_SPECIES}; do \
		echo "	Collecting upstream sequences for $${org}	${NOORF}" ; \
		retrieve-seq ${NOORF} -org $${org} -noov -all >> ${ALLUP_TAXON} ; \
	done
	@echo ${ALLUP_TAXON}


################################################################
## Create a background file for MotifSampler using CreateBackground
##(INCLUSive suite)
BG_OL=`perl -pe 'print ${BG_ORDER}+1'`
##BG_OLIGOS=${RSAT}/data/genomes/${ORG}/
BG_ORDER=5
BG_DIR=bg_models
ms_bg:
	@mkdir -p ${BG_DIR_MS}
	CreateBackgroundModel -f ${ALL_ORTHO_SEQ} -o ${BG_ORDER} -b ${BG_FILE_MS} -n ${TAXON}
	@echo  ${BG_FILE_MS}

################################################################
## Compute dyad frequencies in all promoters
TAXON_DYAD_DIR=bg_models/dyad_frequencies
TAXON_DYADS=${TAXON_DYAD_DIR}/dyads_3nt_sp0-20_upstream${NOORF}_${TAXON}${NOOV}${DYAD_STR}.tab
#TAXON_DYADS=${TAXON_DYAD_DIR}/dyads_${TAXON}_3nt_sp0-20${DYAD_STR}${NOOV}.tab
DYAD_STR=-2str
NOOV=-noov
taxon_bg_dyads:
	@mkdir -p ${TAXON_DYAD_DIR}
	dyad-analysis -v ${V} -timeout 360000 -zeroocc \
		-i  ${ALLUP_TAXON} ${NOOV} ${DYAD_STR} \
		-l 3 -sp 0-20 -return occ,freq \
		-o ${TAXON_DYADS}
	@echo ${TAXON_DYADS}

################################################################
## Compute 6nt frequencies in all promoters
TAXON_OLIGO_DIR=bg_models/oligo_frequencies
TAXON_OLIGOS=${TAXON_OLIGO_DIR}/oligos_${OL}nt_upstream${NOORF}_${TAXON}${NOOV}${OLIGO_STR}.tab
OLIGO_STR=-2str
OL=6
NOOV=-noov
taxon_bg_oligos:
	@mkdir -p ${TAXON_OLIGO_DIR}
	oligo-analysis -v ${V} -zeroocc \
		-i  ${ALLUP_TAXON} ${NOOV} ${OLIGO_STR} \
		-l ${OL} -return occ,freq \
		-o ${TAXON_OLIGOS}
	@echo ${TAXON_OLIGOS}

## Convert dyads with spacing 0 nito a transition table
BG_DIR_TAB=${BG_DIR}/tab
BG_FILE_TAB=${BG_DIR_TAB}/tab_bg_${TAXON}_markov${BG_ORDER}${NOOV}${DYAD_STR}.txt
taxon_bg_dyads2tab:
	@mkdir -p ${BG_DIR_TAB}
	convert-background-model -v 1 -i ${TAXON_DYADS} \
		-from dyads -to tab  \
		-o ${BG_FILE_TAB}
	@echo ${BG_FILE_TAB}

BG_DIR_MEME=${BG_DIR}/meme
BG_FILE_MEME=${BG_DIR_MEME}/meme_bg_${TAXON}_markov${BG_ORDER}${NOOV}${DYAD_STR}.txt
taxon_bg_dyads2meme:
	@mkdir -p ${BG_DIR_MEME}
	convert-background-model -i ${TAXON_DYADS} \
		-from dyads -to meme  \
		-o ${BG_FILE_MEME}
	@echo ${BG_FILE_MEME}

BG_DIR_MS=${BG_DIR}/MotifSampler
BG_FILE_MS=${BG_DIR_MS}/MS_bg_${TAXON}_markov${BG_ORDER}${NOOV}${DYAD_STR}.txt
taxon_bg_dyads2ms:
	@mkdir -p ${BG_DIR_MS}
	convert-background-model -i ${TAXON_DYADS} \
		-from dyads -to MotifSampler  \
		-o ${BG_FILE_MS}
	@echo ${BG_FILE_MS}

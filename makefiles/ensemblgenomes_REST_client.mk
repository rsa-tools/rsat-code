################################################################
## This makefile contains some targets to test the ensemblgenome
## clients developed by Justine Long and Jeanne Cheneby.

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/ensemblgenomes_REST_client.mk

PYTHON=python2.7

## Define parameters
V=2
DATABASE=ensemblgenomes
ORGANISMS_DIR=public_html/data/genomes/ensemblgenome_organisms/
ORGANISMS_ALL=${ORGANISMS_DIR}/organisms_${DATABASE}.tab
ORGANISMS_TAXON=${ORGANISMS_DIR}/organisms_${DATABASE}_${TAXON}.tab
TAXON=83333
SPECIES=escherichia_coli_str_k_12_substr_mg1655


# ################################################################
# ## List of targets
# usage:
# 	@echo "usage: make [-OPT='options'] target"
# 	@echo "implemented targets"
# 	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}


################################################################
## illustration of the different ways to collect organisms

## Get all supported organisms from Ensembl and store them in a file
organisms:
	@echo
	@mkdir -p ${ORGANISMS_DIR}
	@echo "Getting organisms from ${DATABASE}"
	@${PYTHON} python-scripts/ensemblgenomes_REST_client/ensemblgenomes_get_organisms.py -v ${V} \
		--database ${DATABASE} --outfile ${ORGANISMS_ALL} ${OPT}
	@echo "	${ORGANISMS_ALL}"

organisms_for_taxon:
	@echo
	@mkdir -p ${ORGANISMS_DIR}
	@echo "Getting organisms from ${DATABASE}, for taxon ${TAXON}"
	@${PYTHON} python-scripts/ensemblgenomes_REST_client/ensemblgenomes_get_organisms.py -v ${V} \
		--database ${DATABASE} --outfile ${ORGANISMS_TAXON} \
		--taxon ${TAXON}  ${OPT}
	@echo "	${ORGANISMS_TAXON}"

## Collect all substrains of Escherichia coli K12 from ensemblgenomes. 
## Use taxonomic IDs
organisms_k12:
	${MAKE} organisms_for_taxon DATABASE=ensemblgenomes TAXON=83333

## Collect all primates from ensembl
## Use taxon name
organisms_primates:
	${MAKE} organisms_for_taxon DATABASE=ensembl TAXON=Primates


################################################################
## Illustration of the collection of genes to EC relationships.
## Beware: this only works with EnsemblGenomes because the REST
## interface of Ensembl does not support the method LookupGenome.

## Basic request: only return gene-protein-EC (GPE) file for one organism
GENOME_DIR=public_html/data/genomes
SPECIES_DIR=${GENOME_DIR}/${SPECIES}
one_org_gpe:
	@echo
	@echo "Getting genome annotations for species ${SPECIES}"
	@${PYTHON} python-scripts/ensemblgenomes_REST_client/ensemblgenomes_get_annotations.py \
		-v ${V} \
		--outdir ${GENOME_DIR} \
		--species ${SPECIES} ${RETURN} ${OPT}
	@echo "	${GENOME_DIR}/${SPECIES}"



## Full request: export all the available information
one_org_annotations:
	${MAKE} one_org_gpe RETURN='--export_genes --export_names'

################################################################
## Targets for genome management on http://bacteria.rsat.eu
organisms_bact:
	${MAKE} organisms_for_taxon  DATABASE=ensemblgenomes TAXON=Bacteria

## Collect gene annotations for Escherichia coli
annotations_coli:
	${MAKE} DATABASE=ensemblgenomes SPECIES=escherichia_coli_str_k_12_substr_mg1655 one_org_annotations


################################################################
## Targets for genome management on http://fungi.rsat.eu
organisms_fungi:
	${MAKE} organisms_for_taxon  DATABASE=ensemblgenomes TAXON=Fungi

## Collect gene annotations for Brachypodium distachyon
annotations_yeast:
	${MAKE} DATABASE=ensemblgenomes SPECIES=saccharomyces_cerevisiae one_org_annotations

################################################################
## Targets for genome management on http://metazoa.rsat.eu
organisms_metazoa:
	${MAKE} organisms_for_taxon  DATABASE=ensemblgenomes TAXON=Metazoa
	${MAKE} organisms_for_taxon  DATABASE=ensembl TAXON=Metazoa

## Collect gene annotations for Brachypodium distachyon
annotations_elegans:
	${MAKE} DATABASE=ensemblgenomes SPECIES=caenorhabditis_elegans one_org_annotations

################################################################
## Targets for genome management on http://plants.rsat.eu
organisms_plants:
	${MAKE} organisms_for_taxon  DATABASE=ensemblgenomes TAXON=Viridiplantae

## Collect gene annotations for Arabidopsis thaliana
annotations_ara:
	${MAKE} DATABASE=ensemblgenomes SPECIES=arabidopsis_thaliana one_org_annotations

## Collect gene annotations for Brachypodium distachyon
annotations_brachy:
	${MAKE} DATABASE=ensemblgenomes SPECIES=brachypodium_distachyon one_org_annotations


################################################################
## Run consecutively ensemblgenome_get_organisms.py and ensemblgenomes_get_annotations.py for escherichia coli K12
consecutively_test:
	@echo
	${MAKE} organisms_for_taxon  
	@${PYTHON} python-scripts/ensemblgenomes_REST_client/ensemblgenomes_get_annotations.py -v ${V} \
		--outdir ${GENOME_DIR} \
		--species_file ${ORGANISMS_TAXON} ${OPT}
	@echo "	Consecutively test successful "



################################################################
## Re-generate the HTML-formatted index of organisms from the markdown
## file.
index_html:
	@echo
	@echo "Generating HTML index of organisms"
	${MAKE} ${GENOME_DIR}/index.html
	@echo "Index file	${GENOME_DIR}/index.html"


## Define a rule to convert markdown into html
%.html: %.md
	@pandoc --from markdown --to html $< -o $@

%.docx: %.md
	@pandoc --from markdown --to docx $< -o $@




################################################################
## Download GTF files from ensemblgenomes
GTF_GROUP=plants
GTF_SPECIES=chlamydomonas_reinhardtii
GTF_URL=ftp.ensemblgenomes.org/pub/${GTF_GROUP}/release-28/gtf/${GTF_SPECIES}/
GTF_PATH=${RSAT}/downloads/${GTF_URL}
download_gtf:
	@echo
	@echo "Downloading GTF file from	ftp://${GTF_URL}"
	cd ${RSAT}/downloads; \
	wget  -rNL ftp://${GTF_URL}
	@echo "	${GTF_PATH}"

GTF_GZ=`ls -1 ${GTF_PATH}/*.gtf.gz`
ORG_ID=${GTF_SPECIES}
PARSE_DIR=${RSAT}/data/ensemblgenomes/${ORG_ID}/genome/
parse_gtf:
	@echo
	@echo "Parsing GTF file	${GTF_GZ}"
	@echo "ORG_ID	${ORG_ID}"
	parse-gtf -v ${V} -i ${GTF_GZ} -o ${PARSE_DIR}
	@echo "	${PARSE_DIR}"


GTF_TASK=download_gtf parse_gtf
gtf_ara:
	make -f makefiles/ensemblgenomes_REST_client.mk   GTF_SPECIES=arabidopsis_thaliana GTF_GROUP=plants ${GTF_TASK}

gtf_worm:
	make -f makefiles/ensemblgenomes_REST_client.mk   GTF_SPECIES=caenorhabditis_elegans GTF_GROUP=metazoa ${GTF_TASK}

gtf_fly:
	make -f makefiles/ensemblgenomes_REST_client.mk   GTF_SPECIES=drosophila_melanogaster GTF_GROUP=metazoa ${GTF_TASK}

gtf_yeast:
	make -f makefiles/ensemblgenomes_REST_client.mk   GTF_SPECIES=saccharomyces_cerevisiae GTF_GROUP=fungi ${GTF_TASK}

gtf_ecoli:
	make -f makefiles/ensemblgenomes_REST_client.mk   GTF_SPECIES=escherichia_coli_str_k_12_substr_mg1655_gca_000801205_1 GTF_GROUP=bacteria ${GTF_TASK}

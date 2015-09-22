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


################################################################
## List of targets
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}


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
OUTDIR=results
SPECIES_DIR=${OUTDIR}/${SPECIES}
one_org_gpe:
	@echo
	@echo "Getting genes-proteins-ECs (GPE) for species ${SPECIES}"
	@${PYTHON} python-scripts/ensemblgenomes_REST_client/ensemblgenomes_get_annotations.py -v ${V} \
		--outdir ${OUTDIR} \
		--species ${SPECIES} ${RETURN} ${OPT}
	@echo "	${OUTDIR}/${SPECIES}"



## Full request: export all the available information
one_org_annotations:
	${MAKE} one_org_gpe RETURN=--export_all


################################################################
## Targets for genome management on http://plants.rsat.eu
organisms_plants:
	${MAKE} organisms_for_taxon  DATABASE=ensemblgenomes TAXON=Viridiplantae

## Install one illustrative plant genome
get_brachypodium:
	${MAKE} SPECIES=brachypodium_distachyon


################################################################
## Run consecutively ensemblgenome_get_organisms.py and ensemblgenomes_get_annotations.py for escherichia coli K12
consecutively_test:
	@echo
	${MAKE} organisms_for_taxon  
	@${PYTHON} python-scripts/ensemblgenomes_REST_client/ensemblgenomes_get_annotations.py -v ${V} \
		--outdir ${OUTDIR} \
		--species_file ${ORGANISMS_TAXON} ${OPT}
	@echo "	Consecutively test successful "



################################################################
## Re-generate the HTML-formatted index of organisms from the markdown
## file.
index_html:
	@echo
	@echo "Generating HTML index of organisms"
	${MAKE} ${OUTDIR}/index.html
	@echo "Index file	${OUTDIR}/index.html"


## Define a rule to convert markdown into html
%.html: %.md
	@pandoc --from markdown --to html $< -o $@

%.docx: %.md
	@pandoc --from markdown --to docx $< -o $@





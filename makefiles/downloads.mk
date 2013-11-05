############################################################
#
# $Id: downloads.mk,v 1.44 2013/06/12 15:49:07 rsat Exp $
#
# Time-stamp: <2003-10-09 14:02:21 jvanheld>
#
############################################################

DOWNLOAD_DIR=${RSAT}/downloads

################################################################
#### programs
MAKEFILE=${RSAT}/makefiles/downloads.mk
MAKE = make -sk -f ${MAKEFILE}

DATE = `date +%Y%m%d_%H%M%S`
LOGFILE=-o logs/wget_${DATE}_log.txt
WGET=wget --passive-ftp -np -rNL ${LOGFILE}
# #WGET = wget -rNL -o logs/wget_${DATE}_log.txt
RSYNC = rsync -ruptvl -e ssh

### target
usage:
	@echo "usage: make [-OPT='options'] target"
	@echo "implemented targets"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}

update: genbank kegg expasy ebi_genomes taxonomy ensembl go prosite bind dbtbs

################################################################
#
# Complete genomes at EBI
#
GENOME_REVIEWS=ftp://ftp.ebi.ac.uk/pub/databases/genome_reviews/
ebi_genome_reviews:
	@echo importing genomereviews from EBI
	${WGET} ${GENOME_REVIEWS}

EBI_GENOMES=ftp://ftp.ebi.ac.uk/pub/databases/genomes/
ebi_genomes:
	@echo importing genomes from EBI
	${WGET} --accept=.embl.Z --accept .con --accept .prot ${EBI_GENOMES}

################################################################
#
# Genbank genome repository
#
GENBANK_DIRS =					\
	genomes					\
	genbank/genomes				\
	refseq 

NCBI_DIR=Fungi/Saccharomyces_cerevisiae_uid128
NCBI_EXCLUDE=	\
		--exclude '*_alt_*'							\
		--exclude 'lproks*_*'							\
		--exclude Assembled_chromosomes						\
		--exclude IDS								\
		--exclude SNP								\
		--exclude maps								\
		--exclude CLONEEND							\
		--exclude CLUSTERS							\
		--exclude FOSMIDS							\
		--exclude SARS								\
		--exclude TOOLS								\
		--exclude WGS_BACTERIA							\
		--exclude ARCHIVE							\
		--exclude BACENDS							\
		--exclude CLONEEND							\
		--exclude Bacteria.OLD							\
		--exclude '*.tar.gz'						


one_ncbi_dir_from_mirror:
	@echo "DOWNLOAD_DIR	${DOWNLOAD_DIR}"
	@mkdir -p ${DOWNLOAD_DIR}/ftp.ncbi.nih.gov/genomes/${NCBI_DIR}
	rsync ${NCBI_EXCLUDE}						\
		-av ${OPT} rsync://bio-mirror.net/biomirror/ncbigenomes/${NCBI_DIR}/*	\
		${DOWNLOAD_DIR}/ftp.ncbi.nih.gov/genomes/${NCBI_DIR}/
	@echo "NCBI directory downloaded to ${DOWNLOAD_DIR}/ftp.ncbi.nih.gov/genomes/${NCBI_DIR}/"

one_ncbi_dir:
	${MAKE} one_ncbi_dir_from_mirror

ncbi:
	rsync --delete	${NCBI_EXCLUDE}						\
		-avz rsync://bio-mirror.net/biomirror/ncbigenomes/*	\
		ftp.ncbi.nih.gov/genomes/

NCBI_GENOMES_FTP=ftp://ftp.ncbi.nih.gov/genomes
#NCBI_DIR=Fungi/Saccharomyces_cerevisiae_uid128
one_ncbi_dir_wget:
	@mkdir -p logs
	@echo "${DATE}	updating dir	$${NCBI_DIR}" >> wget_updates.txt
	${WGET}							\
		--exclude-directories 'Bacteria.OLD'		\
		--exclude-directories ARCHIVE			\
		--exclude-directories BACENDS			\
		--accept=gpff --accept=gpff.gz			\
		--accept=gbk --accept=README --accept=gbff	\
		--accept=gaa --accept=faa --accept=gbk.gz	\
		--accept=README.gz --accept=gbff.gz		\
		--accept=gaa.gz --accept=faa.gz			\
		"${NCBI_GENOMES_FTP}/${NCBI_DIR}" 
	@echo "${DATE}	updated dir	$${NCBI_DIR}" >> wget_updates.txt

#one_genbank_dir:
#	${MAKE} one_ncbi_dir_from_mirror

#genbank_ori:
#	@echo "${DATE}	starting to update	${GENBANK_DIRS}" >> wget_updates.txt; 
#	@for dir in ${GENBANK_DIRS} ; do		\
#		make one_genbank_dir_ori NCBI_DIR=$${dir} ;	\
#	done
#	@${WGET} ${NCBI_GENOMES_FTP}/acc
#	@${WGET} ${NCBI_GENOMES_FTP}/README
#	@${WGET} ${NCBI_GENOMES_FTP}/Bacteria/README
#	@${WGET} ${NCBI_GENOMES_FTP}/Bacteria/accessions
#	@echo "${DATE}	finished to update	${GENBANK_DIRS}" >> wget_updates.txt; 




################################################################
#
# Download KEGG databases
#
# does not take the full genome sequences

kegg: kegg_ligand kegg_genomes kegg_pathways


KEGG_FTP=ftp://ftp.genome.ad.jp/pub/kegg/

KEGG_GENOMES=${KEGG_FTP}/genomes/
kegg_genomes:
	@mkdir -p logs
	${WGET} -X sequences.old,sequences,sequences.weekly.last.tar.Z,genes.weekly.last.tar.Z ${KEGG_GENOMES}


KEGG_LIGAND=${KEGG_FTP}/ligand/
kegg_ligand:
	@mkdir -p logs
	${WGET} ${KEGG_LIGAND}

KEGG_PATHWAYS=${KEGG_FTP}/pathways/
kegg_pathways:
	@mkdir -p logs
	${WGET} ${KEGG_PATHWAYS}


################################################################
# Databases at Expasy :
# - swiss-prot
# - trembl
# - ENZYME
################################################################
EXPASY=ftp://ftp.expasy.org/databases
EXPASY_DIRS=						\
	README						\
	enzyme						\
	sp_tr_nrdb					\
	swiss-prot/release_compressed		
expasy:
	@mkdir -p logs
	for dir in ${EXPASY_DIRS} ; do		\
		${WGET} ${EXPASY}/$${dir} ;	\
	done

UNIPROT=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/
uniprot:
	${WGET} ${UNIPROT}/uniprot_sprot.*.gz

################################################################
#
# Bacilus subtilis transcription factor database

DBTBS=http://elmo.ims.u-tokyo.ac.jp/dbtbs/tfac/
dbtbs:
	${WGET} ${DBTBS}


################################################################
#
#   SGD database

SGD_FTP=ftp://genome-ftp.stanford.edu
SGD_dir=/pub/yeast/data_download
SGD_excl=${SGD_dir}/obsolete_files/,${SGD_dir}/sequence/,${SGD_dir}/chromosomal_feature/archive/,${SGD_dir}/gene_registry/archive/,${SGD_dir}/literature_curation/archive/,${SGD_dir}/oracle_schema/archive/,${SGD_dir}/protein_info/archive/,${SGD_dir}/sequence_similarity/archive/,${SGD_dir}/systematic_results/archive/,${SGD_dir}/systematic_results/SAGE/archive/
sgd:
	@mkdir -p logs
	${WGET} --accept '*.gbf' ftp://genome-ftp.stanford.edu/pub/yeast/data_download/sequence/NCBI_genome_source/
	${WGET} ftp://genome-ftp.stanford.edu/pub/yeast/data_download/sequence/genomic_sequence/chromosomes/fasta/
	${WGET} -X ${SGD_excl} ${SGD_FTP}${SGD_dir}/

################################################################
#
# Alternative yeast genomes
#
yeasts: yeast_mti

## Fungal genomes at Stanford
SGD_FUNGAL_GENOMES=ftp://genome-ftp.stanford.edu/pub/yeast/data_download/sequence/fungal_genomes/
fungal_genomes:
	${WGET} -X '*archive' -X '*.gcg'  ${SGD_FUNGAL_GENOMES}

SPECIES=saccharomyces_bayanus
one_fungal_genome_duke:
	@echo "Getting genome from duke	${SPECIES}"
	${WGET} http://fungal.genome.duke.edu/annotations/${SPECIES}/gff/
	${WGET} http://fungal.genome.duke.edu/annotations/${SPECIES}/nt/${SPECIES}'*'.nt.gz
	${WGET} http://fungal.genome.duke.edu/annotations/${SPECIES}/gff/${SPECIES}.SGD.gff3.gz
	${WGET} http://fungal.genome.duke.edu/annotations/${SPECIES}/gff/${SPECIES}.SNAP.gff3.gz

DUKE_SPECIES= \
	ashbya_gossypii \
	aspergillus_fumigatus \
	aspergillus_nidulans \
	aspergillus_oryzae \
	aspergillus_terreus \
	botrytis_cinerea \
	candida_albicans \
	candida_dubliniensis \
	candida_glabrata \
	candida_guilliermondii \
	candida_lusitaniae \
	candida_tropicalis \
	chaetomium_globosum \
	coccidioides_immitis \
	coprinus_cinereus \
	cryptococcus_neoformans_H99 \
	cryptococcus_neoformans_JEC21 \
	cryptococcus_neoformans_R265 \
	cryptococcus_neoformans_WM276 \
	debaryomyces_hansenii \
	fusarium_graminearum \
	fusarium_verticillioides \
	histoplasma_capsulatum_186R \
	kluyveromyces_lactis \
	kluyveromyces_waltii \
	magnaporthe_grisea \
	neurospora_crassa \
	phanerochaete_chrysosporium \
	pneumocystis_carnii \
	podospora_anserina \
	rhizopus_oryzae \
	saccharomyces_bayanus \
	saccharomyces_castellii \
	saccharomyces_cerevisiae_rm11-1a_1 \
	saccharomyces_cerevisiae_s288c \
	saccharomyces_cerevisiae_yjm789 \
	saccharomyces_kluyveri \
	saccharomyces_kudriavzevii \
	saccharomyces_mikatae_MIT \
	saccharomyces_paradoxus \
	schizosaccharomyces_pombe \
	sclerotinia_sclerotiorum \
	stagonospora_nodorum \
	trichoderma_reesei \
	uncinocarpus_reesii \
	ustilago_maydis \
	yarrowia_lipolytica
fungal_genomes_duke:
	@for sp in ${DUKE_SPECIES} ; do \
		${MAKE} one_fungal_genome_duke SPECIES=$${sp} ; \
	done

################################################################
# Several Saccharomyces genomes from the MIT 
#

#MIT_YEAST=http://www-genome.wi.mit.edu/personal/manoli/yeasts/
MIT_YEAST_DIRS=											\
	http://www-genome.wi.mit.edu/ftp/pub/annotation/fungi/comp_yeasts/			\
	http://www-genome.wi.mit.edu/annotation/fungi/comp_yeasts/downloads.html	
yeast_mit:
	@for dir in ${MIT_YEAST_DIRS} ; do	\
		echo "getting $${dir}";		\
		${WGET} $${dir} ;		\
	done


################################################################
## Genomes from Genolevure
#GENOLEVURE=http://cbi.labri.fr/Genolevures/download.php
GENOLEVURE=http://cbi.labri.fr/Genolevures/raw/seq/
genolevure:
	${WGET}  ${GENOLEVURE}

################################################################
#
# Yeast genome from the MIPS
#
################################################################
MIPS_YEAST = ftp://ftpmips.gsf.de/yeast/
DIRS =						\
	CYGD					\
	catalogues				\
	chri					\
	chrii					\
	chriii					\
	chriv					\
	chrix					\
	chrv					\
	chrvi					\
	chrvii					\
	chrviii					\
	chrx					\
	chrxi					\
	chrxii					\
	chrxiii					\
	chrxiv					\
	chrxv					\
	chrxvi					\
	eurofan					\
	mito

yeast_from_mips:
	for d in ${DIRS}; do					\
		${MAKE} DIR=$${d} yeast_from_mips_one_dir ;	\
	done

DIR=CYGD
yeast_from_mips_one_dir:
	${WGET} ${MIPS_YEAST}/${DIR}

yeast_from_mips_all:
	${WGET} ${MIPS_YEAST}

################################################################
#
# Caenorhabditis genome from euGenes
#
################################################################

WORM_EUGENE=ftp://iubio.bio.indiana.edu/eugenes/worm/features/
worm_from_eugene:
	${WGET} ${WORM_EUGENE}

WORM_WORMBASE =						\
	ftp://ftp.wormbase.org/pub/wormbase/DNA_DUMPS/	\
	ftp://ftp.wormbase.org/pub/wormbase/GENE_DUMPS/
worm_from_wormbase:
	${WGET} ${WORM_WORMBASE}

WORM_SANGER = ftp://ftp.sanger.ac.uk/pub/C.elegans_sequences/ 
worm_from_sanger:
	${WGET} ${WORM_SANGER}

################################################################
#
# Candida albicans genome from Stanford
#
################################################################
candida: candida_pasteur candida_sgd

CANDIDA = ftp://cycle.stanford.edu/pub/projects/candida/
candida_sgd:
	${WGET} ${CANDIDA}


CANDIDA_PASTEUR=ftp://ftp.pasteur.fr/pub/GenomeDB/CandidaDB/FlatFiles/
candida_pasteur:
	${WGET} ${CANDIDA_PASTEUR}

################################################################
#
# Schizosaccharomyces pombe
#
################################################################
POMBE_CONTIGS = ftp://ftp.sanger.ac.uk/pub/yeast/Pombe/CONTIGS
POMBE_GENOME=ftp://ftp.sanger.ac.uk/pub/yeast/sequences/pombe/pombe.dbs.Z
POMBE_PROTEINS=ftp://ftp.sanger.ac.uk/pub/yeast/sequences/pombe/
POMBE_PAULUS = jvanheld@paulus:/permeke/dsk2/genomics/S_pombe/ftp.sanger.ac.uk .
pombe:
	${WGET} ${POMBE_CONTIGS}
	${WGET} ${POMBE_GENOME}
	${WGET} ${POMBE_PROTEINS}
#	${RSYNC} ${POMBE_PAULUS}

################################################################
#
# Ralstonia eutropha
#
################################################################
RALSTONIA = http://genome.ornl.gov/microbial/reut/04dec00/
ralstonia:
	${WGET} ${RALSTONIA}


scpd:
	${WGET} http://cgsigma.cshl.org/jian/


################################################################
#
# Homo sapiens from UCSC

HUMAN_UCSC=ftp://genome.ucsc.edu/goldenPath/14nov2002/
#HUMAN_UCSC=http://genome.cse.ucsc.edu/goldenPath/14nov2002/
HUMAN_UCSC_FILES=				\
	bigZips/chromFa.zip			\
	bigZips/contigFa.zip			\
	bigZips/chromFaMasked.zip		\
	bigZips/contigFaMasked.zip		\
	bigZips/upstream5000.zip		\
	database/genename.sql			\
	database/refGene.txt.gz			\
	database/refGene.sql			\
	database/refLink.sql			\
	database/refLink.txt.gz			\
	database/geneName.sql			\
	database/geneName.txt.gz		\
	database/geneid.sql			\
	database/geneid.txt.gz			\
	database/genieAlt.txt			\
	database/genieAltPep.txt 

human_from_ucsc:
	echo "downloading ${HUMAN_UCSC}"
	${WGET} ${HUMAN_UCSC}

#human_from_ucsc:
#	@for f in ${HUMAN_UCSC_FILES} ; do	\
#		echo "downloading ${HUMAN_UCSC}/$${f}" ;	\
#		${WGET} ${HUMAN_UCSC}/$${f} ;	\
#	done


################################################################
#
# Rice 

RICE_GENOME=ftp://ftp.dna.affrc.go.jp/pub/RiceGAAS/current/*.tar.gz
rice:
	${WGET} ${RICE_GENOME}



################################################################
#
# List new bacteria from a wget file
#

UPDATE=200305
new_bacteria:
	cat logs/wget_${UPDATE}*_log.txt		\
		| grep saved			\
		| grep gbk			\
		| grep Bacteria			\
		| perl -pe 's|Bacteria/|\t|'	\
		| cut -f 2			\
		| perl -pe 's|/|\t|'		\
		| cut -f 1			\
		| sort -u

################################################################
#
# Microarray data
#

CELL_CYCLE=http://genome-www.stanford.edu/cellcycle/data/rawdata/individual.html
spellman_cell_cycle:
	${WGET} --accept=.html,.out		\
		${CELL_CYCLE}


################################################################
#
# NCBI taxonomy
#

TAXONOMY=ftp://ftp.ncbi.nih.gov/pub/taxonomy/
taxonomy:
	${WGET} ${TAXONOMY}

#########################################
#  BIND database
#########################################
BIND_URL=ftp.blueprint.org
BIND_FTP=ftp://${BIND_URL}
BIND_dir=/pub/BIND/current
BIND_excl=${BIND_dir}/MMDBBIND
bind:
	@mkdir -p logs
	${WGET} -X ${BIND_excl} ${BIND_FTP}${BIND_dir}/
	chmod -R g+w ${BIND_URL}

#########################################
#   PROSITE database
#########################################
PROS_URL=ftp.expasy.org
PROS_FTP=ftp://${PROS_URL}
PROS_DIR=/databases/prosite
PROS_excl=""
prosite:
	@mkdir -p logs
	${WGET} -X ${PROS_excl} ${PROS_FTP}${PROS_DIR}/
	chmod -R g+w ${PROS_URL}${PROS_DIR}

#########################################
#   Gene Ontology database
#########################################
GO_URLS= \
	www.geneontology.org/ontology/gene_ontology.obo \
	www.geneontology.org/doc/GO.terms_and_ids

GO_URL=www.geneontology.org/ontology/
GO_COMPLETE=http://archive.godatabase.org/latest/
GO_HTTP=http://${GO_URL}
GO_FTP=ftp://ftp.geneontology.org/pub/go/
go:
	for url in ${GO_URLS} ; do ${MAKE} go_one_url GO_URL=$${url}; done

go_one_url:
	@mkdir -p logs
	wget -np -r -l 1 -N ${GO_HTTP}
	chmod -R g+w ${GO_URL}

## Gene Ontology Associations
goa_ebi:
#	wget -rNL ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/gene_association.goa_human.gz
	wget -rNL http://www.geneontology.org/gene-associations/gene_association.goa_human.gz

################################################################
## JASPAR

#JASPAR=http://jaspar.cgb.ki.se/DOWNLOAD/
JASPAR=http://jaspar.genereg.net/html/DOWNLOAD/
jaspar:
	@mkdir -p logs
	@echo "${DATE}	updating dir	${JASPAR}" 
	${WGET} ${JASPAR}
	@echo "${DATE}	updated dir	${JASPAR}"

################################################################
#
# Homo sapiens from ENSEMBL

HUMAN_ENSEMBL=ftp://ftp.ensembl.org/pub/current_human/data/
human_from_ensembl:
	${WGET} ${HUMAN_ENSEMBL}/mysql
	${WGET} ${HUMAN_ENSEMBL}/flatfiles/genbank

ENSEMBL_DIRS=					\
	current_mouse/data/mysql		\
	current_mouse/data/flatfiles/genbank	\
	current_human/data/mysql		\
	current_human/data/flatfiles/genbank
ensembl:
	@for dir in ${ENSEMBL_DIRS} ; do				\
		echo "downloading directory $${dir} from ENSEMBL";	\
		${MAKE} one_ensembl_dir ENSEMBL_DIR=$${dir} ; 		\
	done

################################################################
#
# ENSEMBL
#
ENSEMBL_BASE=ftp://ftp.ensembl.org/pub/
ENSEMBL_DIR=current_celegans/data/flatfiles/genbank/
one_ensembl_dir:
	@echo "Getting ENSEMBL dir" ${ENSEMBL_DIR} ${LOG_FILE}
	${WGET} ${ENSEMBL_BASE}/${ENSEMBL_DIR}

################################################################
## Anopheles 
anopheles:
	${MAKE} one_ensembl_dir ENSEMBL_DIR=current_mosquito/data/flatfiles/genbank/

MACHIN=worm
current_machin:
	${MAKE} one_ensembl_dir ENSEMBL_DIR=current_${MACHIN}/data/flatfiles/genbank/



################################################################
## Drosophila genomes from UCSC
UCSC_DIR=dp3
UCSC_URL=http://hgdownload.cse.ucsc.edu/goldenPath/${UCSC_DIR}/
UCSC_FILE=bigZips/chromFa.zip 
UCSC_FILES= \
	database/blastDm2FB.sql \
	database/blastDm2FB.txt.gz \
	database/cds.sql \
	database/cds.txt.gz \
	database/all_mrna.sql  \
	database/all_mrna.txt.gz  \
	database/description.sql  \
	database/description.txt.gz  \
	database/organism.sql   \
	database/organism.txt.gz   \
	database/geneName.sql  \
	database/geneName.txt.gz  \
	database/refLink.sql \
	database/refLink.txt.gz \
	database/source.sql \
	database/source.txt.gz \
	database/tableDescriptions.sql \
	database/tableDescriptions.txt.gz \
	bigZips/chromFa.zip \
	bigZips/chromFaMasked.zip
ucsc_one_dir:
	@for f in ${UCSC_FILES} ; do  \
		${MAKE} ucsc_one_file UCSC_FILE=$${f}; \
	done

ucsc_one_file:
	${WGET} ${UCSC_URL}/${UCSC_FILE}

UCSC_DROSO_DIRS=dp3 dm2 droSim1 droSec1 droYak2 droEre1 droAna2 droPer1 droVir2 droGri1 apiMel anoGam1
ucsc_droso:
	@for d in ${UCSC_DROSO_DIRS} ; do \
		${MAKE} ucsc_one_dir UCSC_DIR=$${d}; \
	done

################################################################
## Drosophila genomes from Flybase
## Obsolete: these genomes are distibuted in a better shape at UCSC
pseudoobscura:
	${WGET} ftp://flybase.net/genomes/Drosophila_pseudoobscura/current/gff/
	${WGET} ftp://flybase.net/genomes/Drosophila_pseudoobscura/current/dna/


DROSO=dpse
one_drosophila:
#	${WGET} http://rana.lbl.gov/drosophila/sechellia.html
	wget -rNL http://insects.eugenes.org/species/data/${DROSO}/gff
#	wget -rNL http://insects.eugenes.org/species/data/${DROSO}/fasta/

ALL_DROSO=dpse dmel dmel2 
all_drosophila:
	for d in ${ALL_DROSO} ; do \
		${MAKE} one_drosophila DROSO=$${d} ; \
	done


################################################################
## Plasmodium genome
PLASMO_SANGER=ftp://ftp.sanger.ac.uk/pub/pathogens/Plasmodium/falciparum/3D7/genome.version.2.1.1
plasmo_sanger:
	${WGET} ${PLASMO_SANGER}


################################################################
## Ciona intestinalis genome
ciona:
	${WGET} ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v2.0/
#	${WGET} ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v2.0/ciona050324.unmasked.fasta.gz
#	${WGET} ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v1.0/ciona.fasta.gz
#	${WGET} ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v1.0/ciona_genes.gff.gz
#	${WGET} ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v1.0/ciona_utr.fasta.gz
#	${WGET} ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v1.0/ciona.prot.fasta.gz
#	${WGET} ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v1.0/ciona.mrna.fasta.gz
	${WGET} http://crfb.univ-mrs.fr/aniseed/downloads/ci-GO-annotations.gz
	${WGET} http://crfb.univ-mrs.fr/aniseed/downloads/ci-interpro-annotations.gz
	${WGET} http://crfb.univ-mrs.fr/aniseed/downloads/ci-ortholog.gz

################################################################
## Leishmania

LEISHMANIA_DIR=ftp.sanger.ac.uk/pub/databases/L.major_sequences/DATASETS
leishmania_sanger:
	${WGET} --exclude 'ARCHIVE' \
		--exclude 'EMBL_SUBMISSIONS' \
		--exclude 'T_cruzi' \
		--exclude 'tribemcl_leishgenome.embl' \
		${LEISHMANIA_DIR}

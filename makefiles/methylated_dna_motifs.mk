################################################################
## Testers for the conversion and analysis of motifs with an alphabet
## represnting DNA with methylated cytosine (Cytomod alphabet).
##
## Source data from Coby Viner in MEME format
##
## Reference:
##      http://www.biorxiv.org/content/early/2016/03/15/043794
##
## ALPHABET "DNA with covalent modifications"
## A "Adenine" 8510A8 ~ T "Thymine" A89610
## C "Cytosine" A50026 ~ G "Guanine" 313695
## h "5-Hydroxymethylcytosine" F46D43 ~ 2 "Guanine:5-Hydroxymethylcytosine" 74ADD1
## m "5-Methylcytosine" D73027 ~ 1 "Guanine:5-Methylcytosine" 4575B4
## ? = ACGThm12
## N = ACGT
## X = ACGT
## V = ACG
## H = ACT
## D = AGT
## B = CGT
## z = Chm
## 9 = G12
## M = AC
## R = AG
## W = AT
## S = CG
## Y = CT
## K = GT
## x = hm
## 7 = 12
## END ALPHABET

MATRIX_TYPE=cytomod

include ${RSAT}/makefiles/util.mk
MAKEFILE=${RSAT}/makefiles/methylated_dna_motifs.mk

MOTIF_FOLDER=data/motifs
MOTIFS=${MOTIF_FOLDER}/mod_test_motifs
OUT_FORMAT=transfac

convert_from_meme:
	@echo "Converting methyl motifs from MEME to TRANSFAC format"
	@echo "	Input: ${MOTIFS}.meme"
	@convert-matrix -v 1 -residue_type ${MATRIX_TYPE} -from meme -to ${OUT_FORMAT} \
		-i ${MOTIFS}.meme \
		-o ${MOTIFS}.${OUT_FORMAT}
	@echo "	${MOTIFS}.${OUT_FORMAT}"

#		-return counts,parameters \

convert_from_transfac:
	@echo "Converting methyl motifs from  TRANSFAC to ${OUT_FORMAT} format"
	@echo "	Input: ${MOTIFS}.transfac"
	@convert-matrix -v 1 -residue_type ${MATRIX_TYPE} -from transfac -to ${OUT_FORMAT} \
		-i ${MOTIFS}.transfac \
		-return counts,parameters \
		-o ${MOTIFS}_converted.${OUT_FORMAT}
	@echo "	${MOTIFS}_converted.${OUT_FORMAT}"

#produce_logo:
#	@echo "Generating a logo with the expanded alphabet"
#    weblogo --fin ${MOTIFS}.${OUT_FORMAT} --datatype transfac --format pdf --show-yaxis YES --show-xaxis YES --errorbars YES --fout Methyl_motif_logo.pdf --color '#8510A8' A 'Adenine' --color '#A89610' T 'Thymine' --color '#A50026' C 'Cytosine' --color '#313695' G 'Guanine' --color '#D73027' m '5-Methylcytosine' --color '#4575B4' 1 'Guanine:5-Methylcytosine' --color '#F46D43' h '5-Hydroxymethylcytosine' --color '#74ADD1' 2 'Guanine:5-Hydroxymethylcytosine' --alphabet 'ACGTm1h2' 
#	@echo "	${MOTIFS}.${OUT_FORMAT}"




#PREFIX_CLUSTERING=results/cytomod_motifs/clustering_cytomod
#cluster_motifs:
#	@echo "Clustering motifs with the expanded alphabet"
#	matrix-clustering  -v 2 -matrix_format ${OUT_FORMAT} -matrix 'cytomod' ${MOTIFS}.${OUT_FORMAT} -hclust_method average -calc sum -title "Cytomod_motifs" -metric_build_tree Ncor -lth w 5 -lth cor 0.6 -lth Ncor 0.4 -label_in_tree name -return json,heatmap  -o ${PREFIX_CLUSTERING}
#	@echo "	${PREFIX_CLUSTERING}_SUMMARY.html"


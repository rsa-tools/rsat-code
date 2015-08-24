################################################################
## R-script to complemment matrix-quality analysis
## this script is called by compare-qualities

## Authors
## Samuel Collombet <samuel.collombet@ens.fr>
## Morgane Thomas-Chollier <mthomas@biologie.ens.fr>
## Alejandra Medina-Rivera <amedina@lcg.unam.mx>

library(RColorBrewer)
library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)

## Get the prefix list of files to be analyzed, this list is passed by compare-qualities
mtx.quality.prefix.list.file <- args[1]

## For debugging:
##mtx.quality.prefix.list.file <-"../results/matrix_quality/20150428/zoo_chip_enrichment_quality_prefix_files.txt"

mtx.quality.prefix.list <-as.vector(as.matrix( read.table(file=mtx.quality.prefix.list.file, header=FALSE, stringsAsFactors=FALSE)))
## directory with all the matrix quality results
##DirList <-  list.files("mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/matrixQuality", pattern="mFuzz*", full.names=TRUE)

# Import some lists of selcted motifs
##SelectedMotifs <- read.table("mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/PeakMotifs_SelectedMotifs", header=T, sep ="\t")
##SelectedMotifsStringent <- read.table("mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/PeakMotifs_SelectedMotifs_stringent.tsv", header=T, sep ="\t")

# import all data in one: for each matrix-quality run, import the ...compare-scores.tab file (row=pVal, col=sequence NWD), get the max of each column (MNWD for each sequence) and make a list of all that

ListAll<- lapply(mtx.quality.prefix.list , FUN=function(dir){
    dir <- dirname(dir)
    as.data.frame(t(as.matrix(as.data.frame(lapply(list.files(dir, pattern=".+_nwd_.+_compare-scores.tab", full.names=TRUE), FUN=function(file){
        print (paste("Reading in", file))
        MNWD <- apply(read.table(file, header=T, row.names=1, sep ="\t",na.strings="<NULL>",comment.char="", stringsAsFactors=FALSE), 2, FUN=function(x){
            max(x, na.rm=T)
        })
        ## This is a vector with the max NWD values
        return(MNWD)
        
	})))))
})

# remove part of the motif name to make it more readable
# # second line depens on the name of the file that was given for the output in matrix quality... ALE:Not sure what the rownames should be

ListAll<-lapply(ListAll,function(DF) {rownames(DF) <- paste("",seq(1)); DF})
#ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8_","",colnames(DF)); DF})


ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("_score_distrib","",colnames(DF)); DF})
ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("_nwd.tab","",colnames(DF)); DF})

# reformat the data
## ALE: Since the rownames setting is not working properly this transformation is not correct

TableAll<-do.call("cbind",ListAll)
TableAll$cluster <- rownames(TableAll)
TableAll_melt<-melt(TableAll, id.vars="cluster")
TableAll_melt$cluster<-factor(TableAll_melt$cluster, level=(unique(TableAll_melt$cluster)))

# put a higher threeshold on MNWD for scaling
TableAll_melt$value[TableAll_melt$value>0.8] <- 0.8
# heatmap
heatmap<-ggplot(TableAll_melt, aes(cluster,variable)) +
		geom_tile(aes(fill=value)) +
		scale_fill_gradientn(colours=colorRampPalette(brewer.pal(9, "Greys"))(20)) +
		theme(axis.ticks=element_blank()) + labs(x="clusters",y="") + scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0))

heatmap
#ggsave(file="mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/matrixQuality/MNWD_Heatmap.pdf",plot=heatmap, height=5, width=20)



# same for a subset of selected motifs
## TableAll_selected<-TableAll[,as.vector(SelectedMotifs$MotifId)]
## colnames(TableAll_selected)<-SelectedMotifs$KnownMotif
## TableAll_selected$cluster <- rownames(TableAll_selected)
## TableAll_selected_melt<-melt(TableAll_selected, id.vars="cluster")
## TableAll_selected_melt$cluster<-factor(TableAll_selected_melt$cluster, level=(unique(TableAll_selected_melt$cluster)))

## heatmap<-ggplot(TableAll_selected_melt, aes(cluster,variable)) +
## 		geom_tile(aes(fill=value)) +
## 		scale_fill_gradientn(colours=colorRampPalette(brewer.pal(9, "Greys"))(20)) +
## 		theme(axis.ticks=element_blank()) + labs(x="clusters",y="") + scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0))
## ggsave(file="mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/matrixQuality/MNWD_Heatmap_selected.pdf",plot=heatmap, height=5, width=20)


## # same for a subset of selected motifs
## TableAll_selectedStringent<-TableAll[,as.vector(SelectedMotifsStringent$MotifId)]
## colnames(TableAll_selectedStringent)<-SelectedMotifsStringent$KnownMotif
## TableAll_selectedStringent$cluster <- rownames(TableAll_selectedStringent)
## TableAll_selectedStringent_melt<-melt(TableAll_selectedStringent, id.vars="cluster")
## TableAll_selectedStringent_melt$cluster<-factor(TableAll_selectedStringent_melt$cluster, level=(unique(TableAll_selectedStringent_melt$cluster)))

## heatmap<-ggplot(TableAll_selectedStringent_melt, aes(cluster,variable)) +
## 		geom_tile(aes(fill=value)) +
## 		scale_fill_gradientn(colours=colorRampPalette(brewer.pal(9, "YlGnBu"))(20)) +
## 		theme(text = element_text(size=30),axis.ticks=element_blank()) + labs(x="clusters",y="") + scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0))

## heatmap

##                                         #ggsave(file="mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/matrixQuality/MNWD_Heatmap_selected2.pdf",plot=heatmap, height=10, width=5)


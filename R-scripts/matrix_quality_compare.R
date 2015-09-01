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

mtx.quality.prefix.list.file <-"./results/matrix_quality/20150428/zoo_chip_enrichment_cluster_files.txt"

mtx.quality.prefix.list <-as.list(as.matrix( read.table(file=mtx.quality.prefix.list.file, header=FALSE, stringsAsFactors=FALSE)))
## directory with all the matrix quality results
##DirList <-  list.files("mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/matrixQuality", pattern="mFuzz*", full.names=TRUE)

# Import some lists of selcted motifs
##SelectedMotifs <- read.table("mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/PeakMotifs_SelectedMotifs", header=T, sep ="\t")
##SelectedMotifsStringent <- read.table("mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/PeakMotifs_SelectedMotifs_stringent.tsv", header=T, sep ="\t")

# import all data in one: for each matrix-quality run, import the ...compare-scores.tab file (row=pVal, col=sequence NWD), get the max of each column (MNWD for each sequence) and make a list of all that

################################################################
## TODO ALE
## 1) Read in full files (general matrix-quality run files) to keep using them for other processing options
## 2) Program geting the maxvalue for each column in each nwd file
## 3) Program getin only the maxvalue of significant pvalues
## 4) Program geting the area under the curve of the NWDs
## 5) Weigth the NWD area under the curve 
ListAll<- lapply(mtx.quality.prefix.list , FUN=function(dir){
    dir <- dirname(dir)
    table <- as.data.frame(t(as.matrix(as.data.frame(lapply(list.files(dir, pattern=".+_nwd_.+_compare-scores.tab", full.names=TRUE), FUN=function(file){
        print (paste("Reading in", file))
        MNWD <- apply(read.table(file, header=T, row.names=1, sep ="\t",na.strings="<NULL>",comment.char="", stringsAsFactors=FALSE), 2, FUN=function(x){
            max(x, na.rm=T)
        })
        ## This is a vector with the max NWD values
        print (MNWD)
        return(MNWD)
        
	})))))
    #print (table)
})
################
##

# remove part of the motif name to make it more readable
# # second line depens on the name of the file that was given for the output in matrix quality... ALE:Not sure what the rownames should be

ListAll<-lapply(ListAll,function(DF) {rownames(DF) <- paste("",seq(dim(ListAll[[1]])[1])); DF})

#ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8_","",colnames(DF)); DF})


ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("_score_distrib","",colnames(DF)); DF})
ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("_nwd.tab","",colnames(DF)); DF})
ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("_top_.+","",colnames(DF)); DF})


draw.heatmap <- function (ListAll,heatmap.file){
    ## reformat the data
    
    TableAll<-do.call("cbind",ListAll)
    TableAll$sequence <- rownames(TableAll)
    TableAll_melt<-melt(TableAll, id.vars="sequence")
    TableAll_melt$sequence<-factor(TableAll_melt$sequence, level=(unique(TableAll_melt$sequence)))
    
    ## put a higher threeshold on MNWD for scaling
    
    
    TableAll_melt$value[TableAll_melt$value>0.8] <- 0.8
    heatmap<-ggplot(TableAll_melt, aes(sequence,variable)) +
        geom_tile(aes(fill=value)) +
            scale_fill_gradientn(colours=colorRampPalette(brewer.pal(9, "Greys"))(20)) +
		theme(axis.ticks=element_blank()) + labs(x="sequence",y="") + scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0))
    
    heatmap
    ggsave(file=heatmap.file,plot=heatmap, height=5, width=20)
    
}

draw.heatmap(ListAll=ListAll,heatmap.file="maxNWD.pdf")

#################

ListAll.SigPval<- lapply(mtx.quality.prefix.list , FUN=function(dir){
    dir <- dirname(dir)
    table <- as.data.frame(t(as.matrix(as.data.frame(lapply(list.files(dir, pattern=".+_nwd_.+_compare-scores.tab", full.names=TRUE), FUN=function(file){
        print (paste("Reading in", file))
        nwd.table <- read.table(file, header=T, row.names=1, sep ="\t",na.strings="<NULL>",comment.char="", stringsAsFactors=FALSE)
        print (dim( nwd.table))
        nwd.table <-nwd.table[row.names(nwd.table)<=0.005,]
        #print (head(nwd.table))
        print (dim( nwd.table))
        #stop()
        
        MNWD <- apply( nwd.table, 2, FUN=function(x){            
            max(x, na.rm=T)
        })
        ## This is a vector with the max NWD values
        #print (MNWD)
        return(MNWD)
        
	})))))
    #print (table)
})

# remove part of the motif name to make it more readable
# # second line depens on the name of the file that was given for the output in matrix quality... ALE:Not sure what the rownames should be

ListAll.SigPval<-lapply(ListAll.SigPval,function(DF) {rownames(DF) <- paste("",seq(dim(ListAll[[1]])[1])); DF})

#ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8_","",colnames(DF)); DF})


ListAll.SigPval<-lapply(ListAll.SigPval,function(DF) {colnames(DF)<-gsub("_score_distrib","",colnames(DF)); DF})
ListAll.SigPval<-lapply(ListAll.SigPval,function(DF) {colnames(DF)<-gsub("_nwd.tab","",colnames(DF)); DF})
ListAll.SigPval<-lapply(ListAll.SigPval,function(DF) {colnames(DF)<-gsub("_top_.+","",colnames(DF)); DF})

draw.heatmap(ListAll=ListAll.SigPval,heatmap.file="maxNWD_on_significant_pvalues.pdf")

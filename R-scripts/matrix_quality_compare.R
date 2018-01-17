################################################################
## R-script to complemment matrix-quality analysis
## this script is called by compare-qualities

## Authors
## Samuel Collombet <samuel.collombet@ens.fr>
## Morgane Thomas-Chollier <mthomas@biologie.ens.fr>
## Alejandra Medina-Rivera <amedina@lcg.unam.mx>
## Jaime Castro-Mondragon <j.a.c.mondragon@ncmm.uio.no>

required.libraries <- c("RColorBrewer", "gplots" ,"flux")
for (lib in required.libraries) {
  if (!require(lib, character.only=TRUE)) {
    install.packages(lib)
    library(lib, character.only=TRUE)
  }
}

if(!"flux" %in% x[,1]){
  message("Missing package: flux")
}

args <- commandArgs(TRUE)

## Get the prefix list of files to be analyzed, this list is passed by compare-qualities
mtx.quality.nwds.file <- args[1]
print (mtx.quality.nwds.file)
plot.folder <- args[2]
formats <- args[3]
formats <- unlist(strsplit (formats,split=","))
print (formats)

print.heatmap <- args[4]

## For debugging:

#mtx.quality.nwds.file <-"./results/matrix_quality/20150428/zoo_chip_enrichment/all_nwd_files.txt"
#plot.folder <- "./results/matrix_quality/20150428/zoo_chip_enrichment/"
formats <-c("pdf","png")
#print (plot.folder)
#stop()

# mtx.quality.nwds.file <-"/Users/amedina/work_area/prueba/debug_matrix_quality/test/HOCOMOCO_motifs_CapStarrseq_Active_Prom_common_HeLa_K562_IP_vs_CapStarrseq_InactiveProm_FDR95_All_samples_bg_mkv_2_all_nwd_files.txt"
# plot.folder <- "/Users/amedina/work_area/prueba/debug_matrix_quality/test/HOCOMOCO_motifs_CapStarrseq_Active_Prom_common_HeLa_K562_IP_vs_CapStarrseq_InactiveProm_FDR95_All_samples_bg_mkv_2_all_nwd_plot"

dir.create(plot.folder, showWarnings = TRUE, recursive = TRUE)

################
## Read in table with nwd files
mtx.quality.nwds <- read.table(file=mtx.quality.nwds.file, header=FALSE, stringsAsFactors=FALSE)
colnames(mtx.quality.nwds) <- c("matrix","sequence","file")

################
## Read in each NWD file and store it in a list

matrices <- unique(mtx.quality.nwds$matrix)
nwd.files <- list()

for (i in 1:dim(mtx.quality.nwds)[1]){
    #i<-2
    matrix <- mtx.quality.nwds[i,1]
    seq <- mtx.quality.nwds[i,2]
    print (paste("Reading file for matrix and sequence",matrix,seq))
    one.table <- read.table(mtx.quality.nwds[i,3], header=TRUE, comment.char=";")
    colnames(one.table) <- c("Pvalue", "score_theor", "score_seq", "NWD")
    #print(head(one.table))
    #nwd.files [[paste("mtxID",matrix,"_seqID",seq, sep="")]][["table"]] <-  one.table
    nwd.files [[matrix]][[seq]][["table"]] <-  one.table
}

################
## Get max.nwd, max.sig.nwd, auc.all, auc.sig for each sequence matrix set

print ("Reading in files with nwd data")
lapply(names(nwd.files) ,function(matrix.name){
    print (matrix.name)
    print (names(nwd.files[[matrix.name]]))
    lapply( names(nwd.files[[matrix.name]]), function(seq=x, matrix.name2=matrix.name){
        print (matrix.name2)
        print (seq)
        local.table <- nwd.files[[matrix.name2]][[seq]][["table"]]
        #local.table$NWD <- scale(local.table$NWD)
        #local.table$NWD <- data.Normalization( local.table$NWD,type="n1",normalization="column")
        ## max.nwd
        max.nwd <- max(local.table$NWD)
        print(max.nwd )
        nwd.files[[matrix.name2]][[seq]][["max.nwd"]] <<- max.nwd 
        
        ## ##############
        
        ## get significant values
        local.sig.table <- local.table[local.table$score_theor>=0,]
        
        ## max.nwd.sig maximun taken from significant scores 0
        if ( dim( local.sig.table )[1]<=1){
            nwd.files[[matrix.name2]][[seq]][["max.sig.nwd"]] <<-"NA"
        }
        else{
            max.sig.nwd <- max(local.sig.table[,"NWD"])
            nwd.files[[matrix.name2]][[seq]][["max.sig.nwd"]] <<- max.sig.nwd
        }
        ## ##############
        ## AUC
        #if (min(local.table[,"NWD"])<0){
         #   local.table[,"NWD"] <- local.table[,"NWD"] + (-min(local.table[,"NWD"]))
        #}
        nwd.files[[matrix.name2]][[seq]][["auc.all"]] <<- auc(x=local.table[,"Pvalue"],y=local.table[,"NWD"])

        ## ##############
        ## AUC in significant values
        if ( dim( local.sig.table )[1]<=1){
            nwd.files[[matrix.name2]][[seq]][["auc.sig"]] <<-"NA"
        }
        else{
            nwd.files[[matrix.name2]][[seq]][["auc.sig"]] <<- auc(x=local.sig.table [,"Pvalue"],y=local.sig.table[,"NWD"])
        }
        ##stop()
        
    })
})

################
## Function draw.heatmap
## Takes as input the original list file with all data, the metric to be used in the heatmap
## the output heatmap file 
draw.heatmap <- function (ListAll,
                          metric="max.nwd",
                          heatmap.file, 
                          formats=c("pdf"), 
                          mtx.quality.nwds,
                          ...) {
    metric.table <- matrix(ncol = length(unique(mtx.quality.nwds$sequence)), nrow =length(unique(mtx.quality.nwds$matrix))) 
    colnames(metric.table) <- unique(mtx.quality.nwds$sequence)
    rownames(metric.table) <- unique(mtx.quality.nwds$matrix)
    
    for (mtx in unique(mtx.quality.nwds$matrix) ) {
        for (seq in unique(mtx.quality.nwds$sequence) ){
            metric.table[mtx,seq] <-  ListAll[[mtx]][[seq]][[metric]]
        }
    }

    write.table(file=paste(heatmap.file, "txt", sep="" ),metric.table,quote = FALSE, sep = "\t", na = "NA"  )
    if (print.heatmap==1){
        
        for (format in formats){
            if (format=="pdf"){
                pdf(file=paste(heatmap.file, format, sep="" ))
            } else if (format =="png"){
                png(file=paste(heatmap.file, format, sep="" ))
            } else {
                stop ("Format not available")
            }
                                        # metric.table <- scale(metric.table)
            
            heatmap.2(as.matrix(metric.table), col=colorRampPalette(brewer.pal(11,"RdBu"))(100),
                      trace="none", 
                      margins=c(6,10), 
                      cexRow = 0.75, 
                      cexCol = 0.75,
                                        #                                         , Colv=FALSE
                                        #                                         , breaks=breaks.hm
                                        #                                         , dendrogram= "none"
                      main = metric,
                      xlab = "Sequences",
                      ylab = "Motifs",
                      key = TRUE, 
                      key.title = "",
                      key.xlab = paste(metric, "value"), 
                      key.ylab = "",
                      scale="none",
                   density.info = "none", ...
                      )
            
            dev.off()
        }
    }
}

#####
## Draw heatmap using max.nwd

## If required draw heatmap

max.nwd.heatmap.file <- paste(plot.folder,"maxNWD_heatmap_compare.", sep="/")
draw.heatmap(ListAll=nwd.files,metric="max.nwd",heatmap.file=max.nwd.heatmap.file , formats=formats, mtx.quality.nwds=mtx.quality.nwds)

## ###
## Draw heatmap using max.nwd
max.sig.nwd.heatmap.file <- paste(plot.folder,"maxNWDsignificantScore_heatmap_compare.", sep="/")
draw.heatmap(ListAll=nwd.files,metric="max.sig.nwd",heatmap.file=max.sig.nwd.heatmap.file  , formats=formats, mtx.quality.nwds=mtx.quality.nwds)

## ###
## Draw heatmap using auc.all
auc.all.heatmap.file <- paste(plot.folder,"AUC_NWD_heatmap_compare.", sep="/")
draw.heatmap(ListAll=nwd.files,
             metric="auc.all",heatmap.file=auc.all.heatmap.file, 
             formats=formats, mtx.quality.nwds=mtx.quality.nwds,
             zlim=c(0,1))

## ###
## Draw heatmap using auc.all
auc.sig.heatmap.file <- paste(plot.folder,"AUC_NWDsignificantScore_heatmap_compare.", sep="/")
draw.heatmap(ListAll=nwd.files,metric="auc.sig",heatmap.file=auc.sig.heatmap.file  , formats=formats, mtx.quality.nwds=mtx.quality.nwds)


## ################################################################
## ## TODO ALE
## ## 1) Read in full files (general matrix-quality run files) to keep using them for other processing options
## ## 2) Program geting the maxvalue for each column in each nwd file
## ## 3) Program getin only the maxvalue of significant pvalues
## ## 4) Program geting the area under the curve of the NWDs
## ## 5) Weigth the NWD area under the curve

## ## List, one entry is one matrix, and inside a list with tables for the files per sequence
## ## Then we can try everything we want.



################################################################
### OLD CODE 


## lapply(matrices, function(matrix){
##     print (paste("Matrix", matrix))
##     sequences <- unique(mtx.quality.nwds[mtx.quality.nwds$matrix==matrix,"sequence"])
##     lapply(sequences, function(seq){
##         seq

##     })
    
    
## })
## mtx.quality.prefix.list.file <-"./results/matrix_quality/20150428/zoo_chip_enrichment/all_nwd_files.txt"



## mtx.quality.prefix.list <-as.list(as.matrix( read.table(file=mtx.quality.prefix.list.file, header=FALSE, stringsAsFactors=FALSE)))
## ## directory with all the matrix quality results
## ##DirList <-  list.files("mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/matrixQuality", pattern="mFuzz*", full.names=TRUE)

## # Import some lists of selcted motifs
## ##SelectedMotifs <- read.table("mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/PeakMotifs_SelectedMotifs", header=T, sep ="\t")
## ##SelectedMotifsStringent <- read.table("mFuzz/mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8/rsat/PeakMotifs_SelectedMotifs_stringent.tsv", header=T, sep ="\t")

## # import all data in one: for each matrix-quality run, import the ...compare-scores.tab file (row=pVal, col=sequence NWD), get the max of each column (MNWD for each sequence) and make a list of all that




## ListAll<- lapply(mtx.quality.prefix.list , FUN=function(dir){
##     dir <- dirname(dir)
##     table <- as.data.frame(t(as.matrix(as.data.frame(lapply(list.files(dir, pattern=".+_nwd_.+_compare-scores.tab", full.names=TRUE), FUN=function(file){
##         print (paste("Reading in", file))
##         MNWD <- apply(read.table(file, header=T, row.names=1, sep ="\t",na.strings="<NULL>",comment.char="", stringsAsFactors=FALSE), 2, FUN=function(x){
##             max(x, na.rm=T)
##         })
##         ## This is a vector with the max NWD values
##         print (MNWD)
##         return(MNWD)
        
## 	})))))
##     #print (table)
## })
## ################
## ##

## # remove part of the motif name to make it more readable
## # # second line depens on the name of the file that was given for the output in matrix quality... ALE:Not sure what the rownames should be

## ListAll<-lapply(ListAll,function(DF) {rownames(DF) <- paste("",seq(dim(ListAll[[1]])[1])); DF})

## #ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8_","",colnames(DF)); DF})


## ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("_score_distrib","",colnames(DF)); DF})
## ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("_nwd.tab","",colnames(DF)); DF})
## ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("_top_.+","",colnames(DF)); DF})


## draw.heatmap <- function (ListAll,heatmap.file){
##     ## reformat the data
    
##     TableAll<-do.call("cbind",ListAll)
##     TableAll$sequence <- rownames(TableAll)
##     TableAll_melt<-melt(TableAll, id.vars="sequence")
##     TableAll_melt$sequence<-factor(TableAll_melt$sequence, level=(unique(TableAll_melt$sequence)))
    
##     ## put a higher threeshold on MNWD for scaling
    
    
##     TableAll_melt$value[TableAll_melt$value>0.8] <- 0.8
##     heatmap<-ggplot(TableAll_melt, aes(sequence,variable)) +
##         geom_tile(aes(fill=value)) +
##             scale_fill_gradientn(colours=colorRampPalette(brewer.pal(9, "Blues"))(20)) +
## 		theme(axis.ticks=element_blank()) + labs(x="sequence",y="") + scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0))
    
##     heatmap
##     ggsave(file=heatmap.file,plot=heatmap, height=5, width=20)
    
## }

## draw.heatmap(ListAll=ListAll,heatmap.file="maxNWD.pdf")

## #################

## ListAll.SigPval<- lapply(mtx.quality.prefix.list , FUN=function(dir){
##     dir <- dirname(dir)
##     table <- as.data.frame(t(as.matrix(as.data.frame(lapply(list.files(dir, pattern=".+_nwd_.+_compare-scores.tab", full.names=TRUE), FUN=function(file){
##         print (paste("Reading in", file))
##         nwd.table <- read.table(file, header=T, row.names=1, sep ="\t",na.strings="<NULL>",comment.char="", stringsAsFactors=FALSE)
##         print (dim( nwd.table))
##         nwd.table <-nwd.table[row.names(nwd.table)<=0.005,]
##         #print (head(nwd.table))
##         print (dim( nwd.table))
##         #stop()
        
##         MNWD <- apply( nwd.table, 2, FUN=function(x){            
##             max(x, na.rm=T)
##         })
##         ## This is a vector with the max NWD values
##         #print (MNWD)
##         return(MNWD)
        
## 	})))))
##     #print (table)
## })

## # remove part of the motif name to make it more readable
## # # second line depens on the name of the file that was given for the output in matrix quality... ALE:Not sure what the rownames should be

## ListAll.SigPval<-lapply(ListAll.SigPval,function(DF) {rownames(DF) <- paste("",seq(dim(ListAll[[1]])[1])); DF})

## #ListAll<-lapply(ListAll,function(DF) {colnames(DF)<-gsub("mFuzz_counts_rlog_table_enh_Fd2_PulseFd05_ClustStdrz_k8_","",colnames(DF)); DF})


## ListAll.SigPval<-lapply(ListAll.SigPval,function(DF) {colnames(DF)<-gsub("_score_distrib","",colnames(DF)); DF})
## ListAll.SigPval<-lapply(ListAll.SigPval,function(DF) {colnames(DF)<-gsub("_nwd.tab","",colnames(DF)); DF})
## ListAll.SigPval<-lapply(ListAll.SigPval,function(DF) {colnames(DF)<-gsub("_top_.+","",colnames(DF)); DF})

## draw.heatmap(ListAll=ListAll.SigPval,heatmap.file="maxNWD_on_significant_pvalues.pdf")

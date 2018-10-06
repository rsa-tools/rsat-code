#!/usr/bin/env Rscript 
#
# Computes and plots the distribution of SNPs in upstream regions
# It is called by makefile/upstream-variation.mk
# Author: Chesco Montardit Tarda 2018 (edited by BContreras)

## Define the local directory for R libs
dir.rsat <- Sys.getenv("RSAT")
if (dir.rsat == "") {
  stop(paste("The environment variable RSAT is not defined. Command: ", commandArgs()))
}

dir.rsat.rscripts <- file.path(dir.rsat, "R-scripts")
dir.rsat.rlib <- file.path(dir.rsat.rscripts, "Rpackages")

## Load required libraries or warn user
required.packages = c("BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", # required by GenomicRanges
                      "GenomicRanges",
                      "zoo", # required by changepoint
                      "changepoint",
                      "crayon","pillar","withr","labeling", # required by ggplot2
                      "ggplot2",
                      "dplyr"
)

for (pkg in required.packages) {
  suppressPackageStartupMessages(library(pkg, warn.conflicts=TRUE,
    character.only = TRUE, lib.loc=c(dir.rsat.rlib, .libPaths())))
}



### Load data, delete genes without upstream sequence and separate D/R strands:
snp_vs_upstream <- function(){
  inputs <- commandArgs(trailingOnly = TRUE)
  snps <- read.table(inputs[1])
  colnames(snps) <- c("chr", "pos")
  snps$chr <- as.factor(snps$chr)
  snps$pos <- as.numeric(snps$pos)
  upstream <- read.table(inputs[2], sep = ":")
  colnames(upstream) <- c("chr", "low_pos", "up_pos", "strand")
  upstream$chr <- as.factor(upstream$chr)
  upstream$low_pos <- as.numeric(upstream$low_pos)
  upstream$up_pos <- as.numeric(upstream$up_pos)
  upstream$end <- upstream$low_pos - upstream$up_pos + 500
  upstream <- upstream[upstream$end < 0,]
  upstream_D <- upstream[upstream$strand=="D",]
  upstream_R <- upstream[upstream$strand=="R",]

  # Overlap polymorphisms in D/R upstream sequences, separately:
  position <- GRanges(snps$chr, IRanges(start=snps$pos, end=snps$pos))
  ranges_D <- GRanges(upstream_D$chr, IRanges(start=upstream_D$low_pos, end=upstream_D$up_pos))
  snp_overlaps_D <- findOverlaps(position,ranges_D, select = "last", type="within")
  snps$ID <- snp_overlaps_D
  upstream_D$ID <- seq(1,nrow(upstream_D),1)
  upstream_snps_D <- inner_join(snps, upstream_D, by="ID")
  upstream_snps_D$snp_pos <- -(upstream_snps_D$up_pos - upstream_snps_D$pos - 500)
  snp_counts_D <- as.data.frame(table(upstream_snps_D$snp_pos))
  colnames(snp_counts_D) <- c("seq", "freq")
  snp_counts_D$seq <- as.numeric(as.character(snp_counts_D$seq))
  rm(ranges_D)
  rm(snp_overlaps_D)
  ranges_R <- GRanges(upstream_R$chr, IRanges(start=upstream_R$low_pos, end=upstream_R$up_pos))
  snp_overlaps_R <- findOverlaps(position,ranges_R, select = "last", type="within")
  snps$ID <- snp_overlaps_R
  upstream_R$ID <- seq(1,nrow(upstream_R),1)
  upstream_snps_R <- inner_join(snps, upstream_R, by="ID")
  upstream_snps_R$snp_pos <- upstream_snps_R$low_pos - upstream_snps_R$pos + 500
  snp_counts_R <- as.data.frame(table(upstream_snps_R$snp_pos))
  colnames(snp_counts_R) <- c("seq", "freq")
  snp_counts_R$seq <- as.numeric(as.character(snp_counts_R$seq))
  rm(ranges_R, snp_overlaps_R)

  # Compute the number of genes per position:
  limits <- as.data.frame(table(upstream$end))
  colnames(limits) <- c("seq", "Freq")
  seq5000 <- seq(-5000,0,1)
  tmpseq <- data.frame(seq5000, rep(0, length(seq5000)))
  colnames(tmpseq) <- c("seq", "rep")
  tmp_limits <- merge(tmpseq, limits, all=TRUE)
  tmp_limits$Freq[is.na(tmp_limits$Freq)] <- 0
  genes <- cumsum(tmp_limits$Freq)
  genes_total <- c(genes,rep(nrow(upstream),500))
  gene_per_nt <- data.frame(seq(-5000,500,1), genes_total)
  colnames(gene_per_nt) <- c("seq", "genes")

  # Merge occurrences of polymorphisms per position from D/R strand and calculate normalized frequency (norm) and conservation
  snp_joined <- inner_join(snp_counts_R, snp_counts_D, by="seq")
  snp_joined$freq <- snp_joined$freq.x + snp_joined$freq.y
  snp_counts_tmp <- merge(tmpseq, snp_joined, by="seq", all=TRUE)
  snp_counts_tmp$freq[is.na(snp_counts_tmp$freq)] <- 0
  snp_counts <- inner_join(snp_counts_tmp, gene_per_nt, by="seq")
  rm(snp_counts_R, snp_counts_D, upstream_snps_D, upstream_snps_R, upstream_R, upstream_D, tmp_limits, tmpseq)
  snp_counts[2:4] <- list(NULL)
  colnames(snp_counts) <- c("seq", "freq", "genes")
  snp_counts$norm <- snp_counts$freq / snp_counts$genes
  snp_counts$seq <- as.numeric(as.character(snp_counts$seq))
  snp_counts$conservation <- (snp_counts$genes - snp_counts$freq) / snp_counts$genes
  
  # Detect points of variance change in distribution (library(changepoint))
  var_object <- cpt.var(snp_counts$conservation, method="PELT", penalty="AIC", minseglen = 100)
  var_points <- var_object@cpts - 5001
  utr_points <- head(var_points[var_points > 50], -1)
  prom_points <- var_points[var_points < - 50 & var_points > - 1000]
  recommended_point <- mean(prom_points)
  
  # Save data from occurrences and limits
  write.table(var_points, "tmp_points.tab", sep="\t", quote = FALSE, row.names = FALSE)
  output_list <- list(var_points, utr_points, prom_points, mean(prom_points))  
  names(output_list) <- c("Detected change points in variance", "5'UTR points", 
                            "Proximal promoter region points", "Recommended promoter delimitation point")
  capture.output(output_list, file="Changepoints.txt")
  colnames(snp_counts) <- c("Sequence position", "Count", "Number of genes", "Normalized frequency", "Normalized conservation")
  write.table(snp_counts, file="Result_dataframe_snps.tab", quote = FALSE, row.names = FALSE, sep="\t")
}

snp_vs_upstream()

### Plot saved data in two plots (total and zoomed)
snp_counts <- read.table("Result_dataframe_snps.tab", header=TRUE, sep="\t")
colnames(snp_counts) <- c("seq", "freq", "genes", "norm", "conservation")
var_points <- read.table("tmp_points.tab", sep="\t", header=FALSE, skip=1)
utr_points <- head(var_points[var_points > 50], -1)
prom_points <- var_points[var_points < - 50 & var_points > - 1000]
recommended_point <- mean(prom_points)

plot_filename1 = "snp_plot.png"

png(plot_filename1, width = 1000, height= 843)
ggplot(snp_counts, aes(seq, conservation)) + geom_point(alpha=0.8) + geom_line(alpha=0.3) + geom_smooth(method="loess", data=snp_counts, aes(seq, conservation)) + 
  xlab("Sequence position (nt)") + ylab("Normalized conservation frequency") + ggtitle("Normalized conservation by genetic variant") +
  scale_x_continuous(breaks=seq(-5000,500,250), limits=c(-2000,500)) + geom_vline(xintercept = 0, color="black") +
  geom_vline(xintercept= utr_points, color="red", linetype="dashed") + 
  geom_vline(xintercept = prom_points, color="black", linetype="dotted") + geom_vline(xintercept = recommended_point, color="red")
dev.off()

message("# output file (zoomed): ", plot_filename1)

plot_filename2 = "snp_plot_total.png"

png(plot_filename2, width = 1000, height= 843)
ggplot(snp_counts, aes(seq, conservation)) + geom_point(alpha=0.7) + geom_line(alpha=0.3) + geom_smooth(method="loess", data=snp_counts, aes(seq, conservation)) + 
  xlab("Sequence position (nt)") + ylab("Normalized conservation frequency") + ggtitle("Normalized conservation by genetic variant") +
  scale_x_continuous(breaks=seq(-5000,500,250)) + geom_vline(xintercept = 0, color="black") +
  geom_vline(xintercept= utr_points, color="red", linetype="dashed") + 
  geom_vline(xintercept = prom_points, color="black", linetype="dotted") + geom_vline(xintercept = recommended_point, color="red")
dev.off()

message("# output file (total): ", plot_filename2, "\n\n")

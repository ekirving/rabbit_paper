#!/usr/bin/env Rscript
suppressWarnings(library(ggplot2))
suppressWarnings(library(scales))

# get the command line arguments
args = commandArgs(trailingOnly=TRUE)
pcaname <- args[1]

v=c(15, 17, 15, 20, 18, 15, 17, 15, 16, 20, 18, 15)
t1 = read.table(paste(pcaname, "data", sep = "."))
col=c("#524FA1", "#FDB913", "lightgreen", "#00ADDC", "Darkgreen", "#ED1C24", "Black", "Pink", "Brown", "Cyan", "midnightblue", "palevioletred3", "lightcoral", "yellow4", "wheat4")
names(t1)[1]<-paste("Population")

alpha=c(1, 1, 1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1)
png(filename=paste(pcaname, "png", sep = "."), width = 1024, height = 768, pointsize=16)
pdf(file=paste(pcaname, "pdf", sep = "."), width = 10, height = 7)
ggplot(t1, aes(V3, V4)) + aes(shape=factor(Population)) + scale_shape_manual(values=v) + geom_point(aes(colour = factor(Population)), size=3, alpha=1) +  xlab("PC1") + ylab("PC2") + scale_colour_manual(values=col) + theme_bw() + theme(legend.title=element_blank(), legend.key = element_blank()) + guides(colour = guide_legend(override.aes = list(size=4)))
dev.off()
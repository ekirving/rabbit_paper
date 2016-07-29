#!/usr/bin/env Rscript
suppressWarnings(library(ggplot2))
suppressWarnings(library(scales))
library(stringr)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
pca_file=args[1]
pve_file=args[2]
pdf_file=args[3]
comp1=strtoi(args[4])
comp2=strtoi(args[5])
labeled=strtoi(args[6])

# TODO remove when done testing
#setwd("/Users/Evan/Dropbox/Code/rabbits")
#pca_file = "flashpca/pca_no-outgroup.data"
#pve_file = "flashpca/pve_no-outgroup.txt"
#comp1=1
#comp2=2
#labeled=strtoi("0")

# get the percentage of variance each component explains
pve <- round(read.table(pve_file)[,1]*100, 1)

v=c(15, 17, 15, 20, 18, 15, 17, 15, 16, 20, 18, 15)
t1 = read.table(pca_file)
col=c("#524FA1", "#FDB913", "lightgreen", "#00ADDC", "Darkgreen", "#ED1C24", "Black", "Pink", "Brown", "Cyan", "midnightblue", "palevioletred3", "lightcoral", "yellow4", "wheat4")
names(t1)[1]<-paste("Population")

alpha=c(1, 1, 1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1)
pdf(file=pdf_file, width = 10, height = 7)
gg <- ggplot(t1, aes(t1[[paste('V', comp1+2, sep='')]], t1[[paste('V', comp2+2, sep='')]])) + aes(shape=factor(Population)) + scale_shape_manual(values=v) + geom_point(aes(colour = factor(Population)), size=3, alpha=1) +  xlab(paste("PC", comp1, " (", pve[comp1], "%)", sep='')) + ylab(paste("PC", comp2, " (", pve[comp2], "%)", sep='')) + scale_colour_manual(values=col) + theme_bw() + theme(legend.title=element_blank(), legend.key = element_blank()) + guides(colour = guide_legend(override.aes = list(size=4)))
if (labeled) {
  # label all the points
  gg <- gg + geom_text(aes(label=str_sub(t1$V2, -3)),hjust=-.3, vjust=0)
}  
# display the plot
gg
dev.off()
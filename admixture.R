#!/usr/bin/env Rscript
suppressWarnings(library(ggplot2))
library(reshape2)

# get the command line arguments
args = commandArgs(trailingOnly=TRUE)
groupname <- args[1]

#groupname <- "admix/no-outgroup.pruned.3"

dat = read.table(paste(groupname, "data", sep = "."), header=TRUE)
dat.m <- melt(dat,id=c("Samples"))

# TODO how to assign the colour scheme?
# TODO sort the results so there is a smoother transition within the populations
# TODO show population names

pdf(file=paste(groupname, "pdf", sep = "."), width = 10, height = 7)
ggplot(dat.m, aes(x = Samples, y = value, fill=variable)) + geom_bar(stat='identity') + theme(legend.title=element_blank(), legend.key = element_blank(), axis.text.x  = element_text(angle=45, vjust=0.9, hjust=1), panel.background = element_blank()) +  xlab("Samples") + ylab("Ancestry")
dev.off()

#!/usr/bin/env Rscript
suppressWarnings(library(ggplot2))

# get the command line arguments
args = commandArgs(trailingOnly=TRUE)
groupname <- args[1]

#groupname <- "admix/no-outgroup.CV"

dat = read.table(paste(groupname, "data", sep = "."), header=TRUE)

# TODO highlight the 3 lowest points

pdf(file=paste(groupname, "pdf", sep = "."), width = 10, height = 7)
ggplot(dat, aes(x = K, y = CV, group=1)) + geom_point() + geom_line(stat='identity') + theme(legend.title=element_blank(), legend.key = element_blank(), panel.background = element_blank()) +  xlab("Ancestral populations (K)") + ylab("Cross-validation Error") + scale_x_continuous(breaks=c(0:10)) + scale_y_continuous(breaks=seq(0, 3, 0.1))
dev.off()
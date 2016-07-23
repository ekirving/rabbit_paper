#!/usr/bin/env Rscript

# get the command line arguments
args = commandArgs(trailingOnly=TRUE)
admixname <- args[1]

tbl = read.table(paste(admixname, "Q", sep = "."))

pdf(file=paste(admixname, "pdf", sep = "."), width = 10, height = 7)
barplot(t(as.matrix(tbl)), col=rainbow(3), xlab="Samples", ylab="Ancestry", border=NA)
dev.off()
#!/usr/bin/env Rscript
suppressWarnings(library("ape"))

args <- commandArgs(trailingOnly = TRUE)
data=args[1]
out=args[2]
m<-as.matrix(read.table(data, head=T, row.names=1))
tr=bionj(m)
write.tree(tr, file=out)
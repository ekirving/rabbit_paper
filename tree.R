#!/usr/bin/env Rscript
suppressWarnings(library("ape"))

args <- commandArgs(trailingOnly = TRUE)
data=args[1]
out=args[2]
setwd("/Users/Evan/Dropbox/Code/rabbits")

data <- "tree/all-pops.data"

m<-as.matrix(read.table(data, head=T, row.names=1))
tr=bionj(m)
tr<-root(tr, outgroup = "OUT-842", resolve.root = TRUE)

# sort the tree
tr <- ladderize(tr)

# colour the tips of the tree
col <- substring(tr$tip, 0, 3)
col[col=="DOM"] = "#524FA1"
col[col=="FRE"] = "#FDB913"
col[col=="IB1"] = "lightgreen"
col[col=="IB2"] = "#00ADDC"
col[col=="OUT"] = "Darkgreen"

# save the tree data
write.tree(tr, file=out)

# plot the tree
pdf(file=paste(out, "pdf", sep = "."), width = 10, height = 7)
plot(tr, type='phylogram', tip.color=col)
dev.off()
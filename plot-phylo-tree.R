#!/usr/bin/env Rscript
suppressWarnings(library("ape"))

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file=args[1]
tree_file=args[2]
pdf_file=args[3]

# TODO remove when done testing
#setwd("/Users/Evan/Dropbox/Code/rabbits")
#data_file <- "tree/all-pops.data"

m<-as.matrix(read.table(data_file, head=T, row.names=1))
tr=bionj(m)
tr<-root(tr, outgroup = "OUT-842", resolve.root = TRUE)

# sort the tree
tr <- ladderize(tr)

collist=c("#524FA1", "#FDB913", "lightgreen", "#00ADDC", "Darkgreen", "#ED1C24", "Black", "Pink", "Brown", "Cyan", "midnightblue", "palevioletred3", "lightcoral", "yellow4", "wheat4")

# colour the tips of the tree
col <- substring(tr$tip, 0, 3)
col[col=="DOM"] = collist[1]
col[col=="FRE"] = collist[2]
col[col=="IB1"] = collist[3]
col[col=="IB2"] = collist[4]
col[col=="OUT"] = collist[5]

# save the tree data
write.tree(tr, file=tree_file)

# plot the tree
pdf(file=pdf_file, width = 10, height = 7)
plot(tr, type='phylogram', tip.color=col)
dev.off()
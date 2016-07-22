setwd("~/Dropbox/Code/rabbits")
#dat <- read.table("./flashpca/pcs_all-pops.txt",h=F)
#dat <- read.table("./flashpca/pcs_no-outgroup.txt",h=F)
#dat <- read.table("./flashpca/pcs_no-out-ib1.txt",h=F)

#head(dat)
#plot(dat[,1],dat[,2])


library(ggplot2)
library(scales)
v=c(15, 17, 15, 20, 18, 15, 17, 15, 16, 20, 18, 15)
t1 = read.table("./flashpca/pcs_no-outgroup.labeled.txt")
#col=c("#FDB913", "#ED1C24", "#00ADDC", "#524FA1", "Gray", "Darkgreen", "Orange", "Black");
col=c("#524FA1", "#FDB913", "lightgreen", "#00ADDC", "Darkgreen", "#ED1C24", "Black", "Pink", "Brown", "Cyan", "midnightblue", "palevioletred3", "lightcoral", "yellow4", "wheat4")
names(t1)[1]<-paste("Population")
alpha=c(1, 1, 1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1)
#scale_alpha_manual(values=c(1, 1, 1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1))
#pdf("aIrish_full_pca_PC12.pdf")
ggplot(t1, aes(V2, V3)) + aes(shape=factor(Population)) + scale_shape_manual(values=v) + geom_point(aes(colour = factor(Population)), size=3, alpha=1) +  xlab("PC1") + ylab("PC2") + scale_colour_manual(values=col) + theme_bw() + theme(legend.title=element_blank(), legend.key = element_blank()) + guides(colour = guide_legend(override.aes = list(size=4)))
dev.off()

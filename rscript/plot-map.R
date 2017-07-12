require('maps')
require('mapdata')

setwd("/Users/Evan/Dropbox/Documents/Oxford/DPhil/Papers/TREE paper/Figs")

pts <- read.table("map-coordinates2.txt", header=T, sep="\t", comment.char='')
pts$color <- sapply(pts$color, as.character)

# setup the legend data
legend <- unique(pts[c('date','color')])
legend <- legend[with(legend, order(date)), ]

par(xpd=F)

# native range color
col1 <- '#fdbf6f'

# fetch the French regions
data(franceMapEnv)

# get the regions of France
regs <- map("france", namesonly=TRUE, plot=FALSE)

# draw the map
map('france', xlim=c(-11,19.5),ylim=c(36,60), regions=regs[c(74,79,80,86:88,91:94, 96:110,113)], fill=TRUE, col=col1, border=col1)
map("world", region=c('Portugal', "Spain(?!:.+)"), fill=TRUE, col=col1, border=col1, add=TRUE)
map("worldHires", xlim=c(-11,19.5),ylim=c(36,60), boundary=T, interior=T, add=TRUE)

# add the sample locations
points(pts$longitude, pts$latitude, pch=19, col=pts$color, cex=2)  #plot my sample sites

par(xpd=T)
# add the legend
legend('topright', inset=c(-0.3,-0.05), legend=legend$date, col=legend$color, cex=1, pch=19, title="Dates")

# -------------
# nams1 <- map("france", namesonly=TRUE, plot=FALSE)
# map("france")
# map('france', regions=nams1[c(74,79,80,86:88,91:94, 96:110,113)], fill=TRUE, col=col1, border=col1, add=T)
# map('france', regions=nams1[c(74)], fill=TRUE, col='orange', add=T,border='orange')
# map('france', regions=nams1[c(79)], fill=TRUE, col='orange', add=T,border='orange')


# ------------
# library(splancs)
# xy <- cbind(x=runif(100), y=runif(100))
# plot(xy)
# poly <- getpoly() # this gets the polygon on-screen
# plot(xy)
# polygon(poly)
# io <- inout(xy, poly)
# # this returns a logical vector for points in the polygon
# points(xy[io,], pch=16, col="blue")

# -----------

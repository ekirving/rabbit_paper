require('maps')
require('mapdata')
# install.packages(c('maps','mapdata', 'scales'))

setwd("/Users/Evan/Dropbox/Documents/Oxford/DPhil/Papers/TREE paper/Figs")

pts <- read.table("map-coordinates2.txt", header=T, sep="\t", comment.char='')
pts$color <- sapply(pts$color, as.character)
# cols <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f')

# setup the legend data
legend <- unique(pts[c('date','color')])
legend <- legend[with(legend, order(date)), ]

par(xpd=TRUE)

# native range color
col1 <- 'red'

# draw the map
map("worldHires", xlim=c(-11,19.5),ylim=c(38,55))

# color the native range
map("world", region=c('Portugal', "Spain(?!:.+)"), fill=TRUE, add=TRUE, col=col1)

points(pts$longitude, pts$latitude, pch=19, col=pts$color, cex=2)  #plot my sample sites
legend('topright', inset=c(-0.2,-0.05), legend=legend$date, col=legend$color, cex=1.3, pch=19, title="Century")


# ------------
library(splancs)
xy <- cbind(x=runif(100), y=runif(100))
plot(xy)
poly <- getpoly() # this gets the polygon on-screen
plot(xy)
polygon(poly)
io <- inout(xy, poly)
# this returns a logical vector for points in the polygon
points(xy[io,], pch=16, col="blue")

# -----------

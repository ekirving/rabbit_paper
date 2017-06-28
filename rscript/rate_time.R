library(ggplot2)
library(scales)

setwd("/Users/Evan/Dropbox/Documents/Oxford/DPhil/Papers/TREE paper/Figs")
pdf_file='rate-curve.pdf'

rates <- c(1.62e-9, 1.74e-9, 2.02e-9, 2.35e-9, 1.18e-6)
# times <- c(26324, 24539, 21137, 18169, 36)

# times <- c(17726,      16523,   14233,   12234,      24)
# cols  <- c('black',    'black', 'black', 'black',    'black')
# dates <- c('17.7 kya', '',      '',      '12.2 kya', '24 ya')

rates <- c(1.62e-9,    1.74e-9, 2.02e-9, 2.35e-9,    1.37e-8,    2.05e-8,    9.07e-8,  1.18e-6)
times <- c(17726,      16523,   14233,   12234,      2117,       1417,       317,       24)
col1 <- '#2166ac' # blue
col2 <- '#b2182b' # red
cols  <- c(col1,    col1, col1, col1,    col2,  col2,  col2, col1)
dates <- c('17.7 kya', '',      '',      '12.2 kya', '100 BC',   '600 AD',   '1700 AD', '24 ya')
# dates <- c('17.7 kya', '',      '',      '12.2 kya', '1st century BC',   '6st century AD',   '18th century AD', '24 ya')

# rates2 <- c(1.62e-9, 1.74e-9, 2.02e-9, 2.35e-9)
# times2 <- c(26324, 24539, 21137, 18169)

# dates <- c('26.3 kya', '', '', '18.2 kya', '36 ya')


dat = data.frame(x=times, y=rates)
# dat2 = data.frame(x=times2, y=rates2)

# the time estimate in units of mu
# t <- 4.269762246E-05
t <- 2.875081736E-05

breaks <- 10**(-9:-6)

fn1 <- function(x) {log10(t/x)}
fn2 <- function(x) {t/x}

pdf(file=pdf_file, width = 3.5, height = 2.6)

ggplot(data=dat, aes(x=x, y=y)) +
    theme_bw() +
    # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black"),
          text = element_text(size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    xlab("Time (years ago)") +
    ylab("Mutation rate (µ)") +
    geom_point(size=4, shape=19, colour=cols) +
    geom_text(aes(label=dates), nudge_y=log10(1.3), nudge_x=1300, colour=cols, cex=7) +
    scale_y_log10(breaks=breaks, labels=breaks) +
    scale_x_continuous(breaks=seq(0, 15000, by = 5000), limits=c(0, 20000), labels = comma, expand = c(0,500)) +
    stat_function(fun=fn1, geom="line")
    # stat_function(fun=fn1, geom="area", xlim = c(18169, 26324), fill = 'red', alpha = 0.2)
    # geom_area(stat = "function", fun = fn1, xlim = c(18169, 26324, fill = 'red', alpha = 0.2)

dev.off()

ggplot(data=dat2, aes(x=x, y=y)) +
    xlab("Time (t)") +
    ylab("Mutation rate (μ)") +
    geom_point(colour="red", size=4, shape=21, fill="white") +
    scale_x_continuous(limits=c(0, 30000), labels = comma) +
    stat_function(fun=fn2, geom="line") +
    stat_function(fun=fn2, geom="area", xlim = c(18169, 26324), fill = 'red', alpha = 0.2)


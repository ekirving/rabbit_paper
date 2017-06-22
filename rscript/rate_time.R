library(ggplot2)
library(scales)

rates <- c(1.62e-9, 1.74e-9, 2.02e-9, 2.35e-9, 1.18e-6)
times <- c(26324, 24539, 21137, 18169, 36)

rates2 <- c(1.62e-9, 1.74e-9, 2.02e-9, 2.35e-9)
times2 <- c(26324, 24539, 21137, 18169)

dates <- c('26.3 kya', '', '', '18.2 kya', '36 ya')

dat = data.frame(x=times, y=rates)
dat2 = data.frame(x=times2, y=rates2)

# the time estimate in units of mu
t <- 4.269762246E-05

breaks <- 10**(-9:-6)

fn1 <- function(x) {log10(t/x)}
fn2 <- function(x) {t/x}

ggplot(data=dat, aes(x=x, y=y)) +
    theme_bw() +
    xlab("Time (years)") +
    ylab("Mutation rate (μ)") +
    geom_point(colour="red", size=4, shape=21, fill="white") +
    geom_text(aes(label=dates), nudge_y=log10(1.2), nudge_x=1600) +
    scale_y_log10(breaks=breaks, labels=breaks) +
    scale_x_continuous(breaks=seq(0, 30000, by = 5000), limits=c(0, 30000), labels = comma) +
    stat_function(fun=fn1, geom="line")
    # stat_function(fun=fn1, geom="area", xlim = c(18169, 26324), fill = 'red', alpha = 0.2)
    # geom_area(stat = "function", fun = fn1, xlim = c(18169, 26324, fill = 'red', alpha = 0.2)

ggplot(data=dat2, aes(x=x, y=y)) +
    xlab("Time (t)") +
    ylab("Mutation rate (μ)") +
    geom_point(colour="red", size=4, shape=21, fill="white") +
    scale_x_continuous(limits=c(0, 30000), labels = comma) +
    stat_function(fun=fn2, geom="line") +
    stat_function(fun=fn2, geom="area", xlim = c(18169, 26324), fill = 'red', alpha = 0.2)


library(ggplot2)
library(scales)

rates <- c(1.62e-9, 1.74e-9, 2.02e-9, 2.35e-9, 1.18e-6)
times <- c(26324, 24539, 21137, 18169, 36)
dates <- c('26 Ka', '25 Ka', '21 Ka', '18 Ka', '36 yr')

# the time estimate in units of mu
t <- 4.269762246E-05

breaks <- 10**(-9:-6)

fn1 <- function(x) {log10(t/x)}
fn2 <- function(x) {t/x}

ggplot(data=dat, aes(x=x, y=y)) +
    theme_bw() +
    xlab("Time (t)") +
    ylab("Mutation rate (μ)") +
    geom_point(colour="red", size=4, shape=21, fill="white") +
    # geom_text(aes(label=dates), nudge_y=log10(1.2)) +
    scale_y_log10(breaks=breaks, labels=breaks) +
    scale_x_continuous(limits=c(0, 30000), labels = comma) +
    stat_function(fun=fn1, geom="line")
    # stat_function(fun=fn1, geom="area", xlim = c(2000, 3000), fill = 'red', alpha = 0.2) +
    # geom_area(stat = "function", fun = fn1, xlim = c(2000, 3000), fill = 'red', alpha = 0.2)

ggplot(data=dat, aes(x=x, y=y)) +
    xlab("Time (t)") +
    ylab("Mutation rate (μ)") +
    geom_point(colour="red", size=4, shape=21, fill="white") +
    scale_x_continuous(limits=c(0, 30000), labels = comma) +
    stat_function(fun=fn2, geom="line") +
    stat_function(fun=fn2, geom="area", xlim = c(2000, 3000), fill = 'red', alpha = 0.2)


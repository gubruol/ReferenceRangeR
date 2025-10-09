scale.fact <- 1
x.lo.limit <- NA
x.hi.limit <- NA
use.oh <- NA
use.dev <- NA
ana.hour.min <- NA
ana.minu.min <- NA
ana.hour.max <- NA
ana.minu.max <- NA
print.log.message <- TRUE
smooth.hist1 <- TRUE
smooth.hist2 <- FALSE
lambda.min <- 0
lambda.max <- 1
fastnull     <- 1.e-10
fastnull.chi <- 0.1
x.tr.prop.min <- 0.60
x.tr.prop.max <- 0.95
x.tr.prop.limits <- seq(x.tr.prop.min, x.tr.prop.max, length.out=5)
x.tr.prop.ints.n <- length(x.tr.prop.limits) - 1
l.fact <- 0.0
p.fact <- 0.0 
s.fact <- 2/3
w.fact <- 1.0 
path.tab.stra <- NA
path.fig.stra <- NA
figtype <- "bmp"
gtab.list <- c("tmc") 
gtab.list.n <- length(gtab.list)
meth.list <- gtab.list
meth.list.n <- length(meth.list)
age.limits <- c(NA, NA)
RL1.p <- 0.025
RL2.p <- 0.975
x.clip.type <- "as.specified"
x.clip.by1 <- NA
x.clip.by2 <- NA
par.las <- 1
par.tcl <- 0.5
histcol    <- "lightgray"
bordercol <- "black"
tmccol      <- "darkgreen"
difbordercol <- "red"
difhistcol <- adjustcolor("red", alpha.f = 0.10)
kdecol <- "black"
nbs <-  0
user <- "user"
info.file <- NA
r.start <- NA
r.ende <- NA
nrep <- r.ende - r.start + 1
eval.rep <- !is.na(nrep)
idx.fix <- NULL
idx.est <- c(1, 2, 3)
RunId <- NA
n.per.bin.min  <- 10
df.est <- 3
df.con <- 1
x.tr.bins.min <- df.est + df.con + 2
kernel.n <- 4096
bins.n.max <- 100
bins.n.min <- df.est + df.con + 4
red.limit <- 1000
time.start <- Sys.time()
lambda.seq <- c(0.0, 0.1, 1.0)
prev.acc.lo <- -0.05
prev.acc.hi <-  1.05
alpha <- 0.05
r.fact <- 100
p.fit.min <- 0.20
crit123456 <- TRUE
crit6.p <- 0.05
outname <- NA
outname.stra <- NA
sexlabel <- "sex"
agelabel <- "age"
x.unaffected.lo <- NA
x.unaffected.hi <- NA
plot.fig100.001 <- FALSE
plot.fig100.002 <- FALSE
plot.fig100.003 <- FALSE
plot.fig100.004 <- FALSE
plot.fig100.005 <- FALSE
plot.fig100.006 <- FALSE
plot.fig100.007 <- FALSE
plot.fig100.008 <- FALSE
plot.fig100.009 <- FALSE
plot.fig100.010 <- FALSE
plot.fig100.011 <- FALSE
plot.fig100.012 <- FALSE
plot.fig100.013 <- FALSE
plot.fig100.014 <- FALSE
plot.fig100.015 <- FALSE
plot.fig100.016 <- TRUE
plot.fig100.017 <- FALSE
plot.fig108.010 <- FALSE
irep <- 0
iage <- 0
igtab <- 0
subset.type <- 0
round.unit <- 0.01
xlabel <- NULL
detect.limits.max <- NA
spalte.g <- NA
lambda.gen <- NA
x.n1 <- NA
x.n14 <- NA
iblock <- 0
datafile <- NA
infile <- NA
tab.names <- c(NA,NA)
options(encoding = "latin1")



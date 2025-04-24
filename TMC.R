##################################################################
# TMC - Reference limit estimator                              ###
# Author: W. Wosniok                                           ###
# https://www.degruyter.com/document/doi/10.1515/cclm-2018-1341/html
# 08.05.2024 Slight modifications by G. Brandhorst             ###
##################################################################


### CHANGES FROM ORIGINAL FUNCTIONS ###
#588: # win.graph()
#1453: x.hist, xlabel, NULL, NULL, NULL,
#1498: #savePlot(file=fig100.016.file, type=figtype)
#1499: myplot <- recordPlot()
#2181: if !eval.rep) { sink(file=outname.stra.txt,split=TRUE) }
#6842: main=NULL,
#6843: #
#6844: # 
#6845: sub=NULL, cex.sub=0.7)
#7031: # savePlot(file=figA.file,type=figA.type)


scale.fact <- 1
x.lo.limit <- NA
x.hi.limit <- NA
use.oh <- NA
use.dev <- NA
ana.hour.min <- NA
ana.minu.min <- NA
ana.hour.max <- NA
ana.minu.max <- NA
print.log.message <- FALSE
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
round.unit <- 1
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


tmc <- function (x) {

q10 <- quantile(x, probs = 0.1)
q90 <- quantile(x, probs = 0.9)
x.clip.min <- q10 - (q90 - q10) / 1.3
x.clip.max <- q90 + (q90 - q10) / 1.3

#  TMC_seg100_Analysis.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Segment doing the analysis per stratum

#  #  CHANGE HISTORY

#  02.05.2021 Result plots per method moved to method-specific segments
#  25.04.2021 Tukey approach added
#  03.03.2021 Basic information per stratum always in output table
#  22.02.2021 Calculation of x.tr.prop.limits moved to 040_ProcessSettings 
#  15.02.2021 Subsegments outsourced (100_Analysis_QQW, ...)
#  26.01.2021 meth.list made effective
#  13.01.2021 r.fact introduced
#  11.01.2021 Search sequence for TI changed to large -> small with 
#             intermediate checks for acceptability
#  04.01.2021 Smoothing strategy changed, applied to uncollapsed equidistant 
#             data 
#  17.12.2020 Graph windows opened by seg046
#  07.12.2002 Safety actions if QQW produces no result
#  01.12.2020 Organisation of histogram calculation, smoothing and reduction
#             changed, now done by FindHisto()
#  16.11.2020 x.kde added to call QQW()
#  14.11.2020 bins.n.min added to call QQW()
#  09.11.2020 New approach of analysis started. TMC is now a function, 
#             qqw and Bhattacharya added. New approach for kde bandwidth 
#             and histogram bandwidth.

#  21.08.2020 Save plot100.040 only for regular analysis, not simulation 
#             Unaffected range taken from DataInfoFile, not calculated here
#  18.08.2020 x.clip.by1 == NA or x.clip.by2 == NA correctly processed 
#  21.07.2020 Pr.lo, Pr.hi set to 0.001, 0.999 (previously 0.01, 0.95)
#             Setting moved to TMC_seg015_DefaultSettings.R
#  20.07.2020 Detail output for startum in txt file re-organized
#  16.07.2020 Common par parameters for plotting
#  01.07.2020 Check generated data unified
#  30.05.2020 Logarithmic x axis in fig100.010
#  25.05.2020 RB3[2, ] expanded: rule for y[1]==0 added
#     05.2020 Plotting of density with new approach for density bw added
#  11.05.2020 RB3[2, ] geändert
#  29.04.2020 Plots newly arranged. kde plot still excluded.
#  21.04.2020 x.kde calculated from only an interior range of the data
#  29.03.2020 opt.crit als einheitliche Bezeichnung
#  11.03.2020 Strategie vorübergehend zurückgesetzt, neuen Namen beibehalten
#  10.03.2020 Strategie für Bandbreitenberechnung geändert 
#  19.02.2020 tmc.master4 replaced by tmc.master5 (hierarchical search of
#             truncation interval)
#  24.01.2020 New structure for stopping execution if not enough data
#  23.01.2020 Break values for histogram calculated by find.breaks,
#             construction keeps n.per.bin.min and x.bins.min, if possible
#  18.01.2020 QQ plot based calculation of initial values
#  05.01.2020 Plot numbers changed
#  04.01.2020 Notation unified: always a . after x, y in x.hist, x.breaks, ...
#             No x or y with x.counts: only 'counts'
#  03.01.2020 New organisation: only tmc is done, new strategy for finding
#             the optimal truncation interval, initial value selection 
#             simplified, results are stored in 'tab' only  
#  05.11.2019 Plots for bootstrapping updated
#  27.10.2019 More information for BS result in output table
#  23.10.2019 TMC estimation applied to smoothed histogram
#  12.10.2019 BS calculation of CI activated
#  10.10.2019 TMC part of the analysis moved to TMC_seg101_Analysis_TMC.R
#  28.08.2019 prop.tmc.seq calculation revised to chi2 contribution
#             p.fit.min introduced, see Expert settings
#  27.08.2019 prop.tmc.seq calculation revised
#  14.08.2019 prop.tmc.seq shifted to interval above detection limit
#  26.07.2019 In RL.est: xc.mode, xc.sig, xc.P50 added, also for
#             methods other than IFCC 
#             In RL.est: x.mode,  y.mode,  y.sig replaced by
#                        xc.mode, yc.mode, yc.sig 
#  24.07.2019 Calculation prop.qqw changed
#             New file naming convention, see seg050
#  22.07.2019 Typing error in optimality criterion removed
#  13.07.2019 Sequence for prop.tmc candidates defined in 
#             TMC_seg060_ExpertSettings.R
#  02.07.2019 names of TMC result variables corrected
#  01.07.2019 RL1.p, RL2.p  instead of fixed probabilities
#  08.06.2019 Empirische Quantile grundsätzlich berechnet
#  12.05.2019 errcode für tmc noch einmal geändert
#             prop.qqw and prop.tmc are dynamcally calculated, depending
#             on the number of values < detection limit
#  08.05.2019 errcode für tmc geändert
#  07.05.2019 qqw operates on x with ties replaced by increasing sequence
#  05.04.2019 y.mode hier berechnet
#  12.03.2019 Start from WW17_seg40  (directory ProgWW17_RH)

# ===========================================================================

if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Start\n") }

# ===========================================================================

# ......................................................................
#  Clean up. From here on processing of a new data subset.

if (exists("x.kde"))  { rm(x.kde) }
if (exists("y.kde"))  { rm(y.kde) }
if (exists("x.hist")) { rm(x.hist) }
if (exists("y.hist")) { rm(y.hist) }

# ......................................................................
#  Initialise error messages

ErrTxt   <- c("Not enough values in stratum",
              "Not enough histogram bins of requested size in stratum",
              "Requested number of bins reduced")

ErrTxt.n <- length(ErrTxt)

ErrMsg <- data.frame(ErrNum=1:ErrTxt.n,
                     Occurred=rep(FALSE, times=ErrTxt.n),
                     ErrTxt,
                     stringsAsFactors=FALSE)

#  Set back execution control. go.on == FALSE means no further execution
#  in this stratum
go.on <- TRUE

# ......................................................................
#  Define stratum specific output file names

#  TMC_seg110_NamesPerStratum.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Segment defining output file names that are stratum-specific.
#  See head of seg050 for the naming convention
#  File-specific names are defined in segment 080
#  Scenario-specific names are defined in segment 211

#  TODO       
#  -  

#  #  CHANGE HISTORY
#  20.07.2020 Sex code and age class added to outfile.stra 
#  24.07.2019 File names changed according to new convention
#             outfile.base -> outname.stra
#             outfile.txt  <- outname.stra.txt
#  12.03.2019 Start from WWW17_seg40  (directory ProgWW17_RH)

# ===========================================================================

# ......................................................................
outname.stra <- paste(outname,"_",sexlabel,"_",agelabel,sep="") 

#  Text file with results per stratum (= separated by sex and age), 
#  no smoothing
outname.stra.txt  <- paste(path.tab.stra,outname.stra, ".txt",sep="")
if (file.exists(outname.stra.txt)) { file.remove(outname.stra.txt) }

#  Histogram data from Fig02, equidistant bins 
outhist.rdata <- paste(path.tab.stra,outname.stra,".RData",sep="")
if (file.exists(outhist.rdata)) { file.remove(outhist.rdata) }

#  Histogram data for estimated pathological data, as in fig 16.1 
outhistpath.rdata <- paste(path.tab.stra,outname.stra,"_path.RData",sep="")
if (file.exists(outhistpath.rdata)) { file.remove(outhistpath.rdata) }

# ......................................................................

# ......................................................................
#  Graph windows

#  TMC_seg099_OpenWindows.R
#  (c) wwosniok@math.uni-bremen.de

#  Open graph windows and define plot file names for plots requested in
#  TMC_seg045_PlotRequest.R

#  17.12.2020 Start
# ============================================================================

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg099_OpenWindows.R  Start\n") }
# ============================================================================

# ============================================================================
# Open.100.001.R
# ==============================================================================
  if (plot.fig100.001)
  { #  Plot is wanted
    if (!exists("fig100.001"))
    { #  Open graph new window
      win.graph()
      fig100.001 <- dev.cur()
    }
    fig100.001.file <- paste(path.fig.stra,outname.stra,"-F100.001.",
                             figtype,sep="")
    if (file.exists(fig100.001.file)) { file.remove(fig100.001.file) }
  } else
  { #  Plot is not wanted
    fig100.001 <- NA
  }

# ============================================================================
# Open.100.002.R
# ==============================================================================
  if (plot.fig100.002)
  { #  Plot is wanted
    if (!exists("fig100.002"))
    { #  Open graph new window
      win.graph()
      fig100.002 <- dev.cur()
    }
    fig100.002.file <- paste(path.fig.stra,outname.stra,"-F100.002.",
                             figtype,sep="")
    if (file.exists(fig100.002.file)) { file.remove(fig100.002.file) }
  } else
  { #  Plot is not wanted
    fig100.002 <- NA
  }

# ============================================================================
# Open.100.003.R
# ==============================================================================
  if (plot.fig100.003)
  { #  Plot is wanted
    if (!exists("fig100.003"))
    { #  Open graph new window
      win.graph()
      fig100.003 <- dev.cur()
    }
    fig100.003.file <- paste(path.fig.stra,outname.stra,"-F100.003.",
                             figtype,sep="")
    if (file.exists(fig100.003.file)) { file.remove(fig100.003.file) }
  } else
  { #  Plot is not wanted
    fig100.003 <- NA
  }

# ============================================================================
# Open.100.004.R
# ==============================================================================
  if (plot.fig100.004)
  { #  Plot is wanted
    if (!exists("fig100.004"))
    { #  Open graph new window
      win.graph()
      fig100.004 <- dev.cur()
    }
    fig100.004.file <- paste(path.fig.stra,outname.stra,"-F100.004.",
                             figtype,sep="")
    if (file.exists(fig100.004.file)) { file.remove(fig100.004.file) }
  } else
  { #  Plot is not wanted
    fig100.004 <- NA
  }


# ============================================================================
# Open.100.005.R
# ============================================================================
  if (plot.fig100.005)
  { #  Plot is wanted
    if (!exists("fig100.005"))
    { #  Open graph new window
      win.graph()
      fig100.005 <- dev.cur()
    }
    fig100.005.file <- paste(path.fig.stra,outname.stra,"-F100.005.",
                             figtype,sep="")
    if (file.exists(fig100.005.file)) { file.remove(fig100.005.file) }
  } else
  { #  Plot is not wanted
    fig100.005 <- NA
  }

# ============================================================================
# Open.100.006.R
# ============================================================================
  if (plot.fig100.006)
  { #  Plot is wanted
    if (!exists("fig100.006"))
    { #  Open graph new window
      win.graph()
      fig100.006 <- dev.cur()
    }
    fig100.006.file <- paste(path.fig.stra,outname.stra,"-F100.006.",
                             figtype,sep="")
    if (file.exists(fig100.006.file)) { file.remove(fig100.006.file) }
  } else
  { #  Plot is not wanted
    fig100.006 <- NA
  }


# ============================================================================
# Open.100.007.R
# ============================================================================
  if (plot.fig100.007)
  { #  Plot is wanted
    if (!exists("fig100.007"))
    { #  Open graph new window
      win.graph()
      fig100.007 <- dev.cur()
    }
    fig100.007.file <- paste(path.fig.stra,outname.stra,"-F100.007.",
                             figtype,sep="")
    if (file.exists(fig100.007.file)) { file.remove(fig100.007.file) }
  } else
  { #  Plot is not wanted
    fig100.007 <- NA
  }


# ============================================================================
# Open.100.008.R
# ============================================================================
  if (plot.fig100.008)
  { #  Plot is wanted
    if (!exists("fig100.008"))
    { #  Open graph new window
      win.graph()
      fig100.008 <- dev.cur()
    }
    fig100.008.file <- paste(path.fig.stra,outname.stra,"-F100.008.",
                             figtype,sep="")
    if (file.exists(fig100.008.file)) { file.remove(fig100.008.file) }
  } else
  { #  Plot is not wanted
    fig100.008 <- NA
  }


# ============================================================================
# Open.100.009.R
# ============================================================================
  if (plot.fig100.009)
  { #  Plot is wanted
    if (!exists("fig100.009"))
    { #  Open graph new window
      win.graph()
      fig100.009 <- dev.cur()
    }
    fig100.009.file <- paste(path.fig.stra,outname.stra,"-F100.009.",
                             figtype,sep="")
    if (file.exists(fig100.009.file)) { file.remove(fig100.009.file) }
  } else
  { #  Plot is not wanted
    fig100.009 <- NA
  }


# ============================================================================
# Open.100.010.R
# ============================================================================
  if (plot.fig100.010)
  { #  Plot is wanted
    if (!exists("fig100.010"))
    { #  Open graph new window
      win.graph()
      fig100.010 <- dev.cur()
    }
    fig100.010.file <- paste(path.fig.stra,outname.stra,"-F100.010.",
                             figtype,sep="")
    if (file.exists(fig100.010.file)) { file.remove(fig100.010.file) }
  } else
  { #  Plot is not wanted
    fig100.010 <- NA
  }


# ============================================================================
# Open.100.011.R
# ============================================================================
  if (plot.fig100.011)
  { #  Plot is wanted
    if (!exists("fig100.011"))
    { #  Open graph new window
      win.graph()
      fig100.011 <- dev.cur()
    }
    fig100.011.file <- paste(path.fig.stra,outname.stra,"-F100.011.",
                             figtype,sep="")
    if (file.exists(fig100.011.file)) { file.remove(fig100.011.file) }
  } else
  { #  Plot is not wanted
    fig100.011 <- NA
  }


# ============================================================================
# Open.100.012.R
# ============================================================================
  if (plot.fig100.012)
  { #  Plot is wanted
    if (!exists("fig100.012"))
    { #  Open graph new window
      win.graph()
      fig100.012 <- dev.cur()
    }
    fig100.012.file <- paste(path.fig.stra,outname.stra,"-F100.012.",
                             figtype,sep="")
    if (file.exists(fig100.012.file)) { file.remove(fig100.012.file) }
  } else
  { #  Plot is not wanted
    fig100.012 <- NA
  }


# ============================================================================
# Open.100.013.R
# ============================================================================
  if (plot.fig100.013)
  { #  Plot is wanted
    if (!exists("fig100.013"))
    { #  Open graph new window
      win.graph()
      fig100.013 <- dev.cur()
    }
    fig100.013.file <- paste(path.fig.stra,outname.stra,"-F100.013.",
                             figtype,sep="")
    if (file.exists(fig100.013.file)) { file.remove(fig100.013.file) }
  } else
  { #  Plot is not wanted
    fig100.013 <- NA
  }


# ============================================================================
# Open.100.014.R
# ============================================================================
  if (plot.fig100.014)
  { #  Plot is wanted
    if (!exists("fig100.014"))
    { #  Open graph new window
      win.graph()
      fig100.014 <- dev.cur()
    }
    fig100.014.file <- paste(path.fig.stra,outname.stra,"-F100.014.",
                             figtype,sep="")
    if (file.exists(fig100.014.file)) { file.remove(fig100.014.file) }
  } else
  { #  Plot is not wanted
    fig100.014 <- NA
  }


# ============================================================================
# Open.100.015.R
# ============================================================================
  if (plot.fig100.015)
  { #  Plot is wanted
    if (!exists("fig100.015"))
    { #  Open graph new window
      win.graph()
      fig100.015 <- dev.cur()
    }
    fig100.015.file <- paste(path.fig.stra,outname.stra,"-F100.015.",
                             figtype,sep="")
    if (file.exists(fig100.015.file)) { file.remove(fig100.015.file) }
  } else
  { #  Plot is not wanted
    fig100.015 <- NA
  }


# ============================================================================
# Open.100.016.R
# ============================================================================
  if (plot.fig100.016)
  { #  Plot is wanted
    if (!exists("fig100.016"))
    { #  Open graph new window
#      win.graph()
      fig100.016 <- dev.cur()
    }
    fig100.016.file <- paste(path.fig.stra,outname.stra,"-F100.016.",
                             figtype,sep="")
    if (file.exists(fig100.016.file)) { file.remove(fig100.016.file) }
  } else
  { #  Plot is not wanted
    fig100.016 <- NA
  }


# ============================================================================
# Open.100.017.R
# ============================================================================
  if (plot.fig100.017)
  { #  Plot is wanted
    if (!exists("fig100.017"))
    { #  Open graph new window
      win.graph()
      fig100.017 <- dev.cur()
    }
    fig100.017.file <- paste(path.fig.stra,outname.stra,"-F100.017.",
                             figtype,sep="")
    if (file.exists(fig100.017.file)) { file.remove(fig100.017.file) }
  } else
  { #  Plot is not wanted
    fig100.017 <- NA
  }

# ============================================================================
# Open.108.010.R
# ============================================================================
  if (plot.fig108.010)
  { #  Plot is wanted
    if (!exists("fig108.010"))
    { #  Open graph new window
      win.graph()
      fig108.010 <- dev.cur()
    }
    fig108.010.file <- paste(path.fig.stra,outname.stra,"-F108.010.",
                             figtype,sep="")
    if (file.exists(fig108.010.file)) { file.remove(fig108.010.file) }
  } else
  { #  Plot is not wanted
    fig108.010 <- NA
  }

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg099_OpenWindows.R  End\n") }
# ============================================================================


# ......................................................................
#  Create results tables
#  All results from this segment, which refer to one stratum,
#  go into 'tab.stra'. Column names are given in tab.names, defined
#  in seg050.
#  Information on data generation was read by read.table and is kept in 
#  the relevant variables - see below

#  Create output table for the actual stratum
tab.stra           <- matrix(NA,nrow=meth.list.n, ncol=length(tab.names))
dimnames(tab.stra) <- list(meth.list, tab.names)
tab.stra           <- data.frame(tab.stra)

# ......................................................................
#  Separator lines for output tables
sep.p <- str_dup(".", 78)
sep.m <- str_dup("-", 78)
sep.e <- str_dup("=", 78)

# ============================================================================
#  Decribe the actual stratum

if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 050\n") }

x.n <- length(x)
subtitle <- paste(RunId," / ",outname.stra," / n =",x.n)

#  Wieviele Werte unter NWG?
if (!is.na(detect.limits.max))
{
  x.lt.detect.n <- sum(x <= detect.limits.max)
} else
{
  x.lt.detect.n <- 0
}
x.lt.detect.pct <- 100*x.lt.detect.n/x.n

#  @@@@@ check if  x.lt.detect.pct < allowed maximum

if (exists("age")) 
{ Age.mea <- mean(age) } else
{ Age.mea <- NA }

# ----------------------------------------------------------------------------
#  Put basic information in tab.stra, also for subsets that will not be 
#  analysed below due to too few cases

if ("bha" %in% meth.list)
{
  tab.stra["bha", "irep"]    <- irep
  tab.stra["bha", "method"]  <- "bha"

  tab.stra["bha", "sexlabel"]    <- sexlabel
  tab.stra["bha", "agelabel"]    <- agelabel
  tab.stra["bha", "Age.num"]     <- iage
  tab.stra["bha", "Age.mea"]     <- Age.mea

  tab.stra["bha","subset.type"] <- subset.type
  tab.stra["bha", "n"]          <- x.n
  tab.stra["bha", "pct.lt.DL"]  <- x.lt.detect.pct
}

if ("npa" %in% meth.list)
{
  tab.stra["npa", "irep"]    <- irep
  tab.stra["npa", "method"]  <- "npa"

  tab.stra["npa", "sexlabel"]    <- sexlabel
  tab.stra["npa", "agelabel"]    <- agelabel
  tab.stra["npa", "Age.num"]     <- iage
  tab.stra["npa", "Age.mea"]     <- Age.mea

  tab.stra["npa","subset.type"] <- subset.type
  tab.stra["npa", "n"]          <- x.n
  tab.stra["npa", "pct.lt.DL"]  <- x.lt.detect.pct
}

if ("qqw" %in% meth.list)
{
  tab.stra["qqw", "irep"]    <- irep
  tab.stra["qqw", "method"]  <- "qqw"

  tab.stra["qqw", "sexlabel"]    <- sexlabel
  tab.stra["qqw", "agelabel"]    <- agelabel
  tab.stra["qqw", "Age.num"]     <- iage
  tab.stra["qqw", "Age.mea"]     <- Age.mea

  tab.stra["qqw","subset.type"] <- subset.type
  tab.stra["qqw", "n"]          <- x.n
  tab.stra["qqw", "pct.lt.DL"]  <- x.lt.detect.pct
}

if ("tmc" %in% meth.list)
{
  tab.stra["tmc", "irep"]    <- irep
  tab.stra["tmc", "method"]  <- "tmc"

  tab.stra["tmc", "sexlabel"]    <- sexlabel
  tab.stra["tmc", "agelabel"]    <- agelabel
  tab.stra["tmc", "Age.num"]     <- iage
  tab.stra["tmc", "Age.mea"]     <- Age.mea

  tab.stra["tmc","subset.type"] <- subset.type
  tab.stra["tmc", "n"]          <- x.n
  tab.stra["tmc", "pct.lt.DL"]  <- x.lt.detect.pct
}

if ("tmu" %in% meth.list)
{
  tab.stra["tmu", "irep"]    <- irep
  tab.stra["tmu", "method"]  <- "tmu"

  tab.stra["tmu", "sexlabel"]    <- sexlabel
  tab.stra["tmu", "agelabel"]    <- agelabel
  tab.stra["tmu", "Age.num"]     <- iage
  tab.stra["tmu", "Age.mea"]     <- Age.mea

  tab.stra["tmu","subset.type"] <- subset.type
  tab.stra["tmu", "n"]          <- x.n
  tab.stra["tmu", "pct.lt.DL"]  <- x.lt.detect.pct
}

if ("tuk" %in% meth.list)
{
  tab.stra["tuk", "irep"]    <- irep
  tab.stra["tuk", "method"]  <- "tuk"

  tab.stra["tuk", "sexlabel"]    <- sexlabel
  tab.stra["tuk", "agelabel"]    <- agelabel
  tab.stra["tuk", "Age.num"]     <- iage
  tab.stra["tuk", "Age.mea"]     <- Age.mea

  tab.stra["tuk","subset.type"] <- subset.type
  tab.stra["tuk", "n"]          <- x.n
  tab.stra["tuk", "pct.lt.DL"]  <- x.lt.detect.pct
}

# -------------------------------------------------------------------------
#  Analysis only if data subset has sufficient cases (absolute minimum)
#  Data distribution may additionally cause too few bins of required size - 
#  is checked below. Check here is not finally, values below detection limit
#  may play an additional role.

#  Required minimum number of observations in a stratum
x.n.min <- n.per.bin.min*x.tr.bins.min

if (x.n < x.n.min)
{ #  Not enough observations - no analysis
  cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
      "\nData subset for    ", sexlabel,"  ",agelabel,
      "\n  has only         ", x.n," values",
      "\n  Required minimum:", x.n.min,
      "\n                   ==> no analysis",
      "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
 
  #  Set switch and error message 
  ErrMsg[1,"Occurred"] <- TRUE
  go.on <- FALSE

}  else
{ 
  # ===========================================================================
  #  Required minimum number of observations in a stratum exists: do analysis

  if (print.log.message) { cat("%%%   TMC_seg100_Analysis  100\n") }

  #  Save test data  
  #if (outname.stra == "Krebs_GGT_2_M_40-49")
  #if (outname.stra == "ALT_Bochum_7_stationaer_M_18-29")
  #{ 
  #  save(x, file=paste("../temp/", outname.stra, ".RData", sep="")) 
  #}
 
  #  Calculate descriptive quantities
  x.min <- min(x)
  x.max <- max(x)

  #  Nonparametric quantiles of the total dataset. Subsequent RL estimates 
  #  must lie in [x.Q1, x.Q2]

  x.Q1 <- Quantile(x,probs=RL1.p)
  x.Q2 <- Quantile(x,probs=RL2.p)

  #  Wie differenziert sind die Daten (wieviel verschiedene Werte) ?
  x.table <- table(x)

  x.val   <- as.numeric(names(x.table))
  x.val.n <- length(x.table)

  # .........................................................................
  #  Empirical distribution function

  x.emp.cdf <- c(0, cumsum(x.table)/x.n)
  x.emp.x   <- c(0, x.val)

  # ==========================================================================
  if (print.log.message) { cat("%%%   TMC_seg100_Analysis  200\n") }

  #  Calculate kde
  KDE <- FindKDE(x, round.unit, kernel.n, 
                 fig100.001, subtitle, kdecol,
                 Q.lo=0.025, Q.hi=0.975, b.fact=1.05)
  x.kde      <- KDE[["x.kde"]]
  x.kde.mode <- KDE[["x.kde.mode"]]

  xsupp <- x.kde$x

  x.lt.mode.n <- sum(x <  x.kde.mode)
  x.ge.mode.n <- sum(x >= x.kde.mode)

  # ==========================================================================
  if (print.log.message) { cat("%%%   TMC_seg100_Analysis  300\n") }

  #  Produce a reasonable histogram on the original scale, together with a 
  #  smoothed version

  HIST <- FindHisto(x, round.unit, detect.limits.max, x.table, x.val, x.kde, 
                    n.per.bin.min, bins.n.min, bins.n.max,
                    smooth.hist1, smooth.hist2, red.limit, xlabel, subtitle, 
                    bordercol1, histcol1, bordercol2, histcol2, 
                    bordercol3, histcol3, kdecol, fig100.002)

  x.hist.eqd     <- HIST[["x.hist.eqd"]]      # equidistant, raw
  x.hist.eqd.smo <- HIST[["x.hist.eqd.smo"]]  # equidistant, smoothed
  x.hist         <- HIST[["x.hist"]]          # collapsed, possibly smoothed
  x.reso.red     <- HIST[["x.reso.red"]]      # reduced data set

  #  Histogram construction may have left too few bins (==> x.hist == NA).
  x.hist.exists <- (length(x.hist) > 1)

  if (x.hist.exists)
  {
    bins.n.act <-  length(x.hist$counts)

    #  Which bin in x.hist contains x.kde.mode?
    x.kde.mode.idx <- FindModeIndex(x.hist$breaks, x.hist$counts, x.kde.mode)

    #  counts.range refers to x.hist
    counts.range <- paste(formatC(min(x.hist$counts), format="f", digits=0),
                        "-", 
                        formatC(max(x.hist$counts), format="f", digits=0), 
                        sep="")

    # =========================================================================
    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  400\n") }

    cat("\n\n", sep.p,
    "\nStart of analysis                      ", format(time.start),
    "\nData file                              ",datafile, 
    "\nAnalyte                                ",xlabel, 
    "\nAnalysis of the stratum defined by",
    "\n  Age                                  ",agelabel,
    "\n  Sex                                  ",sexlabel,    
    "\nRecords for analysis in this stratum   ",x.n,
    "\nValues lie in the range                ",x.min, x.max,
    "\nNumber of distinct values              ",x.val.n,
    "\nSmoothing parameter for x              ",x.kde$bw,
    "\nmode(x)                                ",x.kde.mode,
    "\n", sep.p, 
    "\n\n")

    # ---------------------------------------------------------------------------
    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  500\n") }

    #  Is there indication of problematic rounding (eg "go to the even value"
    #  or inconsistent rounding)?  Do diagnostic plot.

    if (plot.fig100.004)
    {
      dev.set(fig100.004)

      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(as.numeric(names(x.table)), as.numeric(x.table),type="h",
         log="x",
         ylim=c(0,max(x.table)),
         xlab=xlabel, ylab="Count",
         main="Frequencies of observed values",
         sub=subtitle,cex.sub=0.7)
      # if (!eval.rep) savePlot(file=fig100.004.file,type=figtype)
    }

  }  else    #  x.hist.exists
  {
    #  x.hist does not exist
    bins.n.act <- 0
  } 

  # =========================================================================
  #  Start RL calculation

  if (print.log.message) { cat("%%%   TMC_seg100_Analysis  600\n") }

  #  Initialize output tables 

#  TMC_seg101_Ini_Analysis.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Initialize output tables for TMC_seg100

#  CHANGE HISTORY

#  15.02.2021 Start 
# ===========================================================================


tab.npa <- data.frame(method="npa",  
                      x.tr.lo.npa=NA, x.tr.hi.npa=NA,
                      x.tr.n.npa=NA, x.tr.n.npa=NA,
                      xl.n.npa=NA, xc.n.npa=NA, xr.n.npa=NA,
                      prev.l.npa=NA,prev.c.npa=NA, prev.r.npa=NA,
                      x.RL1.npa=NA, x.RL2.npa=NA)   

tab.qqw.names <- c(
 "lambda.qqw" ,     "bin.start"       ,"bin.end"    ,     "bins.n"    ,     
 "prop.qqw"   ,     "x.tr.lo"         ,"x.tr.hi"    ,     "y.tr.lo"   ,     
 "y.tr.hi"    ,     "mue.qqw"         ,"sigma.qqw"  ,     "r2"        ,     
 "opt.crit"   ,     "x.tr.n"          ,"rank.r2"    ,     "y.RL1.qqw" ,     
 "y.RL2.qqw"  ,     "x.RL1.qqw"       ,"x.RL2.qqw"  ,     "n.l.qqw"   ,     
 "prev.l.qqw" ,     "n.c.qqw"         ,"prev.c.qqw" ,     "n.r.qqw"   ,     
 "prev.r.qqw" ,     "neg.prev.sum.qqw","x.RL1.cilo" ,     "x.RL1.cihi",     
 "x.RL2.cilo" ,     "x.RL2.cihi")

tab.qqw <- matrix(NA, nrow=1, ncol=length(tab.qqw.names))
colnames(tab.qqw) <- tab.qqw.names

tab.tmc.names <- c(
 "ilo"           ,"ihi"         ,  "x.tr.lo"    ,   "x.tr.hi"       ,
 "x.tr.n"        ,"x.tr.prop"   ,  "rc"         ,   "iter"          ,
 "lambda.tmc"    ,"mue.tmc"     ,  "sigma.tmc"  ,   "x.RL1.tmc"     ,
 "x.RL2.tmc"     ,"x.RL1.cilo"  ,  "x.RL1.cihi" ,   "x.RL2.cilo"    ,
 "x.RL2.cihi"    ,"reldist.tmc" ,  "p.fit"      ,   "opt.crit"      ,
 "prev.l.tmc"    ,"prev.c.tmc"  ,  "prev.r.tmc" ,   "prev.l.tmc.pen",
 "prev.r.tmc.pen","xc.n.tmc"    ,  "chi2.total" ,   "chi2.total.df" ,
 "chi2.trun"     ,"chi2.trun.df",  "chi2.path"  ,   "sol.score" )

tab.tmc <- matrix(NA, nrow=1, ncol=length(tab.tmc.names))
colnames(tab.tmc) <- tab.tmc.names

crit1 <- NA
crit2 <- NA
crit3 <- NA
crit4 <- NA
crit5 <- NA
crit6 <- NA

# colnames(tab.tmu)
tab.tmu.names <- c(
 "mue.tmu"    , "sigma.tmu"   ,"distance.tmu","rc.tmu"     , "iter.tmu"  , 
 "x.RL1.tmu"  , "x.RL2.tmu"   ,"x.RL1.cilo"  ,"x.RL1.cihi" , "x.RL2.cilo", 
 "x.RL2.cihi" , "opt.crit.tmu","p.fit.tmu"   ,"prev.l.tmu" , "prev.c.tmu", 
 "prev.r.tmu") 

tab.tmu <- matrix(NA, nrow=1, ncol=length(tab.tmu.names))
colnames(tab.tmu) <- tab.tmu.names



  # =========================================================================
  #  If group membership is known (typically for simulated data):
  #  Calculate the ideal direct result  
  #  Ideal direct result: all nonpathological persons are correctly identified
  #  No pathological cases lie in the truncation interval (though no 
  #  truncation interval is needed for RLx.npa)
  # 
  
  if (!is.na(spalte.g))
  {
    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  700\n") }

# TMC_seg100_Analysis_npa.R
#

#  15.02.2021
# ==============================================================================

    # xl.n.npa  <- sum(grp=="L") # also under perfect identification unknown
    xc.n.npa  <- sum(grp=="C")   # this number would be known under perfect
                                 # identification of non-diseased subject
    # xr.n.npa  <- sum(grp=="R") # also under perfect identification unknown
    prev.c.npa <- xc.n.npa/x.n

    #  Calculation of prev.l.npa and prev.r.npa here is different from
    #  the other methods, because prev.c.n is exactly known
    xl.n.npa <- sum(x < x.kde.mode)  - sum((x < x.kde.mode) & (grp=="C")) 
    prev.l.npa <- xl.n.npa/x.n
    xr.n.npa <- sum(x >= x.kde.mode) - sum((x >= x.kde.mode) & (grp=="C")) 
    prev.r.npa <- xr.n.npa/x.n

    #  Optimal truncation interval (left limits are contained, right not)
    if (xl.n.npa == 0)
    { x.tr.lo.npa <- min(x[grp=="C"]) - round.unit/2} else
    { x.tr.lo.npa <- max(x[grp=="L"]) + round.unit/2}
    if (xr.n.npa == 0)
    { x.tr.hi.npa <- max(x[grp=="C"]) + round.unit/2} else
    { x.tr.hi.npa <- min(x[grp=="R"]) - round.unit/2}
    x.tr.n.npa    <- length(x[(x.tr.lo.npa <= x) & (x < x.tr.hi.npa)])
    x.tr.prop.npa <- x.tr.n.npa / x.n  

    x.RL1.npa <- Quantile(x[grp=="C"], probs=RL1.p)
    x.RL2.npa <- Quantile(x[grp=="C"], probs=RL2.p)

    tab.npa <- data.frame(method="npa",
                          x.tr.lo.npa=x.tr.lo.npa, x.tr.hi.npa=x.tr.hi.npa,
                          x.tr.n.npa=x.tr.n.npa, x.tr.prop.npa=x.tr.prop.npa, 
                          xl.n.npa=xl.n.npa, xc.n.npa=xc.n.npa,
                          xr.n.npa=xr.n.npa,
                          prev.l.npa=prev.l.npa, prev.c.npa=prev.c.npa,
                          prev.r.npa=prev.r.npa,
                          x.RL1.npa=x.RL1.npa, x.RL2.npa=x.RL2.npa)   
    row.names(tab.npa) <- NULL


  }  else
  {
    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  800\n") }

    #  Output table tab.npa initialized in seg_1010 
  }

  # ==========================================================================
  if (print.log.message) { cat("%%%   TMC_seg100_Analysis  900\n") }

  if (x.hist.exists)
  { 
    #  There is a reasonable number of bins
    #  Find TI candidates with required properties (proportion of values)
    #  bins.n.min in next call replaced by x.tr.bins.min 23.02.2021

    TI.can <- proportions(x.hist$counts, x.kde.mode.idx, 
                          x.tr.bins.min, x.tr.prop.min, x.tr.prop.max)

    #  Go through the candidate list in groups of x.tr.prop
    NoQqwSolution <- TRUE   
    NoTmcSolution <- TRUE
    first.cycle   <- TRUE   
    iti <- x.tr.prop.ints.n + 1     #  set in seg040

    while (NoTmcSolution & iti >= 2 )
    {
      iti <- iti - 1     # index for prop range

      x.tr.prop.range <- paste(formatC(100*x.tr.prop.limits[iti], width=2,
                                       digits=0, format="f"), "-",
                               formatC(100*x.tr.prop.limits[iti+1], width=3,
                                       digits=0, format="f"), "%", sep="")

      TI.subset <- (x.tr.prop.limits[iti] <= TI.can[ ,"prop"]) &
                   (TI.can[ ,"prop"] < x.tr.prop.limits[iti+1] )
      TI.subset.n <- sum(TI.subset, na.rm=TRUE)
          
      if (print.log.message) { cat("%%%   TMC_seg100_Analysis  910\n") }

      if (TI.subset.n > 0)
      {   
        #  Non-empty subset, use the prop.limits as preliminary 
        #  limits for QQW, TMC, TMU estimation

        #  Calculate QQW solution
        #  Use reduced and resolved dataset if the dataset is large

        if (print.log.message) 
        { cat("%%%   TMC_seg100_Analysis  920", date(), "\n") }

#  TMC_seg103_Analysis_qqw.R
#
#  15.02.2021
# ==============================================================================

if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw  100\n") }

tab.qqw <- QQW(x, x.reso.red, round.unit, lambda.seq, 
                   x.hist, x.kde, x.kde.mode, x.kde.mode.idx, x.val, 
                   x.unaffected.lo, x.unaffected.hi,  
                   n.per.bin.min, bins.n.min, x.tr.bins.min, 
                   counts.range, x.tr.prop.range,
                   x.tr.prop.limits[iti], x.tr.prop.limits[iti+1], 
                   prev.acc.lo, prev.acc.hi, 
                   xsupp,
                   x.Q1, x.Q2, RL1.p, RL2.p,
                   lambda.gen, yc.mode.gen, yc.sig.gen,
                   df.est,df.con,
                   l.fact, p.fact, r.fact, w.fact, 
                   xlabel, subtitle,
                   fig100.006, fig100.007,
                   fig100.008, fig100.009,
                   fig100.010, 
                   fig100.011, fig100.011.file,
                   fig100.012, fig100.012.file, figtype,
                   x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                   x.clip.type, par.las, par.tcl,
                   bordercol2, histcol2, kdecol, qqrawcol, polycol, gencol, 
                   difcol, denlty,
                   print.log.message)

# --------------------------------------------------------------------
if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 200\n") }

#  QQW may have produced no result. Then go to next range.
       
if (length(tab.qqw) == 1)
{ #  No result
  if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 300\n") }

  #cat("\n +++ No QQW result for proportion range ",
  #            x.tr.prop.range,
  #    "\n +++", subtitle,
  #    "\n\n")
}  else
{
  # QQW solution exists, add more details 

  if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 400\n") }

  NoQqwSolution <- FALSE

  temp.qqw <- chi2trunc(x.hist, tab.qqw[1, "x.tr.lo"],
                        tab.qqw[1, "x.tr.hi"],
                        sum(x < tab.qqw[1, "x.tr.lo"]), 
                        sum(x >= tab.qqw[1, "x.tr.hi"]),
                        x.Q1, x.Q2, RL1.p, RL2.p,
                        tab.qqw[1, "lambda.qqw"],
                        tab.qqw[1, "mue.qqw"],
                        tab.qqw[1, "sigma.qqw"],
                        df.est, df.con,
                        l.fact, p.fact, r.fact, w.fact, 
                        opt.crit.only=FALSE,fastnull=fastnull)

  # -------------------------------------------------------------------------
  #  Plot the QQW result
  if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 420\n") }

  if (!is.na(fig100.012))
  { dev.set(fig100.012)
    #bringToTop(fig100.012)

    #  Plot result - basic part x.kde set to NA
    fig100.012.lim <- PlotHistFit(fig100.012, fig100.012.file, figtype, 
                  x.hist, xlabel, "QQW result", subtitle, n.per.bin.min,
                  NA, x.val, counts.range, x.tr.prop.range,
                  x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                  x.clip.type, NA, par.las, par.tcl,
                  bordercol2, histcol2, difcol, kdecol, denlty)

    if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 440\n") }

    #  Plot result - qqw specific part 
    #  Add estimated QQW solution, RLs, TI
    #  Calculate estimated QQW density (not yet scaled by prev.c.qqw)
    xc.pdf.qqw <- pdf.PN(xsupp, unlist(tab.qqw[1, "lambda.qqw"]), 
                                unlist(tab.qqw[1, "mue.qqw"]), 
                                unlist(tab.qqw[1, "sigma.qqw"])) 

    #  Calculate difference between x.kde and xc.pdf.qqw
    #  Truncate difference below to yplotmin
    #  Switched off 01.05.2021
    #x.kde.qqw <- x.kde$y -  tab.qqw["prev.c.qqw"] * xc.pdf.qqw
    #x.kde.qqw[x.kde.qqw < fig100.012.lim["yplotmin"]] <- 
    #                                              fig100.012.lim["yplotmin"]

    if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 460\n") }

    #  No plot file name, save is done farther below
    PlotMetRes(fig100.012, NA, figtype, 
               fig100.012.lim["yplotmin"], 
               xsupp, tab.qqw[1, "prev.c.qqw"], xc.pdf.qqw, qqwcol,
               x.kde.mode, kdecol, NA, difcol, difbordercol, difhistcol,
               unlist(tab.qqw[1, "x.RL1.qqw"]), unlist(tab.qqw[1, "x.RL2.qqw"]), 
               0.85*fig100.012.lim["yplotmax"], 
               unlist(tab.qqw[1, "x.tr.lo"]), unlist(tab.qqw[1, "x.tr.hi"]), 
               0.90*fig100.012.lim["yplotmin"], qqwcol,
               temp.qqw)

    # --------------------------------------------------------------------------
    #  If the data is generated : show unaffected range
    if (!is.na(x.unaffected.lo) & !is.na(x.unaffected.hi))
    { 
      lines(c(x.unaffected.lo, x.unaffected.hi), 
      0.60*fig100.012.lim["yplotmin"]*c(1,1), col="chartreuse", lwd=2)
    }

    # -------------------------------------------------------------
    # Save qqw result figure

    savePlot(file=fig100.012.file, type=figtype) 
  }  
} 

if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 500\n") }



        # ==================================================================
        if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1100\n") }

        #  Calculate Bhattacharya solution, if requested. Needs equidistant 
        #  histogram.
          
        if ("bha" %in% meth.list)
        {  
          data.text <- xlabel
          tab.bha   <- Bhatta(x.hist.eqd.smo, x.kde, data.text, xlabel, subtitle,
                            RL1.p, RL2.p, 
                            fig100.013, fig100.014,
                            kdecol, bhacol, bordercol1, histcol1)
        }

        # ==================================================================
        if (print.log.message) 
        { cat("%%%   TMC_seg100_Analysis 1200", date(), "\n") }

        # If a QQW solution exists: TMC

        QqwSolution <- length(tab.qqw) > 1

        if ( QqwSolution )         
        {
          #  QQW solution exists.
          #  Calculate TMC solution. x.hist may be smoothed,
          #  depending on smooth.hist1 /2, but further processing is the same
          #  in all cases. tab.qqw contains truncation intervals and initial
          #  parameter estimates.
    
          #  Do TMC analysis for one TI subset group, find optimal among
          #  all actual size groups
          
#  TMC_seg105_Analysis_tmc.R
#
#  07.04.2021 crit1 not used if values < DL exist 
#  15.02.2021 Start
# ==============================================================================
  tab.tmc0 <- tmc.master(x, x.hist, tab.qqw, 
                         theta.fix, theta.est, idx.fix, idx.est,
                         lambda.min, lambda.max, 
                         df.est, df.con, p.fit.min,
                         l.fact, p.fact, r.fact, w.fact, 
                         x.Q1, x.Q2, RL1.p, RL2.p, fastnull, fastnull.chi,                     
                         print.log.message,
                         unlist(tab.npa["x.RL1.npa"]), 
                         unlist(tab.npa["x.RL2.npa"]), 
                         xsupp, gencol, kdecol, tmccol, fig100.015)

  #  Assessment of TMC result
  #  Ideally, the solution fulfills all 7 criteria (fit + crit 1-6).
  #  If there is no such solution, the final solution is found by 
  #  as follows:
  #  solution 1 is better than solution 2, if it has more criteria
  #  fulfilled,
  #  if two solutions have the same number of fulfulled criteria,
  #  the one with higher p.fit is better

  #  Calculate assessment criteria
  #  Estimated RLs in [x.Q1, x.Q2]?
  #  Consider crits 1 only if no values < DL, because x.Q1 is
  #  very unreliable if values < DL exist
  if (is.na(detect.limits.max))
  { crit1 <- unname(x.Q1 <= tab.tmc0["x.RL1.tmc"]) } else
  { crit1 <- TRUE } 
    
  crit2 <- unname(x.Q2 >= tab.tmc0["x.RL2.tmc"])

  #  Estimated prevalences in [prev.acc.lo, prev.acc.hi]?
  crit3 <- unname((prev.acc.lo <= tab.tmc0["prev.l.tmc"])) &
                   unname((tab.tmc0["prev.l.tmc"] <= prev.acc.hi))  
  crit4 <- unname((prev.acc.lo <= tab.tmc0["prev.c.tmc"])) &
                   unname((tab.tmc0["prev.c.tmc"] <= prev.acc.hi))  
  crit5 <- unname((prev.acc.lo <= tab.tmc0["prev.r.tmc"])) &
                   unname((tab.tmc0["prev.r.tmc"] <= prev.acc.hi))
  crit6 <- NA
                    
  #  Basic solution score
  solution.score <- crit1 + crit2 + crit3 + crit4 + crit5

  NoTmcSolution <- (unname(tab.tmc0["p.fit"]) < p.fit.min) |
                           (solution.score < 5) 

  #  More sophisticated solution score, if desired
  
  if (crit123456)
  { #  Solution score is calculated including crit6 (random 
    #  fluctuation of residuals)
    temp.tmc0 <- chi2trunc(x.hist,tab.tmc0["x.tr.lo"],
                             tab.tmc0["x.tr.hi"],
                             sum(x <  tab.tmc0["x.tr.lo"]), 
                             sum(x >= tab.tmc0["x.tr.hi"]),
                             x.Q1, x.Q2, RL1.p, RL2.p,
                             tab.tmc0["lambda.tmc"],
                             tab.tmc0["mue.tmc"],
                             tab.tmc0["sigma.tmc"],
                             df.est, df.con,
                             l.fact, p.fact, r.fact, w.fact, 
                             opt.crit.only=FALSE,fastnull=fastnull)
          
    sign <- 1.0 * (temp.tmc0$tab["diff"] <= 0)          
    runstest <- runs.test(sign, exact=FALSE, alternative="two.sided")

    #  p level for significance reduced because observations are not
    #  independent 
    crit6 <- unname((runstest["p.value"] > crit6.p))
  
    solution.score <- solution.score + crit6  
    NoTmcSolution <- (unname(tab.tmc0["p.fit"]) < p.fit.min) |
                             (solution.score < 6)
  } 

  tab.tmc0 <- c(tab.tmc0, sol.score=solution.score)
          
  if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1210\n") }

  #  If this is the first loop cycle, store actual 
  #  result  (to get a result at all)
  if (first.cycle)
  { 

    if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1220\n") }
  
    tab.tmc     <- tab.tmc0
    temp.tmc    <- temp.tmc0  
  }

  #  If this is not the first loop cycle (which means that previous
  #  solutions were not satisfactory), compare actual 
  #  result with previous and keep the better
  if (!first.cycle)
  { 
    if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1230\n") }

    #cat("\n [Analysis_tmc] Comparison of tab.tmc0 and tab.tmc\n") 

    #temp <- matrix(c(tab.tmc0["x.RL1.tmc"], tab.tmc0["x.RL2.tmc"],
    #tab.tmc0["p.fit"], tab.tmc0["opt.crit"], 
    #tab.tmc0["x.tr.prop"], tab.tmc0["sol.score"],
    #tab.tmc0["RL1.pen"], tab.tmc0["RL2.pen"], 
    #tab.tmc["x.RL1.tmc"], tab.tmc["x.RL2.tmc"],
    #tab.tmc["p.fit"], tab.tmc["opt.crit"], 
    #tab.tmc["x.tr.prop"], tab.tmc["sol.score"],
    #tab.tmc["RL1.pen"], tab.tmc["RL2.pen"]), 
    #byrow=TRUE, nrow=2)
    #dimnames(temp) <- list(c("tmc0", "tmc"), c("RL1", "RL2", 
    #                               "p.fit", "opt.crit", "x.tr.prop",
    #                               "sol.score", "RL1.pen", "RL2.pen"))
    # print(temp)

    #  If scores are identical, **opt.crit** decides 
    #  (with bad fit, p.fit is unreliable)
    if (tab.tmc0["sol.score"] == tab.tmc["sol.score"])
    { # Both scores are identical
      if (tab.tmc0["opt.crit"] < tab.tmc["opt.crit"])
      { # tmc0 is better
        tab.tmc  <- tab.tmc0
        temp.tmc <- temp.tmc0
      }
    } 

    #  Scores differ, solution with higher score is better
    #  @@ To discuss: what if solution with higher score has worse fit?
    if (tab.tmc0["sol.score"] > tab.tmc["sol.score"])
    { # tab.tmc0 is better
      tab.tmc  <- tab.tmc0
      temp.tmc <- temp.tmc0
    }
    # No action needed, if score(tab.tmc) > score(tab.tmc0)
  }

  first.cycle <- FALSE             

  #  Check: which decision taken?
  #cat("\n [Analysis] p.fit of actual optimum ", tab.tmc["p.fit"],
  #    "\n            score of actual solution", tab.tmc["sol.score"],
  #    "\n            iti                     ", iti,
  #    "\n")

  #  Report solution status
  #cat("\n [105_tmc] p.fit of actual optimum ", tab.tmc["p.fit"],
  #    "\n           score of actual solution", tab.tmc["sol.score"],
  #    "\n           iti                     ", iti,
  #    "\n           NoTmcSolution           ", NoTmcSolution,
  #    "\n")


          #  No further x.tr.prop ranges if !NoTmcSolution 
     
          cat("%%%   TMC_seg100_Analysis 1235", date(), "\n")
        
        } #  QQW solution exists
      }   #  TI subset has values
    }     #  loop over TI subsets 

    # ------------------------------------------------------------------------
    #  Is there a TMC solution?
    
    if (print.log.message) 
    { cat("%%%   TMC_seg100_Analysis 1240", date(), "\n") }

    TmcSolution <- (length(tab.tmc) > 1) & (!is.na(tab.tmc[1]))

    if ( TmcSolution )
    {
      #  Replace pseudo-zero estimate by zero for nicer printout
      if (abs(tab.tmc["lambda.tmc"]) < fastnull) { tab.tmc["lambda.tmc"] <- 0 }

      # TMC Solution exists, plot result   

      if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1700\n") }

      if (plot.fig100.016)
      { 
        #  Plot result - basic part
        #  No plot file name, save is done here below
        fig100.016.lim <- PlotHistFit(fig100.016, NA, figtype, 
                                      x.hist, xlabel, NULL, NULL, NULL,
                NA, x.val, counts.range,  
                paste(formatC(100*tab.tmc["x.tr.prop"], format="f", width=5, 
                      digits=1), "%", sep=""),
                x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                x.clip.type, NA, par.las, par.tcl,
                bordercol, histcol, difcol, kdecol, denlty)

        #  Plot result - tmc specific part 
        #  Add estimated tmc solution, RLs, TI
        #  Calculate estimated TMC density (not yet scaled by prev.c.tmc)
        xc.pdf.tmc <- pdf.PN(xsupp, tab.tmc["lambda.tmc"], 
                                   tab.tmc["mue.tmc"], tab.tmc["sigma.tmc"])    

        #  Calculate difference between x.kde and xc.pdf.tmc

        #  Truncate difference below to yplotmin
        #  Switched off 01.05.2021
        #x.kde.tmc <- x.kde$y -  tab.tmc["prev.c.tmc"] * xc.pdf.tmc
        #x.kde.tmc[x.kde.tmc < fig100.016.lim["yplotmin"]] <- 
        #                                           fig100.016.lim["yplotmin"]
  
        PlotMetRes(fig100.016, fig100.016.file, figtype, 
                   fig100.016.lim["yplotmin"],
                   xsupp, tab.tmc["prev.c.tmc"], xc.pdf.tmc, tmccol,
                   x.kde.mode, kdecol, NA, difcol, difbordercol, difhistcol, 
                   tab.tmc["x.RL1.tmc"], tab.tmc["x.RL2.tmc"], 
                   0.85*fig100.016.lim["yplotmax"], 
                   tab.tmc["x.tr.lo"], tab.tmc["x.tr.hi"], 
                   0.90*fig100.016.lim["yplotmin"], tmccol,
                   temp.tmc)

        # -------------------------------------------------------------
        if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1800\n") }

        #  If the data is generated : show unaffected range
        if (!is.na(x.unaffected.lo) & !is.na(x.unaffected.hi))
        { 
          lines(c(x.unaffected.lo, x.unaffected.hi), 
                    0.60*fig100.016.lim["yplotmin"]*c(1,1), col="chartreuse", 
                    lwd=2)
        }

        # -------------------------------------------------------------
        # Save tmc result figure
        #savePlot(file=fig100.016.file, type=figtype)
        myplot <- recordPlot()
      }      #  if (plot.fig100.016)

      # ================================================================
      if (print.log.message) 
      { cat("%%%   TMC_seg100_Analysis 1900", date(), "\n") }

      if ("tmu" %in% meth.list)
      { #  Calculate TMU solution

#  TMC_seg106_Analysis_tmu.R
#
#  02.05.2021  Plot of residual histogram instead of residual density
#  15.02.2021
# ==============================================================================

y.tr <- BoxCox(x[(tab.tmc["x.tr.lo"] <= x) &
                 (x < tab.tmc["x.tr.hi"]) ],
               tab.tmc["lambda.tmc"] )
y.tr.lo <- BoxCox(tab.tmc["x.tr.lo"], 
                  tab.tmc["lambda.tmc"] )
y.tr.hi <- BoxCox(tab.tmc["x.tr.hi"], 
                  tab.tmc["lambda.tmc"] )

tab.tmu <- EstParentDist(y.tr.lo, y.tr.hi, 
                         y.tr, c(tab.tmc["mue.tmc"], tab.tmc["sigma.tmc"]),
                         alpha, tab.tmc["lambda.tmc"], RL1.p, RL2.p)

temp.tmu <- chi2trunc(x.hist,tab.tmc["x.tr.lo"],
                      tab.tmc["x.tr.hi"],
                      sum(x <  tab.tmc["x.tr.lo"]), 
                      sum(x >= tab.tmc["x.tr.hi"]),
                      x.Q1, x.Q2, RL1.p, RL2.p,
                      tab.tmc["lambda.tmc"],
                      tab.tmu["mue.tmu"],
                      tab.tmu["sigma.tmu"],
                      df.est, df.con,
                      l.fact, p.fact, r.fact, w.fact, 
                      opt.crit.only=FALSE,fastnull=fastnull)
tab.tmu <- c(tab.tmu, 
             opt.crit.tmu=unname(temp.tmu[["res"]]["opt.crit"]),
             p.fit.tmu=unname(temp.tmu[["res"]]["chi2.trun.p"]),
             prev.l.tmu=unname(temp.tmu[["res"]]["prev.l.tmc"]),
             prev.c.tmu=unname(temp.tmu[["res"]]["prev.c.tmc"]),
             prev.r.tmu=unname(temp.tmu[["res"]]["prev.r.tmc"]) )

# ------------------------------------------------------------
if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2000\n") }
  
if (plot.fig100.017)
{ 
  #  Plot result - basic part 
  #  No plot file name, save is done here below
  fig100.017.lim <- PlotHistFit(fig100.017, NA, figtype, 
                       x.hist, xlabel, "TMU result", subtitle, n.per.bin.min,
                       NA, x.val, counts.range, 
                       paste(formatC(100*tab.tmc["x.tr.prop"], format="f", width=5, 
                             digits=1), "%", sep=""),
                       x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                       x.clip.type, NA, par.las, par.tcl,
                       bordercol, histcol, difcol, kdecol, denlty)

  #  Plot result - tmu specific part 
  #  Add estimated tmu solution, RLs, TI
  #  Calculate estimated TMC density (not yet scaled by prev.c.tmu)
  xc.pdf.tmu <- pdf.PN(xsupp, tab.tmc["lambda.tmc"], 
                              tab.tmu["mue.tmu"],
                              tab.tmu["sigma.tmu"]) 

  #  Calculate difference between x.kde and xc.pdf.tmu
  #  Truncate difference below to yplotmin
  #  Switched off 01.05.2021
  # x.kde.tmu <- x.kde$y -  tab.tmu["prev.c.tmu"] * xc.pdf.tmu
  # x.kde.tmu[x.kde.tmu < fig100.017.lim["yplotmin"]] <- 
  #                                                 fig100.017.lim["yplotmin"]
             
  PlotMetRes(fig100.017, fig100.017.file, figtype, fig100.017.lim["yplotmin"],
             xsupp, tab.tmu["prev.c.tmu"], xc.pdf.tmu, tmucol,
             x.kde.mode, kdecol, NA, difcol, difbordercol, difhistcol,
             tab.tmu["x.RL1.tmu"], tab.tmu["x.RL2.tmu"], 
             0.85*fig100.017.lim["yplotmax"], 
             tab.tmc["x.tr.lo"], tab.tmc["x.tr.hi"], 
             0.90*fig100.017.lim["yplotmin"], tmccol,
             temp.tmu)

  # -----------------------------------------------------------
  #  If the data is generated : show unaffected range
  if (!is.na(x.unaffected.lo) & !is.na(x.unaffected.hi))
  { 
    lines(c(x.unaffected.lo, x.unaffected.hi), 
            0.60*fig100.017.lim["yplotmin"]*c(1,1), col="chartreuse", lwd=2)
  }

  # -------------------------------------------------------------
  # Save tmu result figure
  savePlot(file=fig100.017.file, type=figtype)
}


 
      }  

      # ================================================================
      if (print.log.message) 
      { cat("%%%   TMC_seg100_Analysis 2000", date(), "\n") }           

      #  Calculate Tukey solution, using lambda from tmc

      if ("tuk" %in% meth.list)
      { #  Calculate Tukey solution

#  TMC_seg108_Analysis_tuk.R 
#
#  Tukey analysis
#
#  25.04.2021
# ==============================================================================

if (print.log.message) { cat("%%%   TMC_seg108_Analysis  100\n") }

#  Transform data according to tmc result
y <- BoxCox(x, tab.tmc["lambda.tmc"])

#  Run Tukey analysis

#tab.tuk <- TukeyWW(y, xlabel, IQR.fact=IQR.fact, RL1.p=RL1.p, RL2.p=RL2.p, 
#                    figA=NA)
tab.tuk <- TukeyWW(y, xlabel, IQR.fact=IQR.fact, iter.max=iter.max.tuk, 
                    x.tr.prop.min=x.tr.prop.min, RL1.p=RL1.p, RL2.p=RL2.p, 
                    print.details=FALSE, figA=NA)

print(tab.tuk)

#return(c(mean=mean.hat, sd.raw=sigma.raw, sd.adj=sigma.adj, 
#           x.tr.lo=x.lo, x.tr.hi=x.hi, 
#           out.lo=out.lo, out.hi=out.hi, x.RL1=x.RL1, x.RL2=x.RL2, xc.n=xc.n)) 

#  Back transformation to x scale
x.tr.lo.tuk <- BoxCoxInv(tab.tuk["x.tr.lo"], tab.tmc["lambda.tmc"]) 
x.tr.hi.tuk <- BoxCoxInv(tab.tuk["x.tr.hi"], tab.tmc["lambda.tmc"]) 
x.RL1.tuk <- BoxCoxInv(tab.tuk["x.RL1"], tab.tmc["lambda.tmc"]) 
x.RL2.tuk <- BoxCoxInv(tab.tuk["x.RL2"], tab.tmc["lambda.tmc"]) 

tab.tuk["x.tr.lo"] <- x.tr.lo.tuk
tab.tuk["x.tr.hi"] <- x.tr.hi.tuk
tab.tuk["x.RL1.tuk"] <- x.RL1.tuk
tab.tuk["x.RL2.tuk"] <- x.RL2.tuk

temp.tuk <- chi2trunc(x.hist, x.tr.lo.tuk, x.tr.hi.tuk,
                      sum(x < x.tr.lo.tuk), 
                      sum(x >= x.tr.hi.tuk),
                      x.Q1, x.Q2, RL1.p, RL2.p,
                      tab.tmc["lambda.tmc"],
                      tab.tuk["mean"],
                      tab.tuk["sd.adj"],
                      df.est, df.con,
                      l.fact, p.fact, r.fact, w.fact, 
                      opt.crit.only=FALSE,fastnull=fastnull)

#save(temp.tuk, file="../dev/temp.tuk")   # used for development of HistoResid


#  Estimated prevalence cannot be taken from chi2trunc because the Tukey 
#  approach decomposes the data into non-overlapping components, while tmc, 
#  tmu assume overlaps.
              
tab.tuk <- c(tab.tuk, 
             opt.crit.tuk=unname(temp.tuk[["res"]]["opt.crit"]),
             p.fit.tuk=unname(temp.tuk[["res"]]["chi2.trun.p"]),
             prev.l.tuk=unname(tab.tuk["out.lo"]/x.n),
             prev.c.tuk=unname((x.n-tab.tuk["out.lo"]-tab.tuk["out.hi"])/
                               x.n), 
             prev.r.tuk=unname(tab.tuk["out.hi"]/x.n) )

# ------------------------------------------------------------
if (print.log.message) { cat("%%%   TMC_seg108_Analysis  200\n") }
  
if (plot.fig108.010)
{ 
  #  Plot result - basic part 
  #  No plot file name, save is done here below
  fig108.010.lim <- PlotHistFit(fig108.010, NA, figtype,
                    x.hist, xlabel, paste("Tukey result, IQR.fact", IQR.fact), 
                    subtitle, n.per.bin.min,
                    NA, x.val, counts.range, 
                 paste(formatC(100*tab.tuk["xc.n"]/x.n, format="f", width=5,
                   digits=1), "%", sep=""),
                   x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                   x.clip.type, NA, par.las, par.tcl,
                   bordercol, histcol, difcol, kdecol, denlty)

  #  Plot result - tuk specific part 
  #  Add estimated tuk solution, RLs, TI
  #  Calculate estimated tuk density (not yet scaled by prev.c.tuk)
  xc.pdf.tuk <- pdf.PN(xsupp, tab.tmc["lambda.tmc"], 
                              tab.tuk["mean"],
                              tab.tuk["sd.adj"]) 

  #  Calculate difference between x.kde and xc.pdf.tuk
  #  Truncate difference below to yplotmin
  #  Switched off 01.05.2021
  #x.kde.tuk <- x.kde$y -  tab.tuk["prev.c.tuk"] * xc.pdf.tuk
  #x.kde.tuk[x.kde.tuk < fig108.010.lim["yplotmin"]] <- 
  #                                             fig108.010.lim["yplotmin"]
           
  PlotMetRes(fig108.010, fig108.010.file, figtype, fig108.010.lim["yplotmin"], 
             xsupp, tab.tuk["prev.c.tuk"], xc.pdf.tuk, tukcol,
             x.kde.mode, kdecol, NA, difcol, difbordercol, difhistcol,
             tab.tuk["x.RL1.tuk"], tab.tuk["x.RL2.tuk"], 
             0.85*fig108.010.lim["yplotmax"], 
             tab.tuk["x.tr.lo"], tab.tuk["x.tr.hi"], 
             0.90*fig108.010.lim["yplotmin"], tukcol,
             temp.tuk)

  # -----------------------------------------------------------
  #  If the data is generated : show unaffected range
  if (!is.na(x.unaffected.lo) & !is.na(x.unaffected.hi))
  { 
    lines(c(x.unaffected.lo, x.unaffected.hi), 
          0.60*fig108.010.lim["yplotmin"]*c(1,1), col="chartreuse", lwd=2)
  }

  # -----------------------------------------------------------
  #  Save the tuk result plot

  savePlot(file=fig108.010.file, type=figtype) 
}

      }  

      if (print.log.message) 
      { cat("%%%   TMC_seg100_Analysis 2100", date(), "\n") }           

      # ================================================================

    }              #  TMC solution exists
  }                #  x.hist exists

  # ==========================================================================
  if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2200\n") }
  # ==========================================================================

  #  Combine essential results from all methods for this stratum
  #  tab.stra defined in the beginning

  if ("npa" %in% meth.list)
  {
    tab.stra["npa", "x.tr.lo"] <- tab.npa["x.tr.lo.npa"]
    tab.stra["npa", "x.tr.hi"] <- tab.npa["x.tr.hi.npa"]
    tab.stra["npa", "RL1"]     <- tab.npa["x.RL1.npa"]
    tab.stra["npa", "RL2"]     <- tab.npa["x.RL2.npa"]
    tab.stra["npa", "prop"]    <- tab.npa["x.tr.prop.npa"]
    tab.stra["npa", "prev.l"]  <- tab.npa["prev.l.npa"]
    tab.stra["npa", "prev.c"]  <- tab.npa["prev.c.npa"]
    tab.stra["npa", "prev.r"]  <- tab.npa["prev.r.npa"]
  }

  if ( ("qqw" %in% meth.list) & !NoQqwSolution )
  {
    tab.stra["qqw", "x.tr.lo"] <- tab.qqw[1, "x.tr.lo"]
    tab.stra["qqw", "x.tr.hi"] <- tab.qqw[1, "x.tr.hi"]
    tab.stra["qqw", "bins.n"]  <- tab.qqw[1, "bins.n"]
    tab.stra["qqw", "prop"]    <- tab.qqw[1, "prop.qqw"]
    tab.stra["qqw", "lambda"]  <- tab.qqw[1, "lambda.qqw"]
    tab.stra["qqw", "mue"]     <- tab.qqw[1, "mue.qqw"]
    tab.stra["qqw", "sigma"]   <- tab.qqw[1, "sigma.qqw"]
    tab.stra["qqw", "RL1"]     <- tab.qqw[1, "x.RL1.qqw"]
    tab.stra["qqw", "RL2"]     <- tab.qqw[1, "x.RL2.qqw"]
    tab.stra["qqw", "RL1.cilo"] <- tab.qqw[1, "x.RL1.cilo"]
    tab.stra["qqw", "RL1.cihi"] <- tab.qqw[1, "x.RL1.cihi"]
    tab.stra["qqw", "RL2.cilo"] <- tab.qqw[1, "x.RL2.cilo"]
    tab.stra["qqw", "RL2.cihi"] <- tab.qqw[1, "x.RL2.cihi"]
    tab.stra["qqw", "bins.n"]  <- tab.qqw[1, "bins.n"]
    tab.stra["qqw", "x.tr.n"]  <- tab.qqw[1, "x.tr.n"]   # may refer to reduced
                                                         # histogram 
    tab.stra["qqw", "prop"]    <- tab.qqw[1, "prop.qqw"] # may refer to reduced
                                                         # histogram 
    tab.stra["qqw", "opt.crit"] <- unname(temp.qqw[["res"]]["opt.crit"])
    tab.stra["qqw", "p.fit"]   <- unname(temp.qqw[["res"]]["chi2.trun.p"])
    tab.stra["qqw", "prev.l"]  <- tab.qqw[1, "prev.l.qqw"]
    tab.stra["qqw", "prev.c"]  <- tab.qqw[1, "prev.c.qqw"]
    tab.stra["qqw", "prev.r"]  <- tab.qqw[1, "prev.r.qqw"]
  }

  if ("bha" %in% meth.list)
  {
    #  @@@ subpopulation auswählen
    tab.stra["bha", "RL1"]    <- tab.bha[["subpop"]]["comp01", "RL1.2"]
    tab.stra["bha", "RL2"]     <- tab.bha[["subpop"]]["comp01", "RL2.2"]
  }

  if (("tmc" %in% meth.list) & x.hist.exists & TmcSolution)
  {
    tab.stra["tmc", "x.tr.lo"] <- tab.tmc["x.tr.lo"]
    tab.stra["tmc", "x.tr.hi"] <- tab.tmc["x.tr.hi"]
    tab.stra["tmc", "x.tr.n"]  <- tab.tmc["x.tr.n"]
    tab.stra["tmc", "bins.n"]  <- tab.tmc["ihi"] - tab.tmc["ilo"] + 1
    tab.stra["tmc", "prop"]    <- tab.tmc["x.tr.prop"]  
    tab.stra["tmc", "rc"]      <- tab.tmc["rc"]  
    tab.stra["tmc", "iter"]    <- tab.tmc["iter"]  
    tab.stra["tmc", "lambda"]  <- tab.tmc["lambda.tmc"]
    tab.stra["tmc", "mue"]     <- tab.tmc["mue.tmc"]
    tab.stra["tmc", "sigma"]   <- tab.tmc["sigma.tmc"]
    tab.stra["tmc", "RL1"]     <- tab.tmc["x.RL1.tmc"]
    tab.stra["tmc", "RL2"]     <- tab.tmc["x.RL2.tmc"]

    tab.stra["tmc", "RL1.cilo"] <- tab.tmc["x.RL1.cilo"]
    tab.stra["tmc", "RL1.cihi"] <- tab.tmc["x.RL1.cihi"]
    tab.stra["tmc", "RL2.cilo"] <- tab.tmc["x.RL2.cilo"]
    tab.stra["tmc", "RL2.cihi"] <- tab.tmc["x.RL2.cihi"]

    tab.stra["tmc", "opt.crit"] <- tab.tmc["opt.crit"]
    tab.stra["tmc", "p.fit"]   <- tab.tmc["p.fit"]
    tab.stra["tmc", "bins.n"]  <- tab.tmc["ihi"] - 
                                  tab.tmc["ilo"] + 1
    tab.stra["tmc", "prop"]    <- tab.tmc["x.tr.prop"]
    tab.stra["tmc", "prev.l"]  <- tab.tmc["prev.l.tmc"]
    tab.stra["tmc", "prev.c"]  <- tab.tmc["prev.c.tmc"]
    tab.stra["tmc", "prev.r"]  <- tab.tmc["prev.r.tmc"]
  }

  if ("tmu" %in% meth.list & x.hist.exists & TmcSolution)
  {
    tab.stra["tmu", "x.tr.lo"] <- tab.tmc["x.tr.lo"]
    tab.stra["tmu", "x.tr.hi"] <- tab.tmc["x.tr.hi"]
    tab.stra["tmu", "x.tr.n"]  <- tab.tmc["x.tr.n"]
    tab.stra["tmu", "bins.n"]  <- tab.tmc["ihi"] - tab.tmc["ilo"] + 1
    tab.stra["tmu", "prop"]    <- tab.tmc["x.tr.prop"]  
    tab.stra["tmu", "rc"]      <- tab.tmu["rc.tmu"]  
    tab.stra["tmu", "iter"]    <- tab.tmu["iter.tmu"]  
    tab.stra["tmu", "lambda"]  <- tab.tmc["lambda.tmc"]
    tab.stra["tmu", "mue"]     <- tab.tmu["mue.tmu"]
    tab.stra["tmu", "sigma"]   <- tab.tmu["sigma.tmu"]
    tab.stra["tmu", "RL1"]     <- tab.tmu["x.RL1.tmu"]
    tab.stra["tmu", "RL2"]     <- tab.tmu["x.RL2.tmu"]

    tab.stra["tmu", "RL1.cilo"] <- tab.tmu["x.RL1.cilo"]
    tab.stra["tmu", "RL1.cihi"] <- tab.tmu["x.RL1.cihi"]
    tab.stra["tmu", "RL2.cilo"] <- tab.tmu["x.RL2.cilo"]
    tab.stra["tmu", "RL2.cihi"] <- tab.tmu["x.RL2.cihi"]
    tab.stra["tmu", "opt.crit"] <- tab.tmu["opt.crit.tmu"]
    tab.stra["tmu", "p.fit"]   <- tab.tmu["p.fit.tmu"]
    tab.stra["tmu", "prev.l"]  <- tab.tmu["prev.l.tmu"]
    tab.stra["tmu", "prev.c"]  <- tab.tmu["prev.c.tmu"]
    tab.stra["tmu", "prev.r"]  <- tab.tmu["prev.r.tmu"]
  }                #  do tmu

  if ("tuk" %in% meth.list & x.hist.exists & TmcSolution)
  {
    tab.stra["tuk", "x.tr.lo"] <- tab.tuk["x.tr.lo"]
    tab.stra["tuk", "x.tr.hi"] <- tab.tuk["x.tr.hi"]
    tab.stra["tuk", "x.tr.n"]  <- tab.tuk["x.tr.n"]
    tab.stra["tuk", "bins.n"]  <- NA
    tab.stra["tuk", "prop"]    <- (x.n-tab.tuk["out.lo"]-tab.tuk["out.hi"])/ 
                                  x.n  
    tab.stra["tuk", "rc"]      <- NA 
    tab.stra["tuk", "iter"]    <- NA  
    tab.stra["tuk", "lambda"]  <- tab.tmc["lambda.tmc"]
    tab.stra["tuk", "mue"]     <- tab.tuk["mean"]
    tab.stra["tuk", "sigma"]   <- tab.tuk["sd.adj"]
    tab.stra["tuk", "RL1"]     <- tab.tuk["x.RL1.tuk"]
    tab.stra["tuk", "RL2"]     <- tab.tuk["x.RL2.tuk"]

    tab.stra["tuk", "RL1.cilo"] <- NA
    tab.stra["tuk", "RL1.cihi"] <- NA
    tab.stra["tuk", "RL2.cilo"] <- NA
    tab.stra["tuk", "RL2.cihi"] <- NA
    tab.stra["tuk", "opt.crit"] <- tab.tuk["opt.crit.tuk"]
    tab.stra["tuk", "p.fit"]    <- tab.tuk["p.fit.tuk"]
    tab.stra["tuk", "prev.l"]   <- tab.tuk["prev.l.tuk"]
    tab.stra["tuk", "prev.c"]   <- tab.tuk["prev.c.tuk"]
    tab.stra["tuk", "prev.r"]   <- tab.tuk["prev.r.tuk"]
  }                #  do tuk

  # --------------------------------------------------------------------------
  #  Error code for result. Considers only goodness of fit. Possible 
  #  extension: check estimated prevalence for reasonable values
  tab.stra[ , "errcode"] <- " ?"
  ok <- (tab.stra[ , "p.fit"] > p.fit.min)
  ok[is.na(ok)] <- FALSE
  tab.stra[ok, "errcode"] <- "ok"

  ok <- (tab.stra[ , "p.fit"] <= p.fit.min)
  ok[is.na(ok)] <- FALSE
  tab.stra[ok, "errcode"] <- "!!"

  # --------------------------------------------------------------------------
  #  Toleranzintervalle für die tmc-RLs per bootstrapping berechnen
  #  NOT YET ADAPTED TO NEW CALL OF tmc.master

#  TMC_seg101_TMC_BS.R

  if (go.on & (nbs > 0))
  {  
    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 910\n") }

    #  Generate plot windows, if requested
    if (plot.fig100.115)          # Criteria for prop.tmc selection when 
                                  # bootstrapping
    { 
      if (!exists("fig100.115"))
      { #  Open graph new window
        win.graph()
        fig100.115 <- dev.cur()
      }
      fig100.115.file <- paste(path.fig.stra,outname.stra,"-F100.115.",
                             figtype,sep="")
      if (file.exists(fig100.115.file)) { file.remove(fig100.115.file) }
    }

    if (plot.fig100.125)          # Estimated RLs vs prop.tmc when bootstrapping
    { if (!exists("fig100.125"))
      { #  Open graph new window
        win.graph()
        fig100.125 <- dev.cur()
      }
      fig100.125.file <- paste(path.fig.stra,outname.stra,"-F100.125.",
                               figtype,sep="")
      if (file.exists(fig100.125.file)) { file.remove(fig100.125.file) }
    }

    if (plot.fig100.130)          #  Convergence of RL estimates when bootstrapping 
    { if (!exists("fig100.130"))
      { #  Open graph new window
        win.graph()
        fig100.130 <- dev.cur()
      }
      fig100.130.file <- paste(path.fig.stra,outname.stra,"-F100.130.",
                               figtype,sep="")
      if (file.exists(fig100.130.file)) { file.remove(fig100.130.file) }
    }
     
    if (plot.fig100.140)          #  Bootstrapped histograms and densities

    { if (!exists("fig100.140"))
      { #  Open graph new window
        win.graph()
        fig100.140 <- dev.cur()
      }
      fig100.140.file <- paste(path.fig.stra,outname.stra,"-F100.140.",
                               figtype,sep="")
      if (file.exists(fig100.140.file)) { file.remove(fig100.140.file) }
    }
     
    if (plot.fig100.145)          #  Densities of bootstrapped results
    { if (!exists("fig100.145"))
      { #  Open graph new window
        win.graph()
        fig100.145 <- dev.cur()
      }
      fig100.145.file <- paste(path.fig.stra,outname.stra,"-F100.145.",
                               figtype,sep="")
      if (file.exists(fig100.145.file)) { file.remove(fig100.145.file) }
    }
     
    #  Berechnung von Toleranzintervallen per Bootstrap ist erwünscht
    #  Tabelle für Einzelergebnisse anlegen
    tab.bs.names     <- c("lambda","mue","sigma","RL1","RL2","prev",
                          "prop", "x.tr.lo", "x.tr.hi", "p.fit")
    tab.bs           <- matrix(NA,nrow=nbs,ncol=length(tab.bs.names))
    colnames(tab.bs) <- tab.bs.names

    theta.ini <- c(lambda.tmc, mue.tmc, sigma.tmc)

    if (plot.fig100.140)
    { dev.set(fig100.140)
      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(density(x), col="black", lwd=3,
           xlim=c(0,40),
           main="Data and BS densities",
           sub=subtitle, cex.sub=0.7)
    }

    for (ibs in 1:nbs)
    {
        t1 <- Sys.time() 

        #  Bootstrap sample herstellen
  
        xbs <- sample(x, x.n, replace=TRUE)   #  empirical BS

        #xbs <- rxref(x.kde.pdf$x,x.cdf,x.n) # nichtparametrischer bootstrap,
        #                                # beruhend auf der Dichteschätzung der 
        #                                # Daten

        xbs <- sort(xbs)
        #save(xbs, file="SmoothHistTest.RData")
        #save(x.breaks, file="SmoothHistBreaks.RData")

        #cat("   ibs:",ibs, "mean(x) ", mean(x), " mean(xbs) ", mean(xbs), "\n")
        cat("   ibs:",ibs, "mean(x) ", mean(x), " mean(xbs) ", mean(xbs), "\x0D")
        #cat("   ibs:",ibs,  "\x0D")
         
        #  Generate histogram and density estimate for the xbs sample
        #  Parameter 'freq' not used here
        xbs.hist <- hist(xbs,breaks=x.breaks, right=FALSE,
                         plot=FALSE,
                         warn.unused=FALSE)
        xbs.kde  <- density(xbs,n=kernel.n,adjust=x.adjust,from=x.from,to=x.to)

        #  Smooth the histogram
        xbs.hist.smo <- SmoothHist(xbs.hist, fig=FALSE)

        save(xbs.hist, file="xbs.hist.RData")
        save(xbs.hist.smo, file="xbs.hist.smo.RData")

        if (plot.fig100.140)
        { dev.set(fig100.140)       
          lines(xbs.kde, col="green3", lwd=1)
          par(las=par.las)    
          par(tcl=par.tcl)    
          plot(xbs.hist, add=TRUE, border="cyan")
        } 

        #  Find mode from kde
        xbs.kde.mode.idx <- which.max(xbs.kde$y)

        #  If there is more than 1 maximum: message
        if (length(xbs.kde.mode.idx) > 1)
        {
          cat(
               "\n\n +++ [Analysis] BS density estimate has > 1 maximum - first is taken\n\n" )
          xbs.kde.mode.idx <- xbs.kde.mode.idx[1] 
        } 
        xbs.kde.mode     <- xbs.kde$x[xbs.kde.mode.idx]

        xbs.breaks.mode.idx <- which((x.breaks[1:(x.breaks.n-1)] <= xbs.kde.mode) &
                                     (xbs.kde.mode < x.breaks[2:x.breaks.n] ) ) 

        #  Plot history: show used time
        #t2 <- Sys.time() 
        #delta21 <- difftime(t2,t1)
        ##  Verbrauchte Zeit darstellen
        #if (ibs==1)
        #{ xplotmin11 <- 1
        #  xplotmax11 <- nbs
        #  yplotmin11 <- 0      
        #  yplotmax11 <- 10 * as.numeric(delta21)

        #  plot(ibs,delta21,type="h",col="blue",
        #         xlim=c(xplotmin11,xplotmax11),ylim=c(yplotmin11,yplotmax11),
        #         main="Bootstrap progress",
        #         xlab="# of BS sample",ylab="Time used")
        #}  else
        #{  
        #  lines(ibs,delta21,type="h",col="blue")
        #}

        #  Free choice of prop.tmc produces extreme fluctuation of the 
        #  truncation intervel and thus in the estimates.   
        TMC.bs <- tmc.master(df.est, df.con, ErrMsg, figtype, 
                       outname.stra, p.fact, path.fig.stra,
                       pathol.position, 
                       plot.fig100.115, plot.fig100.125,
                       fig100.115, fig100.125,
                       fig100.115.file, fig100.125.file,   
                       prop.tmc.seq0, prop.tmc.seq.n,
                       RL1.p, RL2.p, theta.ini, w.fact,
                       xbs.breaks.mode.idx, xbs.hist.smo, xbs.kde,  
                       x.n, x.n.min, x.tr.bins.min,
                       print.log.message, fastnull, fastnull.chi)

        tab.bs[ibs,"lambda"] <- TMC.bs[["est"]]["lambda"]
        tab.bs[ibs,"mue"]    <- TMC.bs[["est"]]["mue"]
        tab.bs[ibs,"sigma"]  <- TMC.bs[["est"]]["sigma"]
        tab.bs[ibs,"RL1"]    <- TMC.bs[["est"]]["x.RL1"]
        tab.bs[ibs,"RL2"]    <- TMC.bs[["est"]]["x.RL2"]
        tab.bs[ibs,"prev"]   <- TMC.bs[["fit"]]["prev.tmc"]

        tab.bs[ibs,"prop"]      <- TMC.bs[["est"]]["prop.tmc"]
        tab.bs[ibs,"x.tr.lo"]   <- TMC.bs[["est"]]["x.tr.lo"]
        tab.bs[ibs,"x.tr.hi"]   <- TMC.bs[["est"]]["x.tr.hi"]
        tab.bs[ibs,"p.fit"]     <- TMC.bs[["est"]]["p.fit"]

        # show estimated density 
        lines(x.kde.pdf$x, (1-tab.bs[ibs,"prev"])*
                       pdf.PN(x.kde.pdf$x, tab.bs[ibs,"lambda"],
                                       tab.bs[ibs,"mue"],
                                       tab.bs[ibs,"sigma"]), 
              col=boscol)

        #t3 <- Sys.time()
        #delta32 <- difftime(t3,t2)
        #print(delta32)
        #lines(c(ibs,ibs),c(delta21,delta21+delta32),col="green3")
        #cat("   ibs:",ibs,  difftime(t3,t1),"\x0D")

        # Alternative to elapsed time: show convergence of RLs
        if (plot.fig100.130)
        { dev.set(fig100.130)
          if (ibs==1)
          { xplotmin11 <- 1
            xplotmax11 <- nbs
            yplotmin11 <- 0      
            yplotmax11 <- 2 * tab.bs[ibs,"RL2"]

            par(las=par.las)    
            par(tcl=par.tcl)    
            plot(ibs,tab.bs[ibs,"RL1"],type="p",col=boscol, pch=3, 
               xlim=c(xplotmin11,xplotmax11),ylim=c(yplotmin11,yplotmax11),
               main="Bootstrap progress",
               xlab="BS sample #",ylab="Mean RIL",
               sub=subtitle, cex.sub=0.7)
            points(ibs,tab.bs[ibs,"RL2"],type="p",col=boscol, pch=3)  
          }  else
          {  
            points(ibs,tab.bs[ibs, "RL1"],type="p",col="gray", pch=20)
            points(ibs,tab.bs[ibs, "RL2"],type="p",col="gray", pch=20)
            points(ibs,median(tab.bs[1:ibs, "RL1"]),type="p",col=boscol, pch=3)
            points(ibs,median(tab.bs[1:ibs, "RL2"]),type="p",col=boscol, pch=3)
          }
        }

    }  # if (plot ...
    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 920\n") }

    if ((nbs > 0) & plot.fig100.130)
    { #  Add tmc estimate
      dev.set(fig100.130)
      abline(h=x.RL1.tmc, col=tmccol)
      abline(h=x.RL2.tmc, col=tmccol)
      abline(h=TMC[["ci"]][3], col=tmccol, lty=3)
      abline(h=TMC[["ci"]][4], col=tmccol, lty=3)

      abline(h=median(tab.bs[1:nbs, "RL1"]), col=boscol)
      abline(h=median(tab.bs[1:nbs, "RL2"]), col=boscol)
      abline(h=Quantile(tab.bs[1:nbs, "RL2"], probs=0.025), col=boscol, lty=3)
      abline(h=Quantile(tab.bs[1:nbs, "RL2"], probs=0.975), col=boscol, lty=3)
    }
    savePlot(file=fig100.130.file, type=figtype)

    #cat("Single nonparametric bootstrap results for tolerance interval\n") 
    #print(tab.bs[order(tab.bs[ ,"RL1"]), ])

    if (plot.fig100.145)
    { dev.set(fig100.145)
      par(mfrow=c(2,2))
      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(density(x),main="Estimated density of x")
      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(density(tab.bs[ ,"RL1"]),main="BS density of RL1",
         sub=paste("nbs =",nbs))
      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(density(tab.bs[ ,"RL2"]),main="BS density of RL2",
         sub=paste("nbs =",nbs))
      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(density(tab.bs[ ,"prev"]),main="BS density of prev",
         sub=paste("nbs =",nbs))
    }  # if (plot...
    
    x.RL.est["tmc.tmc","Prevtilo"] <- Quantile(tab.bs[ ,"prev"],probs=0.025) 
    x.RL.est["tmc.tmc","Prevtihi"] <- Quantile(tab.bs[ ,"prev"],probs=0.975) 
    x.RL.est["tmc.tmc","x.RL1bsmea"] <- mean(tab.bs[ ,"RL1"]) 
    x.RL.est["tmc.tmc","x.RL1tilo"] <- Quantile(tab.bs[ ,"RL1"],probs=0.025) 
    x.RL.est["tmc.tmc","x.RL1tihi"] <- Quantile(tab.bs[ ,"RL1"],probs=0.975) 
    x.RL.est["tmc.tmc","x.RL2bsmea"] <- mean(tab.bs[ ,"RL2"]) 
    x.RL.est["tmc.tmc","x.RL2tilo"] <- Quantile(tab.bs[ ,"RL2"],probs=0.025) 
    x.RL.est["tmc.tmc","x.RL2tihi"] <- Quantile(tab.bs[ ,"RL2"],probs=0.975) 
 
    #cat("[seg100] Bootstrap-Ergebnis für tmc\n") 
    #print(x.RL.est["tmc.tmc", ])
  
    #  Nutzung der CIs in SplineVerlauf()

    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 930\n") }
  }     #   if (go.on & nbs > 0) 



  # --------------------------------------------------------------------------

  if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 5000\n") }

 
  #  ==================================================================
  #  Print results tables per stratum (not if replicates are analysed) 

  # if (!eval.rep) { sink(file=outname.stra.txt,split=TRUE) }

  cat("\n\n", sep.e, 
      "\n\nRL ESTIMATES, ALL REQUESTED METHODS, ONE age-sex STRATUM",
      "\nUser                                   ", user,
      "\nData information taken from            ", info.file,
      "\nAnalysis date                          ", date(),
      "\nCopy of this table in                  ",
      "\n    ",                                   outname.stra.txt,
      "\nData file name                         ",
      "\n    ",                                   infile,
      "\nFilter  outpatient/hospitalized        ", use.oh,
      "\nFilter  device                         ", use.dev,
      "\nFilter  analysis time of day           ", paste(ana.hour.min, ":",
                                                         ana.minu.min, " - ",
                                                         ana.hour.max, ":",
                                                         ana.minu.max, sep=""),
      "\nSex group                              ", sexlabel,    
      "\nAge interval                           ", agelabel,
      "\nInput values are multiplied by         ", scale.fact,
      "\nScaled values are limited by           ", x.lo.limit,x.hi.limit,
      "\nScaled values are rounded here to      ",round.unit,
      "\nTotal number of all values in file     ",x.n1,
      "\nTotal number of numerical values       ",x.n14,
      "\nNumber of values in this stratum       ",x.n,
      "\nNumber of distinct values in stratum   ",x.val.n,
      "\nMaximal detection/determination limit  ",detect.limits.max,
      "\n# of values / % below detection limit  ",x.lt.detect.n, " / ",
                                                  format(x.lt.detect.pct,digits=3),
                                                  "%",
      "\nValue range                            ",x.min, x.max,
      "\nEmpirical P025, P975, complete stratum ",x.Q1, x.Q2,
      "\n",
      "\nSmooth empirical histogram?            ",smooth.hist1, smooth.hist2,
      "\nMin. # of values per bin               ",n.per.bin.min,
      "\nMin. # of bins in stratum              ",bins.n.min,     
      "\nMin. # of bins in truncation interval  ",x.tr.bins.min,     
      "\nMin. # of values in stratum            ",x.n.min,     
      "\nActual # of bins in stratum            ",bins.n.act,
      "\nUser-specified min. prop. in TI        ", x.tr.prop.min, 
      "\nUser-specified max. prop. in TI        ", x.tr.prop.max, 
      "\nWeight of truncation interval length   ",l.fact,
      "\nPenalty factor 'wrong prevalence'      ",p.fact,
      "\nPenalty factor 'wrong prediction'      ",w.fact,
      "\nMinimal lambda                         ",lambda.min,
      "\nMaximal lambda                         ",lambda.max,
      "\n",
      "\nTMC: Truncation interval               ",tab.tmc["x.tr.lo"], 
                                                  tab.tmc["x.tr.hi"],   
      "\nTMC: TI contains proportion ...        ",tab.tmc["x.tr.prop"],
      "\nTMC: estimated lambda                  ",tab.tmc["lambda.tmc"],
      "\nTMC: estimated mu                      ",tab.tmc["mue.tmc"],
      "\nTMC: estimated sigma                   ",tab.tmc["sigma.tmc"],
      "\nTMC: chi^2 from truncation interval    ",format(
                                                    tab.tmc["chi2.trun"],
                                                    digits=4),
      "\nTMC: chi^2 from outside TI             ",format(
                                                   tab.tmc["chi2.total"]-
                                                   tab.tmc["chi2.trun"],
                                                   digits=4),
      "\n",
      "\nFalse negative probability  <RL1 / <RL2", RL1.p,"/ ",RL2.p,
      "\nTMC: estimated RLs                       RL1 = ", 
                                             format(tab.tmc["x.RL1.tmc"],
                                                    digits=4),      
      "\n                                         RL2 = ", 
                                             format(tab.tmc["x.RL2.tmc"],
                                                    digits=4), 
      "\nAssessment of TMC solution",
      "\np for chi^2 goodness of fit            ",format(
                                                    tab.tmc["p.fit"],
                                                    digits=4), 
      "\nRL1 >= x.Q1, RL2 <= x.Q2               ", crit1, crit2,
      "\nprev.l ok / prev.c ok /prev.r ok?      ", crit3, crit4, crit5,
      "\nError structure ok?                    ", crit6,
      "\nTotal score for solution (max =",5+crit123456,")    ", 
                                                tab.tmc["sol.score"],    
      "\n", sep.e, 
      "\n\n")


  if (any(ErrMsg[ ,"Occurred"]))
  { cat("\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    print(ErrMsg[ErrMsg[ ,"Occurred"], "ErrTxt"])
    cat(   "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  } 

  if (user == "WW")
  { #  Ausführung für WW

    #cat("\n All details \n")
    #print(tab.stra[meth.list, ])
    #cat("\n")
  }

  if (!eval.rep) { sink() }
 
  # ==========================================================================
  if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 6000\n") }
  # ==========================================================================

}  #  if (x.n < x.n.min)

# ---------------------------------------------------------------------------
#  Copy results in tab.stra (one stratum) to tab.f (one file)
#  Preliminary: gtab instead of tab.f

# ==========================================================================
if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 6100\n") }
# ==========================================================================

#  block.start controls the identification by sexlabel and agelabel in the 
#  printout

#  Put only rows with information in gtab

if (!is.na(tab.stra[1, "irep"]) )
{
  block.start          <- TRUE
  iblock               <- iblock + 1
  for (meth in gtab.list)
  { igtab              <- igtab + 1

    meth.idx <- tab.stra[ ,"method"] == meth
  
    #  Append tab.stra to gtab
    if (iblock == 1 & igtab == 1)
    {
      #  First block at all, establish gtab
      gtab <- tab.stra[meth.idx, ]
    }  else
    {
      #  Append block to gtab
      gtab <- rbind(gtab, tab.stra[meth.idx, ])
    }

    if (block.start)
    { 
      gtab[igtab,"Sex"]     <- sexlabel
      gtab[igtab,"Age"]     <- agelabel
    } else 
    { 
      gtab[igtab,"Sex"]     <- " "
      gtab[igtab,"Age"]     <- " "
    }
    block.start          <- FALSE

    # ==========================================================================
    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 6400\n") }
    # ==========================================================================

  }    #  information for stratum is available
  cat("\n", sep.e, 
    "\n   End analysis by seg100 for", 
    "\n   ", outname.stra,
    "\n", sep.e, 
    "\n", sep="")
   
}      #  x.n > x.n.min

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg100_Analysis  End\n") }
# ============================================================================




#############################################################################
return(list(myplot = myplot, RL1 = tab.tmc["x.RL1.tmc"], RL2 = tab.tmc["x.RL2.tmc"]))
}




#  TMC_seg022_Functions.R 
#  (c) wwosniok@math.uni-bremen.de 
 
#  Truncated minimum chi-square estimation 
 
#  Compilation of all functions required by TMC 
#  see ../func/CompileFunctions.R for the generation of this file 
 
#  #  CHANGE HISTORY 
#  10.04.2021 All actual experimental functions included here 
#  16.02.2021 Actual functions collected in TMC_seg022_Functions.R 
#  26.01.2021 Start as update of TMC_seg020_Functions.R 
# ============================================================================== 

#  Functions compiled in this file:
# [1] "F_Bhatta_V2.R"       "F_BoxCox.R"          "F_BoxCoxInv.R"      
# [4] "F_CalcChi2.R"        "F_CalcHistLoop.R"    "F_CalcParDistance.R"
# [7] "F_CalcPrev_V2.R"     "F_CalcR2InSubInt.R"  "F_cdf.PN.R"         
#[10] "F_Ceiling.R"         "F_CheckRounding.R"   "F_chi2.PNV.tr.lms.R"
#[13] "F_chi2trunc_V4.R"    "F_chi2trunc0_V3.R"   "F_CIQuant.LNV.R"    
#[16] "F_CIQuant.PNV.R"     "F_CollapseHist.R"    "F_Compose.R"        
#[19] "F_Drift.R"           "F_EstParentDist.R"   "F_Expand.R"         
#[22] "F_Extract.R"         "F_FC.R"              "F_FilterWDay.R"     
#[25] "F_FindBw_V4.R"       "F_FindHisto_V2.R"    "F_FindKDE.R"        
#[28] "F_FindModeIndex.R"   "F_Floor.R"           "F_GenDatCalcChi2.R" 
#[31] "F_hist.table.R"      "F_Infoblock.R"       "F_IniMatrix.R"      
#[34] "F_Interpolx.R"       "F_Interpoly.R"       "F_o2R.R"            
#[37] "F_pdf.PN.R"          "F_PermissibleDiff.R" "F_plot.res.R"       
#[40] "F_PlotConfElli.R"    "F_PlotDS.R"          "F_PlotHistFit.R"    
#[43] "F_PlotHistResid.R"   "F_PlotMetRes.R"      "F_PlotRL.R"         
#[46] "F_ProcessDL.R"       "F_proportions.R"     "F_q.PN.R"           
#[49] "F_QQW.R"             "F_Quantile.R"        "F_R2o.R"            
#[52] "F_reldist.R"         "F_RL.age.spline.R"   "F_Round.R"          
#[55] "F_RowToMatrix.R"     "F_SmoothHist1_V5.R"  "F_SmoothHist2.R"    
#[58] "F_Splineverlauf4.R"  "F_TIQQ.R"            "F_tmc.master.R"     
#[61] "F_tmc0_V3.R"         "F_tmc1.R"            "F_TransformKDE.R"   
#[64] "F_TukeyWW.R"         "F_Untie.R"           "F_x.mean.PN.R"      
#[67] "F_x.mode.PN.R"       "F_x.pdf.PN.R"        "F_x.Var.PN.R"       
#[70] "F_xx.pdf.PN.R"      


# ****************************************************************************
#  F_Bhatta_V2.R
#
#  Bhattacharya-Methode nach CG Bhattacharya, Biometrics 23 (1967) 115-135
#
#  (c) wwosniok@math.uni-bremen.de
#
#  To do
#  - Warum die getroffene Auswahl unter den 3 Lösungen für s - siehe @@@@ 
#    unten?
#  - Bhatta-Transformation glätten, siehe Oosterhuis 1990
#    Könnte leicht geschehen durch Verwendung von x.hist.smo

#  11.10.2020 xhist changed to x.hist in accordance with TMC use
#             x.hist must be provided, condition h not < round.unit must be 
#             ensured during histogram generation
#  10.10.2020 Histogram generation as in tmc
#  06.10.2020 Input revised
#  28.11.2018 Start
# 
# ============================================================================

Bhatta <- function(x.hist, x.kde, data.text, xlabel, subtitle, RL1.p,RL2.p, 
                   figD, figE,
                   kdecol, bhacol, bordercol1, histcol1)
{ #
  #  INPUT
  #
  # x.hist     data as equidistant histogram, bin width >= round.unit 
  #            in the original raw data
  # x.kde      kernel density estimate, for plotting only, may  be NA
  # data.text  description of the data
  # RL1.p      level for RL1
  # RL2.p      level for RL2
  # figD       figure for data histogram + kde + fitted distributions 
  # figE       figure for Bhattacharya plot

  #  OUTPUT
  #
  # =========================================================================

  #  Extract necessary information
  xcounts <- x.hist$counts
  xbreaks <- x.hist$breaks
  xmids   <- x.hist$mids
  xcounts.n <- length(xcounts)
  xbreaks.n <- length(xbreaks)
  x.n       <- sum(xcounts)
  h         <- xmids[2] - xmids[1]

  xsupp <- seq(x.kde$x[1], tail(x.kde$x, 1), length.out=101) 

  # ..........................................................................
  # Plot histogram

  if (!is.na(figD))
  {
    dev.set(figD)
    plot(x.hist,freq=FALSE,
         main="Histogram for Bhattacharya method",
         xlab=xlabel, 
         border=bordercol1, col=histcol1, sub=subtitle,cex.sub=0.7)
    if (length(x.kde) > 0) { lines(x.kde,col=kdecol) }
  }

  # ..........................................................................
  #  Calculate delta.log.y, the Bhattacharya transformation. 
  #  To avoid problems with empty cells: replace 0 by fastnull

  fastnull <- 0.5/x.n
  xcounts.1 <- xcounts
  xcounts.1[xcounts==0] <- fastnull

  delta.log.y   <- log(xcounts.1[2:xcounts.n]) - 
                   log(xcounts.1[1:(xcounts.n-1)])

  # ..........................................................................
  # Plot delta.log.y

  if (!is.na(figE))
  {  
    dev.set(figE)
    plot(xmids[1:(xcounts.n-1)], delta.log.y, type="o", col=bhacol,
         main="Bhattacharya transformation",
         xlab=xlabel, ylab="log(f_i+1 / f_i)",
         sub=subtitle,cex.sub=0.7)  
    abline(h=0)
  }

  # ..........................................................................
  #  Select subsets (truncation areas, at least 2 distributions are suspected)
  #  Automatic selection
  #  Find minima / maxima
  delta.log.y.n <- length(delta.log.y)
  M.names <- c("xc","xcounts","yl","yc","yr","is.min","is.max","is.ext")
  M <- data.frame(matrix(NA,nrow=delta.log.y.n,ncol=length(M.names)))
  colnames(M) <- M.names

  M[ ,"xc"]      <- xmids[1:delta.log.y.n]
  M[ ,"xcounts"] <- xcounts[1:delta.log.y.n]
  M[2:delta.log.y.n, "yl"] <- delta.log.y[1:(delta.log.y.n-1)]
  M[ ,"yc"] <- delta.log.y
  M[1:(delta.log.y.n-1), "yr"] <- delta.log.y[2:(delta.log.y.n)]

  innen <- 2:(delta.log.y.n-1)

  M[1,"is.min"]     <- (M[1,"yc"] < M[1,"yr"])
  M[innen,"is.min"] <- (M[innen,"yl"] > M[innen,"yc"]) & 
                       (M[innen,"yc"] < M[innen,"yr"]) 
  M[delta.log.y.n,"is.min"] <- (M[delta.log.y.n,"yl"] > M[delta.log.y.n,"yc"])

  M[1,"is.max"]     <- (M[1,"yc"] > M[1,"yr"])
  M[innen,"is.max"] <- (M[innen,"yl"] < M[innen,"yc"]) & 
                       (M[innen,"yc"] > M[innen,"yr"])
  M[delta.log.y.n,"is.max"] <- (M[delta.log.y.n,"yl"] < M[delta.log.y.n,"yc"])

  # is.min and is.max may be NA due to zeroes in the data. Set to FALSE
  M[is.na(M[ ,"is.min"]), "is.min"] <- FALSE
  M[is.na(M[ ,"is.max"]), "is.max"] <- FALSE

  #  Identical values of delta.log.y are possible. Check for regular min/max
  #  alternating in the extrema sequence.
  is.ext <- rep(NA,times=delta.log.y.n)
  M[M[ ,"is.min"],"is.ext"]  <-  "min"
  M[M[ ,"is.max"],"is.ext"]  <-  "max"

  #cat("\n[Bhatta] Matrix M\n")
  #print(M)

  #  @@@@ hier fehlt noch was (min/max-Reihenfolge überpüfen)

  #  Interesting pieces start with max
  #  Indices are indices in xmids
  max.idx <- which(M[ ,"is.max"])
  max.idx.n <- length(max.idx)

  min.idx <- which(M[ ,"is.min"])

  res.names <- c("istart","iend","nobs","adjrsq","mue.hut","sig.hut",
                 "RL1","RL2","RL1.1","RL2.1","RL1.2","RL2.2","RL1.3","RL2.3",
                 "p.hut","n.hut","prop.hut")
  res <- matrix(NA,nrow=length(max.idx),ncol=length(res.names))
  colnames(res) <- res.names

  i <- 0
  for (istart in max.idx)
  { i <- i + 1
    
    #cat("\n Starting search for component ", i, "------------------------ \n")

    res[i,"istart"] <- istart
    #  Find next minimum
    #iend <- min(min.idx[min.idx>istart])  #### @@@ iend könnte inf sein!
    
    # Avoid warnings:
    min.idx.up <- min.idx>istart
    if (sum(min.idx.up)== 0 )
    { # no more minima
      iend <- -Inf
    } else
    { # more minima
      iend <- min(min.idx[min.idx.up])
    }

    res[i,"iend"] <- iend

    #cat("\n[Bhatta] istart, iend ",istart, iend,"\n")

    if (is.finite(iend) & ((iend-istart)>1) )
    {
      if (!is.na(figE))
      {
        dev.set(figE)
        abline(v=M[istart,"xc"],col="red")
        abline(v=M[iend,"xc"],col="green4")
      }

      # subset: indices of points in the truncation interval
      subset <- istart:iend
      x.sub  <- M[subset,"xc"]
      y.sub  <- M[subset,"yc"]

      delta.log.y.lm <- lm(y.sub ~ x.sub)
      delta.log.y.lm.sum <- summary(delta.log.y.lm)   
      #print(delta.log.y.lm.sum)

      if (!is.na(figE))
      {
        dev.set(figE)
        abline(delta.log.y.lm$coefficients,col="green3",lty=2)
      } 

      beta0 <- delta.log.y.lm$coefficients[1]
      beta1 <- delta.log.y.lm$coefficients[2]
      res[i,"adjrsq"] <- unlist(delta.log.y.lm.sum["adj.r.squared"])

      lambda    <- -beta0 / beta1   #  this is not the PND lambda!
      mue.hut   <-  lambda + h/2
      s.hut.1   <-  sqrt(-h/beta1) - h^2/12   #  Bhattacharya, (2)
      s.hut.2   <-  sqrt(-h/beta1)            #  Sparre & Venema

      #  Lösung der quadratischen Gleichung für sigma^2
      discrim <- 1/(4*beta1^2) + h/(12*beta1)

      sigma2.1 <- NA
      sigma2.2 <- NA
      s.hut.3  <- 0
      if (discrim >= 0)
      {
        sigma2.1 <- -h/(2*beta1) + h*sqrt(discrim)
        sigma2.2 <- -h/(2*beta1) - h*sqrt(discrim)
        s.hut.3  <- sqrt(sigma2.1)
      }

      res[i,"mue.hut"] <- mue.hut
      res[i,"sig.hut"] <- s.hut.2  # @@@@@ warum diese Lösung als Resultat?

      #  Estimate the size of this subpopulation
      #  Probability in this truncation interval 
      p.tr <- pnorm(tail(x.sub,1),mean=mue.hut,sd=s.hut.2) -
              pnorm(x.sub[1],mean=mue.hut,sd=s.hut.2)
      res[i,"p.hut"] <- p.tr

      #  Counts in the truncation interval
      xcounts.tr.n  <- sum(M[subset,"xcounts"])
      res[i,"nobs"] <- xcounts.tr.n

      #  Bezeichnung xc ist nicht gut
      xc.hut <- xcounts.tr.n/p.tr
      res[i,"n.hut"] <- xc.hut
      res[i,"prop.hut"]  <- xc.hut/x.n
      res[i,"RL1"]   <- qnorm(RL1.p,mean=res[i,"mue.hut"], sd=res[i,"sig.hut"] )
      res[i,"RL2"]   <- qnorm(RL2.p,mean=res[i,"mue.hut"], sd=res[i,"sig.hut"] )
      res[i,"RL1.1"] <- qnorm(RL1.p,mean=res[i,"mue.hut"], sd=s.hut.1 )
      res[i,"RL2.1"] <- qnorm(RL2.p,mean=res[i,"mue.hut"], sd=s.hut.1 )
      res[i,"RL1.2"] <- qnorm(RL1.p,mean=res[i,"mue.hut"], sd=s.hut.2 )
      res[i,"RL2.2"] <- qnorm(RL2.p,mean=res[i,"mue.hut"], sd=s.hut.2 )
      res[i,"RL1.3"] <- qnorm(RL1.p,mean=res[i,"mue.hut"], sd=s.hut.3 )
      res[i,"RL2.3"] <- qnorm(RL2.p,mean=res[i,"mue.hut"], sd=s.hut.3 )

      #cat("\nBHATTACHARYA method",
      #"\nVersion 1   Mean  ", mue.hut,"  sd  ",s.hut.1,
      #"\nVersion 2   Mean  ", mue.hut,"  sd  ",s.hut.2,
      #"\nVersion 3   Mean  ", mue.hut,"  sd  ",s.hut.3,
      #"\nEstimated subpopulation size            ",round(xc.hut),
      #"\nCorresponding proportion of total data  ",xc.hut/x.n,
      #"\n")

      if (!is.na(figD))
      {
        dev.set(figD)
        lines(xsupp,(xc.hut/x.n)*dnorm(xsupp,mean=mue.hut,sd=s.hut.1),col="blue")
        lines(xsupp,(xc.hut/x.n)*dnorm(xsupp,mean=mue.hut,sd=s.hut.2),col="red",
            lwd=2)
        lines(xsupp,(xc.hut/x.n)*dnorm(xsupp,mean=mue.hut,sd=s.hut.3),
            col="orange")

        legend("topleft",c("Bhattacharya (2)","Sparre & Venema","Quad. equation"),
             col=c("blue","red","orange"),lwd=c(1,2,1),cex=0.7)
      } 

    }   # is.finite ...

    #cat("\n Search for component ", i, "  finished ---------------------- \n")

  }     # for (istart ...

  # Make sure that res is a matrix
  if (!is.matrix(res))
  { res <- matrix(res,nrow=1) 
    #  no sorting need
  } else
  { 
    res <- res[!is.na(res[ ,"n.hut"]), ]
    if (!is.matrix(res))
    { res <- matrix(res,nrow=1) 
      #  no sorting need
    } else
    {   
      res <- res[order(-res[ ,"n.hut"]), ]
    }
  }
  colnames(res) <- res.names
  rownames(res) <- paste("comp", 
                         formatC(1:nrow(res), format="f", digits=0, width=2, 
                                 flag="0"), sep="")

  return(list(x.hist.bha=x.hist,subpop=res))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_BoxCox.R
#
#  Box-Cox-Transformation mit Parameter lambda
#
#  (c) wwosniok@math.uni-bremen.de
#
#  CHANGE HISTORY
#
#  18.07.2018 WW   xgt0   <- (!is.na(x)) & (     x > fastnull)
#                  geändert in 
#                  xgt0   <- (!is.na(x)) & (     x >= fastnull)
#
# --------------------------------------------------------------------------
# Eingabe
# x         zu transformierender Vektor ("Originalskala" = beobachtet)
# lambda    Box-Cox-Parameter,
#           lambda = 0: entspricht log(x)
#           lambda = 1: entspricht x-1
# fastnull  Absicherung gegen kleine Abbruchfehler

# Ausgabe
# Tx        transformierter Vektor 
#                           x < 0  Tx = NA   (obwohl sinnvolle Werte für 
#                                             ganzzahliges lambda möglich)
#           lambda <= 0 und x = 0  Tx = -Inf
#           lambda =  0 und x > 0  Tx = log(x)
#           lambda >  0 und x = 0  Tx = 0
#           lambda <> 0 und x > 0  Tx = (x^lambda-1)/lambda
# ============================================================================
             
BoxCox <- function(x,lambda,fastnull=1.e-10)
{ 
  xgt0   <- (!is.na(x)) & (     x >= fastnull)

  xeq0   <- (!is.na(x)) & (abs(x) < fastnull)

  Tx     <- rep(NA,times=length(x))

  if (lambda < -fastnull)
  { Tx[xeq0] <- -Inf
    Tx[xgt0] <- (x[xgt0]^lambda-1)/lambda 
  }
  if (abs(lambda) <= fastnull) 
  { 
    Tx[xeq0] <- -Inf
    Tx[xgt0] <- log(x[xgt0])
  }
  if (lambda > fastnull)
  { 
    # Tx[xeq0] <- (x[xeq0]^lambda-1)/lambda 
    # Tx[xgt0] <- (x[xgt0]^lambda-1)/lambda 
    Tx         <- (x^lambda-1)/lambda 
  }

  return(Tx)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test
#BCx <- BoxCox(10.5,1.29247e-25)
#BoxCox(0,1.001) 

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# ****************************************************************************
#  F_BoxCoxInv.R

#  Inverse Box-Cox-Transformation withrameter lambda
#
#  (c) wwosniok@math.uni-bremen.de
#
#  CHANGE HISTORY

#  27.06.2016 Start
# --------------------------------------------------------------------------

# INPUT
# Tx       vector to transform 
# lambda   Box Cox parameter
#          lambda = 0: same as exp(x)
#          lambda = 1: produces x+1
# fastnull value replacing zero 

# Ausgabe
# x        transformed vector, x = (1 + lambda*Tx)^(1/lambda)
# ============================================================================

BoxCoxInv <- function(Tx,lambda,fastnull=1.e-10)
{ 
  x <- rep(NA,times=length(Tx))
 
  if (lambda < -fastnull)
  { x <- (Tx*lambda+1)^(1/lambda) }

  if (lambda > fastnull)
  { x <- (Tx*lambda+1)^(1/lambda) }
  
  if (abs(lambda) <= fastnull)
  { xneInf       <- is.finite(Tx)
    x[xneInf]    <- exp(Tx[xneInf])
    if (any(!xneInf)) 
    { x[!xneInf] <- 0  }
  }
  
  return(x)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test
#x <- BoxCoxInv(5, 1)   # 6
#x <- BoxCoxInv(5, 0)   # 148.4132 = exp(5)
#x <- 5
#y <- BoxCox(x, 0)      # 1.609438 = log(5)
#xx <- BoxCoxInv(y, 0)  # 5
#x <- 5
#y <- BoxCox(x, 1)      # 4
#xx <- BoxCoxInv(y, 1)  # 5

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_CalcChi2.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculate chi squared from multinomial data, given the parameters
#  of a lognormal distribution 
 
#  #  CHANGE HISTORY
#  28.08.2020 Start
#
# ==============================================================================

CalcChi2 <- function(breaks, counts, mue, sig)
{
  #  INPUT 
  #  breaks   
  #  counts
  #  mue
  #  sig

  #  OUTPUT 
  #  
  # ---------------------------------------------------------------------------

  x.n        <- sum(counts)
  counts.n   <- length(counts)
  prop.emp   <- counts/x.n
  cdf        <- plnorm(breaks, meanlog=mue, sdlog=sig)
  prop.the   <- cdf[2:(counts.n+1)]-cdf[1:counts.n]
  counts.the <- x.n * prop.the
  tab <- data.frame(x.lo=breaks[1:counts.n],
                    x.hi=breaks[2:(counts.n+1)],
                    count=counts,
                    prop.emp=prop.emp,
                    prop.the=prop.the,
                    chi2.bin=((counts-counts.the)^2)/counts.the )
  return(tab)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#  F_CalcHistLoop.R
#
#  - Checks an equidistant set of breaks for required properties
#  - generates an equidistant histogram 
#  - generates a collapsed version of the (smoothed) equidistant histogram 
#    which has all cell frequencies >= n.per.bin.min
 
#  04.01.2021  Smoothing of histogram done before collapsing
#  28.12.2020  Start
# ==============================================================================

CalcHistoLoop <- function(x, x.kde, x.min, x.max, eqd.breaks0, h, n.per.bin.min,
                          round.unit, x.kde.mode, smooth.hist1, smooth.hist2, 
                          xlabel, subtitle, 
                          histcol2, bordercol2, histcol3, bordercol3, kdecol)
{
  eqd.breaks <- seq(eqd.breaks0[1], tail(eqd.breaks0, 1), by=h) 

  #  Make sure that the maximum is included
  if (tail(eqd.breaks, 1) < x.max)
  { eqd.breaks <- c(eqd.breaks, tail(eqd.breaks, 1) + h) }
  
  # The smallest break must be > 0 and < min(data)
  # Zeroes in the data were set to 0.1*round.unit by the reading procedure
  eqd.breaks[1] <- max(0.05*round.unit, eqd.breaks[1])    

  #  Make sure again that the maximum is included
  if (tail(eqd.breaks, 1) < x.max)
  { eqd.breaks <- c(eqd.breaks, tail(eqd.breaks, 1) + round.unit) }

  #  Generate an equidistant histogram with the (modified) breaks
  x.hist.eqd <- hist(x, breaks=eqd.breaks, right=FALSE, plot=FALSE)

  # -------------------------------------------------------------------------  #  Collapse histogram, if necessary to achieve a minimum count per bin 
  #  Smooth the equidistant histogram, if requested
  if ((!smooth.hist1) & (!smooth.hist2))
  { # No smoothing
    x.hist.eqd.smo <- x.hist.eqd
  }

  if (( smooth.hist1) & (!smooth.hist2))
  { # Smoothing by moving average
    x.hist.eqd.smo <- SmoothHist1(x.hist.eqd, figA=NA )
  } 

  if ((!smooth.hist1) & ( smooth.hist2))
  { # Smoothing by kde fit
    x.hist.eqd.smo <- SmoothHist2(x.kde, x.hist.eqd, NA, 
                                  histcol2, bordercol2, 
                                  histcol3, bordercol3, kdecol, subtitle)
  } 

  # -------------------------------------------------------------------------  #  Collapse histogram, if necessary to achieve a minimum count per bin 
  #  Collapse the (possibly) smoothed histogram
  #  Preparation: Find the index of the interval containing the mode
  x.hist.eqd.smo.mode.idx <- FindModeIndex(x.hist.eqd.smo$breaks, 
                                           x.hist.eqd.smo$counts, 
                                           x.kde.mode)

  x.hist <- CollapseHist(x, x.hist.eqd.smo, x.hist.eqd.smo.mode.idx, 
                         n.per.bin.min, 
                         xlabel, subtitle, histcol2, bordercol2, NA)

  return(H=list(x.hist.eqd=x.hist.eqd, 
                x.hist.eqd.smo=x.hist.eqd.smo,
                x.hist=x.hist))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


# ****************************************************************************
#  F_CalcParDistance.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculates distance between observed truncated parameters
#  and estimated truncated parameters resulting
#  from the presumed parent distribution, given the truncation limits
 
#  #  CHANGE HISTORY
#  09.09.2020 Operation restricted to normally distributed data
#  19.08.2020 Start
#
# ==============================================================================

CalcParDistance <- function(theta.parent, theta.trunc.obs, a, b)
{
  #  INPUT
  #  theta.parent      Vector with components 
  #          mue.parent   (presumed) mean of the parent normal distribution
  #          sig.parent   (presumed) sd of the parent normal distribution
  #  theta.trunc.obs   Vector with components 
  #          mue.trunc    observed mean of the truncated distribution 
  #          sig.trunc    observed sd of the truncated distribution 
  #  a             left truncation limit
  #  b             right truncation limit

  #  OUTPUT 
  #  SS     squared distance between 
  #         (mue.trunc.presumed, sig.trunc.presumed) and (mue, sig)
  # ---------------------------------------------------------------------------

  #  Auxiliary terms
  mue.parent <- theta.parent[1] 
  sig.parent <- theta.parent[2] 
  mue.trunc.obs  <- theta.trunc.obs[1]
  sig.trunc.obs  <- theta.trunc.obs[2]

  alpha <- (a - mue.parent) / sig.parent
  beta  <- (b - mue.parent) / sig.parent

  A  <-  dnorm(alpha, mean=0, sd=1)
  B  <-  dnorm(beta, mean=0, sd=1)
  AA <-  pnorm(alpha, mean=0, sd=1)
  BB <-  pnorm(beta, mean=0, sd=1)
  C  <-  (B - A) / (BB - AA)

  #  Calculate the mean of the truncated distribution, given theta.parent
  mue.trunc.est <- mue.parent - sig.parent * C

  #  Calculate the sd of the truncated distribution, given theta.parent
  sig.trunc.est <- sig.parent * sqrt((1 - (beta*B - alpha*A) / (BB - AA) - C^2 ))

  R <- matrix(c(a, b,
                mue.parent, sig.parent,
                mue.trunc.est, sig.trunc.est,
                mue.trunc.obs, sig.trunc.obs), byrow=TRUE, ncol=2)
  d  <- c(mue.trunc.obs - mue.trunc.est, sig.trunc.obs-sig.trunc.est)
  #cat("\n [CarlcParDist] R\n")
  #print(R)              
  #cat("\n [CarlcParDist] d\n")
  #print(d)

  SS <- sqrt(d %*% d)

  return(SS)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#theta.parent <- c(140.011879, 2.671882)
#theta.trunc  <- c(139.933013, 1.326838)
#a <- 138
#b <- 142

#HW <- CalcParDistance(theta.parent, theta.trunc, a, b)
#print(HW)
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_CalcPrev_V2.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculate the sizes of xl, xc, xr using the truncation limits and the 
#  estimated parameters of xc. Values inside the truncation interval 
#  are assumed to be completely from xc.

#  Calculation on the untransformed basis using cdf.PN(x, lambda, mue, sigma)
 
#  #  CHANGE HISTORY
#
# 06.12.2020 Start
# ============================================================================

CalcPrev <- function(x.n, x.tr.n, x.lt.tr.n, x.ge.tr.n, 
                     x.tr.lo, x.tr.hi, lambda, mue, sigma)
{
  #  Central component of the mixture

  xc.tr.p <- (cdf.PN(x.tr.hi, lambda, mue, sigma) - 
              cdf.PN(x.tr.lo, lambda, mue, sigma))

  # May get zero if provided parameters are too wrong. Prevent overflow.
  xc.tr.p   <- max(1.e-6, xc.tr.p) 

  c.n    <- x.tr.n/xc.tr.p
  prev.c <- c.n/x.n

  #  Size of the left component of the mixture. This is assumed to exist only 
  #  < x.tr.lo, because the truncation interval is assumed to contain only
  #  values from xc.
 
  l.n <- x.lt.tr.n - c.n * cdf.PN(x.tr.lo, lambda, mue, sigma)
  prev.l <- l.n/x.n

  #  Size of the rightt component of the mixture. This is assumed to exist only 
  #  >= x.tr.hi, because the truncation interval is assumed to contain only
  #  values from xc.
  r.n <- x.ge.tr.n - c.n * (1-cdf.PN(x.tr.hi, lambda, mue, sigma))
  prev.r <- r.n/x.n

  neg.prev.sum <- sum(min(prev.l, 0), min(prev.r, 0)) 

  l.n <- unname(l.n)
  c.n <- unname(c.n)
  r.n <- unname(r.n)
  prev.l <- unname(prev.l)
  prev.c <- unname(prev.c)
  prev.r <- unname(prev.r)
  neg.prev.sum <- unname(neg.prev.sum)

  return(c(l.n=l.n, prev.l=prev.l, c.n=c.n, prev.c=prev.c,  
           r.n=r.n, prev.r=prev.r, neg.prev.sum=neg.prev.sum))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
#
#  Test for version 2. Data used is FileNo 10021, M only, x.n = 42400

#x.n <- 42400
#lambda.c <- 0
#mue.c    <- 2.928666
#sigma.c  <- 0.3201523

#x.kde.mode  <- 17.04995
#x.lt.kde.mode.n <- 16874    # sum(x < x.kde.mode) 
#x.ge.kde.mode.n <- 25526    # sum(x >= x.kde.mode) 

#x.tmc.mode <- exp(mue.c-sigma.c^2) # 16.88066
#x.lt.tmc.mode.n <- 13971    # sum(x < x.tmc.mode) 
#x.ge.tmc.mode.n <- 28429    # sum(x >= x.tmc.mode) 

#x.tr.lo <-  7.5
#x.tr.hi <- 27.5

#x.tr.n <- 35615     # sum((x.tr.lo <= x) & (x < x.tr.hi))
#x.lt.tr.n <- 85     # sum(x < x.tr.lo)
#x.ge.tr.n <- 6700   # sum(x >= x.tr.hi)

#res <- CalcPrev(x.n, x.tr.n, x.lt.tr.n, x.ge.tr.n, 
#                x.tr.lo, x.tr.hi, lambda.c, mue.c, sigma.c)
#res
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_CalcR2InSubInt.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculate r^2 for subintervals of the QQ plot with sufficiently many bins
#  and cases, and also contining the mode (if requested)
#  Subintervals are defined by the breaks of the collapsed histogram.
 
#  #  CHANGE HISTORY
#  02.12.2020 x.tr.n added to output
#  24.11.2020 Estimation of xc.n, xl.n, xr.n removed
#  16.11.2020 data.frames constructed with stringsAsFactors=FALSE
#  26.09.2020 Call changed
#  31.08.2020 More criteria for subinterval selection added
#  29.08.2020 Start
#
# ==============================================================================

CalcR2InSubInt <- function(x.hist, y.hist, y.poly.qq, lambda, 
                           bins.n.min.act,  x.tr.prop.min, x.tr.prop.max, 
                           y.hist.mode.idx,
                           df.est,df.con,
                           x.Q1, x.Q2, RL1.p, RL2.p,
                           l.fact, p.fact, r.fact, w.fact,
                           print.log.message)
{
  #  INPUT 
  #  y.hist
  #  y.poly.qq
  #  xlabel   label for data
  #  subtitle
  #  bins.n.min.act
  #  figA

  #  OUTPUT 
  #  
  # ===========================================================================

  # --------------------------------------------------------------------------
  #if (print.log.message) { cat("%%%   CalcR2InSubInt   Start\n") }
  # --------------------------------------------------------------------------

  y0.tab <- hist.table(y.hist$breaks, y.hist$counts) 
  
  #  Due to possibly dealing with a reduced dataset, some lines in y.tab might 
  #  be empty 

  bins.n <- nrow(y0.tab)
  y.tab <- matrix(NA, nrow=bins.n, ncol=ncol(y0.tab))
  colnames(y.tab) <- colnames(y0.tab)

  y.tab[1, "x.lo"] <- y0.tab[1, "x.lo"]

  j <- 0
  for (i in 1:bins.n)
  {  
    if (y0.tab[i, "count"] > 0 )
    {
      j <- j + 1

      y.tab[j, "x.hi"]  <-  y0.tab[i, "x.hi"]
      y.tab[j, "count"] <-  y0.tab[i, "count"]
      y.tab[j, "prop"]  <-  y0.tab[i, "prop"]

      if (j > 1) { y.tab[j, "x.lo"] <- y.tab[j-1, "x.hi"] } 
    }
  }  

  #  Remove empty lines
  y.tab  <- y.tab[!is.na(y.tab[ ,"count"]), ]
  bins.n <- nrow(y.tab)
  
  #  Calculate partial r2 subintervals having
  #  - at least bins.n.min.act bins
  #  - at least a proportion of x.tr.prop.min values
  #  - at most  a proportion of x.tr.prop.max values
  #  - contains the mode
  #  Note: y.hist may be a reduced version of the true histogram.
  #  Proportions like x.tr.prop.min are used here, but applied
  #  possibly to the reduced histogram. Absolute sizes of xl, xc, xr
  #  are not computed here. Must be done in the calling programme.

  bins.n <- nrow(y.tab)
  x.n    <- sum(y.hist$counts)

  Pr2.names <- c("bin.start", "bin.end", "bins.n", "prop",  
                 "x.tr.lo", "x.tr.hi", "y.tr.lo", "y.tr.hi", 
                 "mue", "sig", "r2", "opt.crit")
  Pr2 <- matrix(NA, nrow=bins.n * (bins.n-1) / 2, 
               ncol=length(Pr2.names))
  colnames(Pr2) <- Pr2.names
  Pr2 <-  data.frame(Pr2, stringsAsFactors=FALSE)

  Pr2.n <- 0
  for (bin.start in 1:(bins.n-bins.n.min.act+1) )
  { y.tr.lo <- y.hist$breaks[bin.start]

    for (bin.end in (bin.start+bins.n.min.act-1):bins.n)
    { 
      # Check the subinterval for having required properties
      bins.act <- bin.end - bin.start + 1
      bins.n.ok <- (bins.act >= bins.n.min.act)
      ok <- (y.tab[bin.start, "x.lo"] <= y.poly.qq$y) &
            (y.poly.qq$y < y.tab[bin.end, "x.hi"] )

      x.tr.n <- sum(ok)
      prop   <- x.tr.n/x.n

      n.ok   <- (x.tr.prop.min <= prop) & (prop <= x.tr.prop.max)
       
      mode.ok <- TRUE

      #  TI must contain the mode 

      if (!is.na(y.hist.mode.idx))
      { mode.ok <- ( (bin.start <= y.hist.mode.idx) &
                     (y.hist.mode.idx <= bin.end )   )
      }

      # Calculations for this interval only if all conditions hold

      if (bins.n.ok & n.ok & mode.ok)
      {
        Pr2.n <- Pr2.n + 1

        y.tr.hi <- y.hist$breaks[bin.end+1]
        Pr2[Pr2.n, "bin.start"] <- bin.start 
        Pr2[Pr2.n, "bin.end"]   <- bin.end
        Pr2[Pr2.n, "y.tr.lo"]   <- y.tr.lo 
        Pr2[Pr2.n, "y.tr.hi"]   <- y.tr.hi
        Pr2[Pr2.n, "x.tr.n"]    <- x.tr.n   # may refer to a reduced histogram
        Pr2[Pr2.n, "prop"]      <- prop     # may refer to a reduced histogram

        #  Calculate regression per accepted interval

        y.lm <- lm(y.poly.qq$y[ok] ~ y.poly.qq$x[ok])
        y.lm.sum <- summary(y.lm)
        Pr2[Pr2.n, "mue"]      <- y.lm$coefficients[1]
        Pr2[Pr2.n, "sig"]      <- y.lm$coefficients[2]
        Pr2[Pr2.n, "r2"]       <- y.lm.sum$adj.r.squared

        #  Calculate the optimality criterion that is used by tmc, tmu
        #  No, no more used
        #if (bin.start == 1)

        #{ x.lt.tr.n <- 0 } else
        #{ 
        #  idx <- 1:(bin.start-1)
        #  x.lt.tr.n <- sum(x.hist$counts[idx])
        #}

        #if (bin.end == bins.n)
        #{ x.ge.tr.n <- 0 } else
        #{ 
        #  idx <- (bin.end+1):bins.n
        #  x.ge.tr.n <- sum(x.hist$counts[idx])
        #}

        ##  C2T <- chi2trunc(x.hist, 
        #opt.crit <- chi2trunc(x.hist, 
        #                 BoxCoxInv(y.tr.lo, lambda), 
        #                 BoxCoxInv(y.tr.hi, lambda), 
        #                 x.lt.tr.n, x.ge.tr.n,
        #                 x.Q1, x.Q2, RL1.p, RL2.p,
        #                 lambda, Pr2[Pr2.n, "mue"], Pr2[Pr2.n, "sig"],
        #                 df.est,df.con,
        #                 l.fact, p.fact, r.fact, w.fact, 
        #                 opt.crit.only=TRUE)
        # Pr2[Pr2.n, "opt.crit"]  <- C2T$res["opt.crit"]      
        #Pr2[Pr2.n, "opt.crit"]  <-  opt.crit      
      }
    }
  }
  Pr2 <- Pr2[1:Pr2.n, ]

  #  Check for reasonable regression lines
  #  @ tbd

  #  Sort by descending r2
  Pr2 <- Pr2[order(-Pr2[ ,"r2"]), ]

  #  Sort by ascending opt.crit
  #Pr2 <- Pr2[order(Pr2[ ,"opt.crit"]), ]

  Pr2[ , "bins.n"] <- Pr2[ , "bin.end"] - Pr2[ , "bin.start"] + 1
  Pr2 <- data.frame(Pr2, rank.r2=1:nrow(Pr2), stringsAsFactors=FALSE)
  Pr2 <- data.frame(Pr2, y.RL1=qnorm(RL1.p, mean=Pr2[ ,"mue"], sd=Pr2[ ,"sig"]),
                    stringsAsFactors=FALSE)
  Pr2 <- data.frame(Pr2, y.RL2=qnorm(RL2.p, mean=Pr2[ ,"mue"], sd=Pr2[ ,"sig"]),
                    stringsAsFactors=FALSE)

  # --------------------------------------------------------------------------
  #if (print.log.message) { cat("%%%   CalcR2InSubInt   End\n") }
  # --------------------------------------------------------------------------

  return(list(tab=y.tab, Pr2=Pr2))
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#
# F_cdf.PN.R

# ----------------------------------------------------------------------       
cdf.PN <- function(x,lambda,mue,sigma,fastnull=1.e-10)
{ 
  #  Verteilungsfunktion der Power-Normalverteilung. 
  #  mue, sigma gelten auf der Y-Skala
  #  siehe Freeman, S. 767

  #  23Jul2018 WW K = 0 abgefangen

  cdf  <- rep(0, times=length(x))
  xgt0 <- x > fastnull
  T <- 1/(lambda*sigma) + mue/sigma
  Z <- (BoxCox(x[xgt0],lambda)-mue)/sigma
 
  if (lambda > fastnull)
  { 
    K <- pnorm(T,mean=0,sd=1)

    #  Overflow protection
    if (K > fastnull)
    {
      cdf[xgt0]  <- (1/K) * (pnorm(Z) - pnorm(-T))
    } else
    {
      cdf[xgt0]  <- 1
    }
  }
  if (abs(lambda) <= fastnull)
  { # LNV
    
    cdf <- plnorm(x,meanlog=mue,sdlog=sigma)
  }
  if (lambda < -fastnull)
  { 
    K <- pnorm(-T,mean=0,sd=1)

    #  Overflow protection
    if (K > fastnull)
    {
      cdf[xgt0]  <- (1/K) * pnorm(Z)
    } else
    {
      cdf[xgt0] <- 0 
    }
  }
  return(cdf)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#lambda <- 0
#lambda <- -0.09989172
#mue    <- 11.25949
#sigma  <-  0.01440878
#x.P50  <- q.PN(0.50, lambda, mue, sigma)
#x      <-  seq(1,20)
#x      <-  seq(77600, 77700, by=5)
#x.pdf  <- pdf.PN(x, lambda, mue, sigma)
#x.cdf  <- cdf.PN(x, lambda, mue, sigma)
#cat("\n x, x.pdf, x.cdf")
#print(cbind(x, x.pdf, x.cdf))

#cdf.PN(135.8146, 0.9, 140.0014, 2.646370)  
#cdf.PN(135.8146, 0.99, 140.0014, 2.646370)  
#cdf.PN(135.8146, 0.999, 140.0014, 2.646370)  
#cdf.PN(135.8146, 0.9999, 140.0014, 2.646370)  
#cdf.PN(135.8146, 1.00, 140.0014, 2.646370)  
#cdf.PN(135.8146, 1.0001, 140.0014, 2.646370)  
#cdf.PN(135.8146, 1.001, 140.0014, 2.646370)  
#cdf.PN(135.8146, 1.01, 140.0014, 2.646370)  
#cdf.PN(135.8146, 1.1, 140.0014, 2.646370)  

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_Ceiling.R

# 08.01.2020

# =============================================================================
Ceiling <- function(x,unit)
{ # Rundet nach oben auf die nächsthöhere unit, sofern x nicht 
  # ganzzahlig Vielfaches der unit ist
  vielfaches <- trunc(x/unit)
  x.ceiling <- (vielfaches+(abs(x-vielfaches*unit) > 1e-15)) * unit
  return(x.ceiling)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_CheckRounding.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Check data for inconsistent rounding. 
#  Do all possible values of the last digit after rounding have similar
#  frequencies? Produce a frequency table of the last digit. The rounding unit
#  is determined from the smallest difference between different values.
#  Rounding to powers of ten is assumed (not e.g. to units of 0.2)
 
#  #  CHANGE HISTORY
#  15.10.2020 Problems due to small numerical errors removed
#  08.06.2020 Start
# ============================================================================

CheckRounding <- function(x, nearlyzero=1.e-10)
{
  #  INPUT 
  #  x                input data (vector)
  #  nearlyzero       absolute smaller values are treated as zero

  #  OUTPUT 
  #  last.dig.table   frequency table of the last digit
  # ==========================================================================

  #  Find rounding unit

  x             <- sort(x)
  x.diff        <- diff(x, lag=1)
  x.diff        <- x.diff[x.diff > nearlyzero]
  x.diff.min    <- min(x.diff)
  x.diff.min.lg <- log10(x.diff.min)
  round.uni.data    <- 10^Round(x.diff.min.lg,1)

  #  Determine last digit 

  last.dig0      <- x - (round.uni.data*10) * 
                    Floor((x+nearlyzero)/(round.uni.data*10), 1)

  #  last.dig values may be different due to representation problems
  #  (differences after the tenth position)
  #  Reduce to a reasonable number of digits
  last.dig1 <- signif(last.dig0, digits=3)
  last.dig  <- Round(last.dig0, round.uni.data)

  D <- data.frame(x, last.dig0, last.dig1, last.dig)

  D <- D[order(D[ ,"last.dig"]), ]

  last.dig.table.values <- (0:9) * round.uni.data 

  temp.exp              <- data.frame(val=last.dig.table.values,
                                      exp=rep(length(x)/10, times=10))

  last.dig.table.obs    <- table(last.dig)

  #  There may be no rounding in the data
  if (all(last.dig.table.obs == 1))
  { 
    cat("\n  +++ No rounding detected  +++ \n") 
    last.dig.table <- NA
  } else
  {
    val.obs               <- as.numeric(names(last.dig.table.obs))
    val.obs               <- Round(val.obs, round.uni.data)
    temp.obs              <- data.frame(val=val.obs,
                                        obs=c(last.dig.table.obs))

    #  Merge observed and expected counts
    last.dig.table        <- merge(temp.exp, temp.obs, all.x=TRUE)
    last.dig.table[is.na(last.dig.table[ ,"obs"]) ,"obs"] <- 0

    #  Calculate chi square
    last.dig.table <- data.frame(last.dig.table, 
                                 chi2=(last.dig.table[ ,"obs"] - 
                                       last.dig.table[ ,"exp"])^2 / 
                                       last.dig.table[ ,"exp"])

    colnames(last.dig.table) <- c("Value", "Expected count", "Observed count", 
                                  "chi2")
    last.dig.table <- last.dig.table[ ,
                                c("Value", "Observed count", "Expected count",
                                  "chi2")]

    return(last.dig.table)
  } 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test

#x <- c(21, 23, 25, 25, 28)
#x <- x / 100

#ldt <- CheckRounding(x) 
#ldt
#ldt <- CheckRounding(seq(79.75, 79.85, by=0.01)) 

#  roughly glucose
#x0   <- rlnorm(1000, meanlog=4.94, sdlog=0.34)
#x0   <- seq(79.75, 79.85, by=0.01)
#x0   <- sort(x0)

#x1   <- round(x0, digits=2)
#x2   <- Round(x0, 0.01)
#xxx  <- data.frame(x0, x1, x2, diff=abs(x1-x2)) 
#print(xxx)
#print(max(xxx[ ,"diff"]))

#ldt <- CheckRounding(x1) 
#ldt
#1-pchisq(sum(ldt[ ,"chi2"]), df=9)

#ldt <- CheckRounding(x2) 
#ldt
#1-pchisq(sum(ldt[ ,"chi2"]), df=9)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ***************************************************************************
#  F_chi2.PNV.tr.lms.R

# ===========================================================================
#  Zielfunktion für gestutzte gruppierte Schätzung von lambda, MW und SD
#  einer PNV 

#  14.05.2021 Call changed, processing of fixed theta components added
#  13.01.2012 Call changed, parameters for RL@.pen added, q.fact now 
#             alphabetically
#  06.12.2020 Call changed, x.kde.mode removed, x.lt.mode, x.ge.mode changed 
#  05.12.2020 Call changed, x.kde.mode added
#  24.11.2020 Call changed, x.lt.mode, x.ge.mode added
#  17.09.2019 Faktor für Strafterm "negative Prävalenz" eingeführt
#  28.01.2018
# ===========================================================================

chi2.PNV.tr.lms <- function(ttheta.ini.est, RB, 
                            ttheta.ini.fix, idx.fix, idx.est,
                            x.hist,x.tr.lo,x.tr.hi,
                            x.lt.tr.n, x.ge.tr.n,
                            l.fact, p.fact, r.fact, w.fact,
                            x.Q1, x.Q2, RL1.p, RL2.p,
                            df.est,df.con, opt.crit.only, fastnull 
                            )
{
  #  LF auf Grundlage der bedingten Dichte einer gestutzten PNV
  #  Diskrete Version
  #  Berechnet wird -log-LF (wegen Nutzung durch nlm() )

  #  ttheta    auf R transformierte Parameter der PNV
  #  RB        Matrix der Randbedingungen an theta
  #  idx.fix   contains indices of fixed components in theta
  #  idx.est   contains indices of theta components tat are being estimated
  #  x.hist    Histogramm mit Komponenten $breaks, $counts, lineare Skala
  #  x.tr.lo   unterer truncation point, lineare Skala
  #  x.tr.hi   oberer truncation point, lineare Skala
  #            truncation points müssen breaks sein!
  #  df.est    Freiheitsgrade für geschätzte Parameter
  #  opt.crit.only Steuert die Ausgabe, nur opt.crit oder mehr, siehe chi2trunc 
  #  fastnull  Ersatzwert für kleine Erwartungswerte
  #  l.fact    Faktor für Gewichtung "Länge des truncation Intervalls"
  #  p.fact    Faktor für Strafterm "negative Prävalenz" in chi2-Funktion
  #  r.fact    Faktor für Strafterm ...
  #  w.fact    Faktor für Strafterm "zu hohe Schätzung" in chi2-Funktion

  # Compose the full parameter
  ttheta <- Compose(ttheta.ini.fix, ttheta.ini.est, idx.fix, idx.est)

  theta  <- R2o(ttheta,RB)
  lambda <- theta[1] 
  mue    <- theta[2]
  sigma  <- theta[3]

  #  chi^2-Anteil aus der gestutzten Verteilung berechnen
  
  res <- chi2trunc(x.hist,x.tr.lo,x.tr.hi,
                   x.lt.tr.n, x.ge.tr.n,
                   x.Q1, x.Q2, RL1.p, RL2.p,
                   lambda,mue,sigma,df.est, df.con, 
                   l.fact, p.fact, r.fact, w.fact, 
                   opt.crit.only=opt.crit.only,fastnull=fastnull)

  return(res)  
}

# ****************************************************************************
#  F_chi2trunc_V4.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  CHANGE HISTORY
#  13.01.2021  x.Q1, x.Q2, RL1.p, RL2.p, r.fact added to call for penalty 
#              calculation. Sequence of @.fact now alphabetically.
#  06.12.2020  Calculation of xl.n, xr.n and prevalence changed back to old
#              approach using the TI limits, because usage of mode, though 
#              theoretically correct, gives sometimes nonsense, because
#              of rounding effects when counting empirical frequencies
#              below / above the mode
#  04.12.2020  New calculation of xl.n, xr.n and prevalence new, now separation
#              by mode 
#              Strategy x.tr.eff changed
#              'check' removed
#  19.11.2020  Optimality criterion set back to chi2/df, because standard test
#              data creates p values = 0 for initial values from qqw,
#              therefore no progress in nlm. Attempt to "normalize" chi2 to get
#              more independent from df? 
#  29.03.2020  Optimality criterion renamed to opt.crit
#              chi2.per.df removed from output
#  06.01.2020  Prevalences l, c, r introduced
#  01.01.2020  prev.tmc.pen new defined
#              chi2.per.df  new defined
# ============================================================================

chi2trunc <- function(x.hist,x.tr.lo,x.tr.hi, 
                      x.lt.tr.n, x.ge.tr.n,
                      x.Q1, x.Q2, RL1.p, RL2.p,
                      lambda,mue,sigma,
                      df.est,df.con,
                      l.fact=1.0, p.fact=1.0, r.fact=1, w.fact=0.15, 
                      opt.crit.only=FALSE,fastnull=1.e-10)
{  
   #  Berechnet modifizierte chi2-Anpassung für eine gestutzte 
   #  Power-Normalverteilung
   #
   #  EINGABE
   #  x.hist     Histogramm von x (gesamter Datensatz)
   #  x.tr.lo    unterer Kappungspunkt (truncation point)
   #  x.tr.hi    oberer Kappungspunkt (truncation point)
   #  x.lt.tr.n  number of values <   lower TI limit
   #  x.ge.tr.n  number of values >=  upper TI limit
   #  lambda )
   #  mue    )   Parameter der angepassten PNV 
   #  sigma  )
   #  df.est     Freiheitsgrade für geschätzte Parameter (sind nicht 
   #             zwangsläufig 3!)
   #  l.fact     Gewichtungsfaktor "Länge des truncation intervals"
   #  w.fact     Gewichtungsfaktor für Strafterm zu hohe Vorhersage
   #  p.fact     Gewichtungsfaktor für Strafterm negative Prävalenz
   #  opt.crit.only  TRUE: nur der chi2-Wert für Optimierung wird zurückgegeben
   #             FALSE: vollständige Information, siehe AUSGABE
   #
   #  AUSGABE
   #  If opt.crit.only=TRUE
   #  Wert des Optimierungskriteriums (kleiner ist besser)

   #  If opt.crit.only=FALSE
   #  Full information, see return statement at end of function

  # ========================================================================= 

  #  Gesamtzahl Intervalle
  xcounts.n <- length(x.hist$counts)

  #  Gesamtzahl breaks
  xbreaks.n <- length(x.hist$breaks)

  #  Gesamtzahl Fälle
  x.n <- sum(x.hist$counts)

  #  Zu nutzender Bereich von x (weil als nichtpathologisch erkannt)
  #  Werte nicht notwendig Klassengrenzen des Histogramms
  #  Angegebene Werte gehören zum Intervall

  #  Die dazugehörigen Intervall-Untergrenzen und ihre Indices (in breaks) 
  #  finden
  #  Version before 24.11.2020

  x.tr.eff.lo <- min(x.hist$breaks[x.hist$breaks >= x.tr.lo]) 
  x.tr.eff.hi <- max(x.hist$breaks[x.hist$breaks <= x.tr.hi])

  # Index in breaks der linken Grenze des am weitesten links liegenden 
  # im Truncation-Bereich (in breaks!)
  x.tr.eff.lo.idx <- which(x.hist$breaks == x.tr.eff.lo)

  # Index der rechten Grenze des am weitesten rechts liegenden Intervalls 
  # im Truncation-Bereich (in breaks!)

  x.tr.eff.hi.idx <- which(x.hist$breaks == x.tr.eff.hi)

  # sink("Kontrolle.txt", append=TRUE)
  #cat("\n  chi2trunc_V3 ", subtitle,  
  #    "\nx.tr.lo        ",x.tr.lo,        "  x.tr.hi        ",x.tr.hi,
  #    "\nx.tr.eff.lo    ",x.tr.eff.lo,    "  x.tr.eff.hi    ",x.tr.eff.hi,
  #    "\nx.tr.eff.lo.idx",x.tr.eff.lo.idx,"  x.tr.eff.hi.idx",x.tr.eff.hi.idx,
  #    "\n")
  #print(x.hist$breaks)
  # sink()
  
  if (length(x.tr.eff.hi.idx) != 1)
  { cat("\n  [chi2trunc_V4] xbreaks\n")
    print(x.hist$breaks)

    cat("\n  [chi2trunc_V4] x.tr.eff.hi\n")
    print(x.tr.eff.hi)

    cat("\n  [chi2trunc_V4] x.tr.eff.hi.idx\n")
    print(x.tr.eff.hi.idx)
  }

  #  Mögliches Problem: x.tr.lo und x.tr.hi definieren zusammen weniger 
  #  als df.est + df.con + 1  Intervalle (bins). 
  x.tr.hist.counts.n <- x.tr.eff.hi.idx - x.tr.eff.lo.idx 

  if (x.tr.hist.counts.n < df.est + df.con + 1)
  { cat("+++ [chi2trunc]: Angegebene Grenzen definieren zu wenige Intervalle +++ \n")
    print(data.frame(UG=x.hist$breaks[1:xcounts.n],
                     OG=x.hist$breaks[2:(xcounts.n+1)],
                     Anzahl=x.hist$counts))
  
    #cat("\n [chi2trunc]",
    #    "\nx.tr.lo        ",x.tr.lo,        "  x.tr.hi        ",x.tr.hi,
    #    "\nx.tr.eff.lo    ",x.tr.eff.lo,    "  x.tr.eff.hi    ",x.tr.eff.hi,
    #    "\nx.tr.eff.lo.idx",x.tr.eff.lo.idx,"  x.tr.eff.hi.idx",x.tr.eff.hi.idx,
    #    "\n")

    if (opt.crit.only)
    { return(NA) }  else
    { return(list(res=c(opt.crit=NA,
                        chi2.total=NA, chi2.total.df=NA, 
                        chi2.trun=NA, chi2.trun.df=NA,
                        chi2.trun.p=NA,
                        chi2.path=NA,
                        chi2.df=NA,
                        prev.l.tmc=NA,
                        prev.c.tmc=NA,
                        prev.r.tmc=NA,
                        prev.l.tmc.pen=NA,
                        prev.r.tmc.pen=NA,
                        xc.n.tmc=NA),
                  tab=NA,
                  tab.lo=NA,
                  tab.hi=NA))
    }  

  } else
  { #  Kein Problem mit Anzahl Intervalle
    #  Nicht pathologischer Teil des Histogramms
   
    #  breaks des nichtpathologischen Teils
    x.tr.hist.breaks <- x.hist$breaks[x.tr.eff.lo.idx:(x.tr.eff.hi.idx)]

    #  counts des nichtpathologischen Teils
    x.tr.hist.counts <- x.hist$counts[x.tr.eff.lo.idx:(x.tr.eff.hi.idx-1)]

    #  Anzahl im nichtpathologischen Teil
    x.tr.n <- sum(x.tr.hist.counts)

    #  Anzahl bins im nichtpathologischen Teil ist
    #  x.tr.hist.counts.n, schon oben bestimmt
    #  Estimate prevalence per component

    cp <- CalcPrev(x.n, x.tr.n, x.lt.tr.n, x.ge.tr.n, 
                   x.tr.lo, x.tr.hi, lambda, mue, sigma)
    prev.l.tmc <- cp["prev.l"]
    prev.c.tmc <- cp["prev.c"]
    prev.r.tmc <- cp["prev.r"]
    xl.n.tmc   <- cp["l.n"]
    xc.n.tmc   <- cp["c.n"]
    xr.n.tmc   <- cp["r.n"]

    #  Unbedingte Wahrscheinlichkeiten der nichtpathologischen Intervalle

    x.tr.hist.cdf <- cdf.PN(x.tr.hist.breaks,lambda,mue,sigma)
    x.tr.hist.pdf <- diff(x.tr.hist.cdf) 

    x.tr.hist.p   <- sum(x.tr.hist.pdf)  # proportion of all nonpath in TI

    # Letzteres wird Null, wenn die Parameter zu falsch. Überlauf verhindern.
    x.tr.hist.p   <- max(1.e-5,x.tr.hist.p) 

    #  Bedingte Wahrscheinlichkeiten der nichtpathologischen Intervalle
    x.tr.hist.pdfc <- x.tr.hist.pdf / x.tr.hist.p
    nerw           <- x.tr.hist.pdfc * x.tr.n
    nerw[nerw < fastnull]  <- fastnull
    diff           <- x.tr.hist.counts - nerw  
    chi2.i         <- (x.tr.hist.counts - nerw)^2 / nerw

    #  Überlauf?
    if (any(!is.finite(chi2.i)))
    { cat("\n Invalid chi2.i value",
          "\nlambda =",lambda, "  mue =",mue, "  sigma =",sigma,"\n")
      cat("\nProbabilities  x.tr.hist.p =",x.tr.hist.p,"\n")
      print(data.frame(breaks=x.tr.hist.breaks[1:x.tr.hist.counts.n],
                       cdf=x.tr.hist.cdf[1:x.tr.hist.counts.n],  
                       pdfc=x.tr.hist.pdfc))
      print(data.frame(nobs=x.tr.hist.counts,nerw=nerw,chi2.i=chi2.i)) 
    }

    #  Chi-Quadrat-Beiträge des nicht-pathologischen Bereichs,
    #  für Anpassungstest benutzt  
    chi2.trun      <- sum(chi2.i)
    chi2.trun.df   <- x.tr.hist.counts.n - df.est - 1
    chi2.trun.p    <- 1-pchisq(chi2.trun,chi2.trun.df) 

    #  Beobachtete, erwartete Anzahlen und chi2-Beiträge
    tab <-  data.frame(UG=x.tr.hist.breaks[1:x.tr.hist.counts.n],
                       OG=x.tr.hist.breaks[2:(x.tr.hist.counts.n+1)],
                       nobs=x.tr.hist.counts,
                       nerw=nerw,
                       diff=diff,
                       chi2=chi2.i)

    #  Meldung, wenn Überlauf 
    if (!is.numeric(chi2.trun))
    { cat("\n+++ [chi2trunc]: Invalid chi-square, non pathological subset",
          "\n                 Sum chi^2           =", chi2.trun,
          "\n                 Degrees of freeedom =", chi2.trun.df,
          "\n")
      print(tab)
    } 

    #  Bis hier wurde nur die Anpassung in den vermutet nichtpathologischen
    #  Intervallen bewertet.
    #  Zusätzlich die zweite Information nutzen: in den pathologischen 
    #  Bereichen darf die angepasste Kurve nicht über den Daten liegen.
    #  Entsprechenden chi2-Beitrag berechnen. Nur Abweichungen nach oben
    #  spielen eine Rolle.  

    #  Intervalle unterhalb des nichtpathologischen Bereichs
    if (1 < x.tr.eff.lo.idx)
    { 
      x.lo.hist.breaks <- x.hist$breaks[1:x.tr.eff.lo.idx]
      x.lo.hist.counts <- x.hist$counts[1:(x.tr.eff.lo.idx-1)]

      #  Unbedingte Wahrscheinlichkeiten der Intervalle unterhalb
      x.lo.hist.cdf <- cdf.PN(x.lo.hist.breaks,lambda,mue,sigma)
      x.lo.hist.pdf <- x.lo.hist.cdf[2:x.tr.eff.lo.idx] - 
                       x.lo.hist.cdf[1:(x.tr.eff.lo.idx-1)] 

      #  Erwartete Anzahlen in den Intervallen unterhalb des truncation 
      #  intervals(unter Verwendung der Schätzung für xc.n.tmc)
      nerw.lo       <- x.lo.hist.pdf * xc.n.tmc
      nerw.lo[nerw.lo < fastnull] <- fastnull
      diff.lo       <- x.lo.hist.counts - nerw.lo  # positiv, wenn erwarteter
                                                   # Wert < beobachtet (ok)

      #  Ungewichtete chi^2-Beiträge
      chi2.lo.i     <- (diff.lo)^2 / nerw.lo

      #  Weights
      w.lo.i             <- w.fact * pchisq(chi2.lo.i,1)
      w.lo.i[diff.lo>=0] <- 0

      #  Gewichtete chi^2-Beiträge
      chi2.lo.w.i     <- w.lo.i * chi2.lo.i

      #  Summen
      chi2.lo       <- sum(chi2.lo.i)
      chi2.lo.w     <- sum(chi2.lo.w.i)
      df.lo         <- length(x.lo.hist.counts)

      #  Beobachtete, erwartete Anzahlen und chi2-Beiträge
      tab.lo <-  data.frame(UG=x.lo.hist.breaks[1:(x.tr.eff.lo.idx-1)],
                            OG=x.lo.hist.breaks[2:(x.tr.eff.lo.idx)],
                       nobs=x.lo.hist.counts,
                       nerw=nerw.lo,
                       diff=diff.lo,
                       w=w.lo.i,
                       chi2=chi2.lo.i,
                       chi2.w=chi2.lo.w.i)

    } else
    { # Keine pathologischen Intervalle links 
      tab.lo    <- NA
      chi2.lo.w <-  0 
    }

    #  Intervalle oberhalb des nichtpathologischen Bereichs
    if (x.tr.eff.hi.idx < xbreaks.n)
    {
      x.hi.hist.breaks <- x.hist$breaks[x.tr.eff.hi.idx:xbreaks.n]
      x.hi.hist.counts <- x.hist$counts[x.tr.eff.hi.idx:xcounts.n]

      #  Unbedingte Wahrscheinlichkeiten der Intervalle
      x.hi.n        <- length(x.hi.hist.breaks)
      x.hi.hist.cdf <- cdf.PN(x.hi.hist.breaks,lambda,mue,sigma)
      x.hi.hist.pdf <- x.hi.hist.cdf[2:x.hi.n] - 
                       x.hi.hist.cdf[1:(x.hi.n-1)] 

      #  Erwartete Anzahlen in den Intervallen (unter der Annahme,
      #  alles sei nichtpathologisch)
      nerw.hi         <- x.hi.hist.pdf * xc.n.tmc
      nerw.hi[nerw.hi < fastnull] <- fastnull
      diff.hi         <- x.hi.hist.counts - nerw.hi

      #  Ungewichtete chi^2-Beiträge
      chi2.hi.i     <- (diff.hi)^2 / nerw.hi

      #  Penalty weights for terms outside truncation interval
      w.hi.i             <- w.fact * pchisq(chi2.hi.i,1)
      w.hi.i[diff.hi>=0] <- 0

      #  Gewichtete chi^2-Beiträge
      chi2.hi.w.i     <- w.hi.i * chi2.hi.i

      #  Summen
      chi2.hi       <- sum(chi2.hi.i)
      chi2.hi.w     <- sum(chi2.hi.w.i)
      df.hi         <- length(x.hi.hist.counts)

      #  Beobachtete, erwartete Anzahlen und chi2-Beiträge
      tab.hi <-  data.frame(UG=x.hi.hist.breaks[1:(x.hi.n-1)],
                            OG=x.hi.hist.breaks[2:x.hi.n],
                            nobs=x.hi.hist.counts,
                            nerw=nerw.hi,
                       diff=diff.hi,
                       w=w.hi.i,
                       chi2=chi2.hi.i,
                       chi2.w=chi2.hi.w.i)

    } else
    { # Keine pathologischen Intervalle rechts 
      tab.hi    <- NA
      chi2.hi.w <-  0 
    }

    #  (Nur die gewichteten) chi2-Beiträge zusammenfassen
    #  Alte Taktik bis 28.08.2018 
    #chi2    <- chi2 + chi2.lo.w + chi2.hi.w
    #df      <- df + df.lo + df.hi - df.est - 1
    #chi2.df <- chi2 / df 

    #  Neue Taktik ab 28.08.2018 
    #  Alle chi2-Beiträge zählen zum Optimierungskriterium
    chi2.path     <- chi2.lo.w + chi2.hi.w
    chi2.total    <- chi2.trun + chi2.path

    #  Penalty for estimated pathological prevalence < 0
    if (prev.l.tmc >=  0)
    { prev.l.tmc.pen <- 0 } else
    { #  Negative prevalence. Requires penalty. Idea: Treat predicted 
      #  prevalence like a bin with expectation zero and observed count
      #  abs(prev) * x.n. Replace zero by 1.e-2 and calculate chi2 
      #  contribution as penalty 
      # prev.tmc.pen <- p.fact * (abs(prev.tmc)*x.n-0.01)^2/0.01

      #  Previous idea is to radical or needs a very small p.fact.
      #  Approach 01.01.2020
      prev.l.tmc.pen <- p.fact * (prev.l.tmc < 0) * 
                        abs(prev.l.tmc) * chi2.trun
    }
    #cat(" [chi2trunc] chi2, prev.l, pen ", 
    #    chi2.total, prev.l.tmc, prev.l.tmc.pen, "\n")  
   
    if (prev.r.tmc >=  0)
    { prev.r.tmc.pen <- 0 } else
    { #  Negative prevalence. Requires penalty. Idea: Treat predicted 
      #  prevalence like a bin with expectation zero and observed count
      #  abs(prev) * x.n. Replace zero by 1.e-2 and calculate chi2 
      #  contribution as penalty 
      # prev.tmc.pen <- p.fact * (abs(prev.tmc)*x.n-0.01)^2/0.01

      #  Previous idea is to radical or needs a very small p.fact.
      #  Approach 01.01.2020
      prev.r.tmc.pen <- p.fact * (prev.r.tmc < 0) * 
                        abs(prev.r.tmc) * chi2.trun
    }
    #cat(" [chi2trunc] chi2, prev.r, pen ", 
    #    chi2.total, prev.r.tmc, prev.r.tmc.pen, "\n")  

    #  Penalty for RLs outside [x.Q1, x.Q2], where x.Q@ are the RL1.p and
    #  RL2.p quantiles of the complete stratum
    #  Candidates for 2 additional penalty terms
    x.RL1 <- q.PN(RL1.p, lambda, mue, sigma)
    x.RL2 <- q.PN(RL2.p, lambda, mue, sigma)
    RL.diff  <- x.RL2 - x.RL1
    RL1.pen  <- max(0, (x.Q1-x.RL1)/RL.diff )
    RL2.pen  <- max(0, (x.RL2-x.Q2)/RL.diff )

    chi2.total    <- chi2.total + 
                     prev.l.tmc.pen + prev.r.tmc.pen +
                     r.fact*RL1.pen +  r.fact*RL2.pen

    #  Als Freiheitsgrade zählen nur die Intervalle im truncation-Bereich,
    #  weil nur hier beobachtet-erwartet sinnvoll ist. Bei korrekter 
    #  Anpassung entsteht kein chi2-Beitrag außerhalb.

    chi2.total.df <- chi2.trun.df
    # Die Normierung auf die Zahl der DF reicht nicht, weil der Effekt der 
    # DF auf den chi^2-Wert nicht linear ist. 
    # chi2.per.df   <- chi2.total / chi2.total.df

    #  p.fit im Prinzip gut, aber schlecht, wenn Startwerte zu schlecht. Dann
    #  Gefahr des Stehenbleibens am Startpunkt. 
    #  Gleich guter Fit für ein Intervall mit mehr Fällen ist besser.
    #  Entsprechende Gewichtung durch l.fact.
    #  30.03.2020
opt.crit <- (chi2.total / chi2.total.df) * (1 + l.fact *(x.n/x.tr.n)^2) 

    #  11.09.2020  
    #  18.11.2020 Diese Lösung führt schon bei den Routine-Testdaten zum
    #             Stillstand 
####opt.crit <- (1-chi2.trun.p) * (1 + l.fact *(x.n/x.tr.n)^2) 

    #  Ergebnis mitteilen
    if (opt.crit.only)
    { 
      return(opt.crit) 
    }  else
    { 
      opt.crit <- unname(opt.crit) 
      chi2.total <- unname(chi2.total) 
      chi2.total.df <- unname(chi2.total.df) 
      chi2.trun <- unname(chi2.trun) 
      chi2.trun.df <- unname(chi2.trun.df) 
      chi2.trun.p <- unname(chi2.trun.p) 
      chi2.path <- unname(chi2.path) 
      prev.l.tmc <- unname(prev.l.tmc) 
      prev.c.tmc <- unname(prev.c.tmc) 
      prev.r.tmc <- unname(prev.r.tmc) 
      prev.l.tmc.pen <- unname(prev.l.tmc.pen) 
      prev.r.tmc.pen <- unname(prev.r.tmc.pen) 
      xl.n.tmc <- unname(xl.n.tmc) 
      xc.n.tmc <- unname(xc.n.tmc) 
      xr.n.tmc <- unname(xr.n.tmc) 

      return(list(res=c(opt.crit=opt.crit,
                        chi2.total=chi2.total, chi2.total.df=chi2.total.df, 
                        chi2.trun=chi2.trun, chi2.trun.df=chi2.trun.df,
                        chi2.trun.p=chi2.trun.p,
                        chi2.path=chi2.path,
                        prev.l.tmc=prev.l.tmc,
                        prev.c.tmc=prev.c.tmc,
                        prev.r.tmc=prev.r.tmc,
                        prev.l.tmc.pen=prev.l.tmc.pen,
                        prev.r.tmc.pen=prev.r.tmc.pen,
                        xl.n.tmc=xl.n.tmc,
                        xc.n.tmc=xc.n.tmc,
                        xr.n.tmc=xr.n.tmc
                        ),
                  tab=tab,
                  tab.lo=tab.lo,
                  tab.hi=tab.hi))
    }  
  }
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


# ****************************************************************************
#  F_chi2trunc0_V3.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  CHANGE HISTORY
#  13.01.2021  backup version of F_chi2trunc_V3.R
#  06.12.2020  Calculation of xl.n, xr.n and prevalence changed back to old
#              approach using the TI limits, because usage of mode, though 
#              theoretically correct, gives sometimes nonsense, because
#              of rounding effects when counting empirical frequencies
#              below / above the mode
#  04.12.2020  New calculation of xl.n, xr.n and prevalence new, now separation
#              by mode 
#              Strategy x.tr.eff changed
#              'check' removed
#  19.11.2020  Optimality criterion set back to chi2/df, because standard test
#              data creates p values = 0 for initial values from qqw,
#              therefore no progress in nlm. Attempt to "normalize" chi2 to get
#              more independent from df? 
#  29.03.2020  Optimality criterion renamed to opt.crit
#              chi2.per.df removed from output
#  06.01.2020  Prevalences l, c, r introduced
#  01.01.2020  prev.tmc.pen new defined
#              chi2.per.df  new defined
# ============================================================================

chi2trunc0 <- function(x.hist,x.tr.lo,x.tr.hi, 
                      x.lt.tr.n, x.ge.tr.n,
                      lambda,mue,sigma,df.est,df.con,
                      l.fact=1.0, w.fact=0.15, p.fact=1.0, 
                      opt.crit.only=FALSE,fastnull=1.e-10)
{  
   #  Berechnet modifizierte chi2-Anpassung für eine gestutzte 
   #  Power-Normalverteilung
   #
   #  EINGABE
   #  x.hist     Histogramm von x (gesamter Datensatz)
   #  x.tr.lo    unterer Kappungspunkt (truncation point)
   #  x.tr.hi    oberer Kappungspunkt (truncation point)
   #  x.lt.tr.n  number of values <   lower TI limit
   #  x.ge.tr.n  number of values >=  upper TI limit
   #  lambda )
   #  mue    )   Parameter der angepassten PNV 
   #  sigma  )
   #  df.est     Freiheitsgrade für geschätzte Parameter (sind nicht 
   #             zwangsläufig 3!)
   #  l.fact     Gewichtungsfaktor "Länge des truncation intervals"
   #  w.fact     Gewichtungsfaktor für Strafterm zu hohe Vorhersage
   #  p.fact     Gewichtungsfaktor für Strafterm negative Prävalenz
   #  opt.crit.only  TRUE: nur der chi2-Wert für Optimierung wird zurückgegeben
   #             FALSE: vollständige Information, siehe AUSGABE
   #
   #  AUSGABE
   #  If opt.crit.only=TRUE
   #  Wert des Optimierungskriteriums (kleiner ist besser)

   #  If opt.crit.only=FALSE
   #  Full information, see return statement at end of function

  # ========================================================================= 

  #  Gesamtzahl Intervalle
  xcounts.n <- length(x.hist$counts)

  #  Gesamtzahl breaks
  xbreaks.n <- length(x.hist$breaks)

  #  Gesamtzahl Fälle
  x.n <- sum(x.hist$counts)

  #  Zu nutzender Bereich von x (weil als nichtpathologisch erkannt)
  #  Werte nicht notwendig Klassengrenzen des Histogramms
  #  Angegebene Werte gehören zum Intervall

  #  Die dazugehörigen Intervall-Untergrenzen und ihre Indices (in breaks) 
  #  finden
  #  Version before 24.11.2020

  ###cat("\n vergleich breaks / x.tr.lo, x.tr.lo =", x.tr.lo, "\n")
  ###print(x.hist$breaks >= x.tr.lo)

  x.tr.eff.lo <- min(x.hist$breaks[x.hist$breaks >= x.tr.lo]) 
  x.tr.eff.hi <- max(x.hist$breaks[x.hist$breaks <= x.tr.hi])

  #  Version starting 24.11.2020
  ###?x.tr.eff.lo <- max(x.hist$breaks[x.hist$breaks <= x.tr.lo]) 
  ###?x.tr.eff.hi <- min(x.hist$breaks[x.hist$breaks >= x.tr.hi])

  # Index in breaks der linken Grenze des am weitesten links liegenden 
  # im Truncation-Bereich (in breaks!)
  x.tr.eff.lo.idx <- which(x.hist$breaks == x.tr.eff.lo)

  # Index der rechten Grenze des am weitesten rechts liegenden Intervalls 
  # im Truncation-Bereich (in breaks!)

  x.tr.eff.hi.idx <- which(x.hist$breaks == x.tr.eff.hi)

  # sink("Kontrolle.txt", append=TRUE)
  #cat("\n  chi2trunc_V3 ", subtitle,  
  #    "\nx.tr.lo        ",x.tr.lo,        "  x.tr.hi        ",x.tr.hi,
  #    "\nx.tr.eff.lo    ",x.tr.eff.lo,    "  x.tr.eff.hi    ",x.tr.eff.hi,
  #    "\nx.tr.eff.lo.idx",x.tr.eff.lo.idx,"  x.tr.eff.hi.idx",x.tr.eff.hi.idx,
  #    "\n")
  #print(x.hist$breaks)
  # sink()
  
  if (length(x.tr.eff.hi.idx) != 1)
  { cat("\n  [chi2trunc_V3] xbreaks\n")
    print(x.hist$breaks)

    cat("\n  [chi2trunc_V3] x.tr.eff.hi\n")
    print(x.tr.eff.hi)

    cat("\n  [chi2trunc_V3] x.tr.eff.hi.idx\n")
    print(x.tr.eff.hi.idx)
  }

  #  Mögliches Problem: x.tr.lo und x.tr.hi definieren zusammen weniger 
  #  als df.est + df.con + 1  Intervalle (bins). 
  x.tr.hist.counts.n <- x.tr.eff.hi.idx - x.tr.eff.lo.idx 

  if (x.tr.hist.counts.n < df.est + df.con + 1)
  { cat("+++ [chi2trunc]: Angegebene Grenzen definieren zu wenige Intervalle +++ \n")
    print(data.frame(UG=x.hist$breaks[1:xcounts.n],
                     OG=x.hist$breaks[2:(xcounts.n+1)],
                     Anzahl=x.hist$counts))
  
    #cat("\n [chi2trunc]",
    #    "\nx.tr.lo        ",x.tr.lo,        "  x.tr.hi        ",x.tr.hi,
    #    "\nx.tr.eff.lo    ",x.tr.eff.lo,    "  x.tr.eff.hi    ",x.tr.eff.hi,
    #    "\nx.tr.eff.lo.idx",x.tr.eff.lo.idx,"  x.tr.eff.hi.idx",x.tr.eff.hi.idx,
    #    "\n")

    if (opt.crit.only)
    { return(NA) }  else
    { return(list(res=c(opt.crit=NA,
                        chi2.total=NA, chi2.total.df=NA, 
                        chi2.trun=NA, chi2.trun.df=NA,
                        chi2.trun.p=NA,
                        chi2.path=NA,
                        chi2.df=NA,
                        prev.l.tmc=NA,
                        prev.c.tmc=NA,
                        prev.r.tmc=NA,
                        prev.l.tmc.pen=NA,
                        prev.r.tmc.pen=NA,
                        xc.n.tmc=NA),
                  tab=NA,
                  tab.lo=NA,
                  tab.hi=NA))
    }  

  } else
  { #  Kein Problem mit Anzahl Intervalle
    #  Nicht pathologischer Teil des Histogramms
   
    #  breaks des nichtpathologischen Teils
    x.tr.hist.breaks <- x.hist$breaks[x.tr.eff.lo.idx:(x.tr.eff.hi.idx)]

    #  counts des nichtpathologischen Teils
    x.tr.hist.counts <- x.hist$counts[x.tr.eff.lo.idx:(x.tr.eff.hi.idx-1)]

    #  Anzahl im nichtpathologischen Teil
    x.tr.n <- sum(x.tr.hist.counts)

    #  Anzahl bins im nichtpathologischen Teil ist
    #  x.tr.hist.counts.n, schon oben bestimmt
    #  Estimate prevalence per component

    cp <- CalcPrev(x.n, x.tr.n, x.lt.tr.n, x.ge.tr.n, 
                   x.tr.lo, x.tr.hi, lambda, mue, sigma)
    prev.l.tmc <- cp["prev.l"]
    prev.c.tmc <- cp["prev.c"]
    prev.r.tmc <- cp["prev.r"]
    xl.n.tmc   <- cp["l.n"]
    xc.n.tmc   <- cp["c.n"]
    xr.n.tmc   <- cp["r.n"]

    #  Unbedingte Wahrscheinlichkeiten der nichtpathologischen Intervalle

    x.tr.hist.cdf <- cdf.PN(x.tr.hist.breaks,lambda,mue,sigma)
    x.tr.hist.pdf <- diff(x.tr.hist.cdf) 

    x.tr.hist.p   <- sum(x.tr.hist.pdf)  # proportion of all nonpath in TI

    # Letzteres wird Null, wenn die Parameter zu falsch. Überlauf verhindern.
    x.tr.hist.p   <- max(1.e-5,x.tr.hist.p) 

    #  Bedingte Wahrscheinlichkeiten der nichtpathologischen Intervalle
    x.tr.hist.pdfc <- x.tr.hist.pdf / x.tr.hist.p
    nerw           <- x.tr.hist.pdfc * x.tr.n
    nerw[nerw < fastnull]  <- fastnull
    diff           <- x.tr.hist.counts - nerw  
    chi2.i         <- (x.tr.hist.counts - nerw)^2 / nerw

    #  Überlauf?
    if (any(!is.finite(chi2.i)))
    { cat("\n Invalid chi2.i value",
          "\nlambda =",lambda, "  mue =",mue, "  sigma =",sigma,"\n")
      cat("\nProbabilities  x.tr.hist.p =",x.tr.hist.p,"\n")
      print(data.frame(breaks=x.tr.hist.breaks[1:x.tr.hist.counts.n],
                       cdf=x.tr.hist.cdf[1:x.tr.hist.counts.n],  
                       pdfc=x.tr.hist.pdfc))
      print(data.frame(nobs=x.tr.hist.counts,nerw=nerw,chi2.i=chi2.i)) 
    }

    #  Chi-Quadrat-Beiträge des nicht-pathologischen Bereichs,
    #  für Anpassungstest benutzt  
    chi2.trun      <- sum(chi2.i)
    chi2.trun.df   <- x.tr.hist.counts.n - df.est - 1
    chi2.trun.p    <- 1-pchisq(chi2.trun,chi2.trun.df) 

    #  Beobachtete, erwartete Anzahlen und chi2-Beiträge
    tab <-  data.frame(UG=x.tr.hist.breaks[1:x.tr.hist.counts.n],
                       OG=x.tr.hist.breaks[2:(x.tr.hist.counts.n+1)],
                       nobs=x.tr.hist.counts,
                       nerw=nerw,
                       diff=diff,
                       chi2=chi2.i)

    #  Meldung, wenn Überlauf 
    if (!is.numeric(chi2.trun))
    { cat("\n+++ [chi2trunc]: Invalid chi-square, non pathological subset",
          "\n                 Sum chi^2           =", chi2.trun,
          "\n                 Degrees of freeedom =", chi2.trun.df,
          "\n")
      print(tab)
    } 

    #  Bis hier wurde nur die Anpassung in den vermutet nichtpathologischen
    #  Intervallen bewertet.
    #  Zusätzlich die zweite Information nutzen: in den pathologischen 
    #  Bereichen darf die angepasste Kurve nicht über den Daten liegen.
    #  Entsprechenden chi2-Beitrag berechnen. Nur Abweichungen nach oben
    #  spielen eine Rolle.  

    #  Intervalle unterhalb des nichtpathologischen Bereichs
    if (1 < x.tr.eff.lo.idx)
    { 
      x.lo.hist.breaks <- x.hist$breaks[1:x.tr.eff.lo.idx]
      x.lo.hist.counts <- x.hist$counts[1:(x.tr.eff.lo.idx-1)]

      #  Unbedingte Wahrscheinlichkeiten der Intervalle unterhalb
      x.lo.hist.cdf <- cdf.PN(x.lo.hist.breaks,lambda,mue,sigma)
      x.lo.hist.pdf <- x.lo.hist.cdf[2:x.tr.eff.lo.idx] - 
                       x.lo.hist.cdf[1:(x.tr.eff.lo.idx-1)] 

      #  Erwartete Anzahlen in den Intervallen unterhalb des truncation 
      #  intervals(unter Verwendung der Schätzung für xc.n.tmc)
      nerw.lo       <- x.lo.hist.pdf * xc.n.tmc
      nerw.lo[nerw.lo < fastnull] <- fastnull
      diff.lo       <- x.lo.hist.counts - nerw.lo  # positiv, wenn erwarteter
                                                   # Wert < beobachtet (ok)

      #  Ungewichtete chi^2-Beiträge
      chi2.lo.i     <- (diff.lo)^2 / nerw.lo

      #  Weights
      w.lo.i             <- w.fact * pchisq(chi2.lo.i,1)
      w.lo.i[diff.lo>=0] <- 0

      #  Gewichtete chi^2-Beiträge
      chi2.lo.w.i     <- w.lo.i * chi2.lo.i

      #  Summen
      chi2.lo       <- sum(chi2.lo.i)
      chi2.lo.w     <- sum(chi2.lo.w.i)
      df.lo         <- length(x.lo.hist.counts)

      #  Beobachtete, erwartete Anzahlen und chi2-Beiträge
      tab.lo <-  data.frame(UG=x.lo.hist.breaks[1:(x.tr.eff.lo.idx-1)],
                            OG=x.lo.hist.breaks[2:(x.tr.eff.lo.idx)],
                       nobs=x.lo.hist.counts,
                       nerw=nerw.lo,
                       diff=diff.lo,
                       w=w.lo.i,
                       chi2=chi2.lo.i,
                       chi2.w=chi2.lo.w.i)

    } else
    { # Keine pathologischen Intervalle links 
      tab.lo    <- NA
      chi2.lo.w <-  0 
    }

    #  Intervalle oberhalb des nichtpathologischen Bereichs
    if (x.tr.eff.hi.idx < xbreaks.n)
    {
      x.hi.hist.breaks <- x.hist$breaks[x.tr.eff.hi.idx:xbreaks.n]
      x.hi.hist.counts <- x.hist$counts[x.tr.eff.hi.idx:xcounts.n]

      #  Unbedingte Wahrscheinlichkeiten der Intervalle
      x.hi.n        <- length(x.hi.hist.breaks)
      x.hi.hist.cdf <- cdf.PN(x.hi.hist.breaks,lambda,mue,sigma)
      x.hi.hist.pdf <- x.hi.hist.cdf[2:x.hi.n] - 
                       x.hi.hist.cdf[1:(x.hi.n-1)] 

      #  Erwartete Anzahlen in den Intervallen (unter der Annahme,
      #  alles sei nichtpathologisch)
      nerw.hi         <- x.hi.hist.pdf * xc.n.tmc
      nerw.hi[nerw.hi < fastnull] <- fastnull
      diff.hi         <- x.hi.hist.counts - nerw.hi

      #  Ungewichtete chi^2-Beiträge
      chi2.hi.i     <- (diff.hi)^2 / nerw.hi

      #  Penalty weights for terms outside truncation interval
      w.hi.i             <- w.fact * pchisq(chi2.hi.i,1)
      w.hi.i[diff.hi>=0] <- 0

      #  Gewichtete chi^2-Beiträge
      chi2.hi.w.i     <- w.hi.i * chi2.hi.i

      #  Summen
      chi2.hi       <- sum(chi2.hi.i)
      chi2.hi.w     <- sum(chi2.hi.w.i)
      df.hi         <- length(x.hi.hist.counts)

      #  Beobachtete, erwartete Anzahlen und chi2-Beiträge
      tab.hi <-  data.frame(UG=x.hi.hist.breaks[1:(x.hi.n-1)],
                            OG=x.hi.hist.breaks[2:x.hi.n],
                            nobs=x.hi.hist.counts,
                            nerw=nerw.hi,
                       diff=diff.hi,
                       w=w.hi.i,
                       chi2=chi2.hi.i,
                       chi2.w=chi2.hi.w.i)

    } else
    { # Keine pathologischen Intervalle rechts 
      tab.hi    <- NA
      chi2.hi.w <-  0 
    }

    #  (Nur die gewichteten) chi2-Beiträge zusammenfassen
    #  Alte Taktik bis 28.08.2018 
    #chi2    <- chi2 + chi2.lo.w + chi2.hi.w
    #df      <- df + df.lo + df.hi - df.est - 1
    #chi2.df <- chi2 / df 

    #  Neue Taktik ab 28.08.2018 
    #  Alle chi2-Beiträge zählen zum Optimierungskriterium
    chi2.path     <- chi2.lo.w + chi2.hi.w
    chi2.total    <- chi2.trun + chi2.path

    #  Penalty for estimated pathological prevalence < 0
    if (prev.l.tmc >=  0)
    { prev.l.tmc.pen <- 0 } else
    { #  Negative prevalence. Requires penalty. Idea: Treat predicted 
      #  prevalence like a bin with expectation zero and observed count
      #  abs(prev) * x.n. Replace zero by 1.e-2 and calculate chi2 
      #  contribution as penalty 
      # prev.tmc.pen <- p.fact * (abs(prev.tmc)*x.n-0.01)^2/0.01

      #  Previous idea is to radical or needs a very small p.fact.
      #  Approach 01.01.2020
      prev.l.tmc.pen <- p.fact * (prev.l.tmc < 0) * 
                        abs(prev.l.tmc) * chi2.trun
    }
    #cat(" [chi2trunc] chi2, prev.l, pen ", 
    #    chi2.total, prev.l.tmc, prev.l.tmc.pen, "\n")  
   
    if (prev.r.tmc >=  0)
    { prev.r.tmc.pen <- 0 } else
    { #  Negative prevalence. Requires penalty. Idea: Treat predicted 
      #  prevalence like a bin with expectation zero and observed count
      #  abs(prev) * x.n. Replace zero by 1.e-2 and calculate chi2 
      #  contribution as penalty 
      # prev.tmc.pen <- p.fact * (abs(prev.tmc)*x.n-0.01)^2/0.01

      #  Previous idea is to radical or needs a very small p.fact.
      #  Approach 01.01.2020
      prev.r.tmc.pen <- p.fact * (prev.r.tmc < 0) * 
                        abs(prev.r.tmc) * chi2.trun
    }
    #cat(" [chi2trunc] chi2, prev.r, pen ", 
    #    chi2.total, prev.r.tmc, prev.r.tmc.pen, "\n")  

    #  Penalty for RLs outside [x.Q1, x.Q2], where x.Q@ are the RL1.p and
    #  RL2.p quantiles of the complete stratum

    chi2.total    <- chi2.total + prev.l.tmc.pen + prev.r.tmc.pen   

    #  Als Freiheitsgrade zählen nur die Intervalle im truncation-Bereich,
    #  weil nur hier beobachtet-erwartet sinnvoll ist. Bei korrekter 
    #  Anpassung entsteht kein chi2-Beitrag außerhalb.

    chi2.total.df <- chi2.trun.df
    # Die Normierung auf die Zahl der DF reicht nicht, weil der Effekt der 
    # DF auf den chi^2-Wert nicht linear ist. 
    # chi2.per.df   <- chi2.total / chi2.total.df

    #  p.fit im Prinzip gut, aber schlecht, wenn Startwerte zu schlecht. Dann
    #  Gefahr des Stehenbleibens am Startpunkt. 
    #  Gleich guter Fit für ein Intervall mit mehr Fällen ist besser.
    #  Entsprechende Gewichtung durch l.fact.
    #  30.03.2020
opt.crit <- (chi2.total / chi2.total.df) * (1 + l.fact *(x.n/x.tr.n)^2) 

    #  11.09.2020  
    #  18.11.2020 Diese Lösung führt schon bei den Routine-Testdaten zum
    #             Stillstand 
####opt.crit <- (1-chi2.trun.p) * (1 + l.fact *(x.n/x.tr.n)^2) 

    #  Ergebnis mitteilen
    if (opt.crit.only)
    { 
      return(opt.crit) 
    }  else
    { 
      opt.crit <- unname(opt.crit) 
      chi2.total <- unname(chi2.total) 
      chi2.total.df <- unname(chi2.total.df) 
      chi2.trun <- unname(chi2.trun) 
      chi2.trun.df <- unname(chi2.trun.df) 
      chi2.trun.p <- unname(chi2.trun.p) 
      chi2.path <- unname(chi2.path) 
      prev.l.tmc <- unname(prev.l.tmc) 
      prev.c.tmc <- unname(prev.c.tmc) 
      prev.r.tmc <- unname(prev.r.tmc) 
      prev.l.tmc.pen <- unname(prev.l.tmc.pen) 
      prev.r.tmc.pen <- unname(prev.r.tmc.pen) 
      xl.n.tmc <- unname(xl.n.tmc) 
      xc.n.tmc <- unname(xc.n.tmc) 
      xr.n.tmc <- unname(xr.n.tmc) 

      return(list(res=c(opt.crit=opt.crit,
                        chi2.total=chi2.total, chi2.total.df=chi2.total.df, 
                        chi2.trun=chi2.trun, chi2.trun.df=chi2.trun.df,
                        chi2.trun.p=chi2.trun.p,
                        chi2.path=chi2.path,
                        prev.l.tmc=prev.l.tmc,
                        prev.c.tmc=prev.c.tmc,
                        prev.r.tmc=prev.r.tmc,
                        prev.l.tmc.pen=prev.l.tmc.pen,
                        prev.r.tmc.pen=prev.r.tmc.pen,
                        xl.n.tmc=xl.n.tmc,
                        xc.n.tmc=xc.n.tmc,
                        xr.n.tmc=xr.n.tmc
                        ),
                  tab=tab,
                  tab.lo=tab.lo,
                  tab.hi=tab.hi))
    }  
  }
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


# ****************************************************************************
#  F_CIQuant.LNV.R

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Asymptotics 1-alpha confidence interval for quantiles of the 
#  log-normal distribution

#  See R.J.Serfling, Approximation theorems of mathematical statistics,
#  Wiley 1980, S. 77, corollary B
 
#  #  CHANGE HISTORY
#  15.05.2017 Start
# ============================================================================

CIQuant.LNV <- function(n,p,Q,mue,sigma,alpha=0.05)
{ 
  #  INPUT
  #  n       Size of the sample from which to calculate the quantiles
  #  p       Probability associated with the quantile
  #  Q       Quantile
  #  mue     )  Parameters of the log-normal distribution
  #  sigma   )
  #  alpha   1 - confidence probability

  #  OUTPUT
  #  Q.ci.lo )
  #  Q.ci.hi ) confidence interval limits
  # ==========================================================================

  #  pdf at the quantile
  f.Q <- dlnorm(Q,meanlog=mue,sdlog=sigma)

  #  Asymptotic variance of the quantile
  Q.var <- p*(1-p) / (f.Q^2 * n)

  #  Asymptotic standard deviation of the quantile
  Q.sd  <- sqrt(Q.var)

  #  Asymptotic confidene interval
  z       <- qnorm(1-alpha/2)
  Q.ci.lo <- Q - z*Q.sd
  Q.ci.hi <- Q + z*Q.sd

  return(c(Q.ci.lo,Q.ci.hi))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_CIQuant.PNV.R

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Asymptotics 1-alpha confidence interval for quantiles of the 
#  Power normal distribution

#  See R.J.Serfling, Approximation theorems of mathematical statistics,
#  Wiley 1980, S. 77, corollary B
 
#  #  CHANGE HISTORY
#  09.08.2017 Start
# ============================================================================

CIQuant.PNV <- function(n,p,Q,lambda,mue,sigma,alpha=0.05, fastnull=1.e-10)
{ 
  #  INPUT
  #  n       Size of the sample from which to calculate the quantiles
  #  p       Probability associated with the quantile
  #  Q       Quantile
  #  lambda  )
  #  mue     )  Parameters of the power normal distribution
  #  sigma   )
  #  alpha   1 - confidence probability

  #  OUTPUT
  #  Q.ci.lo )
  #  Q.ci.hi ) confidence interval limits
  # ==========================================================================

  #  pdf at the quantile
  f.Q <- pdf.PN(Q,lambda,mue,sigma,fastnull=fastnull)

  #  Asymptotic variance of the quantile
  Q.var <- p*(1-p) / (f.Q^2 * n)

  #  Asymptotic standard deviation of the quantile
  Q.sd  <- sqrt(Q.var)

  #  Asymptotic confidene interval
  z       <- qnorm(1-alpha/2)
  Q.ci.lo <- Q - z*Q.sd
  Q.ci.hi <- Q + z*Q.sd

  return(c(Q.ci.lo,Q.ci.hi))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_CollapseHist.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Collapse a histogram auch that each bin has at least n.per.bin.min
#  members
 
#  #  CHANGE HISTORY
#  04.01.2021 Compressed histogram constructed from components, not via hist()
#  08.10.2020 Call new arranged
#  31.08.2020 Calculation of x.counts.mode.idx moved to call
#  29.08.2020 Call changed
#  28.08.2020 Start
#
# ==============================================================================

CollapseHist <- function(x, x.hist, x.counts.mode.idx, n.per.bin.min,  
                         xlabel, subtitle, histcol, bordercol, dev.figA)
{
  x.hist.counts.n <- length(x.hist$counts) 
  x.n <- sum(x.hist$counts)

  #  Create collapsed / compressed histogram
  #  x.breaks.comp    Bin limits of compressed histogram
  #  x.counts.comp    Counts in the compressed histogram
  #  x.counts.comp.n  Number of bins in the compressed histogram

  x.counts.comp <- c()
  x.breaks.comp <- c()

  #  In x.hist$counts, go down from the mode and collapse bins with 
  #  too few cases

  x.counts.comp.n <- 1
  idx <- x.counts.mode.idx  
  while (idx >= 1)
  { 
    x.counts.comp[x.counts.comp.n] <- x.hist$counts[idx]
    jdx <- idx
    while ( (x.counts.comp[x.counts.comp.n] < n.per.bin.min) & 
            (jdx > 1) )
    { # Count too small, add next bin
      jdx <- jdx - 1
      x.counts.comp[x.counts.comp.n] <- x.counts.comp[x.counts.comp.n] + 
                                        x.hist$counts[jdx]
      # cat(idx,jdx,x.counts.comp[x.counts.comp.n],"\n")
    }
    #  The interval [idx,jdx] either has enough cases or there are not 
    #  sufficiently many intervals
    x.breaks.comp[x.counts.comp.n] <- x.hist$breaks[jdx]  

    #cat("Intervall ",x.counts.comp.n, 
    #    "  Indices",idx,jdx,
    #    "  break li",x.breaks.comp[x.counts.comp.n],
    #    "  Anzahl ",x.counts.comp[x.counts.comp.n], "\n")     

    #  Consider next bin
    x.counts.comp.n   <- x.counts.comp.n + 1
    idx <- jdx - 1
  }
  #  Collapsing went in opposite direction, reverse vectors
  #cat("\n All intervals left of mode inspected, resulte after reversing\n")

  x.breaks.comp <- rev(x.breaks.comp)
  x.counts.comp <- rev(x.counts.comp)

  #  In x.hist$counts, go from first bin to the mode and collapse too small
  #  bins

  idx <- x.counts.mode.idx + 1  
  while (idx <= x.hist.counts.n)
  { 
    x.counts.comp[x.counts.comp.n] <- x.hist$counts[idx]
    jdx <- idx
    while ( (x.counts.comp[x.counts.comp.n] < n.per.bin.min) & 
            (jdx < x.hist.counts.n) )
    { # Anzahl zu klein, nächsten Wert hinzunehmen
      jdx <- jdx + 1
      x.counts.comp[x.counts.comp.n] <- x.counts.comp[x.counts.comp.n] + 
                                        x.hist$counts[jdx]
    }
    #  The interval [idx,jdx] either has enough cases or there are not 
    #  sufficiently many intervals
    x.breaks.comp[x.counts.comp.n] <- x.hist$breaks[idx]  

    #cat("Intervall ",x.counts.comp.n, 
    #    "  Indices",idx,jdx,
    #    "  break li",x.breaks.comp[x.counts.comp.n],
    #    "  Anzahl ",x.counts.comp[x.counts.comp.n], "\n")     

    #  Consider next interval
    x.counts.comp.n   <- x.counts.comp.n + 1
    idx <- jdx + 1
  }
  x.counts.comp.n   <- x.counts.comp.n - 1 

  #  Margin bins may have too few cases, if so, join with the neighbours
  if (x.counts.comp[1] < n.per.bin.min)
  { x.counts.comp[1] <- x.counts.comp[1] + x.counts.comp[2]
    x.counts.comp    <- x.counts.comp[-2]
    x.breaks.comp    <- x.breaks.comp[-2]
    x.counts.comp.n  <- x.counts.comp.n - 1
  }

  #  Note that last break is *left* limit of the last interval
  if (x.counts.comp[x.counts.comp.n] < n.per.bin.min)
  { x.counts.comp[x.counts.comp.n-1] <- x.counts.comp[x.counts.comp.n] + 
                                        x.counts.comp[x.counts.comp.n-1]
    x.counts.comp <- x.counts.comp[-(x.counts.comp.n)]
    x.breaks.comp <- x.breaks.comp[-(x.counts.comp.n)]
    x.counts.comp.n <- x.counts.comp.n - 1
  } 

  # Right limit of last bin is still missing
  x.breaks.comp[x.counts.comp.n+1] <- tail(x.hist$breaks, 1)

  #  Create the full histogram
  x.breaks.comp.n <- length(x.breaks.comp)
  bin.width <- x.breaks.comp[2:x.breaks.comp.n]-
               x.breaks.comp[1:(x.breaks.comp.n-1)]
  x.hist.comp <- list(breaks=x.breaks.comp, 
                  counts=x.counts.comp, 
                  density=x.counts.comp/(x.n*bin.width),
                  mids=(x.breaks.comp[2:x.breaks.comp.n]+
                       x.breaks.comp[1:(x.breaks.comp.n-1)])/2 )
  class(x.hist.comp) <- "histogram"

  # Produce histogram plot, if requested
  if (!is.na(dev.figA))
  { 
    dev.set(dev.figA)
    plot(x.hist.comp, freq=FALSE, col=histcol,border=bordercol,
                      main=paste("Histogram with at least",n.per.bin.min,
                                 "values per bin"),
                      xlab=xlabel, ylab="pdf(x)", 
                      sub=paste(date(), "   ",subtitle),cex.sub=0.7 )
  }

  return(x.hist.comp)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_Compose.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Composes the full parameter vector theta from components theta.est and
#  theta.fix
#  No checks for reasonable vector size and contents for time reasons!
 
#  CHANGE HISTORY
#  12.05.2020 Start
# ==============================================================================

Compose <- function(theta.fix, theta.est, idx.fix, idx.est)
{
  #  INPUT 
  #  theta.fix   the fixed component of the parameter vector. May be NULL.
  #  theta.est   the estimated component of the parameter vector
  #  idx.fix     the indices of the fix component(s) in theta. May be NULL.

  #  OUTPUT 
  #  theta.est   vector with length as idx.est
  # ---------------------------------------------------------------------------

  if (is.null(theta.fix))
  { # No fixed component
    theta <- theta.est
  }  else
  { # Fixed component(s) exist(s)
    theta.n <- length(theta.fix) + length(theta.est)
    theta   <- rep(NA, times=theta.n)
    theta[idx.fix] <- theta.fix
    theta[idx.est] <- theta.est
  } 
  
  return(theta)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test
#t.fix <- NULL
#t.est <- c(10, 20, 30)
#idx.fix <- NULL
#idx.est <- c(1,2,3)
#print(Compose(t.fix, t.est, idx.fix, idx.est))

#t.fix <- 10
#t.est <- c(20, 30)
#idx.fix <- 1
#idx.est <- c(2,3)
#print(Compose(t.fix, t.est, idx.fix, idx.est))


#t.fix <- 10
#t.est <- c(20, 30)
#idx.fix <- 2
#idx.est <- c(1,3)
#print(Compose(t.fix, t.est, idx.fix, idx.est))



# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_Drift.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Generate a plot for drift identification for the monthly median
 
#  CHANGE HISTORY
#  26.03.2021 Start
# ==============================================================================

Drift <- function(date, x, RL1, RL2, figA, xlabel, titletext, alpha=0.05, 
                  gamcol="lightblue",
   subt="The blue confidence range should overlap the grey range at all times")
{
  #  INPUT
  #  date  date, format  "2021-03-26" generated by
  #        as.Date("26.03.2021", format = "%d.%m.%Y"),
  #  x     time series data
  #  RL1   Lower 95% reference limit ( 2.5%) 
  #  RL2   Upper 95% reference limit (97.5%) 
  #  figA  plot window
  #  xlabel label of x
  #  titletext text for plot title

  #  OUTPUT 
  # ---------------------------------------------------------------------------

  #  Median per unique time point

  x.median <- aggregate(x, by=list(date), median)
  colnames(x.median) <- c("date", "median")
 
  #  Fit a gam to the monthly medians, if enough points
  if (nrow(x.median) >= 12)
  {
    #cat("\n Drift estimation by GAM: median = time + s(time)\n")

    #  To simplify writing:
    time  <- as.numeric(x.median[ , "date"])
    x.med <- x.median[ , "median"]

    x.med.gam <- gam(x.med ~ time + s(time))
    #print(x.med.gam)
    #print(summary(x.med.gam))    
  
    x.med.gam.pred <- predict.gam(x.med.gam,se.fit=TRUE)
    # print(x.med.gam.pred)

    x.med.gam.pred.fit <- unlist(x.med.gam.pred$fit)
    x.med.gam.pred.se  <- unlist(x.med.gam.pred$se.fit)
    zalpha             <- qnorm(1-alpha/2)
    x.med.gam.pred.ciu <- x.med.gam.pred.fit - zalpha*x.med.gam.pred.se
    x.med.gam.pred.cio <- x.med.gam.pred.fit + zalpha*x.med.gam.pred.se

    # print(cbind(x.med,x.med.gam.pred.ciu,x.med.gam.pred.cio))

    #  Total median (of all data)
    #x.median.total <- median(x)

    x.median.total <- median(x.med)
    #  Permissible deviation around x.median.total
    x.median.total.pI <- PermissibleDiff(RL1, RL2, x.median.total)

    #  Simplify writing
    xmt.lo <- x.median.total.pI[1, "xi.lo"]
    xmt.hi <- x.median.total.pI[1, "xi.hi"]

    yplotmin <- min(c(x.med, xmt.lo, x.med.gam.pred.ciu)) 
    yplotmax <- max(c(x.med, xmt.hi, x.med.gam.pred.cio)) 

    # Span the plot window
    dev.set(figA)
    #bringToTop(figA)

    plot(c(x.median[1 ,"date"], tail(x.median[ ,"date"], 1)),
         c(yplotmin, yplotmax), 
         type="p", pch=3, col="transparent",
         xlab="Date", ylab=xlabel, main=titletext,
         sub=subt,
         cex.sub=0.8) 

    #  Permissible range of the total median
    polygon(c(x.median[1, "date"], tail(x.median[ , "date"], 1),
              tail(x.median[ , "date"], 1),  x.median[1, "date"] ), 
            c(xmt.lo, xmt.lo, xmt.hi, xmt.hi),
            col="gray95", border="gray95") 
    
    #  CI of drift
    polygon(c(time, rev(time)),
            c(x.med.gam.pred.ciu, rev(x.med.gam.pred.cio)), 
            col=gamcol, border="lightblue") 
    
    #  Monthly medians   
    points(x.median[ ,"date"], x.median[ ,"median"], type="p", col="black",
           pch=3)

    #   Total median
    abline(h=x.median.total, col="red")

    #   Total median
    abline(h=c(xmt.lo, xmt.hi), col="red", lty=2)

    #  GAM fit to monthly medians, + CI
    lines(time, x.med.gam.pred.fit, col="blue",lty=1)
    lines(time, x.med.gam.pred.ciu, col="blue",lty=2)
    lines(time, x.med.gam.pred.cio, col="blue",lty=2)

    #  Check for drift
    drift.likely <- any(x.med.gam.pred.ciu > xmt.hi) |
                    any(x.med.gam.pred.cio < xmt.lo)

  }   else
  {  #  No data for drift analysis
     drift.likely <- NA
  }

  return(drift.likely)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test: see Test_DriftAnalysis.R

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_EstParentDist.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Estimate the parent normal distribution if only limits, mean and sd of a 
#  truncated subset are known. See @@@ for a method description. 

#  CHANGE HISTORY
#  Changes in functions are documented in the functions!

#  02.12.2020 RLs, CI(RL) addidionally calculated here, call changed
#  28.09.2020 ok range changed to { , )
#  08.09.2020 tmu installed, execution via do.tmu.  
#  30.03.2020 eps removed from parameter list
#  17.09.2019 p.fact as additional parameter
#  15.03.2019 Installation from directory ProgWW17_RH

# ============================================================================

EstParentDist <- function(x.tr.lo, x.tr.hi, x.tr, theta.parent.ini,
                          alpha, lambda, RL1.p, RL2.p)
{ 
  # x.tr.lo            lower truncation limit
  # x.tr.hi            upper truncation limit
  # x.tr               truncated data set
  # theta.parent.ini   initial guess for parent parameters

  # ==========================================================================

  # --------------------------------------------------------------------------
  # if (print.log.message) { cat("%%%   F_EstParentDist  Start\n") }
  # --------------------------------------------------------------------------

  #  Estimates of mue and sigma from the truncated, to ND transformed 
  #  interval
  x.tr.n   <- length(x.tr)
  x.tr.mea <- mean(x.tr)
  x.tr.sd  <- sd(x.tr)
 
  #  Improve parent parameter estimation iteratively  

  TMU.nlm <- nlm(CalcParDistance, theta.parent.ini, 
                 theta.trunc=c(x.tr.mea, x.tr.sd), 
                 a=x.tr.lo, b=x.tr.hi)
    
  #  Modified estimates for mue and sigma
  tab.tmu <- c(mue.tmu=TMU.nlm$estimate[1], sigma.tmu=TMU.nlm$estimate[2], 
              distance.tmu=TMU.nlm$minimum, rc.tmu=TMU.nlm$code, 
              iter.tmu=TMU.nlm$iterations)

  #  Meaning of "code" ("rc"):
  #  1: relative gradient is close to zero, current iterate is 
  #     probably solution.
  #  2: successive iterates within tolerance, current iterate is 
  #     probably solution.
  #  3: last global step failed to locate a point lower than estimate. 
  #     Either estimate is an approximate local minimum of the function 
  #     or steptol is too small.
  #  4: iteration limit exceeded.
  #  5: maximum step size stepmax exceeded five consecutive times. 
  #     Either the function is unbounded below, becomes asymptotic 
  #     to a finite value from above in some direction or stepmax is too small.

  # -------------------------------------------------------------------------
  #  Calculate RLs
  x.RL1.tmu <- unname(q.PN(RL1.p, lambda, tab.tmu["mue.tmu"], 
                          tab.tmu["sigma.tmu"]))
  x.RL2.tmu <- unname(q.PN(RL2.p, lambda, tab.tmu["mue.tmu"], 
                          tab.tmu["sigma.tmu"]))
  tab.tmu <- c(tab.tmu, x.RL1.tmu=x.RL1.tmu, x.RL2.tmu=x.RL2.tmu)

  # -------------------------------------------------------------------------
  #  Add asymptotic confidence intervals
  xc.RL1.tmu.CI <- CIQuant.PNV(x.tr.n, RL1.p, 
                               tab.tmu["x.RL1.tmu"],
                               lambda, 
                               tab.tmu["mue.tmu"], 
                               tab.tmu["sigma.tmu"],
                               alpha=alpha)
  tab.tmu <- c(tab.tmu, x.RL1.cilo=unname(xc.RL1.tmu.CI[1]),
                        x.RL1.cihi=unname(xc.RL1.tmu.CI[2]))

  xc.RL2.tmu.CI <- CIQuant.PNV(x.tr.n, RL2.p, 
                               tab.tmu["x.RL2.tmu"],
                               lambda, 
                               tab.tmu["mue.tmu"], 
                               tab.tmu["sigma.tmu"],
                               alpha=alpha)
  tab.tmu <- c(tab.tmu, x.RL2.cilo=unname(xc.RL2.tmu.CI[1]),
                        x.RL2.cihi=unname(xc.RL2.tmu.CI[2]))

  # --------------------------------------------------------------------------
  # if (print.log.message) { cat("%%%   F_EstParentDist  End\n") }
  # --------------------------------------------------------------------------

  return(tab.tmu)
}   
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_Expand.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Expand tied data (due to rounding) to fictive individual data. Counts per
#  bin are as in input, position in bin is log-equidistant, with half distance
#  to the margins
 
#  #  CHANGE HISTORY
#  16.02.2021 Special treatment for count = 1 added
#  06.01.2021 Step function as 2nd possibility of expansion added
#  30.12.2020 Start
#
# ==============================================================================

Expand <- function(breaks, counts, step=TRUE)
{
  #  INPUT 
  #  breaks          Jump positions of the step functions
  #                  Must be > 0 
  #                  Reasonable number of breaks must be checked outside 
  #                  this function
  #  counts          # of observations per bin
  #  step            TRUE: x.fict in bin uniformly distributed
  #                  FALSE: x.fict in bin logarithmically distributed
 
  #  OUTPUT 
  #  x.fict          fictive data points, counts per bin as in input,
  #                  position in bin depends on step
  # ---------------------------------------------------------------------------

  #  Allocate fictive cases in each (original) bin
  #  x positions for all data
  bins.n <- length(counts)

  x.fict.exists <- FALSE             # in place of exists()

  if (step)
  { 
    # Use a simple step function
    for (ibin in 1:bins.n)
    {
      if (counts[ibin] > 0)
      { 
        #  Avoid numerical trouble due to very small differences if count=1
         
        if (counts[ibin] == 1)
        {
          x.fict0 <-  (breaks[ibin+1] + breaks[ibin]) / 2
        } else
        {   
          dx <- (breaks[ibin+1] - breaks[ibin]) / counts[ibin]
          x.fict0 <- seq(breaks[ibin]+dx/2, breaks[ibin+1]-dx/2, by=dx)
        }

        if (!x.fict.exists)
        { 
          x.fict        <- x.fict0
          x.fict.exists <- TRUE 
        } else
        {
          x.fict <- c(x.fict, x.fict0)   
        }  
      }
    }
  }   else
  { # Locate points at logarithmic positions  

    for (ilo in 1:(bins.n))
    { ihi <- ilo + 1
      if (counts[ilo] > 0)
      { 
        delta.logx <- (log(breaks[ihi]) - log(breaks[ilo]))/ (counts[ilo]+1)  
        logx.fict  <- seq(log(breaks[ilo])+delta.logx/2,
                        log(breaks[ihi])-delta.logx/2, length.out=counts[ilo])

        if (!x.fict.exists)
        { 
          x.fict        <- exp(logx.fict)
          x.fict.exists <- TRUE 
        } else
        {
          x.fict <- c(x.fict, exp(logx.fict))   
        }
      }
    }
  } 

  return(x.fict) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test
# see also Test.Expand.R

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_Extract.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Extract the component theta.est from theta
#  No checks for reasonable vector size and contents for time reasons!
 
#  CHANGE HISTORY
#  12.05.2020 Start
# ==============================================================================

Extract <- function(theta, idx.est)
{
  #  INPUT 
  #  theta   the full parameter vector
  #  idx.est the indices of components that are estimated (not fix)

  #  OUTPUT 
  #  theta.est   vector with length as idx.est
  # ---------------------------------------------------------------------------

  theta.est <- theta[idx.est]

  return(theta.est)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test
#t <- c(10, 20, 30)
#idx.est <- c(1,2,3)
#print(Extract(t, idx.est))
#idx.est <- c(1,3)
#print(Extract(t, idx.est))
#idx.est <- 2
#print(Extract(t, idx.est))
#idx.est <- NA
#print(Extract(t, idx.est))
#idx.est <- 4
#print(Extract(t, idx.est))

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_FC.R

#  Truncated minimum chi square estimation

# (c) wwosniok@math.uni-bremen.de

#  Print values using format specifications as in formatC. 
#  formatC itself does not in all cases deal with NA correctly.

# x  <- c(123, 123.4, 123.45, NA)
# y  <- c(NA, NA, NA)
# formatC(x, width=8, digits=2)  
# " 1.2e+02" " 1.2e+02" " 1.2e+02" "      NA"
# formatC(x, format="f", width=8, digits=2)  
# "  123.00" "  123.40" "  123.45" "      NA"
# formatC(y, format="f", width=8, digits=2)
#Fehler in formatC(y, format = "f", width = 8, digits = 2) : 
#  unsupported type 

# formatC(as.numeric(y), format="f", width=8, digits=2)

#  CHANGE HISTORY
#  11.04.2021 Start
# ============================================================================ 

FC <- function(x, w, d)
{ 
  if (is.na(w)) 
  { Fx <- formatC(as.numeric(x), format="f", digits=d, flag="#") } else
  { Fx <- formatC(as.numeric(x), format="f", width=w, digits=d, flag="#") }

  return(Fx)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# x  <- c(123, 123.4, 123.45, NA)
# y  <- c(NA, NA, NA)
#FC(x,  8, 2)
#FC(x, NA, 2)
#FC(y, NA, 2)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv # *****************************************************************************
#  F_FilterWDay.R
#
#  Extract data from given weekdays
#
#  (c) wwosniok@math.uni-bremen.de
#
#  CHANGE HISTORY
#
#  29.10.2020 Start
#
# --------------------------------------------------------------------------
# INPUT
# wday      vector containing the weekdays of observations as numbers
#          (Monday=0)
# use.wday  vector containing the weekdays to select, 
#           either NA or "all" (==> use all data)
#           or a subset of c("mo", "tu", "we", "th", "fr", "sa", "su")

# OUTPUT
# ok        vector (length as wday) indicating the selected components of wday
# --------------------------------------------------------------------------

FilterWDay <- function(wday, use.wday)
{  
  #  If selection of weekdays is required: do it
  if (!is.na(use.wday) && use.wday!="all")
  { 
    ok <- rep(FALSE, times=length(wday))

    #  Check for valid entries

    for (day in use.wday)
    { 
      # print(day)

      idx <- which(tolower(day) == c("mo", "tu", "we", "th", "fr", "sa", "su"))
      # print(idx) 
      if (length(idx) == 1)
      { 
        ok  <- ok | (wday == (idx-1))
      }  else
      { cat("\n ++++++ Weekday selection ", day, " invalid  ++++++\n") } 
    }
  }  else
  {
    # No selection required
    ok <- rep(TRUE, times=length(wday))   
  }

  #  The calling programme has to handle the situation of no data left
  return(ok)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

#wday <- c(0, 1,1, 2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5,5, 6,6,6,6,6,6)

#use.wday <- c("mo", "su", "WE")

#select <- FilterWDay(wday, use.wday)
#print(data.frame(wday, select))

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_FindBw_V4.R

FindBw <- function(x, kernel.n, round.unit, FigA, subtitle, kdecol,
                    min.bw=0.75, b.fact=1.05, inner.lo=0.025, inner.hi=0.975)
{  
  #  Adjust bandwidth of kernel density estimate such that the
  #  inner part of the density has only 1 maximum
  #
  #  INPUT
  #  x           data
  #  round.unit  rounding applied to the data
  #  kernel.n    number of support points for the kernel 
  #  FigA        window for plotting. Must exist. NA: no plot
  #  b.fact      factor by which bw.adjust is increased per iteration
  #  inner.lo    lower limit of the inner range (quantile of total data)
  #  inner.hi    upper limit of the inner range (quantile of total data)

  #  05.04.2021 temp[ , "lo"] <=  temp[ , "mid"]) allowed for maximum
  #  01.12.2020 Missing new calculation of density() added in loop 
  #  26.09.2020 Colours as function input
  #  15.09.2020 Bandwidth must be >= min.bw*round.unit
  #  13.09.2020 Additional condition: bandwidth must be >= round.unit
  #  22.05.2020 Target changed: smoothing is increased until all maxima
  #             are more than the rounding unit apart (or only 1 maximum left)
  #  21.05.2020 Target changed: up to max.n.acc maxima are accepted, the highest
  #             is taken as mode position
  #  20.05.2020 Start
  # ==========================================================================

  x.kde.mode <- NA
  bw.adjust  <- 1

  #  First try. Used to determine the support and to ensure that 
  #  bw > round.unit
  #  First guess. May have negative support.
  x.kde <- density(x, n=kernel.n, adjust=bw.adjust)

  #  Support points must be > 0 to avoid trouble with BC transformation
  x.from <- min(x.kde$x[x.kde$x>0])
  x.to   <- max(x.kde$x)

  x.kde <- density(x, n=kernel.n, adjust=bw.adjust, from=x.from,to=x.to)

  if (x.kde$bw < min.bw*round.unit)
  { #  Force bandwidth to be >= min.bw*round.unit
    bw.adjust <- min.bw*round.unit / x.kde$bw
    x.kde <- density(x, n=kernel.n, adjust=bw.adjust, from=x.from,to=x.to)
  }

  if (!is.na(FigA))
  { 
    dev.set(FigA)
    yplotmax <- max(x.kde$y)
    plot(x.kde, type="l", col="gray", ylim=c(0, yplotmax),
         main=paste("KDE with bandwidth adjustment, bw =", x.kde$bw, 
                    "  adj =", format(bw.adjust, digits=3)),          
         sub=subtitle, cex.sub=0.7)
    axis(1, at=seq(130, 160, by=round.unit), labels=FALSE)
  }

  while (is.na(x.kde.mode))
  {
    x.kde <- density(x, n=kernel.n, adjust=bw.adjust, from=x.from,to=x.to)

    if (!is.na(FigA))
    { 
      lines(x.kde, type="l", col="gray")
    }

    #  Check only the inner range of the density for maxima
    x.kde.cdf <- cumsum(x.kde$y)
    x.kde.cdf <- x.kde.cdf/tail(x.kde.cdf,1)
    Q1Q3      <- (inner.lo <= x.kde.cdf) & (x.kde.cdf < inner.hi)

    x.kde.Q1Q3.x <- x.kde$x[Q1Q3]
    x.kde.Q1Q3.y <- x.kde$y[Q1Q3]
    Q1Q3.n <- sum(Q1Q3)

    if (!is.na(FigA))
    { 
      abline(v=c(x.kde.Q1Q3.x[1], tail(x.kde.Q1Q3.x, 1)), col="gray")
    }

    #  Find maxima
    temp <- cbind(x=x.kde.Q1Q3.x[2:(Q1Q3.n-1)],
                  lo=x.kde.Q1Q3.y[1:(Q1Q3.n-2)],
                  mid=x.kde.Q1Q3.y[2:(Q1Q3.n-1)],
                  hi=x.kde.Q1Q3.y[3:(Q1Q3.n  )])
    is.max <- (temp[ , "lo"] <=  temp[ , "mid"]) & 
              (temp[ , "mid"]>   temp[ , "hi"])
    temp <- cbind(temp, is.max=is.max)

    max.n <- sum(is.max)
    # cat("\n bw.adjust: ", bw.adjust, " number of maxima: ", max.n)

    if (max.n == 0)
    { cat("\n [Find.bw] No maximum found\n")
      print(temp)
      stop("++++++ NO maximum found ++++++")
    }

    if (!is.na(FigA) & max.n > 0)
    { 
      abline(v=temp[is.max, "x"], col="orange", lty=2) 
    }

    #  If > 1 maximum: may be artificial due to strong rounding. If
    #  any distance between maximum positions ~~ round.unit, increase 
    #  bw.adjust, otherwise the solution is reached.
    #  Calculate pairwise distance between maximum positions. Slow method,
    #  but there should not be too may positions.

    if (max.n > 1)
    {
      temp.max <-  temp[is.max, ] 
      #cat("\n More than 1 maximum, maxima positions\n")
      #print(temp.max)

      x.diff <- abs(diff(temp.max[ , "x"]))

      #cat("\n Absolute maxima position differences \n")
      #print(x.diff)

      # Is there any distance ~~ round.unit outside the diagonal?
      d.eq.ru <- (0.55*round.unit <= x.diff) & (x.diff <= 1.45*round.unit)

      if (any(d.eq.ru))
      {  
        # at least 1 pair of maxima with distance ~~ r.u., increase bandwidth
        bw.adjust <- bw.adjust * b.fact
        #  do another round. x.kde.mode is still NA             
      } else   
      { #  Several maxima left, all distances > round.unit. 
        #  Take the one with max density
        #cat("\n\n")

        mode.idx <- which.max(temp[ ,"mid"])
        x.kde.mode <- temp[mode.idx, "x"]  #  this terminates while()

        #cat("\n Local maxima\n")
        #print(temp[is.max, ])
        #cat("\n Mode: ", x.kde.mode, "\n")
      } 
    }  else    
    { 
      # only 1 maximum left, calculate mode
      x.kde.mode <- temp[is.max, "x"]
      #cat("\n Mode: ", x.kde.mode, "\n")
    }
  }    #  while ...

  if (!is.na(FigA))  
  { 
    lines(x.kde, col=kdecol)
    abline(v=x.kde.mode, col="firebrick") 
    legend("topright", c("kde", "Q1, Q3", "mode"), 
           col=c(kdecol, "gray", "firebrick"), lty=1, lwd=2, cex=0.8)   
    abline(v=x.kde.mode, col="firebrick")
  }

  return(list(bw.adjust = bw.adjust,
              x.kde.mode=x.kde.mode,
              x.kde=x.kde))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  See Test_Bandwidth.R 

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_FindHisto_V2.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Find a histogram containing a specified number of bins with specified 
#  minimum size per bin
#  Generate also a reduced version and two smoothed versions
 
#  #  CHANGE HISTORY
#  09.03.2021 Treatment of values < DL extended to second break position
#  23.02.2021 Treatment of values < DL added
#  15.02.2021 bins.n2 corrected to integer
#  11.02.2021 Derivation of breaks corrected
#  07.02.2021 Start with new approach: eqd is determined by Freedman-Diaconis,
#             x.hist gets its breaks from quantiles of the kde-cdf
#
# ==============================================================================

FindHisto <- function(x, round.unit, detect.limits.max, x.table, x.val, x.kde,
                      n.per.bin.min, bins.n.min, bins.n.max,
                      smooth.hist1, smooth.hist2, red.limit, xlabel, subtitle, 
                      bordercol1, histcol1, bordercol2, histcol2, 
                      bordercol3, histcol3, kdecol, figB)
{
  #  INPUT 
  #  x              data vector
  #  round.unit
  #  x.kde          kde of x
  #  x.kde.mode     mode of x.kde
  #  n.per.bin.min  min number of per bin
  #  red.limit      create a reduced data set if original size is > re.limit 
  #  xlabel         label for analyte 
  #  subtitle       subtitle in FigA 
  #  bordercol1     colour for histogram borders, initial
  #  histcol1       colour for bin area, initial
  #  bordercol2     colour for histogram borders, collapsed
  #  histcol2       colour for bin area, collapsed
  #  bordercol2     colour for histogram borders, collapsed and smoothed
  #  histcol2       colour for bin area, collapsed and smoothed
  #  kdecol         colour for kde
  #  figB           initial and collapsed histogram
  #  figC           collapsed and smoothed histogram

  #  OUTPUT 
  #  
  # ---------------------------------------------------------------------------
  
  x.n     <- length(x)  
  x.val.n <- length(x.val)
  if (is.na(detect.limits.max))
  { x.lt.DL.n <- 0 }  else
  { x.lt.DL.n <- sum(x < detect.limits.max )}

  x.min <- min(x)
  x.max <- max(x)

  # ---------------------------------------------------------------------------
  #  Produce a set of reasonable histograms on the original scale
  # --------------------------------------------------------------------------
  #  Calculate equidistant histogram, standard approach using the 
  #  Freedman-Diaconis approach (which, like all others is not suitable for 
  #  rounded data)
  #  Parameter 'freq' not used here.

  #  All histograms have the same first and nearl the same last breaks. 
  #  Note that the data value zero is set to 0.1*round.unit in seg070. 
  #  first: max(x.val[1] - round.unit/2, 0.05*round.unit)
  #  last: for equidistant breaks: such that max(x) is just included
  #        for non-equidistant breaks: max(x) + round/2
  #  Generate equidistant histogram

  x.hist.eqd <- hist(x, breaks="FD", right=FALSE, plot=FALSE)

  #  'FD' produces 'nice' histograms with coinciding limits and values.
  #  Gives a misleading picture if a density estimate is plotted as well. 
  #  Therefore, shift the equidistant breaks by -round.unit/2.
  FD.breaks <- x.hist.eqd$breaks
  FD.breaks <- FD.breaks - round.unit/2 
  h         <- FD.breaks[2] - FD.breaks[1]

  #  The bin width must not be smaller than round.unit
  if (h < round.unit)
  { 
    h <- round.unit
  }

  #  The bin width must be an integer multiple of round.unit
  if (abs(h %% round.unit) > 1.e-5)
  { 
    #  Use the next lower multiple of round.unit 
    h <- Floor(h, round.unit)
  }

  #  Now h is large enough, but may be too large. 

  # ..........................................................................
  #  Equidistant histogram 
  breaks.start <-  max(x.val[1] - round.unit/2, 0.05*round.unit)
  breaks.eqd.n <-  Ceiling((x.val[x.val.n] - breaks.start)/h, 1)
  breaks.end   <-  breaks.start + breaks.eqd.n * h

  breaks.eqd <- seq(breaks.start, breaks.end, by=h)

  #  Generate an equidistant histogram with the (modified) breaks
  x.hist.eqd <- hist(x, breaks=breaks.eqd, right=FALSE, plot=FALSE)

  # -------------------------------------------------------------------------  
  #  Smooth the equidistant histogram, if requested
  if ((!smooth.hist1) & (!smooth.hist2))
  { # No smoothing
    x.hist.eqd.smo <- x.hist.eqd
  }

  if (( smooth.hist1) & (!smooth.hist2))
  { # Smoothing by moving average
    x.hist.eqd.smo <- SmoothHist1(x.hist.eqd, figA=NA )
  } 

  if ((!smooth.hist1) & ( smooth.hist2))
  { # Smoothing by kde fit
    x.hist.eqd.smo <- SmoothHist2(x.kde, x.hist.eqd, NA, 
                                  histcol2, bordercol2, 
                                  histcol3, bordercol3, kdecol, subtitle)
  } 

  # --------------------------------------------------------------------------
  #  Calculate (non-equidistant) histogram with
  #  bin frequencies > n.per.bin.min and further conditions:
  #  n.per.bin >= n.per.bin.min
  #  bins.n.min <= bins.n <= bins.n.max 
  #  max possible bins.n is x.val.n
  #  First guess for bins.n: 
  x.hist.ok  <- TRUE  
  n.per.bin1 <- n.per.bin.min

  #  Find the maximal possible and accepted number of bins
  bins.n1 <- min(x.val.n, bins.n.max) 
  bins.n2 <- NA 

  #  If bins.n1 < bins.n.min: write a warning and give up
  if (bins.n1 < bins.n.min)
  { cat("\n +++ [FindHisto] Default histogram construction impossible",
        "\n                 because of too few distinct values in the data",
        "\n                 # of values          ", x.n,
        "\n                 # of distinct values ", x.val.n,
        "\n                 Max. detection limit ", detect.limits.max,
        "\n                 % of values < DL     ", 
                              formatC(100*x.lt.DL.n / x.n, width=5, digits=1), 
        "\n                 Min. # of bins       ", bins.n.min,
        "\n                 Max. # of bins       ", bins.n.max,
        "\n                 Min. size of a bin   ", n.per.bin1,
        "\n\n")

    x.hist.ok <- FALSE  
  } 

  if (x.hist.ok)
  {
    #  Find the number of bins resulting from sample size and minimum bin size
    bins.n2 <- Round(x.n/n.per.bin.min, 1)

    if (bins.n2 < bins.n.min)
    { #  + 2 in the denominator to guarantee required total number
      n.per.bin1 <- Floor(x.n / (bins.n.min+2), 1)

      #  Resulting corrected # of bins
      bins.n2 <- Ceiling(x.n / n.per.bin1, 1)
     
      cat("\n +++ [FindHisto] Default histogram construction impossible",
          "\n                 because of too few values in the data",
          "\n                 # of values          ", x.n,
          "\n                 # of distinct values ", x.val.n,
          "\n                 Max. detection limit ", detect.limits.max,
          "\n                 % of values < DL     ", 
                              formatC(100*x.lt.DL.n / x.n, width=5, digits=1), 
          "\n                 Min. # of bins       ", bins.n.min,
          "\n                 Max. # of bins       ", bins.n.max,
          "\n                 Min. size of a bin   ", n.per.bin.min,
          "\n                 reduced to           ", n.per.bin1,
          "\n                 leading to bins.n2 = ", bins.n2,
          "\n\n")

      #  @@@   decrease n.per.bin.min?
    }
  } 

  #  Only the smaller number can be realized

  bins.n <- min(bins.n1, bins.n2)

  #  Check if histogram construction is now acceptable 
  x.hist.ok <- x.hist.ok & (bins.n.min <= bins.n) & (bins.n <= bins.n.max)  

  #  Error handling: too few bins?
  if (x.hist.ok)
  {
    #  Probabilities for the first-guess quantiles
    probs   <- (1:(bins.n-1)) / bins.n
    probs.n <- length(probs)

    #  Desired minimum probability per bin
    p.per.bin <- probs[2] - probs[1]

    #  Basic candidates for breaks are the midpoints between (real) values
    x.hist.breaks0 <- (x.val[1:(x.val.n-1)] + x.val[2:x.val.n]) / 2
    
    #  First position depends on existence of < DL values
    #  detect.limits.max is largest detection limit
    if (is.na(detect.limits.max))
    {
      #  No detection limit
      breaks.start1 <- breaks.start
      x.hist.breaks0 <- c(breaks.start1, x.hist.breaks0, breaks.end)

    }  else
    {
      #  Detection limit present
      #  First break is (the replacement of) zero  
      breaks.start1 <- 0.1 * round.unit  

      # Second break is max. DL. There is no point in having a break < DL 
      # (except breaks.start1 above)
      breaks.start2 <- detect.limits.max 

      # If there are break candidates <= breaks.start2 in x.hist.breaks0:
      #  remove
      x.hist.breaks0 <- x.hist.breaks0[x.hist.breaks0 > breaks.start2]
      x.hist.breaks0 <- c(breaks.start1, breaks.start2, x.hist.breaks0, 
                          breaks.end)
    }

    #  Iterative development of break positions such that size requirements
    #  hold

    #  Find quantiles of the empirical distribution function that reproduce
    #  the wanted probabilities as good as possible. Resulting breaks must
    #  contain bins with at least n.per.bin.min values. 
    #  Calculate the empirical cdf, given the break candidates fixed so far 

    #  This cdf uses unsmoothed data
    x.hist0 <- hist(x, breaks=x.hist.breaks0, right=FALSE, plot=FALSE)

    x.cdf   <- unname(cumsum(x.hist0$counts)) 
    x.cdf   <- x.cdf / tail(x.cdf, 1)
    x.cdf   <- c(0, x.cdf)

    #  Realizable breaks may not lie at theoretically desired positions.
    #  Construct breaks such that the minimum count per bin is warranted 
    x.hist.breaks1 <- rep(NA, times=probs.n)
    x.hist.breaks  <- rep(NA, times=probs.n)

    i <- 0
    p.covered <- 0

    while(p.covered <= probs[probs.n-1] )
    { i <- i + 1
      #  The theoretical position
      p.target <- p.covered + p.per.bin
      x.hist.breaks1[i] <- Interpolx(x.hist.breaks0, x.cdf, 
                                   ytarget=p.target)

      #  The next possible position
      idx <- which(x.hist.breaks0 > x.hist.breaks1[i])[1]
      x.hist.breaks[i] <- x.hist.breaks0[idx]
      p.covered        <- x.cdf[idx] 
    }

    x.hist.breaks <- x.hist.breaks[!is.na(x.hist.breaks)]

    #  Add extreme limits
    x.hist.breaks <- c(breaks.start1, x.hist.breaks, breaks.end)

    #  Due to unfavourable bin frequencies, the # of bins may have become
    #  to small - then no histogram

    bins.n.act <- length(x.hist.breaks) - 1

    if (bins.n.act < bins.n.min)
    {
      #  Not enough bins
      cat("\n +++ [FindHisto] Histogram construction impossible",
          "\n                 because of unfavourable bin frequencies",
          "\n                 # of values          ", x.n,
          "\n                 # of distinct values ", x.val.n,
          "\n                 Max. detection limit ", detect.limits.max,
          "\n                 % of values < DL     ", 
                              formatC(100*x.lt.DL.n / x.n, width=5, digits=1), 
          "\n                 Min. # of bins       ", bins.n.min,
          "\n                 Max. # of bins       ", bins.n.max,
          "\n                 Min. size of a bin   ", n.per.bin.min,
          "\n                 Theoretical # of bins", bins.n2,
          "\n                 Available # of bins  ", bins.n.act,
          "\n\n")
      x.hist <- NA
    }  else    
    {
      # No of bins and bin sizes are ok, generate histogram

      x.hist <- hist(x, breaks=x.hist.breaks, right=FALSE, plot=FALSE) 
  
      # ----------------------------------------------------------------------  
      #  Smooth the (collapsed) non-equidistant histogram, if requested

      if (( smooth.hist1) & (!smooth.hist2))
      { # Smoothing by moving average
        x.hist <- SmoothHist1(x.hist, figA=NA )
      } 

      if ((!smooth.hist1) & ( smooth.hist2))
      { # Smoothing by kde fit
        x.hist <- SmoothHist2(x.kde, x.hist, NA, 
                                  histcol2, bordercol2, 
                                  histcol3, bordercol3, kdecol, subtitle)
      } 

    }   
  }   else
  {
    #  Not enough bins alsready before looking at sizes
    x.hist <- NA
  }

  # --------------------------------------------------------------------------
  #  For the QQ plot analysis:
  #  Create a dataset containing fictive individual values from the original
  #  tied data / the equidistant histogram. 
  #  Reduce this dataset to a smaller set of quantiles, if 
  #  the original size > red.limit

  #  Resolve ties
  x.reso <- Expand(x.hist.eqd$breaks, x.hist.eqd$counts)

  if (x.n <= red.limit) 
  { 
    #  No reduction needed, take the full versions as reduced
    x.reso.red    <- x.reso

  } else
  {
    # ------------------------------------------------------------------------ 
    #  Reduction needed
   
    qqw.quant.seq <- seq(0.010, 0.990, by=0.01)
    x.reso.red    <- Quantile(x.reso, probs=qqw.quant.seq)
  }

  # --------------------------------------------------------------------------
  #  Plot the initial histogram (h is already adjusted) and the collapsed
  #  histogram

  if (!is.na(figB) & x.hist.ok)
  { dev.set(figB)
    #bringToTop(figB)

    par(mfcol=c(3,1))

    plot(x.hist,
         border=bordercol1, col=histcol1, freq=FALSE, lty=1,
         xlim=c(0,60), 
         main="Histogram for TMC",
         xlab=xlabel,
         ylab="pdf(y)",
         sub=subtitle, cex.sub=0.7)
    axis(1, at=seq(0, 69, by=2.5), labels=FALSE)

    plot(x.hist.eqd.smo,
         border=bordercol2, col=histcol2, freq=FALSE, lty=1,
         xlim=c(0,60), 
         main="Equidistant histogram, smoothed",
         xlab=xlabel,
         ylab="pdf(y)",
         sub=subtitle, cex.sub=0.7)
    axis(1, at=seq(0, 69, by=2.5), labels=FALSE)

    plot(x.hist.eqd,
         border=bordercol3, col=histcol3, freq=FALSE, lty=1,
         xlim=c(0,60), 
         main="Equidistant histogram, not smoothed",
         xlab=xlabel,
         ylab="pdf(y)",
         sub=subtitle, cex.sub=0.7)
    axis(1, at=seq(0, 69, by=2.5), labels=FALSE)

    legend("topright", c("x.hist", "x.hist.eqd.smo",  "x.hist.eqd"),
           col=c(bordercol1, bordercol2, bordercol3), lwd=3, cex=0.8)
    legend("right", c("x.hist", "x.hist.eqd.smo",  "x.hist.eqd"),
           fill=c(histcol1, histcol2, histcol3), cex=0.8)
  }

  # --------------------------------------------------------------------------

  # x.hist.eqd                            # equidistant, no smooth
  # x.hist.eqd.smo                        # equidistant, smoothed, if required
  # x.hist         <-                     # collapsed, == x.hist.eqd,
  #                                       # smoothed if required. Type of 
  #                                       # smooth given by
  #                                       # smooth.hist1 and smooth.hist2
  # x.reso.red                            # data, reduced size if needed,
  #                                       # ties resolved 

 
  return(list(x.hist.eqd=x.hist.eqd,
              x.hist.eqd.smo=x.hist.eqd.smo, 
              x.hist=x.hist,  
              x.reso.red=x.reso.red)) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#
#  Test: see Test_FindHisto2.R
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_FindKDE.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Find a reasonable kernel density estimate for rounded data. 
#  Only the inner interval (Q.lo, Q.hi) is approximated by kde, to avoid the 
#  need of extremely many support points due to extreme observations
#  Problematic: Q.lo if many values < DL. Choice of surrogate value is 
#  important here.
 
#  #  CHANGE HISTORY
#  08.10.2020 Start
#
# ============================================================================

FindKDE <- function(x, round.unit, kernel.n, 
                    figA, subtitle, kdecol,
                    Q.lo=0.025, Q.hi=0.975, b.fact=1.05)
{
  #  INPUT 
  #  x           data vector
  #  round.unit  unit used when rounding the data
  #  kernel.n    number of support points for kde
  #  figA        device for figure illustrating the steps of finding a 
  #              bandwidth. NA: no plot.
  #  subtitle    subtitle in FigA 
  #  kdecol      color for kde line in FigA 

  #  OUTPUT 
  #  
  # --------------------------------------------------------------------------

  KDE <- FindBw(x, kernel.n, round.unit, figA, subtitle, kdecol,
                    b.fact=1.05, inner.lo=0.025, inner.hi=0.975)

  x.kde      <- KDE$x.kde
  x.kde.mode <- KDE$x.kde.mode

  return(list(x.kde=x.kde, x.kde.mode=x.kde.mode)) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_FindModeIndex.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  In a histogram, find the index of the bin containing the mode 
 
#  #  CHANGE HISTORY
#  31.08.2020 Start
#
# ==============================================================================

FindModeIndex <- function(breaks, counts, mode)
{
  #  INPUT 
  #  breaks   breaks defining the bins of a histogram
  #           bins are right-open: [ , ) except the last bin, which is[ , ]
  #           hist() in R needs right=FALSE for this convention
  #  counts   number of values per bin
  #  mode     usually from a density estimation 

  #  OUTPUT 
  #  idx      the index of the bin (in count) containing the mode
  # ---------------------------------------------------------------------------

  counts.n <- length(counts)
  idx      <- NA

  #  Which bin contains the mode?
  if (mode < breaks[1])
  { # below the first break
    idx <- 1
  }
  if (mode > tail(breaks, 1))
  { # above the last break
    idx <- counts.n
  }
  if (is.na(idx))
  { # somewhere within breaks
    idx <- which( (breaks[1:counts.n] <= mode) & 
                  (mode < breaks[2:(counts.n+1)] ))
  }
  if (is.na(idx))
  {
    idx <- Round(counts.n/2, 1)
    cat("\ ++++++ [FindModeIndex]  Mode index not found ++++++",
        "\n       Set to ", idx,
        "\n")
  }
  return(idx) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_Floor.R

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Function for rounding down to the next smaller unit 

#  CHANGE HISTORY
#  08.01.2020  Last change
# ============================================================================

Floor <- function(x,unit)
{ 
  # INPUT
  # x    number to round down
  # unit rounding unit.
  x.floor <- trunc(x/unit) * unit

  return(x.floor)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

#  F_GenDatCalcChi2.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Generate histogram data from lognormal distribution, collapse the histogram,
#  calculate chi squared from the histogram, given the parameters
#  of the generating distribution
 
#  #  CHANGE HISTORY
#  29.08.2020 Start
#
# ==============================================================================

GenDatCalcChi2 <- function(x.n, mue, sig, xlabel, subtitle, histcol, bordercol,
                           n.per.bin.min, bins.n.min, dev.figA)
{
  #  INPUT 
  #  x.n      number of values to generate  
  #  mue      mean of the generating LND
  #  sig      sd of the generating LND
  #  xlabel   label for data
  #  subtitle
  #  histcol
  #  bordercol
  #  n.per.bin.min
  #  dev.figA

  #  OUTPUT 
  #  
  # ===========================================================================

  x   <- rlnorm(x.n, meanlog=mue, sdlog=sig)

  #  Show standard histogram
  dev.set(dev.figA)
  x.hist0 <- hist(x, breaks="FD", right=FALSE, freq=FALSE, 
                  border="black", lty=1)

  #  Find mode estimate via KDE
  x.kde <- density(x)
  lines(x.kde, col="black")

  x.kde.mode <- x.kde$x[which.max(x.kde$y)]
  abline(v=x.kde.mode, lty=2)

  #  Collapse histogram, if necessary
  x.hist <- CollapseHist(x, x.hist0, xlabel, subtitle, histcol, bordercol, 
                         n.per.bin.min, dev.figA)

  #  Calculate chi2 contributions per bin
  tab <- CalcChi2(x.hist$breaks, x.hist$counts, mue, sig)

  #  Calculate partial chi2 sums over subintervals having at least bins.n.min 
  #  bins
  bins.n <- nrow(tab)
  PS.names <- c("bin.start", "bin.end", "bins.n", "chi2")
  PS <- matrix(NA, nrow=bins.n * (bins.n-1) / 2, 
               ncol=length(PS.names))
  colnames(PS) <- PS.names
  PS <-  data.frame(PS)

  PS.n <- 0
  for (bin.start in 1:(bins.n-bins.n.min+1) )
  { for (bin.end in (bin.start+bins.n.min-1):bins.n)
    { PS.n <- PS.n + 1
      PS[PS.n, "bin.start"] <- bin.start 
      PS[PS.n, "bin.end"]   <- bin.end
      PS[PS.n, "chi2"]      <- sum(tab[bin.start:bin.end, "chi2.bin"])
    }
  }
  PS[ , "bins.n"] <- PS[ , "bin.end"] - PS[ , "bin.start"] + 1
  PS <- PS[1:PS.n, ]
  PS <- PS[order(PS[ ,"bins.n"], PS[ ,"bin.start"]), ]

  return(list(tab=tab, PS=PS))
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test

# *****************************************************************************
#  F_hist.table.R
#  
#  Print histogram data in a readable form
#
#  09.09.2020 print() replaced by return()
#  04.03.2020
# =============================================================================

hist.table <- function(breaks, counts)
{
  #  Print histogram data
  x.n      <- sum(counts)
  counts.n <- length(counts)
  tab <- data.frame(x.lo=breaks[1:counts.n],
                    x.hi=breaks[2:(counts.n+1)],
                    count=counts,
                    prop=counts/x.n)
  return(tab) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
# F_Infoblock.R

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Writes essential control parameters in results plot
#  Version for TMC with dynamic prop.tmc
#

#  TO DO       
#  -  

#  #  CHANGE HISTORY
#  09.04.2021 oh.sym, dev.sym added
#  23.01.2021 x.bins.min replaced by bins.n.min, bins.n.max added
#  15.01.2021 r.fact added to call
#  10.01.2021 pathol.position removed
#  11.12.2020 Doan.adjust removed from call
#  02.12.2020 second switch for smooth.hist added
#  11.05.2020 More variables added
#  27.02.2020 Changed to function
#  24.01.2017 No details into plot if no plot was generated
# =================================================================================
#
Infoblock <- function(RunId, yplotmin, yplotmax, 
                     datafile, x.n, subsample.n, oh.sym, dev.sym, 
                     round.unit, smooth.hist1, 
                     smooth.hist2, 
                     n.per.bin.min, bins.n.min, bins.n.max,
                     lambda.min, lambda.max,   
                     l.fact, p.fact, r.fact, s.fact, w.fact, 
                     x.tr.prop.min, x.tr.prop.max, p.fit.min)
{ 
  #  Find user coordinates of actual window
  usr.coord <- par("usr")
  xplotmin <- usr.coord[1]
  xplotmax <- usr.coord[2]

  #  Use them to position the text
  text(xplotmin+0.010*(xplotmax-xplotmin),
       yplotmin+0.30*(yplotmax-yplotmin),
       labels=paste("\nRun Id            ", RunId,
                    "\ndata              ", datafile,
                    "\ntotal size        ", x.n,"   subsample.n   ",subsample.n, 
                    "\nout-/inpatient    ", oh.sym,
                    "\ndevice            ", dev.sym,
                    "\nround.unit        ", round.unit,
                    "\nsmooth.hist1      ", smooth.hist1,
                    "\nsmooth.hist2      ", smooth.hist2,
                    "\nn.per.bin.min     ", n.per.bin.min,
                    "\nbins.n.min        ", bins.n.min,
                    "\nbins.n.max        ", bins.n.max,
                    "\nl.fact            ", l.fact,
                    "\np.fact            ", p.fact,
                    "\nr.fact            ", r.fact,
                    "\ns.fact            ", s.fact,
                    "\nw.fact            ", w.fact,
                    "\nx.tr.prop.min/max ", x.tr.prop.min, "/", x.tr.prop.max,
                    "\nlambda.min/max    ", lambda.min, "/", lambda.max,
                    "\np.fit.min         ", p.fit.min
                    ),
     cex=0.8,pos=4,family="mono",col="navy") 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

#  F_IniMatrix.R
#
#  26.01.2021
# ==============================================================================

IniMatrix <- function(rnames, cnames, IniVal=NA)
{
  M <- matrix(IniVal, nrow=length(rnames), ncol=length(cnames))
  dimnames(M) <- list(rnames, cnames)
  return(M)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rn <- c("r1", "r2", "r3")
cn <- c("c1", "c2")

A <- IniMatrix(rn, cn)
A
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_Interpolx.R

# Linear interpolation: x, y, y target are given, xtarget will be determined

# INPUT
# x        x values MUST BE sorted (ascending). Sorting is not checked!. 
# y        y values corresponding to x
# ytarget  y target value, the corresponding x will be calculated

# NAs are automatically removed from the input.
# y must be a strictly monotone sequence, otherwise there is no unique 
# solution. Monotonicity is not checked!

# OUTPUT 
# xtarget  required x, y(xtarget) = ytarget

# CHANGE HISTORY
# 19.12.2020 Superfluous checks removed
# 03.07.2015 Start
# =============================================================================

Interpolx <- function(x,y,ytarget)
{  
  # Remove (x,y) pairs containing NA
  ok <- (!is.na(x)) & (!is.na(y))
  xx <- x[ok]
  yy <- y[ok]

  xtarget <- NA
  arg.ok <- TRUE
 
  if (is.na(ytarget)) 
  { cat("\n++++++ [Interpolx]: ytarget is NA\n")
    arg.ok <- FALSE
  }
  n <- length(xx)

  if (n <= 1) 
  { cat("\n++++++ [Interpolx]: x too short\n")
    arg.ok <- FALSE
  }
  if (length(yy) != n) 
  { cat("\n++++++ [Interpolx]: length(x) != length(y)\n")
    arg.ok <- FALSE
  }
  #  Reasonable value for ytarget?
  if (arg.ok)
  {
    if ( (ytarget < min(yy)) | (ytarget > max(yy)) )
    { cat("\n++++++ [Interpolx]: ytarget not in [min(y), max(y)]\n")
      arg.ok <- FALSE
    }
  }
  # No result if
  # ytarget is NA
  # length(x) <= 1
  # length(x) != length(y)
  # ytarget not in [min(y), max(y)]

  if (arg.ok)
  {
    for (i in 2:n)
    {
      if ( ((yy[i-1] <= ytarget) & (ytarget < yy[i  ])) |
           ((yy[i  ] <= ytarget) & (ytarget < yy[i-1])) )
      {
        xtarget <- xx[i-1] + (xx[i]-xx[i-1])*(ytarget-yy[i-1])/(yy[i]-yy[i-1])
        break
      }
    }
    if (ytarget == yy[n]) xtarget <- xx[n]
  }
  return(xtarget)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#D <- matrix(c(
#000, 135.0, 0.000,
#135, 135.5, 0.005,
#136, 136.5, 0.045,
#137, 137.5, 0.075,
#138, 138.5, 0.160,
#139, 139.5, NA,
#140, 140.5, 0.470,
#141, 141.5, 0.645,
#142, 142.5, 0.785,
#143, 143.5, 0.890,
#144, 144.5, 0.940,
#145, 145.5, 0.970,
#146, 146.0, 1.000), byrow=TRUE, ncol=3, nrow=13)
#colnames(D) <- c("u", "x", "y")
#Interpolx(D[ ,"x"], D[ ,"y"], 0.00)
#Interpolx(D[ ,"x"], D[ ,"y"], 0.30)
#Interpolx(D[ ,"x"], D[ ,"y"], 1.00)
#Interpolx(D[ ,"x"], D[ ,"y"], 2.00)

# 0.295


# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_Interpoly.R

# Linear interpolation: x, y, x target are given, ytarget will be determined

# INPUT
# x        x values MUST BE sorted (ascending). Sorting is not checked!. 
# y        y values corresponding to x
# xtarget  x target value, the corresponding y will be calculated
# LOCF     "Last observation carried forward": use largest x as result if
#          xtarget > max(x) ?

# NAs are automatically removed from the input.
# (x, y) must be a strictly monotone sequence, otherwise there is no unique 
# solution. Monotonicity is not checked!

# OUTPUT 
# ytarget  required y, x(ytarget) = xtarget

# CHANGE HISTORY
# 18.01.2021 Superfluous checks removed
# 03.07.2015 Start
# =============================================================================

Interpoly <- function(x, y, xtarget, LOCF=FALSE)
{  
  # Remove (x,y) pairs containing NA
  ok <- (!is.na(x)) & (!is.na(y))
  xx <- x[ok]
  yy <- y[ok]

  ytarget <- NA

  arg.ok <- TRUE
  if (is.na(xtarget)) 
  { cat("++++++ [Interpoly]: xtarget is NA\n")
    arg.ok <- FALSE
  }
  n <- length(yy)
  if (n <= 1) 
  { cat("\n++++++ [Interpoly]: y too short\n")
    arg.ok <- FALSE
  }
  if (length(xx) != n) 
  { cat("\n++++++ [Interpoly]: length(x) != length(y)\n")
    arg.ok <- FALSE
  }
  #  Reasonable value for ytarget? 
  if (arg.ok & !LOCF)
  {
    if ( !(xx[1] <= xtarget & xtarget <= xx[n]) &
         !(xx[n] <= xtarget & xtarget <= xx[1]) )
    { cat("\n++++++ [Interpoly]: xtarget not in [min(x), max(x)]\n")
      arg.ok <- FALSE
    }
  }
  if (arg.ok & LOCF)
  {
    if ( !(xx[1] <= xtarget))
    { cat("\n++++++ [Interpoly]: xtarget < min(x)\n")
      arg.ok <- FALSE
    }
  }

  if (arg.ok)
  {
    for (i in 2:n)
    {
      if ( (xx[i-1] <= xtarget & xtarget < xx[i  ]) |
           (xx[i  ] <= xtarget & xtarget < xx[i-1]) )
      {
        ytarget <- yy[i-1] + (yy[i]-yy[i-1])*(xtarget-xx[i-1])/(xx[i]-xx[i-1])
        break
      }
    }
    if (xtarget == xx[n]) ytarget <- yy[n]

    if ((xtarget > xx[n]) & LOCF) {ytarget <- yy[n]}
  }
  return(ytarget)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
#x0 <- c(1,2,3)
#y0 <- c(10,20,30)
#print(Interpoly(x0,y0,0.5))
#print(Interpoly(x0,y0,1.0))
#print(Interpoly(x0,y0,1.5))
#print(Interpoly(x0,y0,2.0))
#print(Interpoly(x0,y0,2.5))
#print(Interpoly(x0,y0,3.0))
#print(Interpoly(x0,y0,3.5))
#print(Interpoly(x0,y0,3.5,LOCF=TRUE))

#x0 <- c(1,2)
#y0 <- c(10,20,30)
#print(Interpoly(x0,y0,1.5))
#x0 <- c(1,2,3)
#y0 <- c(30,20,10)
#print(Interpoly(x0,y0,0.5))
#print(Interpoly(x0,y0,1.0))
#print(Interpoly(x0,y0,1.5))
#print(Interpoly(x0,y0,2.0))
#print(Interpoly(x0,y0,2.5))
#print(Interpoly(x0,y0,3.0))
#print(Interpoly(x0,y0,3.5))
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
#  F_o2R.R

o2R <- function(theta,RB,fastnull=1.e-10)
{ #  Parameter der PNV auf Originalskala --> unbeschränkte Skala (ganz R)
  #  Transformation durch inverse logit-Transformation mit Grenzen in RB
  #  (offene Intervalle)
  #  16.12.2017 

  #  theta   - (lambda,mue,sigma) einer Power-Normal-Verteilung
  #            Werte müssen innerhalb der durch RB3 gegebenen Grenzen liegen 
  #  RB[1, ] - c(lambda.min,lambda.max)
  #  RB[2, ] - c(mue.min,mue.max)
  #  RB[3, ] - c(sig.min,sig.max) 
  #   

  theta.n <- length(theta)
  ttheta  <- rep(NA,times=theta.n) 

  #  Logistische Transformation nur möglich, wenn Parameter innerhalb der
  #  Randbedingungen liegen, sonst Meldung

  for (k in 1:theta.n)
  {
    if ( (RB[k,1] < theta[k]) & (theta[k] < RB[k,2]) ) 
    { ttheta[k] <- log((theta[k]-RB[k,1])/(RB[k,2]-theta[k])) } else
    { cat("++++++++  Falsche Eingabe für o2R:",
          "\nParameter-Nr.          ",k,           
          "\nParameter-Wert         ",theta[k],
          "\nZulässige Untergrenze >",RB[k,1],
          "\nZulässige Obergrenze  <",RB[k,2],
          "\n") 
    }
  }

  #cat("o2R: ttheta,theta\n")
  #print(cbind(ttheta,theta))

  return(ttheta)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_pdf.PN.R

#  Truncated minimum chi square estimation

# (c) wwosniok@math.uni-bremen.de

#  CHANGE HISTORY
#  16.01.2021 Comments changed to English
# ============================================================================ 

pdf.PN <- function(x,lambda,mue,sigma,fastnull=1.e-10)
{
  #  Probability density function of the power normal distribution (PND) 
  #  Background: see 
  #  Freeman J, Modarres R. Inverse Box-Cox: the power normal 
  #  distribution. Stat Probabil Letters 2006;76; 764-72, particularly p 766

  #  INPUT 
  #  x       value at which to evaluate the pdf
  #  lambda  shape parameter. Real number. Special cases:
  #          1: PND is normal distribution with mean  and sd sigma
  #          0: PND is a logarithmic normal distribution with parameters
  #             mue and sigma
  #          If abs(lambda) < fastnull, lambda== 0 is assumed and the 
  #          intrinsic function plnorm is used
  #  mue     mean of the PND (of BoxCox(x, lambda))
  #  sigma   sd of the PND   (of BoxCox(x, lambda))
  #  fastnull (="nearly zero") absolute values < fastnull are treated as zero
   
  #  OUTPUT
  #  pdf     Probability density function  

  # ==========================================================================

  pdf  <- rep(fastnull, times=length(x))
  xgt0 <- (x > fastnull)
  T    <- 1/(lambda*sigma) + mue/sigma
 
  if (lambda > fastnull)
  { #  lambda > 0  in principle
    K <- pnorm(T,mean=0,sd=1)
    pdf[xgt0] <- (1/K) * (1/(sigma*sqrt(2*pi))) * (x[xgt0])^(lambda-1) *
                exp(-0.5*((BoxCox(x[xgt0],lambda)-mue)/sigma)^2 )
    pdf[!is.finite(pdf)] <- fastnull   
  } 
  if (abs(lambda) <=  fastnull)  
  { # lambda == 0 (in principle)
    #  K <- 1 
    #pdf[xgt0] <-         (1/(x[xgt0] * sigma*sqrt(2*pi))) * 
    #              exp(-0.5*((BoxCox(x[xgt0],lambda)-mue)/sigma)^2 )

    #  intrinsic function is faster
    pdf[xgt0] <- dlnorm(x[xgt0],meanlog=mue,sdlog=sigma)
  }
  if (lambda < -fastnull)
  { #  lambda < 0   in principle
    #  Potential problem: if T -> Inf , K -> 0, which generates an 
    #  overflow below. wenn T groß, dann aber auch R -> 0
    K <- pnorm(-T,mean=0,sd=1)
    K <- max(K,fastnull)

    #cat("T,K",T,K,"\n")
    #R <- (1/(sigma*sqrt(2*pi))) * (x[xgt0])^(lambda-1) *
    #              exp(-0.5*((BoxCox(x[xgt0],lambda)-mue)/sigma)^2 )
    #cat("exp...",R,"\n")

    pdf[xgt0] <- (1/K) * (1/(sigma*sqrt(2*pi))) * (x[xgt0])^(lambda-1) *
                  exp(-0.5*((BoxCox(x[xgt0],lambda)-mue)/sigma)^2 )
  }

  #  This is a definition. @@@ Check for appropriateness eg if lambda=1 
  pdf[!xgt0] <- 0

  if (any(!is.finite(pdf)))
  { cat("[pdf.PN] Infinite values:\n") 
    cat("lambda,mue,sigma: ",lambda,mue,sigma,"\n")
    cat("x, pdf(x) with non-finite pdf\n")
    control <- data.frame(x,pdf)
    print(control[!is.finite(control[ ,2]), ])   
  }
  return(pdf)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #  lambda = -0.09989172   mue = 11.25949   sigma = 0.01440878 
#  x 14 - 19
#lambda <- 0.09989172
#mue    <- 11.25949
#sigma  <-  0.01440878
#x.P50  <- q.PN(0.50, lambda, mue, sigma)
#x      <-  seq(1,20)
#x.pdf  <- pdf.PN(x, lambda, mue, sigma)
#x.cdf  <- cdf.PN(x, lambda, mue, sigma)
#cat("\n x, x.pdf")
#print(cbind(x, x.pdf, x.cdf))

#  Comparison with R provided function
# install.packages("powdist", dependencies=TRUE)

#library(powdist)
#dpnorm(x, lambda = lambda, mu = mue, sigma = sigma, log = FALSE)

# Kundu, D. and Gupta, R. D. (2013) Power-normal distribution. 
#  Statistics 47, 110–125

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_PermissibleDiff.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculate the permissible difference pD according to 
#  Rainer Haeckel*, Werner Wosniok and Eberhard Gurr
#  Diagnostic efficiency in models for permissible
#  measurement uncertainty
#  J Lab Med 2017; 41(6): 309–315
#  and 
#  2021_01_06_Zula__ssige_Messunsicherheit_DGKL.3.10.20_-1.xlsx
#  (from DGKL web site)

 
#  CHANGE HISTORY
#  26.03.2021 Start
# ==============================================================================

PermissibleDiff <- function(RL1, RL2, xi)
{
  #  INPUT 
  #  RL1   Lower 95% reference limit (2.5%)
  #  RL2   Upper 95% reference limit (97.5%) 
  #  xi    Point at which to calculate the pD

  #  OUTPUT 
  #  xi.lo    lower limit of the permissible interval
  #  xi.hi    upper limit of the permissible interval
  # ---------------------------------------------------------------------------

  y.sd      <- (log(RL2) - log(RL1)) / 3.92
  x.med     <- sqrt(RL1 * RL2)
  CV.E.star <- 100*sqrt(exp(y.sd^2) - 1)
  pCV.A     <-  sqrt(CV.E.star - 0.25)
  ps.A.med  <- pCV.A * x.med * 0.01

  ps.A.xi <- 0.8 * ps.A.med * xi / x.med + 0.2 * ps.A.med
  pD.xi   <- 1.28 * ps.A.xi

  xi.lo   <- xi - pD.xi
  xi.hi   <- xi + pD.xi

  #  Not used
  # pCV.A.RL2 <- 100 * ps.A.RL2 / RL2
  # pU.pct    <- 2.39 * pCV.A.RL2

  return(data.frame(xi.lo=xi.lo, xi.hi=xi.hi))
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test  (compare with  
#         2021_01_06_Zula__ssige_Messunsicherheit_DGKL.3.10.20_-1.xlsx
#         sodium)

#pI <- PermissibleDiff(135, 145, 135)
#pI

#pI <- PermissibleDiff(135, 145, 145)
#pI

#pI <- PermissibleDiff(135.9, 144.9, 144.9)
#pI

#pI <- PermissibleDiff(135, 145, c(135, 145))
#pI


# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# F_plot.res.R

plot.res <- function(xsupp, comp.n, lambda.BC, lambda, mue, sig, col, lty)
{
  for (ic in 1:comp.n)
  { x.pdf <- lambda[ic] * pdf.PN(xsupp, lambda.BC, mue[ic], sig[ic])
    lines(xsupp, x.pdf, col=col, lty=lty)
  }
}

  # *****************************************************************************
#  F_PlotConfElli.R
#
#  (c) wwosniok@math.uni-bremen.de
#
#  Plot a confidence ellipse
#
#  16.12.2018 Einrichtung
# =============================================================================

PlotConfElli <- function(mean,cov,alpha,figA)
{ #
  # mean   mean of the ellipse
  # cov    covariance matrix 
  # alpha  area outside the ellipse
  # figA   graphics window, NA: no plot
  # ---------------------------------------------------------------------------
  
  #  Start with ellipse for standard normal distribution
  r <- qnorm(1-alpha/2)
 
  phi <- seq(0,360,by=5) * pi/180
  x0  <- r * cos(phi)
  y0  <- r * sin(phi)

  #  Jordan decomposition of cov
  cov.decomp <- eigen(cov,TRUE)
  #print(cov.decomp)  

  #  Root
  cov.root <- cov.decomp$vectors %*% diag(sqrt(cov.decomp$values)) %*% 
              cov.decomp$vectors 
  #  Check root
  cov.root %*% cov.root

  cov.elli0 <- cbind(x0,y0) %*% cov.root 
  cov.elli  <- cov.elli0 + matrix(rep(1,times=length(phi)),ncol=1) %*%
                           matrix(mean,nrow=1)  

  if (!is.na(figA))
  { dev.set(figA)
    xplotmin <- min(c(x0,cov.elli))
    xplotmax <- max(c(x0,cov.elli))
    yplotmin <- min(c(y0,cov.elli))
    yplotmax <- max(c(y0,cov.elli))

    plot(x0,y0,type="l",col="gray",
         xlim=c(xplotmin,xplotmax),ylim=c(yplotmin,yplotmax))
    abline(h=0,col="gray")
    abline(v=0,col="gray")
    lines(cov.elli,col="green3")
  }

  return(cov.elli)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#alpha <- 0.05
#mean <- c(1,-2)
#cov <- matrix(c(10, 2,
#                 2, 5),byrow=TRUE,ncol=2) 

#c.e <- PlotConfElli(mean,cov,alpha,NA)

#xplotmin <- min(c(x0,cov.elli))
#xplotmax <- max(c(x0,cov.elli))
#yplotmin <- min(c(y0,cov.elli))
#yplotmax <- max(c(y0,cov.elli))

#plot(cov.elli,type="l",col="green3",
#     xlim=c(xplotmin,xplotmax),ylim=c(yplotmin,yplotmax))
#     abline(h=0,col="gray")
#     abline(v=0,col="gray")
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
#  PlotDS.R
#
#  Plot descriptive statistics into plot
#
#  28.02.2020
# ==============================================================================

PlotDS <- function(xvar)
{
  #  Plot must exist. Get extreme coordinates.
  usr.coord <- par("usr")
  xplotmin <- usr.coord[1]
  xplotmax <- usr.coord[2]
  yplotmin <- usr.coord[3]
  yplotmax <- usr.coord[4]

  #  Use them to position the text
  xvar.n    <- sum(!is.na(xvar))
  xvar.min  <- min(xvar, na.rm=TRUE)
  xvar.P025 <- Quantile(xvar, probs=0.025)
  xvar.med  <- median(xvar, na.rm=TRUE)
  xvar.mea  <- mean(xvar, na.rm=TRUE)
  xvar.P975 <- Quantile(xvar, probs=0.975)
  xvar.max  <- max(xvar, na.rm=TRUE)


  text(xplotmin+0.010*(xplotmax-xplotmin),
       yplotmin+0.72*(yplotmax-yplotmin),
       labels=paste("\n# values ", format(xvar.n,   digits=6),
                    "\nMinimum  ", format(xvar.min, digits=4),
                    "\nP025     ", format(xvar.P025,digits=4),
                    "\nMedian   ", format(xvar.med, digits=4),
                    "\nMean     ", format(xvar.mea, digits=4),
                    "\nP975     ", format(xvar.P975,igits=4),
                    "\nMaximum  ", format(xvar.max, digits=4)
                    ),
     cex=0.7,pos=4,family="mono",col="firebrick") 

  # Show data points
  points(xvar, rep(yplotmin+0.10*(yplotmax-yplotmin), times=length(xvar)),
         type="h", col="gray")  

  # Show quantiles etc.
  abline(v=xvar.mea, col="black")
  abline(v=c(xvar.P025, xvar.med,xvar.P975), col="black", lty=2)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_PlotHistFit.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Plot a histogram of the data plus additional components, see INPUT
 
#  #  CHANGE HISTORY
#  01.05.2021 kde plot, difference kde-predicted removed, residual histogram
#             added 
#  12.01.2021 x.tr.prop.range added to call
#  01.12.2020 Call changed (n.per.bin.mea removed, denlty added)
#  19.11.2020 Start
#
# ==============================================================================

PlotHistFit <- function(figA, figA.file, figA.type, 
                        x.hist, xlabel, title, subtitle, n.per.bin.min,
                        x.kde, x.val, counts.range, x.tr.prop.range,
                        x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                        x.clip.type, yplotmaxA, par.las, par.tcl,
                        bordercol, histcol, difcol, kdecol, denlty)
{
  #  =====================================================================
  #  Plot histogram for the subset of the x range defined by x.clip
  #  x plot limits can be forced to coincide with x.breaks 
  #  (x.clip.type=="coincide") or can be used as specified
  #  (x.clip.type=="as.specified")
  #  Additonal components may be added here, see INPUT

  #  INPUT
  #  figA
  #  figA.file 
  #  figA.type
  #
  #  OUTPUT  

  # ==========================================================================

  dev.set(figA)
  #bringToTop(figA)
  
  #  x limits depend on x.clip settings

  if (is.na(x.clip.min))
  { #  Nothing defined
    xplotminA <- x.hist$breaks[1]
  } else
  { if (x.clip.type=="as.specified")
    { 
      xplotminA <- x.clip.min
    } else
    {
      #  Limits defined by x.clip are shifted such that the limits
      #  coincide with break points
      #  Find closest break <= x.clip.min
      x.breaks.cand   <- x.hist$breaks[x.hist$breaks <= x.clip.min]
      x.breaks.cand.n <- length(x.breaks.cand)
      if (x.breaks.cand.n > 0)
      { 
        xplotminA <- x.breaks.cand[x.breaks.cand.n] 
      } else
      {
        xplotminA <- x.clip.min
      }
    
    }  # x.clip.type=="as.specified" ...   
  }    # is.na(x.clip.min) ...

  if (is.na(x.clip.max))
  { #  Nothing defined
    xplotmaxA <- tail(x.hist$breaks, 1)
  } else
  { if (x.clip.type=="as.specified")
    { 
      xplotmaxA <- x.clip.max
    } else
    { 
      #  find closest break <= x.clip.max
      x.breaks.cand <- x.hist$breaks[x.hist$breaks >= x.clip.max]
      x.breaks.cand.n <- length(x.breaks.cand)
      if (x.breaks.cand.n > 0)
      { 
        xplotmaxA <- x.breaks.cand[1] 
      } else
      {
        xplotmaxA <- x.clip.max
      }
    }    # x.clip.type=="as.specified" ...
  }      # is.na(x.clip.max)

  # Use y maximum if specified, otherwise use histogram density maximum 
  if (is.na(yplotmaxA))
  {
    yplotmaxA <- 1.1 * max(x.hist$density)
  }

  #  Set yplotminA such that there is space for indicating observed values
  yplotminA <- -0.05 * yplotmaxA

  #  Start plot
  par(las=par.las)   
  par(tcl=par.tcl)
   
  if (!is.na(x.clip.by1))
  { par(xaxt="n") } else
  { par(xaxt="s") } 

  #  Histogram
  plot(x.hist, freq=FALSE, 
       xlim=c(xplotminA,xplotmaxA),
       ylim=c(yplotminA,yplotmaxA),
       col=histcol, border=bordercol,
       xlab=xlabel, ylab="Density",
       main=NULL,
       #
       # 
       sub=NULL, cex.sub=0.7)
  par(xaxt="s")
  if (!is.na(x.clip.by1) )
  { 
    axis(1,at=seq(xplotminA, xplotmaxA, by=x.clip.by1), labels=TRUE)  
  }

  if (!is.na(x.clip.by2) )
  { 
    axis(1, at=seq(xplotminA, xplotmaxA,by= x.clip.by2), labels=FALSE) 
  }

  #  Indicate observed values
  x.val.n <- length(x.val)
  if (x.val.n > 1)
  {
    points(x.val, rep(0.5*yplotminA, times=x.val.n), type="h", col="gray10")
  }

  #  Density estimation
  if (length(x.kde) > 1)
  {
    lines(x.kde,col=kdecol,lty=denlty)  
  }

  #  Save plot, if requested
  if (!is.na(figA.file))
  {
    savePlot(file=figA.file,type=figA.type)
  }
  return(c(yplotmin=yplotminA,yplotmax=yplotmaxA))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#  F_PlotHistResid.R
#  former name: F_HistoResid.R
#  (c) wwosniok@math.uni-bremen.de

#  TMC: Truncated minimum chi-square estimation

#  Function calculating and plotting the difference
#  between 2 histograms (usually observed and expected)
 
#  #  CHANGE HISTORY
#  26.04.2021 Start
# ==============================================================================

PlotHistResid <- function(temp.met, figA, yplotminA, bordercol, fillcol)
{
  # INPUT
  #  temp.met  object containing observed and expected value of a histogram,
  #             produced by chi2trunc, see there for details
  #  xlabel     histogram label
  #  figA	    window to put the plot in (as addition)

  # OUTPUT
  # delta.hist    the difference histogram

  # ============================================================================

  histo.resid <- temp.met$tab[ , c("UG", "OG", "nobs", "nerw", "diff")] 
  if (length(temp.met$tab.lo) > 1)
  {
    histo.resid <- rbind(temp.met$tab.lo[ ,c("UG", "OG", "nobs", "nerw", "diff")],
                         histo.resid)
  }

  if (length(temp.met$tab.hi) > 1)
  {
    histo.resid <- rbind(histo.resid, 
                     temp.met$tab.hi[ ,c("UG", "OG", "nobs", "nerw", "diff")] )
  }

  #  Create a histogram object
  x.n <- sum(histo.resid[ , "nobs"])
  histo.resid <- data.frame(histo.resid,
                          obs.den=histo.resid[ , "nobs"]/   
                         ( x.n *(histo.resid[ , "OG"]-histo.resid[ , "UG"])),
                          exp.den=histo.resid[ , "nerw"]/   
                         ( x.n *(histo.resid[ , "OG"]-histo.resid[ , "UG"])),
                          delta.den=histo.resid[ , "diff"]/
                         ( x.n *(histo.resid[ , "OG"]-histo.resid[ , "UG"]) ))

  #  Truncate negative densities for aesthetic reasons 
  histo.resid[histo.resid[ , "delta.den"]<yplotminA, "delta.den"] <- yplotminA

  delta.hist <- list(breaks=c(histo.resid[ , "UG"],
                             tail(histo.resid[ , "OG"], 1)),
                   counts=histo.resid[ , "diff"],
                   density=histo.resid[ , "delta.den"],
                   mids=(histo.resid[ , "UG"]+histo.resid[ , "OG"])/2,
                   xname="delta",
                   equidist=FALSE)
  class(delta.hist) <- "histogram"  

  #  Add the histogram to an existing plot
  if (!is.na(figA))
  {
    #plot(delta.hist, border=bordercol, col=fillcol,  density=20, angle=45,
    #     lwd.border=5, add=TRUE)
    plot(delta.hist, border=bordercol, col=fillcol, add=TRUE)
  }

  return(delta.hist)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test: see Test_HistoResid.R

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv



  
  
# ****************************************************************************
#  F_PlotMetRes.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Plot results of truncated estimation by method met in an existing histogram plot
 
#  #  CHANGE HISTORY
#  01.05.2021 kde residuals replaced by histogram residuals
#  08.01.2021 mode added to plot
#  07.12.2020 xsupp added to call
#  19.11.2020 Start
#
# ==============================================================================

PlotMetRes <- function(figA, figA.file, figA.type, yplotminA,  
                       xsupp, prev.c.met, xc.pdf.met, metcol,
                       x.kde.mode, kdecol, x.kde.met, difcol,
                       difbordercol, difhistcol, 
                       x.RL1.met, x.RL2.met, y.RL, 
                       x.tr.lo, x.tr.hi, y.TI, trcol,
                       temp.met) 
{
  #  =====================================================================
  #  Plot histogram for the subset of the x range defined by x.clip
  #  x plot limits can be forced to coincide with x.breaks 
  #  (x.clip.type=="coincide") or can be used as specified
  #  (x.clip.type=="as.specified")
  #  Additonal components may be added here, see INPUT

  #  INPUT
  #  figA
  #  figA.file 
  #  figA.type
  #
  #  OUTPUT  

  # ==========================================================================

  dev.set(figA)
  #bringToTop(figA)

  #  Density estimated by method met  
  lines(xsupp, prev.c.met * xc.pdf.met, col=metcol, lwd=2)

  if (length(x.kde.met) > 1)
  {
    #  Difference between overall kde and estimated density
    lines(xsupp, x.kde.met, col=difcol, lwd=2) 
  }

  #  Histogram of difference between observed and expected
  delta.hist <- PlotHistResid(temp.met, figA, yplotminA, difbordercol, 
                              difhistcol)

  #  Mode from kde
  abline(v=x.kde.mode, col=kdecol, lty=2)

  #  Estimated RLs  
  plotRL(x.RL1.met, y.RL, metcol, pos=2)
  plotRL(x.RL2.met, y.RL, metcol, pos=4)

  #  Truncation interval
  lines(c(x.tr.lo, x.tr.hi), c(y.TI, y.TI), col=trcol, lwd=2)  

  #  Save plot, if requested
  if (!is.na(figA.file))
  {
    # savePlot(file=figA.file,type=figA.type)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_PlotRL.R
#
#  CHANGE HISTORY
#
#  07.01.2020 Start

#  -----------------------------------------------------------------------------
plotRL <- function(RL,y,farbe,pos=4)
{ 
  #  RL    - x-Position der RL
  #  y     - y-Position der RL
  #  farbe - Farbe für Linie und Text
  #  pos   - Position des Textes neben Linie 
  #          1,2,3,4 = unter, links, über, rechts
 
  lines(c(RL,RL),c(0,y),col=farbe)
  text(RL,y,pos=pos,labels=format(RL,digits=3),cex=0.8,col=farbe)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_ProcessDL.R
#
#  (c) wwosniok@math.uni-bremen.de

#  Process data entries of the form "< DL" (DL = detection/determination limit).
#  No action, if no such entries exist.
#  If >= 1 entries exist, several options available. The surrogate factor
#  0 <= s.fact < DL must be supplied. s.fact affects breaks, kde, 
#  QQ estimation.

#  If 1 unique DL exists:
#  option 1 (default):
#  - replace "< DL"                            by s.fact*DL
#    replace values < DL (should never exist)  by s.fact*DL
#    
#  option 2 (not recommended):
#  - replace "< DL"                            by s.fact*DL
#    leave values < DL (should never exist)    as they are

#  If > 1 unique DL exist:
#  option 1 (default):
#  - replace all "< DL"                        by s.fact*max(DL)
#    replace all values < max(DL)              by s.fact*max(DL)
#    
#  option 2 (not recommended):
#  - replace all "< DL"                        by s.fact*max(DL)
#    leave values < some DL (should never exist)    as they are

#  option 3 (not recommended):
#  - replace "< DL"                               by s.fact*DL
#    replace values < some DL (should never exist)  by s.fact*min_DL(DL>value)

#  option 4 (not recommended):
#  - replace "< DL"                               by s.fact*DL
#    leave values < some DL (should never exist)  as they are

#  Non-numeric entries other than "<DL" are left in the data unchanged.

#  Change History
#  27.04.2020 detect.limits.max added to return argument, is needed by
#             other functions
#             Recommendation: use option 1 generally 
#  14.04.2020 Start from 070_ReadData.R
# =============================================================================

ProcessDL <- function(v, s.fact, option=1)
{
  #
  # INPUT
  #
  # v - data vector
  # s.fact - surrogate factor, \in [0, DL]
  # option - controls how values < DL are processed, see above

  # OUTPUT
  # w <- vector with DLs processed.
  # ===========================================================================

  #  Check for permissible option
  if (option %in% c(1,2,3,4))
  {
    detect.limits.max <- NA
    w <- v

    # Are there non-numerical entries?
    options(warn=-1)           #  suppress warning regarding NAs resulting from
                               #  transformation 
    char.val.table <- table(v[is.na(as.numeric(v))]) 
    options(warn=0)

    if (length(char.val.table) == 0)
    { 
      cat("\n No non-numerical values found\n")
    } else
    {
      cat("\n Table of non-numerical entries in the data\n")
      print(char.val.table)
    
      #  Is "<" in the string?
      char.val <- names(char.val.table)
      has.lt <- str_detect(char.val,"<")

      #  Answer may be NA, means here "no"
      has.lt[is.na(has.lt)] <- FALSE
      has.lt.n <- sum(has.lt)        # # of values in char.val with form "< DL"

      if (has.lt.n > 0)
      {
        #  Remove non-numerical values other than "< ..." from char.val.table
        char.val.table <- char.val.table[has.lt]
        char.val       <- char.val[has.lt]

        cat("\n Table of '< DL' entries in the data\n")
        print(char.val.table)

        #  Determine the value(s) of DL. May look different, but be the same 
        #  numerically ("<5" "<5.0", " 5", ...),
        #  so change to numerical before counting  
        DL <- rep(NA, times=has.lt.n)

        for (i.lt in 1:has.lt.n)
        {
          lt.idx <- str_locate(char.val[i.lt], "<")
          DL[i.lt] <- char.val[i.lt]
          str_sub(DL[i.lt], start = lt.idx[ ,1], end = lt.idx[ ,2]) <- " "
        }
        DL <- as.numeric(DL)
        detect.limits.max <- max(DL)

        char.val.table <- rbind(char.val.table, DL)
        char.val.table <- matrix(as.numeric(char.val.table), nrow=2)
        colnames(char.val.table) <- char.val
        rownames(char.val.table) <- c("count", "DL")

        cat("\n Table of '< DL' entries in the data and their numerical values\n")
        print(char.val.table)

        #  Make v numeric, where possible, otherwise (non-numeric gets NA)
        options(warn=-1)      #  suppress warning regarding NAs resulting from
                              #  transformation
        options(warn=-1)
        v.num <- as.numeric(v)
        options(warn=0)       #  return to usual warning style

        #  In the data, replace "< DL" by their surrogate values according to
        #  the chosen option

        if (option == 1)
        { #  Global surrogate for DL and values < max. DL
          for (i.lt in 1:has.lt.n)
          {
            w[v==char.val[i.lt]] <-  s.fact*detect.limits.max 
          }

          #  Real values < DL are replaced by s.fact*detect.limits.max
          subset <- (!is.na(v.num)) & (v.num < detect.limits.max)
          #print(data.frame(v, subset))

          subset[is.na(subset)] <- FALSE         
          w[subset]   <-  s.fact*detect.limits.max
        }

        if (option == 2)
        { #  Global surrogate only for DL
          for (i.lt in 1:has.lt.n)
          {
            w[v==char.val[i.lt]] <-  s.fact*detect.limits.max 
          }
          #  Real values < DL are unchanged
        }

        if (option == 3)
        { #  Individual surrogate for DL and values < indiv. DL
          for (i.lt in seq(has.lt.n, 1, by=-1))
          {
            w[v==char.val[i.lt]] <-  s.fact*char.val.table[2, i.lt]
            #cat("\n option=3 (1), i.lt=", i.lt) 
            #print(data.frame(v, v.num, w))   
 
            #  Real values < actual DL are replaced by s.fact*actual DL
            # look for values < the actual DL
            subset <-  (!is.na(v.num)) & (v.num < char.val.table[2, i.lt])
            subset[is.na(subset)] <- FALSE
            if (sum(subset) > 0)          
            {
              w[subset] <-  s.fact*char.val.table[2, i.lt]
              #cat("\n option=3 (2)") 
              #print(data.frame(v, v.num, w))   
            }
          }
        }

        if (option == 4)
          { #  Individual surrogate for DL only
          for (i.lt in seq(has.lt.n, 1, by=-1))
          {
            w[v==char.val[i.lt]] <-  s.fact*char.val.table[2, i.lt] 
          }
          #  Real values < DL are unchanged
        }

        options(warn=-1)      #  suppress warning regarding NAs resulting from
                              #  transformation
        w0 <- as.numeric(w)   #  as.numeric generates NA from non-numeric
                              #  data that may still be in w. Will be deleted
                              #  later, but shall stay here unchanged.
        is.na.w0 <- is.na(w0)
        w0[is.na.w0] <- w[is.na.w0]
        w <- w0 

        options(warn=0)       #  return to usual warning

      }   # has.lt.n > 0 ...

    }     # length(char.val.table) ...
  
    return(list(data=w, detect.limits.max=detect.limits.max))

  } else       # wrong option
  { cat("\n ++++++ [ProcessDL] Wrong option: ", option,
        "\n ++++++ Values below DL are not processed",
        "\n")
    return(list(data=v, detect.limits.max=NA))
  }         
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_proportions.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Function determining the proportions of observations contained in an 
#  interval with limits ilo, ihi.
#  ilo and ihi are indices in xcounts (e.g. from a histogram)


#  CHANGE HISTORY

#  20.04.2021 Treatment of case with no proportion with requested properties
#             added
#  11.01.2021 Only intervals with proportion in [prop.min, prop.max)
#             and containing the mode are reported
#  18.02.2020 NA lines removed from output
#  28.12.2019 ilo, ihi range  restricted to ensure df.min bins
#  03.12.2019 Start

# ============================================================================

proportions <- function(xcounts, xbreaks.mode.idx, x.bins.min, prop.min, 
                        prop.max)
{ 
  # xcounts          vector of counts (typically from a histogram)
  # xbreaks.mode.idx index of the interval containing the mode of the 
  #                  distribution = index of the left breakpoint
  # x.bins.min       minimum number of bins in the interval
  # prop.min         minimal proportion included iin the interval
  # prop.max         proportion in the interval must be < prop.max
  #                  If prop.max == 1 is required, set prop.max <- 1.1
  # ==========================================================================     

  # ==========================================================================
  # if (print.log.message) { cat("%%%   proportions   Start\n") }
  # ==========================================================================

  xcounts.n <- length(xcounts)
  x.n       <- sum(xcounts)

  ilo.max <- xcounts.n - x.bins.min + 1

  tab.names     <- c("ilo", "ihi", "mode", "prop") 
  tab           <- matrix(NA, nrow=ilo.max^2, ncol=length(tab.names))
  colnames(tab) <- tab.names

  i <- 0
  for (ilo in 1:ilo.max)
  { ihi.min <- ilo + x.bins.min - 1
    for (ihi in ihi.min:xcounts.n)
    { 
      prop <- sum(xcounts[ilo:ihi]) / x.n
      if ( (prop.min <= prop) & (prop < prop.max) ) 
      { 
        i <- i + 1
        tab[i, "ilo"]  <- ilo
        tab[i, "ihi"]  <- ihi
        tab[i, "prop"] <- prop
      }
    }
  }

  if (i > 0)
  { 
    #  Proportion(s) with requested properties found
    tab <- tab[1:i, ]
    tab <- RowToMatrix(tab)
    tab[ , "mode"] <- (tab[ , "ilo"] <= xbreaks.mode.idx) &
                      (xbreaks.mode.idx <= tab[ , "ihi"]) 
    
    #  Reduce tab to intervals containing the mode
    tab <- tab[tab[ ,"mode"]==1, ]

    #  Make tab a matrix, if possible
    if (length(tab) > 0)
    {
      tab <- RowToMatrix(tab)
    }  else
    {
      #  No interval containing the mode
      tab           <- matrix(NA, nrow=1, ncol=length(tab.names))
      colnames(tab) <- tab.names
    }
  } else
  { 
    #  No proportion with requested properties found  
    tab         <- matrix(NA, nrow=1, ncol=length(tab.names))
    colnames(tab) <- tab.names
  }

  # ==========================================================================
  # if (print.log.message) { cat("%%%   proportions   End\n") }
  # ==========================================================================

  return(tab)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#xcounts <- c(10,10,10,10,10,10,10,10,10,10)
#xbreaks.mode.idx <- 5
#x.bins.min <- 3
#prop.min <- 0
#prop.max <- 1.1   # !!

#prop.min <- 0.35
#prop.max <- 0.85

#tab <- proportions(xcounts, xbreaks.mode.idx, x.bins.min, prop.min, prop.max)
#print(tab)

#prop1 <- p[p[ ,"mode"]==1, ]
#prop1 <- prop1[prop1[ ,"prop"] >= 0.55, ]

#prop2 <- aggregate(prop1[ ,"ihi"], by=list(prop1[ ,"ilo"]), FUN=min)
#colnames(prop2) <- c("ilo", "ihi")
#is.matrix(prop2)
#nrow(prop2)
#length(prop2)
#nrow(NULL) > 0

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# ****************************************************************************
#  F_q.PN.R

#  Truncated minimum chi square estimation

# (c) wwosniok@math.uni-bremen.de

#  CHANGE HISTORY
#  29.01.2021 Treatment of p=1 for lambda < 0 addded
#  28.01.2021 Error in negative branch corrected (Definition of T erroneously
#             deleted)
#  23.01.2021 Error in equation (1) of Freeman, Modarres corrected, see
#             PND_properties.tex / pdf
#             Note that for lambda < 0 the quantile function can attain 
#             indefinite values, because a power of a negative quantity arises. 
#  18.01.2021 Comments changed to English
# ============================================================================ 

q.PN <- function(p,lambda,mue,sigma,fastnull=1.e-10)
{
  #  Quantile of the PND, original scale (x)
  #  mue, sigma gelten innerhalb der PNV
  # siehe Freeman, S. 767
  #  Probability density function of the power normal distribution (PND) 
  #  Background: see 
  #  Freeman J, Modarres R. Inverse Box-Cox: the power normal 
  #  distribution. Stat Probabil Letters 2006;76; 764-72, particularly p 766

  #  INPUT   
  #  p       probability for which to determine the quantile, 0 < p < 1
  #  lambda  shape parameter. Real number. Special cases:
  #          1: PND is normal distribution with mean  and sd sigma
  #          0: PND is a logarithmic normal distribution with parameters
  #             mue and sigma
  #          If abs(lambda) < fastnull, lambda== 0 is assumed and the 
  #          intrinsic function plnorm is used
  #  mue     mean of the PND (of BoxCox(x, lambda))
  #  sigma   sd of the PND   (of BoxCox(x, lambda))
  #  fastnull (="nearly zero") absolute values < fastnull are treated as zero
   
  #  OUTPUT
  #  q.PN    Quantiles  

  # ==========================================================================
 
  q.PN <- rep(NA, times=length(p))
 
  if (lambda > fastnull)
  { 
    # cat("\n Positive lambda\n")
    T <- 1/(lambda*sigma) + mue/sigma
    V <- 1 - (1-p)*pnorm(T)    
    qnormV                  <- qnorm(V)
    
    #  qnormV  might be -Inf or Inf, set replacements to avoid overflow or NaN 
    qnormV[qnormV == -Inf]  <-  -10     #  -40 previously - why?
    qnormV[qnormV ==  Inf]  <-   10
    
    W <- lambda*(sigma*qnormV+mue) + 1     # kann durch numerische Fehler
                                               # negativ werden
    #  W might be < 0, avoid trouble

    W.gt.0        <- W > 0
    q.PN[W.gt.0]  <- (W[W.gt.0])^(1/lambda)
    q.PN[!W.gt.0] <- 0     
  }

  if (abs(lambda) < fastnull)
  { # LNV
    # cat("\n lambda = 0\n")
    q.PN  <- qlnorm(p,meanlog=mue,sdlog=sigma)
  } 

  if (lambda < -fastnull)
  { 
    # cat("\n Negative lambda\n")
    #  Calculation according to  F&M  (is wrong)
    # W <- lambda*(sigma*qnorm(p)+mue) + 1   # can get < 0
    # W.gt.0        <- W > 0
    # q.PN[W.gt.0]  <- (W[W.gt.0])^(1/lambda)
    # q.PN[!W.gt.0] <- 0     

    #  Corrected calculation, see PND_properties.tex / pdf  23.01.2021
    #  Holds for p in [0, 1). p = 1 and negative C need special treatment.
    q.PN <- rep(NA, times=length(p))

    T <- 1/(lambda*sigma) + mue/sigma
    A <- p * pnorm(-T)
    B <- lambda * sigma * qnorm(A) + lambda * mue
    C <- (1 + B)

    p.eq.1 <- (p > (1-fastnull)) 
    C.lt.0 <- (C < 0)    
    #  
    ok <- (!p.eq.1) & (!C.lt.0)
 
    q.PN[ok]  <- C[ok]^(1/lambda)    
    
    #  The final expression tends to zero^(1/lambda), which is 0, if p -> 1. 
    #  However, the sequence of q.PN for p -> 1 tends to Inf.
    if (any(p.eq.1)) { q.PN[p.eq.1] <- Inf } 

    #  Next condition seems to coincide with first condition
    #if (any(C.lt.0))
    #{ # This is an error condition
    #  cat("\n [q.PN] ++++++ No q.PN solution for  following p values",
    #      "\n        ++++++ lambda, mue, sigma",
    #      "\n        ++++++", lambda, mue, sigma,
    #      "\n")
    #  print(p[C.lt.0])
    #}  
  }

  return(unname(q.PN))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test

#source("F_BoxCox.R")
#source("F_BoxCoxInv.R")
#source("F_cdf.PN.R")
#source("F_pdf.PN.R")

#lambda <- -0.09989172
#mue    <- 11.25949
#sigma  <-  0.01440878
#pseq   <- seq(0, 1, by=0.1)
#x      <- q.PN(pseq, lambda, mue, sigma)

#x.pdf  <- pdf.PN(x, lambda, mue, sigma)
#x.cdf  <- cdf.PN(x, lambda, mue, sigma)
#cat("\n x, x.pdf")
#print(cbind(pseq, x, x.pdf, x.cdf))


#lambda <- -0.1
#mue    <- 2.184875
#sigma  <-  0.3562688
#pseq   <- c(seq(0, 1, by=0.1), seq(0.91, 1.00, by=0.01), 
#            0.999, 0.9999, 0.99999, 0.999999)
#x      <- q.PN(pseq, lambda, mue, sigma)

#x.pdf  <- pdf.PN(x, lambda, mue, sigma)
#x.cdf  <- cdf.PN(x, lambda, mue, sigma)
#cat("\n x, x.pdf")
#print(cbind(pseq, x, x.pdf, x.cdf))

#q.PN(0.975, -1.000000e-01, 2.521592, 0.2969054)  # 12.30498  41.08872
#cdf.PN(41.08872, -1.000000e-01, 2.521592, 0.2969054)  # 

#Q <- c()
#lambda.seq <- seq(0.99, 1.01, by=0.001)
#for (lambda in lambda.seq)
#{ 
#  Q <- c(Q, q.PN(0.5, lambda, 140, 2.551027))
#}
#print(cbind(lambda.seq, Q))
#plot(lambda.seq, Q, type="o")
#abline(h=140)
#abline(v=1)

#for (lambda in seq(0.90, 1.10, by=0.01))
#{ cat("\n", lambda, cdf.PN(140, lambda, 140, 2.551027) ) }

# ****************************************************************************
#  F_QQW.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculate a QQ-plot based RL estimate
#  Method similar to the modified Tukey based method by G. Hoffmann (?),
#  but modified by WW in several respects.
 
#  #  CHANGE HISTORY
#  17.12.2020 Emergency action if only nonacceptable solutions changed
#  13.12.2020 prev.acc.lo/hi added to call
#  07.12.2020 bins.n.min replaced by bins.n.min.act
#  05.12.2020 Result plot (J) done by standard function
#  01.12.2020 x.red.hist removed from call
#  28.11.2020 Calculation of estimated prevalences moved to CalcPre()
#  16.11.2020 Call of TIQQ changed
#             Calculation of x.poly changed 
#  14.11.2020 bins.n.min added to call and to TIQQ( )
#  08.10.2020 Start
#
# ==============================================================================

QQW <- function(x, x.red, round.unit, lambda.seq, 
                x.hist, x.kde, x.kde.mode, x.kde.mode.idx, x.val, 
                x.unaffected.lo, x.unaffected.hi,   
                n.per.bin.min, bins.n.min.act, x.tr.bins.min, 
                counts.range, x.tr.prop.range, 
                x.tr.prop.min, x.tr.prop.max, prev.acc.lo, prev.acc.hi, 
                xsupp,
                x.Q1, x.Q2, RL1.p, RL2.p,
                lambda.gen, mue.gen, sig.gen,
                df.est,df.con,
                l.fact, p.fact, r.fact, w.fact, 
                xlabel, subtitle,
                figD, figE, figF, figG, figH, 
                figI, figI.file, 
                figJ, figJ.file, figtype,
                x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                x.clip.type, par.las, par.tcl,
                bordercol2, histcol2, kdecol, qqrawcol, polycol, gencol, 
                difcol, denlty,
                print.log.message)
{
  #  INPUT 
  #  x              data vector
  #  round.unit
  #  lambda.seq     candidate lambda values 
  #  x.hist         histogram that will later be used by tmc, breaks are used
  #                 here as candidate limits of the truncation interval
  #  x.tr.bins.min  min number of bins in the truncation interval (TI), 
  #  x.tr.prop.min  min proportion of values in the TI
  #  x.tr.prop.max  max proportion of values in the TI
  #  lambda.gen     if simulated data: lambda for generation 
  #  mue.gen        if simulated data: mue for generation 
  #  sig.gen        if simulated data: sig for generation
  #  xlabel         label for analyte 
  #  subtitle       subtitle in FigA 
  #  figD           in TIQQ: BC transformed collapsed histogram
  #  figE           in TIQQ: QQ plot for BC transformed data, polygon 
  #                       approximation of the QQ plot
  #  figF           in TIQQ: All QQ plot based parameter estimates
  #  figG           in TIQQ: All QQ plot based RL estimates
  #  figH           in TIQQ: QQ plot based relative errors in RL estimates
  #  figI           Collapsed histogram and QQ plot based fits
  #  bordercol2     colour for histogram borders, collapsed
  #  histcol2       colour for bin area, initial, collapsed
  #  kdecol         colour for kde
  #  qqrawcol       symbol color in raw QQ plot 
  #  polycol        symbol color in untied QQ plot
  #  gencol         colour for plot objects showing values from generation

  #  OUTPUT 
  #  
  # ---------------------------------------------------------------------------

  #cat("\n *** Start QQW \n")

  x.n <- length(x)

  #  The original QQ plot method as in Hoffmann 1963 and others is done
  #  by @@@@ elsewhere (execute this programm on the raw data)   
 
  #  Loop over lambda sequence
  #
  ilambda <- 0
  for (lambda in lambda.seq)
  { ilambda <- ilambda + 1 
    
    #  Transform data, histogram breaks and kde
    y        <- BoxCox(x, lambda)
    y.red    <- BoxCox(x.red, lambda)
    y.breaks <- BoxCox(x.hist$breaks, lambda)

    y.red.hist <- hist(y.red, breaks=y.breaks, right=FALSE, plot=FALSE)
    y.kde      <- TransformKDE(x.kde, BoxCox(x.kde$x, lambda))

    # TIs must contain the transformed value of x.kde.mode, not the mode
    # of y.red. Plots contain this transformed mode.
    y.kde.mode <- BoxCox(x.kde.mode, lambda)
                    
    #  Changed: Find the bin (in y.red.hist$counts) containing the 
    #  the transformation of x.kde. This concentrates TI intervals around
    #  the mode of the x distribution on the original scale. 

    # y.red.hist.mode.idx is the same as x.kde.mode.idx
    
    # Get a proposal for reasonable truncation intervals, given lambda

    tab1.qqw <- TIQQ(x.hist, x.kde.mode, 
                     y.red, y.red.hist, y.kde, y.kde.mode, 
                     x.kde.mode.idx,                       
                     lambda, round.unit, kernel.n, bins.n.min.act, 
                     x.tr.bins.min, x.tr.prop.min, x.tr.prop.max,
                     df.est,df.con,
                     x.Q1, x.Q2, RL1.p, RL2.p,
                     l.fact, p.fact, r.fact, w.fact, 
                     xlabel, subtitle,
                     figD, figE, figF, figG, figH,
                     lambda.gen, mue.gen, sig.gen,
                     histcol2, bordercol2, kdecol, qqrawcol, polycol, gencol,
                     print.log.message)

    #  There may be no proposal
    if (length(tab1.qqw) == 1)
    {
      #  No result from TIQQ
      tab.qqw <- NA
    } else
    {
      #  There is a result
    
      tab1.qqw    <- RowToMatrix(tab1.qqw)
      tab1.qqw.n <- nrow(tab1.qqw)

      #  CalcPrev provides general names, no reference to calculating method
      tab2.qqw.names <- c("n.l", "prev.l", "n.c", "prev.c", 
                        "n.r", "prev.r", "neg.prev.sum")
      tab2.qqw       <- matrix(NA, nrow=tab1.qqw.n, ncol=length(tab2.qqw.names))
      colnames(tab2.qqw) <- tab2.qqw.names

      #  Next calculation by row 
      for (i in 1:tab1.qqw.n)
      {
        x.tr.n <- sum( (tab1.qqw[i, "x.tr.lo"] <= x) & 
                     (x < tab1.qqw[i, "x.tr.hi"])  )
        x.lt.tr.n <- sum(x <   tab1.qqw[i, "x.tr.lo"])
        x.ge.tr.n <- sum(x >=  tab1.qqw[i, "x.tr.hi"])

        tab2.qqw[i, ] <- CalcPrev(x.n, 
                                x.tr.n, 
                                x.lt.tr.n,
                                x.ge.tr.n, 
                                tab1.qqw[i, "x.tr.lo"], tab1.qqw[i, "x.tr.hi"], 
                                lambda, 
                                tab1.qqw[i, "mue.qqw"], 
                                tab1.qqw[i, "sigma.qqw"])
      }

      #  Make names in tab2.qqw method-specific. Check with first definitions a 
      #  few  lines higher.
      colnames(tab2.qqw) <- c("n.l.qqw", "prev.l.qqw", 
                            "n.c.qqw", "prev.c.qqw",
                            "n.r.qqw", "prev.r.qqw", "neg.prev.sum.qqw")
   
      if (ilambda==1)
      {
        tab.qqw <- data.frame(tab1.qqw, tab2.qqw, stringsAsFactors=FALSE)
      }  else
      {
        tab.qqw <- rbind(tab.qqw,
                       data.frame(tab1.qqw, tab2.qqw, stringsAsFactors=FALSE))
      }
    }
  }

  if (length(tab.qqw) > 1)
  {

    # --------------------------------------------------------------------------
    #  Decide for the best (some of the best) TIQQ result

    #cat("\n [QQW] Initial result from TIQQ, sorted by lambda, start, end\n")
    #print(tab.qqw[order(tab.qqw[ ,"lambda.qqw"],
    #                    tab.qqw[ ,"x.tr.lo"],
    #                    tab.qqw[ ,"x.tr.hi"]),  ])

    #cat("\n [QQW] Initial result from TIQQ, sorted by RL1, RL2 \n")
    #print(tab.qqw[order(tab.qqw[ ,"x.RL1.qqw"],
    #                    tab.qqw[ ,"x.RL2.qqw"]),  ])
  
    #  If necessary, turn vector into matrix
    tab.qqw    <- RowToMatrix(tab.qqw)

    #  Needed if no acceptable solutions
    tab.qqw.backup <- tab.qqw

    #  Provide overview over estimates
    #temp <- matrix(c(min(tab.qqw.backup[ , "prev.l.qqw"]),
    #               max(tab.qqw.backup[ , "prev.l.qqw"]),
    #               min(tab.qqw.backup[ , "prev.c.qqw"]),
    #               max(tab.qqw.backup[ , "prev.c.qqw"]),
    #               min(tab.qqw.backup[ , "prev.r.qqw"]),
    #               max(tab.qqw.backup[ , "prev.r.qqw"]) ), byrow=TRUE,
    #               nrow=3, ncol=2)
    #dimnames(temp) <- list(c("l", "c", "r"), c("min", "max"))
    #cat("\n Characteristics of QQW prevalence estimates\n")
    #print(temp)

    #  Remove estimates with implausible prevalence estimates
    #  Check prev.l estimates
    tab.qqw.n0  <- nrow(tab.qqw)
    ok <- (prev.acc.lo <= tab.qqw[ ,"prev.l.qqw"]) & 
          (tab.qqw[ ,"prev.l.qqw"] <= prev.acc.hi)
    tab.qqw.n <- sum(ok)
    tab.qqw    <- tab.qqw[ok, ] 
    tab.qqw    <- RowToMatrix(tab.qqw)

    #  Check for nothing more left  
  
    if (tab.qqw.n > 0)
    { #  Check prev.c estimates
      tab.qqw.n0  <- nrow(tab.qqw)
      ok <- (prev.acc.lo <= tab.qqw[ ,"prev.c.qqw"]) & 
            (tab.qqw[ ,"prev.c.qqw"] <= prev.acc.hi)
      tab.qqw.n <- sum(ok)
      tab.qqw   <- tab.qqw[ok, ] 
      tab.qqw   <- RowToMatrix(tab.qqw)
    }
  
    if (tab.qqw.n > 0)
    { #  Check prev.r estimates
      tab.qqw.n0  <- nrow(tab.qqw)
      ok <- (prev.acc.lo <= tab.qqw[ ,"prev.r.qqw"]) & 
          (tab.qqw[ ,"prev.r.qqw"] <= prev.acc.hi)
      tab.qqw.n <- sum(ok)
      tab.qqw    <- tab.qqw[ok, ] 
      tab.qqw    <- RowToMatrix(tab.qqw)
    }
  
    #  Emergency action, if no results are left
  
    if (tab.qqw.n == 0)
    { cat("\n ++++++  Method QQW found no acceptable solution",
        "\n ++++++  Modify prev.acc.lo and / or prev.acc.hi in",
        "\n ++++++  TMC_seg015_DefaultSettings.R",
        "\n ++++++  Below: the first up to 10 prevalence estimates by QQW,",
        "\n ++++++  sorted by the sum of negative prevalences and rsq",
        "\n\n")
      tab.qqw.backup <- tab.qqw.backup[order(-tab.qqw.backup[ ,"neg.prev.sum.qqw"],
                                           -tab.qqw.backup[ ,"r2"]), ] 
      imax <- min(10, nrow(tab.qqw.backup))
      print(tab.qqw.backup[1:imax, c("x.RL1.qqw", "x.RL2.qqw",
                                 "prev.l.qqw", "prev.c.qqw", 
                                 "prev.r.qqw", "neg.prev.sum.qqw", "r2")])

      cat("\n ++++++  Below: the first up to 10 prevalence estimates by QQW,",
          "\n ++++++  sorted by rsq",
          "\n ++++++  First solution willl be taken",
          "\n\n")
      tab.qqw.backup <- tab.qqw.backup[order(-tab.qqw.backup[ ,"r2"]), ]

      print(tab.qqw.backup[1:imax, c("x.RL1.qqw", "x.RL2.qqw",
                                 "prev.l.qqw", "prev.c.qqw", 
                                 "prev.r.qqw", "neg.prev.sum.qqw", "r2")])

      tab.qqw <- tab.qqw.backup
      tab.qqw.n <- nrow(tab.qqw)

    } else
    { 
      #  No problems so far
      rm(tab.qqw.backup)
    }  

    #  Sort by decreasing r2
    tab.qqw    <- tab.qqw[order(-tab.qqw[ ,"r2"]), ]

    #  Sort by increasing opt.crit
    #tab.qqw    <- tab.qqw[order(tab.qqw[ ,"opt.crit"]), ]
  
    #cat("\n [QQW] tab.qqw before reduction\n")
    #print(tab.qqw[1:30, ])

    #  Reduce to the first at most ... candidates.
    tab.qqw.n <- min(5, tab.qqw.n)
    tab.qqw   <- tab.qqw[1:tab.qqw.n, ]
    tab.qqw   <- RowToMatrix(tab.qqw)

    # -------------------------------------------------------------------------
    #  Add asymptotic confidence intervals. Calculation uses the actually
    #  used numbers, which may come from a reduced histogram
    for (i in 1:tab.qqw.n)
    {
      xc.RL1.qqw.CI <- CIQuant.PNV(tab.qqw[i, "x.tr.n"], RL1.p, 
                                 tab.qqw[i, "x.RL1.qqw"],
                                 tab.qqw[i, "lambda.qqw"], 
                                 tab.qqw[i, "mue.qqw"], 
                                 tab.qqw[i, "sigma.qqw"],
                                 alpha=alpha)
      tab.qqw[i, "x.RL1.cilo"] <- xc.RL1.qqw.CI[1]
      tab.qqw[i, "x.RL1.cihi"] <- xc.RL1.qqw.CI[2]

      xc.RL2.qqw.CI <- CIQuant.PNV(tab.qqw[i, "x.tr.n"], RL1.p, 
                                 tab.qqw[i, "x.RL2.qqw"],
                                 tab.qqw[i, "lambda.qqw"], 
                                 tab.qqw[i, "mue.qqw"], 
                                 tab.qqw[i, "sigma.qqw"],
                                 alpha=alpha)
      tab.qqw[i, "x.RL2.cilo"] <- xc.RL2.qqw.CI[1]
      tab.qqw[i, "x.RL2.cihi"] <- xc.RL2.qqw.CI[2]
    }

    # -------------------------------------------------------------------------
    if (!is.na(figI))
    { dev.set(figI)
      #bringToTop(figI)

      #  Show QQ plot for optimal truncation interval
      lambda.qqw.opt <- tab.qqw[1, "lambda.qqw"]
      qqnorm(BoxCox(x.red, lambda.qqw.opt),
           main=paste("Truncation interval according to QQ plot, lambda =",
                      lambda.qqw.opt),
           sub=subtitle, cex.sub=0.7)
      abline(h=BoxCox(tab.qqw[1, "x.tr.lo"], lambda.qqw.opt), col="blue")   
      abline(h=BoxCox(tab.qqw[1, "x.tr.hi"], lambda.qqw.opt), col="blue")
      abline(c(tab.qqw[1, "mue.qqw"], tab.qqw[1, "sigma.qqw"]), col="red")   
      abline(h=y.kde.mode, col=kdecol, lty=2)   

      #  savePlot(file=figI.file, type=figtype)  #  01.05.2021
    }

    # -------------------------------------------------------------------------
    #  Plot figJ moved to calling programme  01.05.2021
    #  @@@ move to calling program?

  }  

  # cat("\n *** End QQW \n")
  return(tab.qqw=tab.qqw) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_Quantile.R
#  
#  (c) wwosniok@math.uni-bremen.de


#  19.05.2021 WW  Check for empty input vector
#  16.03.2021 WW  Determination of round.unit is always necessary
#  14.03.2021 WW  Linear interpolation instead of parabola interpolation
#                 Continuity connection added
#  19.12.2020 WW  New interpretation of input data: x is considered
#                 as a rounded value ==> for a given x, the true value
#                 lies in the interval [x-ru/2, x+ru/2), where ru is the
#                 rounding unit. Due to rounding, ties are likely, therefore
#                 quantiles are calculated from the polygon connecting
#                 (x[1]-ru/2, 0), (x[1], cdf[1]),  (x[2], cdf[2]), ...
#                 (x[n], cdf[n]), where cdf[i] is the # of x values <= x[i].
#                 x[1]-ru/2 may be negative - no precaution for using negative
#                 results so far.  
#  09.07.2019 WW  Behandlung des ersten und letzten Intervalls geändert.
#                 Interpolation verwendet jetzt die korrekte Definition
#                 der empirischen Verteilungsfunktion. Kein Problem mehr,
#                 wenn p-Werte unterhalb / oberhalb der empirischen ps liegen.
#  17.08.2018 WW  Stetigkeitskorrektur eingeführt. Die interpolierende
#                 Linie verbindet die Mittelpunkte der Treppenstufen,
#                 nicht die Eckpunkte.
#  22.11.2017 WW  Einrichtung

# ===========================================================================
Quantile <- function(x, probs=c(0.25,0.50,0.75), na.rm=FALSE, type=0)
{ #  Berechnet Quantile durch Interpolation der empirischen 
  #  Verteilungsfunktion, auch bei Bindungen (das macht quantile() nicht!)
  #  Mit Stetigkeitskorrektur.
  #  EINGABE
  #  x     - Daten
  #  probs - Wahrscheinlichkeiten \in [0,1] zu den gewünschten Quantilen
  #  na.rm - falls TRUE, werden fehlenden Werte aus x entfernt
  #  type  - ohne Bedeutung, nur zur Kompatitibilität mit Aufruf von quantile()
  #  AUSGABE
  #  Quantile mit den Wahrscheinlichkeiten in % als Namen
  # ==========================================================================

  probs.n <- length(probs)

  if (length(x) > 1)
  { # At least 2 values, try execution 
    
    #  Remove missing vales, if required
    if (na.rm) { x <- x[!is.na(x)] }


    #  If there still are NAs, stop, otherwise execute
    if (any(is.na(x)) & !na.rm ) 
    { 
      cat("In Quantile(): x has NA and na.rm == FALSE \n")
      q <- rep(NA, times=probs.n) 
    } else
    { #  No (more) NAs
      x.table <- table(x)
      x.val   <- as.numeric(names(x.table))
      x.val.n <- length(x.val)
      x.val.cdf <- cumsum(x.table)/length(x)

      q <- rep(NA,times=probs.n)
      if (x.val.n == 1)
      {
        # Only 1 unique value  
        q <- x.val 
      } else
      { #  x.val.n > 1
        #  More than 1 unique value
        #  Calculate coordinates for interpolation

        intx <- x.val
        inty <- x.val.cdf

        #  Continuity correction
        round.unit <- min(diff(intx))
        intx       <- intx + round.unit/2 
      
        # If a quantile for a probability < x.val.cdf[1] is requested,
        # we need to construct x[1] - ru
        if (any(probs < x.val.cdf[1]))
        { # yes, find the rounding unit and  add left limit
          intx <- c(intx[1] - round.unit/2, intx)            
          inty <- c(0, inty)
        }   

        for (i in 1:probs.n)
        { 
          #  Interpolate between the shifted jump positions
         
          q[i] <- Interpolx(intx, inty, probs[i]) 
        }
      }     #  if (x.val.n == 1)
    }       #  if (any(is.na(x)) & !na.rm ) 
  
  } else
  {
    #  Too few data
    q <- rep(NA, times=probs.n) 

  }         # if (length(x) > 1)

  names(q) <- paste(100*probs,"%",sep="")
  return(q)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test see TMC6/dev/Test_Quantile.R 
#
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_R2o.R

R2o <- function(ttheta,RB)
{ #  Auf ganz R  transformierte Parameter --> beschränkte Originalskala 
  #  16.12.2017 

  theta.n <- length(ttheta)
  theta   <- rep(NA,times=theta.n)
  
  for (k in 1:theta.n) 
  { 
    #  Schutz vor Überlauf
    ttheta.eff <- min(ttheta[k],50)
    ttheta.eff <- max(-50,ttheta.eff)
    theta[k] <- RB[k,1] + (RB[k,2]-RB[k,1])/(1+exp(-ttheta.eff))
  }

  #cat("R2o: ttheta,theta\n")
  #print(cbind(ttheta,theta))

  return(theta)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test
#RB <- matrix(c(-1.e-10,1+1.e10),nrow=1,ncol=2)
#theta     <- 1-1.e-5
#ttheta    <- o2R(theta,RB,fastnull=1.e-10)
#theta.rec <- R2o(ttheta,RB)
#cat(theta,ttheta,theta.rec,"\n")

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
#  F_reldist.R
#
#  (c) wwosniok@math.uni-bremen.de  
#
#  Calculate relative distance between xtest and xref
#
#  04.01.2020

# =============================================================================
reldist <- function(xtest,xref)
{ 
  temp <- (xtest-xref)/xref
  reldist <- sqrt(temp %*% temp)
  return(reldist)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_RL.age.spline.R

#
#  Show course of RL1, RL2 estimated by method RL.age.meth, vs 
#  mean of age group
#  Smooth by spline function
#  Separation by sex
#  Show tolerance intervals, if present, otherwise asymptotic CIs
#
#  08.05.2021  Estimation method added as call parameter
#  30.04.2021  x.clip.@@ correctly processed
#              par.tcl and par.las as parameters   
#  02.12.2020  Names x.RL, ..  changed to RL, marker for unsafe results
#              changed to "!!"
#  21.04.2020  ylab changed back to 'Reference limits'
#  24.01.2020  Indicator for successful plotting addded to output
#  10.01.2020  External parameters for the RL axis
#  08.01.2020  Var names changed for TMC4
#  12.12.2019  Mark plots in gray, if solution type != 2
#  26.07.2019  Plotten auch, wenn errcode nicht " "
#  08.05.2019  Plotten nur, wenn errcode nicht " "
#              Definition des errcode in analysis geändert
#  17.04.2019  Plotten auch dann, wenn errrcode nicht " "
#              (Überlegen, ob wackelige Schätzungen irgendwie markiert
#               werden können)  
#  02.03.2019  ylab geändert
#  02.05.2018
#  
#  --------------------------------------------------------------------------

RL.age.spline <- function(fig,gtab,sexcode,Verlauf.Meth,shift,
                           xplotmin,xplotmax,yplotmin,yplotmax,
                           yplotby1,yplotby2,  
                           xlabel, 
                           x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                           par.tcl, par.las, RL.age.meth,
                           outname.stra,methcol,methpch,age.limits,
                           arrow.length=0.03,newplot=TRUE, 
                           unsafe.col="dimgray")
{ #  
  #  fig
  #  figfile
  #  gtab
  #  sexcode       
  #  Verlauf.Meth
  #  xplotmin
  #  xplotmax
  #  yplotmin
  #  yplotmax
  #  outname.stra
  #  methcol
  #  methpch
  #  age.limits
  #  arrow.length  arrow length in inches! Default figure width is 7 in

  #  Additional style parameters
  arrlwd <- 1    #  line width for arrows (confidence bars) 
  pchlwd <- 1    #  line width for point symbols 
  spllwd <- 1    #  line width for spline lines 
  pchcex <- 1    #  point symbol size
     
  dev.set(fig)
  #bringToTop(fig)

  #  Select data. Note that some values are formatted (=> of character type).
  errcode  <- gtab[ ,"errcode"]
  age.mea  <- as.numeric(gtab[ ,"Age.mea"])
  RL1.hat  <- as.numeric(gtab[ ,"RL1"])
  RL2.hat  <- as.numeric(gtab[ ,"RL2"])
  RL1.cilo <- as.numeric(gtab[ ,"RL1.cilo"])
  RL1.cihi <- as.numeric(gtab[ ,"RL1.cihi"])
  RL1.tilo <- as.numeric(gtab[ ,"RL1.tilo"])
  RL1.tihi <- as.numeric(gtab[ ,"RL1.tihi"])
  RL2.cilo <- as.numeric(gtab[ ,"RL2.cilo"])
  RL2.cihi <- as.numeric(gtab[ ,"RL2.cihi"])
  RL2.tilo <- as.numeric(gtab[ ,"RL2.tilo"])
  RL2.tihi <- as.numeric(gtab[ ,"RL2.tihi"])

  use     <- (!is.na(age.mea)) & 
             (!is.na(RL1.hat)) & (!is.na(RL2.hat))
  n       <- sum(use)

  plot.ok <- (n >= 2)

  if (plot.ok)
  { #  Plotting is worthwhile
    #  Select data

    errcode  <- errcode[use]
    age.mea  <- age.mea[use]
    RL1.hat  <- RL1.hat[use]
    RL2.hat  <- RL2.hat[use]
    RL1.cilo <- RL1.cilo[use]
    RL1.cihi <- RL1.cihi[use]
    RL1.tilo <- RL1.tilo[use]
    RL1.tihi <- RL1.tihi[use]
    RL2.cilo <- RL2.cilo[use]
    RL2.cihi <- RL2.cihi[use]
    RL2.tilo <- RL2.tilo[use]
    RL2.tihi <- RL2.tihi[use]

    if (newplot)
    { 
      #  Punkte eintragen
      #  RL1
      #  No automatic x legend

      par(las=par.las)
      par(tcl=par.tcl)
      par(xaxt="n") 
      if (!is.na(x.clip.by1))  
      {  # User values are given 
        par(xaxt="n") 
        xplotmin.eff <- x.clip.min
        xplotmax.eff <- x.clip.max
      } else 
      {  # User values are not given 
        xplotmin.eff <- xplotmin
        xplotmax.eff <- xplotmax
      }
  
      if (!is.na(yplotby1))    { par(yaxt="n") }  # User values are given 
 
      plot(age.mea+shift,RL1.hat,type="p",
          col=methcol,pch=methpch,cex=pchcex,lwd=pchlwd,
          xlim=c(xplotmin.eff,xplotmax.eff),
          ylim=c(yplotmin,yplotmax),
          xlab="Age",ylab=paste("Reference limits, method =", RL.age.meth) )

      # Mark unsafe estimates with errcode == "!!" in extra color (gray) 
      points(age.mea[errcode=="!!"]+shift,RL1.hat[errcode=="!!"],
             type="p",col=unsafe.col,pch=methpch,cex=pchcex,lwd=pchlwd) 

      #  Konfidenz- / Toleranzintervslle
      if (all(is.na(RL1.tilo)))
      { #  Kein Bootstrap-Ergebnis - asymptotische Intervalle nehmen
        #cat("[SplineVerlauf3] (1) Asymptotisches CI für RL1 verwendet\n")
        ci.type <- "CI"
        arrows(age.mea+shift,RL1.cilo,
               age.mea+shift,RL1.cihi,
               code=3,angle=90,length=arrow.length,col=methcol,lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe  
        arrows(age.mea[errcode=="!!"]+shift,RL1.cilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.cihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }  else
      { #  Bootstrap-Ergebnis vorhanden - nehmen
        # cat("[SplineVerlauf3] (2) Bootstrap-TI für RL1 verwendet\n")
        ci.type <- "TI"
        arrows(age.mea+shift,RL1.tilo,
               age.mea+shift,RL1.tihi,
               code=3,angle=90,length=arrow.length,col="magenta",lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe
        arrows(age.mea[errcode=="!!"]+shift,RL1.tilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.tihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }
      
      #  User-defined tickmarks for x
      if (!is.na(x.clip.by1))  
      {
        par(xaxt="s")
        axis(1,at=seq(x.clip.min,x.clip.max,by=x.clip.by1),labels=TRUE)
        axis(1,at=seq(x.clip.min,x.clip.max,by=x.clip.by2),labels=FALSE)
      }

      #  User-defined tickmarks for y
      if (!is.na(yplotby1))  
      { 
        par(yaxt="s") 
        axis(2,at=seq(yplotmin,yplotmax,by=yplotby1),labels=TRUE)
        if (!is.na(yplotby2))
        {
          axis(2,at=seq(yplotmin,yplotmax,by=yplotby2),labels=FALSE)
        }
      }

    }  else
    {
      #  Plot existiert schon
      #  Punkte eintragen
      #  RL1 

      points(age.mea+shift,RL1.hat,type="p",
             col=methcol,pch=methpch,cex=pchcex,lwd=pchlwd)

      # Mark unsafe estimates with errcode == "!!" in extra color (gray) 
      points(age.mea[errcode=="!!"]+shift,RL1.hat[errcode=="!!"],
             type="p",col=unsafe.col,pch=methpch,cex=pchcex,lwd=pchlwd) 

      #  Konfidenz- / Toleranzintervslle
      if (all(is.na(RL1.tilo)))
      { #  Kein Bootstrap-Ergebnis - asymptotische Intervalle nehmen
        #cat("[SplineVerlauf3] (1) Asymptotisches CI für RL1 verwendet\n")
        ci.type <- "CI"
        arrows(age.mea+shift,RL1.cilo,
               age.mea+shift,RL1.cihi,
               code=3,angle=90,length=arrow.length,col=methcol,lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe  
        arrows(age.mea[errcode=="!!"]+shift,RL1.cilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.cihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }  else
      { #  Bootstrap-Ergebnis vorhanden - nehmen
        # cat("[SplineVerlauf3] (2) Bootstrap-TI für RL1 verwendet\n")
        ci.type <- "TI"
        arrows(age.mea+shift,RL1.tilo,
               age.mea+shift,RL1.tihi,
               code=3,angle=90,length=arrow.length,col="magenta",lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe
        arrows(age.mea[errcode=="!!"]+shift,RL1.tilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.tihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }
    }   # newplot

    #  Punkte eintragen
    #  RL2 

    points(age.mea+shift,RL2.hat,type="p",
           col=methcol,pch=methpch,cex=pchcex,lwd=pchlwd)

    # Mark unsafe estimates  in extra color (gray) 
    points(age.mea[errcode=="!!"]+shift,RL2.hat[errcode=="!!"],
           type="p",col=unsafe.col,pch=methpch,cex=pchcex,lwd=pchlwd) 

    #  Konfidenz- / Toleranzintervslle
    if (all(is.na(RL2.tilo)))
    { #  Kein Bootstrap-Ergebnis - asymptotische Intervalle nehmen
      #cat("[SplineVerlauf3] (3) Asymptotisches CI für RL2 verwendet\n")
      arrows(age.mea+shift,RL2.cilo,
             age.mea+shift,RL2.cihi,
             code=3,angle=90,length=arrow.length,col=methcol,lwd=arrlwd)

      #  Mark CI in gray, if solution type != 2  
      arrows(age.mea[errcode=="!!"]+shift,RL2.cilo[errcode=="!!"],
             age.mea[errcode=="!!"]+shift,RL2.cihi[errcode=="!!"],
             code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)

    }  else
    { #  Bootstrap-Ergebnis vorhanden - nehmen
      #cat("[SplineVerlauf3] (4) Bootstrap-TI für RL2 verwendet\n")
      arrows(age.mea+shift,RL2.tilo,
             age.mea+shift,RL2.tihi,
             code=3,angle=90,length=arrow.length,col="magenta",lwd=arrlwd)

      #  Mark CI in gray, if solution type != 2  
      arrows(age.mea[errcode=="!!"]+shift,RL2.tilo[errcode=="!!"],
             age.mea[errcode=="!!"]+shift,RL2.tihi[errcode=="!!"],
             code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
    }

    # Legende am Ende, wenn CI-Typ bekannt ist

  } # if (plot.ok)

  # .........................................................................
  #  Spline-Berechnung

  if (n >= 4 )
  { #  Spline ist möglich  

    age.out <- seq(Floor(age.mea[1],1),Ceiling(age.mea[length(age.mea)],1),
                         by=1)

    # RL1.sspline.gcv <- smooth.spline(age.mea, RL1.hat)

    if (( 4 <= n) & (n <= 10)) { df.user <- Round(0.9 * n, 1) }
    if ((11 <= n) & (n <= 12)) { df.user <- Round(0.8 * n, 1) }
    if ((13 <= n) & (n <= 14)) { df.user <- Round(0.6 * n, 1) }
    # if ((15 <= n)            ) { df.user <- RL1.sspline.gcv$df }
    if ((15 <= n)            ) { df.user <- Round(0.5 * n, 1) }

    RL1.sspline.df      <- smooth.spline(age.mea, RL1.hat,df=df.user)
    RL1.sspline.df.pred <- predict(RL1.sspline.df,age.out)

    #  Spline eintragen
    lines(RL1.sspline.df.pred$x+shift,RL1.sspline.df.pred$y,
          col=methcol,lwd=spllwd)

 
    # .........................................................................
    #  RL2 

    RL2.sspline.df      <- smooth.spline(age.mea, RL2.hat,df=df.user)
    RL2.sspline.df.pred <- predict(RL2.sspline.df,age.out)

    #  Spline eintragen
    lines(RL2.sspline.df.pred$x+shift,RL2.sspline.df.pred$y,
          col=methcol,lwd=spllwd)

    return(list(RL.emp=data.frame(age=age.mea,RL1=RL1.hat,RL2=RL2.hat,
                                  RL1.cilo=RL1.cilo,RL2.cilo=RL2.cilo,
                                  RL1.tilo=RL1.tilo,RL2.tilo=RL2.tilo),
                RL.smo=data.frame(age=age.out,
                                  RL1.hut=RL1.sspline.df.pred$y,
                                  RL2.hut=RL2.sspline.df.pred$y),
                plot.ok=plot.ok))

  } else  #  if n>= 4 ..
  { 
    #  Spline nicht möglich, da zu wenig Stützpunkte
    cat(  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
        "\n  No spline smoothing for sex =",sexcode,
        ", because < 4 support points",
        "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

    return(list(RL.emp=data.frame(age=NA, RL1=NA,RL2=NA,
                                  RL1.cilo=NA, RL2.cilo=NA,
                                  RL1.tilo=NA, RL2.tilo=NA),
                RL.smo=data.frame(age=NA,
                                  RL1.hut=NA,
                                  RL2.hut=NA),
                plot.ok=plot.ok))
  }
} 

# ============================================================================
# *****************************************************************************
#  F_Round.R

Round <- function(x,unit)
{ # Rounds in standard manner (round() proceeds differently!)
  
  x.sign  <- sign(x)
  x.round <- unit * trunc((abs(x)+0.5*unit)/unit + 1.e-15)
  x.round <- x.round * x.sign
  return(x.round)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#xxx <- seq(-0.50,0.50,by=0.01)
#xxxR <- Round(xxx,0.10)
#cbind(xxx,xxxR)
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


# ****************************************************************************
#  F_RowToMatrix.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Make input object a matrix / data frame, if possible. Needed after 
#  extracting a single row from a matrix / data frame, which produces a vector.
#  If input object is NULL, return NULL
#  If input object is a vector, return a 1-row matrix with column names taken from 
#  the input
#  If input object is a matrix / data.frame, return it as it is
 
#  #  CHANGE HISTORY
#  09.09.2020 Treatmetn of data frames added 
#  09.03.2020 Start
#
# ****************************************************************************

RowToMatrix <- function(obj)
{
  #  INPUT 
  #  obj   a scalar, a vector, or a matrix

  #  OUTPUT 
  #  obj.out   a 1-row matrix, if possible, otherwise NA
  # ---------------------------------------------------------------------------

  #cat("\n [RowToMatrix] mode(obj)      ", mode(obj) ) 
  #cat("\n [RowToMatrix] class(obj)     ", class(obj) ) 
  #cat("\n [RowToMatrix] typeof(obj)    ", typeof(obj) ) 
  #cat("\n [RowToMatrix] length(obj)    ", length(obj) ) 
  #cat("\n [RowToMatrix] is.matrix(obj) ", is.matrix(obj), "\n" ) 

  if ( (length(obj) > 0) && (is.matrix(obj) & !is.data.frame(obj)) )  
  { 
    obj.out <- obj
  }

  if ( (length(obj) > 0) && (!is.matrix(obj) & is.data.frame(obj)) )  
  { 
    obj.out <- obj
  }

  if ( (length(obj) > 0) && (!is.matrix(obj) & !is.data.frame(obj)) )
  { 
    cnames  <- names(obj)
    obj.out <- matrix(obj, nrow=1)
    colnames(obj.out) <- cnames
  }

  if ( is.null(obj) )
  { 
    obj.out <- obj
  }

  # cat("\n [RowToMatrix] class(obj.out)    ", class(obj.out), "\n") 

  return(obj.out)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test
#A <- matrix(c(11, 12, 13,
#              21, 22, 23), byrow=TRUE, nrow=2)
#dimnames(A) <- list(c("r1", "r2"), c("c1", "c2", "c3"))
## A <- data.frame(A)
#print(A)

#length(A)
#is.matrix(A)
#is.data.frame(A)

#B <- RowToMatrix(A)
#print(B)

#C <- A[1, ]
#print(C)

#D <- RowToMatrix(C)
#print(D)

#E <- NA
#is.null(E)
#print(E)
#F <- RowToMatrix(E)
#print(F)

#G <- NULL
#print(G)
#H <- RowToMatrix(G)
#print(H)
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_SmoothHist1_V5.R
#
#  Smooth a histogram by 5 point moving average
#
#  23.02.2021 Five point weighted parabola smoothing 
# ============================================================================

SmoothHist1 <- function(x.hist,
                        weights1=      c(1, 0, 0),
                        weights2=   c(0, 1, 0, 0),
                        weights3=c(1, 1, 1, 1, 1),
                        figA=NA )
{
  density1   <- x.hist$density
  density1.n <- length(x.hist$density)

  #  No smoothing if less than 5 bins

  if (density1.n < 5)
  { 
    x.hist2 <- x.hist
  } else
  {
    x.n  <- sum(x.hist$counts) 

    #   Normalize weights
    weights1 <- weights1/sum(weights1)
    weights2 <- weights2/sum(weights2)
    weights3 <- weights3/sum(weights3)

    breaks1   <- x.hist$breaks
    breaks1.n <- length(x.hist$breaks)

    #  Leave first and last 2 points unchanged
    #  Smooth interior counts
    density2 <- rep(NA, times=density1.n)
    density2[1:2] <- x.hist$density[1:2]
    density2[(density1.n-1):density1.n] <- tail(x.hist$density, 2)

    for (i in 3:(density1.n-2))
    {
      subset <- (i-2):(i+2)
      x.sub  <- x.hist$mids[subset]
      x.sub2 <- x.sub^2
      y.sub  <- x.hist$density[subset]
      xy.sub.lm <- lm(y.sub ~ x.sub + x.sub2, weight=weights3)
      xy.sub.lm.pred <- predict(xy.sub.lm)
      density2[i] <- xy.sub.lm.pred[3] 
    } 

    #  Normalize smoothed density
    bin.width    <- breaks1[2:breaks1.n] - breaks1[1:(breaks1.n-1)]
    density2.sum <- sum(density2 * bin.width)
    density2     <- density2 / density2.sum

    density2.sum <- sum(density2 * bin.width)

    #  Calculate smoothed counts
    counts2   <- x.n * density2 * bin.width 

    x.hist2 <- list(breaks=breaks1, 
                    counts=counts2, 
                    density=density2, 
                    mids=(breaks1[2:breaks1.n] + 
                          breaks1[1:(breaks1.n-1)])/2 )
    class(x.hist2) <- "histogram"

    if (!is.na(figA))
    { 
      dev.set(figA)
      plot(x.hist, col="cornsilk")
      plot(x.hist2, border="green3", add=TRUE)
    }
  }
  return(x.hist2)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_SmoothHist2.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Smooth a histogram by using the kernel density estimate as
#  smoothing tool

#  CHANGE HISTORY
#  04.01.2021 Concept is problematic, if many bins / small differences
#             between density values. Stop using it for the moment.
#  20.11.2020 Name changed to SmoothHist2 (TMC has already SmoothHist2) 
#  10.10.2020 Calculation of x.n added
#  26.09.2020 Colours as function input
#  25.09.2020 Start
# ==========================================================================

SmoothHist2 <- function(x.kde, x.hist, FigA, histcol2, bordercol2, 
                       histcol3, bordercol3, kdecol, subtitle)
{  
  #
  #  INPUT
  #  x.kde       kernel density estimate of the raw data
  #  x.hist      the histogram to smooth
  #  FigA        window for plotting. Must exist. NA: no plot
  #  subtitle

  # ==========================================================================
  #

  #  Produce a smoothed version of the input histogram
  x.n <- sum(x.hist$counts)

  x.kde.cdf <- cumsum(x.kde$y)
  x.kde.cdf <- x.kde.cdf / tail(x.kde.cdf, 1)

  x.hist.tab   <- hist.table(x.hist$breaks, x.hist$counts)
  x.hist.tab.n <- nrow(x.hist.tab)
  prop.kde     <- rep(NA, times=x.hist.tab.n)
  
  for (i in 1:x.hist.tab.n)
  {
    prop.kde[i] <- max(x.kde.cdf[x.kde$x < x.hist.tab[i, "x.hi"]])
  }
  prop.kde[2:x.hist.tab.n] <- prop.kde[2:x.hist.tab.n] - 
                              prop.kde[1:(x.hist.tab.n-1)] 

  prop.kde[x.hist.tab.n] <- 1 - sum(prop.kde[1:(x.hist.tab.n-1)])

  x.hist.tab <- cbind(x.hist.tab, prop.kde=prop.kde) 
  x.hist.tab <- cbind(x.hist.tab, count.kde=prop.kde*x.n) 
  x.hist.tab <- cbind(x.hist.tab, density.kde=prop.kde/
                             (x.hist.tab[ , "x.hi"] - x.hist.tab[ , "x.lo"]) )

  x.hist.smo <- x.hist
  x.hist.smo$counts <- x.hist.tab[ , "count.kde"]
  x.hist.smo$density <- x.hist.tab[ , "density.kde"]

  if (!is.na(FigA))
  {
    #  Show original (unsmoothed) histogram
    dev.set(FigA)
    plot(x.hist, freq=FALSE, xlim=c(0, 100), col=histcol1, border=bordercol1,
         main="Original and smoothed histogram", 
         sub=subtitle, cex=0.7)
    lines(x.kde,col=kdecol) 

    plot(x.hist.smo, freq=FALSE, add=TRUE, border=bordercol3, 
         col=histcol3)

   legend("topright", c("kde", "Original histogram", "Smoothed histogram"), 
          col=c("black", bordercol2, bordercol3), 
          lwd=3, lty=1, cex=0.7)
  }

  return(x.hist.smo)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test   

#hs <- SmoothHist2(x.kde, x.hist.eqd, 5, histcol3, bordercol3, 
#                  histcol1, bordercol1, "chartreuse", subtitle)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# *****************************************************************************
#  F_Splineverlauf4.R

#
#  Verlauf von RL1, RL2 gegen Mittelwert der Altersgruppe
#  Spline-Glättung
#  Glättungsparameter bei kleinem n vorgeben
#  Nur Ergebnisse einer Methode
#  Eigene Kurven für F / M 
#  Annahme: wenn es RL1.hat gibt, gibt es auch RL2.hat
#  Darstellung der Toleranzintervalle, wenn vorhanden, sonst asymp. CI 
#
#  30.04.2021  x.clip.@@ correctly processed
#              par.tcl and par.las as parameters   
#  02.12.2020  Names x.RL, ..  changed to RL, marker for unsafe results
#              changed to "!!"
#  21.04.2020  ylab changed back to 'Reference limits'
#  24.01.2020  Indicator for successful plotting addded to output
#  10.01.2020  External parameters for the RL axis
#  08.01.2020  Var names changed for TMC4
#  12.12.2019  Mark plots in gray, if solution type != 2
#  26.07.2019  Plotten auch, wenn errcode nicht " "
#  08.05.2019  Plotten nur, wenn errcode nicht " "
#              Definition des errcode in analysis geändert
#  17.04.2019  Plotten auch dann, wenn errrcode nicht " "
#              (Überlegen, ob wackelige Schätzungen irgendwie markiert
#               werden können)  
#  02.03.2019  ylab geändert
#  02.05.2018
#  
#  --------------------------------------------------------------------------

SplineVerlauf4 <- function(fig,gtab,sexcode,Verlauf.Meth,shift,
                           xplotmin,xplotmax,yplotmin,yplotmax,
                           yplotby1,yplotby2,  
                           xlabel, 
                           x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                           par.tcl, par.las,
                           outname.stra,methcol,methpch,age.limits,
                           arrow.length=0.03,newplot=TRUE, unsafe.col="dimgray")
{ #  
  #  fig
  #  figfile
  #  gtab
  #  sexcode       
  #  Verlauf.Meth
  #  xplotmin
  #  xplotmax
  #  yplotmin
  #  yplotmax
  #  outname.stra
  #  methcol
  #  methpch
  #  age.limits
  #  arrow.length  arrow length in inches! Default figure width is 7 in

  #  Additional style parameters
  arrlwd <- 1    #  line width for arrows (confidence bars) 
  pchlwd <- 1    #  line width for point symbols 
  spllwd <- 1    #  line width for spline lines 
  pchcex <- 1    #  point symbol size
     
  dev.set(fig)
  #bringToTop(fig)

  #  Select data. Note that some values are formatted (=> of character type).
  errcode  <- gtab[ ,"errcode"]
  age.mea  <- as.numeric(gtab[ ,"Age.mea"])
  RL1.hat  <- as.numeric(gtab[ ,"RL1"])
  RL2.hat  <- as.numeric(gtab[ ,"RL2"])
  RL1.cilo <- as.numeric(gtab[ ,"RL1.cilo"])
  RL1.cihi <- as.numeric(gtab[ ,"RL1.cihi"])
  RL1.tilo <- as.numeric(gtab[ ,"RL1.tilo"])
  RL1.tihi <- as.numeric(gtab[ ,"RL1.tihi"])
  RL2.cilo <- as.numeric(gtab[ ,"RL2.cilo"])
  RL2.cihi <- as.numeric(gtab[ ,"RL2.cihi"])
  RL2.tilo <- as.numeric(gtab[ ,"RL2.tilo"])
  RL2.tihi <- as.numeric(gtab[ ,"RL2.tihi"])

  use     <- (!is.na(age.mea)) & 
             (!is.na(RL1.hat)) & (!is.na(RL2.hat))
  n       <- sum(use)

  plot.ok <- (n >= 2)

  if (plot.ok)
  { #  Plotting is worthwhile
    #  Select data

    errcode  <- errcode[use]
    age.mea  <- age.mea[use]
    RL1.hat  <- RL1.hat[use]
    RL2.hat  <- RL2.hat[use]
    RL1.cilo <- RL1.cilo[use]
    RL1.cihi <- RL1.cihi[use]
    RL1.tilo <- RL1.tilo[use]
    RL1.tihi <- RL1.tihi[use]
    RL2.cilo <- RL2.cilo[use]
    RL2.cihi <- RL2.cihi[use]
    RL2.tilo <- RL2.tilo[use]
    RL2.tihi <- RL2.tihi[use]

    if (newplot)
    { 
      #  Punkte eintragen
      #  RL1
      #  No automatic x legend

      par(las=par.las)
      par(tcl=par.tcl)
      par(xaxt="n") 
      if (!is.na(x.clip.by1))  
      {  # User values are given 
        par(xaxt="n") 
        xplotmin.eff <- x.clip.min
        xplotmax.eff <- x.clip.max
      } else 
      {  # User values are not given 
        xplotmin.eff <- xplotmin
        xplotmax.eff <- xplotmax
      }
  
      if (!is.na(yplotby1))    { par(yaxt="n") }  # User values are given 
 
      plot(age.mea+shift,RL1.hat,type="p",
          col=methcol,pch=methpch,cex=pchcex,lwd=pchlwd,
          xlim=c(xplotmin.eff,xplotmax.eff),
          ylim=c(yplotmin,yplotmax),
          xlab="Age",ylab="Reference limits")

      # Mark unsafe estimates with errcode == "!!" in extra color (gray) 
      points(age.mea[errcode=="!!"]+shift,RL1.hat[errcode=="!!"],
             type="p",col=unsafe.col,pch=methpch,cex=pchcex,lwd=pchlwd) 

      #  Konfidenz- / Toleranzintervslle
      if (all(is.na(RL1.tilo)))
      { #  Kein Bootstrap-Ergebnis - asymptotische Intervalle nehmen
        #cat("[SplineVerlauf3] (1) Asymptotisches CI für RL1 verwendet\n")
        ci.type <- "CI"
        arrows(age.mea+shift,RL1.cilo,
               age.mea+shift,RL1.cihi,
               code=3,angle=90,length=arrow.length,col=methcol,lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe  
        arrows(age.mea[errcode=="!!"]+shift,RL1.cilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.cihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }  else
      { #  Bootstrap-Ergebnis vorhanden - nehmen
        # cat("[SplineVerlauf3] (2) Bootstrap-TI für RL1 verwendet\n")
        ci.type <- "TI"
        arrows(age.mea+shift,RL1.tilo,
               age.mea+shift,RL1.tihi,
               code=3,angle=90,length=arrow.length,col="magenta",lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe
        arrows(age.mea[errcode=="!!"]+shift,RL1.tilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.tihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }
      
      #  User-defined tickmarks for x
      if (!is.na(x.clip.by1))  
      {
        par(xaxt="s")
        axis(1,at=seq(x.clip.min,x.clip.max,by=x.clip.by1),labels=TRUE)
        axis(1,at=seq(x.clip.min,x.clip.max,by=x.clip.by2),labels=FALSE)
      }

      #  User-defined tickmarks for y
      if (!is.na(yplotby1))  
      { 
        par(yaxt="s") 
        axis(2,at=seq(yplotmin,yplotmax,by=yplotby1),labels=TRUE)
        if (!is.na(yplotby2))
        {
          axis(2,at=seq(yplotmin,yplotmax,by=yplotby2),labels=FALSE)
        }
      }

    }  else
    {
      #  Plot existiert schon
      #  Punkte eintragen
      #  RL1 

      points(age.mea+shift,RL1.hat,type="p",
             col=methcol,pch=methpch,cex=pchcex,lwd=pchlwd)

      # Mark unsafe estimates with errcode == "!!" in extra color (gray) 
      points(age.mea[errcode=="!!"]+shift,RL1.hat[errcode=="!!"],
             type="p",col=unsafe.col,pch=methpch,cex=pchcex,lwd=pchlwd) 

      #  Konfidenz- / Toleranzintervslle
      if (all(is.na(RL1.tilo)))
      { #  Kein Bootstrap-Ergebnis - asymptotische Intervalle nehmen
        #cat("[SplineVerlauf3] (1) Asymptotisches CI für RL1 verwendet\n")
        ci.type <- "CI"
        arrows(age.mea+shift,RL1.cilo,
               age.mea+shift,RL1.cihi,
               code=3,angle=90,length=arrow.length,col=methcol,lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe  
        arrows(age.mea[errcode=="!!"]+shift,RL1.cilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.cihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }  else
      { #  Bootstrap-Ergebnis vorhanden - nehmen
        # cat("[SplineVerlauf3] (2) Bootstrap-TI für RL1 verwendet\n")
        ci.type <- "TI"
        arrows(age.mea+shift,RL1.tilo,
               age.mea+shift,RL1.tihi,
               code=3,angle=90,length=arrow.length,col="magenta",lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe
        arrows(age.mea[errcode=="!!"]+shift,RL1.tilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.tihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }
    }   # newplot

    #  Punkte eintragen
    #  RL2 

    points(age.mea+shift,RL2.hat,type="p",
           col=methcol,pch=methpch,cex=pchcex,lwd=pchlwd)

    # Mark unsafe estimates  in extra color (gray) 
    points(age.mea[errcode=="!!"]+shift,RL2.hat[errcode=="!!"],
           type="p",col=unsafe.col,pch=methpch,cex=pchcex,lwd=pchlwd) 

    #  Konfidenz- / Toleranzintervslle
    if (all(is.na(RL2.tilo)))
    { #  Kein Bootstrap-Ergebnis - asymptotische Intervalle nehmen
      #cat("[SplineVerlauf3] (3) Asymptotisches CI für RL2 verwendet\n")
      arrows(age.mea+shift,RL2.cilo,
             age.mea+shift,RL2.cihi,
             code=3,angle=90,length=arrow.length,col=methcol,lwd=arrlwd)

      #  Mark CI in gray, if solution type != 2  
      arrows(age.mea[errcode=="!!"]+shift,RL2.cilo[errcode=="!!"],
             age.mea[errcode=="!!"]+shift,RL2.cihi[errcode=="!!"],
             code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)

    }  else
    { #  Bootstrap-Ergebnis vorhanden - nehmen
      #cat("[SplineVerlauf3] (4) Bootstrap-TI für RL2 verwendet\n")
      arrows(age.mea+shift,RL2.tilo,
             age.mea+shift,RL2.tihi,
             code=3,angle=90,length=arrow.length,col="magenta",lwd=arrlwd)

      #  Mark CI in gray, if solution type != 2  
      arrows(age.mea[errcode=="!!"]+shift,RL2.tilo[errcode=="!!"],
             age.mea[errcode=="!!"]+shift,RL2.tihi[errcode=="!!"],
             code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
    }

    # Legende am Ende, wenn CI-Typ bekannt ist

  } # if (plot.ok)

  # .........................................................................
  #  Spline-Berechnung

  if (n >= 4 )
  { #  Spline ist möglich  

    age.out <- seq(Floor(age.mea[1],1),Ceiling(age.mea[length(age.mea)],1),
                         by=1)

    # RL1.sspline.gcv <- smooth.spline(age.mea, RL1.hat)

    if (( 4 <= n) & (n <= 10)) { df.user <- Round(0.9 * n, 1) }
    if ((11 <= n) & (n <= 12)) { df.user <- Round(0.8 * n, 1) }
    if ((13 <= n) & (n <= 14)) { df.user <- Round(0.6 * n, 1) }
    # if ((15 <= n)            ) { df.user <- RL1.sspline.gcv$df }
    if ((15 <= n)            ) { df.user <- Round(0.5 * n, 1) }

    RL1.sspline.df      <- smooth.spline(age.mea, RL1.hat,df=df.user)
    RL1.sspline.df.pred <- predict(RL1.sspline.df,age.out)

    #  Spline eintragen
    lines(RL1.sspline.df.pred$x+shift,RL1.sspline.df.pred$y,
          col=methcol,lwd=spllwd)

 
    # .........................................................................
    #  RL2 

    RL2.sspline.df      <- smooth.spline(age.mea, RL2.hat,df=df.user)
    RL2.sspline.df.pred <- predict(RL2.sspline.df,age.out)

    #  Spline eintragen
    lines(RL2.sspline.df.pred$x+shift,RL2.sspline.df.pred$y,
          col=methcol,lwd=spllwd)

    return(list(RL.emp=data.frame(age=age.mea,RL1=RL1.hat,RL2=RL2.hat,
                                  RL1.cilo=RL1.cilo,RL2.cilo=RL2.cilo,
                                  RL1.tilo=RL1.tilo,RL2.tilo=RL2.tilo),
                RL.smo=data.frame(age=age.out,
                                  RL1.hut=RL1.sspline.df.pred$y,
                                  RL2.hut=RL2.sspline.df.pred$y),
                plot.ok=plot.ok))

  } else  #  if n>= 4 ..
  { 
    #  Spline nicht möglich, da zu wenig Stützpunkte
    cat(  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
        "\n  No spline smoothing for sex =",sexcode,
        ", because < 4 support points",
        "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

    return(list(RL.emp=data.frame(age=NA, RL1=NA,RL2=NA,
                                  RL1.cilo=NA, RL2.cilo=NA,
                                  RL1.tilo=NA, RL2.tilo=NA),
                RL.smo=data.frame(age=NA,
                                  RL1.hut=NA,
                                  RL2.hut=NA),
                plot.ok=plot.ok))
  }
} 

# ============================================================================
# ****************************************************************************
#  F_TIQQ.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Use a QQ plot for finding a truncation interval, given a data vector x
#  and transformation lambda 
 
#  #  CHANGE HISTORY
#
# 25.01.2021 Error in plotting regression lines removed
# 13.01.2021 Call changed, parameters for RL@.pen, q.fact alphabetically
# 08.01.2021 Call changed
# 01.12.2020 y removed from call
# 18.11.2020 Call changed. Some calculations done by QQW 
# 16.11.2020 data.frames constructed with stringsAsFactors=FALSE
#            KDE calculation done by KDE(). New call of TIQQ!
#            All QQ calculations are done on a possibly reduced data set.
#            For a full QQ plot analysis use @@@   
# 14.11.2020 bins.n.min.act added to call 
# 08.10.2020 x.bins.min removed from call 
# 26.09.2020 Plot devices re-organized, plotting made conditional on data
# 09.09.2020 Strategy for calculation of h in hist changed, break calculation
#            moved to calling programme
# 31.08.2020 Start
# ============================================================================


TIQQ <- function(x.hist, x.kde.mode, 
                 y.red, y.red.hist, y.kde, y.kde.mode, 
                 y.red.hist.mode.idx,
                 lambda, round.unit, kernel.n, bins.n.min.act, 
                 x.tr.bins.min, x.tr.prop.min, x.tr.prop.max,
                 df.est,df.con,
                 x.Q1, x.Q2, RL1.p, RL2.p,
                 l.fact, p.fact, r.fact, w.fact, 
                 xlabel, subtitle,
                 figD, figE, figF, figG, figH,
                 lambda.gen, mue.gen, sig.gen,
                 histcol2, bordercol2, kdecol, qqrawcol, polycol, gencol,
                 print.log.message)
{
  #  INPUT 
  #  x
  #  x.tr.bins.min
  #  x.bins.min
  #  x.tr.prop.min
  #  x.tr.prop.max

  #  figD        in TIQQ: BC transformed collapsed histogram
  #  figE        in TIQQ: QQ plot for BC transformed data, polygon 
  #                       approximation of the QQ plot
  #  figF        in TIQQ: All QQ plot based parameter estimates
  #  figG        in TIQQ: All QQ plot based RL estimates
  #  figH        in TIQQ: QQ plot based relative errors in RL estimates

  #  OUTPUT 
  #  tab
  # =========================================================================

  y.red.n <- length(y.red)

  # --------------------------------------------------------------------------
 
  # BC transformed collapsed histogram                 
  if (!is.na(figD))
  { dev.set(figD)
    #bringToTop(figD)
    yplotmax <- max(y.kde$y, y.red.hist$density)

    plot(y.red.hist,
         border=bordercol2, col=histcol2, freq=FALSE, lty=1, 
         main=paste("Histogram of reduced y=BC(x", lambda, ")",sep=""),
         xlab=paste("y = BC(", xlabel, ", lambda=", lambda,")", sep=""), 
         ylab="pdf(y)",
         sub=subtitle, cex.sub=0.8)

    lines(y.kde, col=kdecol)
    abline(v=y.kde.mode, lty=2)
  }

  # -------------------------------------------------------------------------
  #  QQ plot for the reduced data

  if (!is.na(figE))
  { dev.set(figE)
    #bringToTop(figE)
  }
  
  y.red.qq <- qqnorm(y.red, plot.it=!is.na(figE), col=qqrawcol,
                 main=paste("Normal QQ plot for BC(x.red, ", lambda, ")", sep=""),
                 sub=subtitle, cex.sub=0.7)

  # -----------------------------------------------------------------------------
  #  Select a subset of y.red.qq for calculating qq estimates.

  if (!is.na(figE))
  { dev.set(figE)
    points(y.red.qq$x, y.red.qq$y, col=polycol)

    #  Show the breaks from y.red.hist
    abline(h=y.red.hist$breaks, col=bordercol2, lty=2)  
  }

  #  Calculate r^2 for those subintervals that have enough bins and values
  #  CalcR2InSubInt assumes that the data is transformed to normally 
  #  distributed nonpath values.

  C <- CalcR2InSubInt(x.hist, y.red.hist, y.red.qq, lambda,
                      bins.n.min.act, x.tr.prop.min, x.tr.prop.max, 
                      y.red.hist.mode.idx, 
                      df.est,df.con,
                      x.Q1, x.Q2, RL1.p, RL2.p,
                      l.fact, p.fact, r.fact, w.fact,
                      print.log.message)

  tab.qqw <- C[["Pr2"]]

  # [1] "bin.start" "bin.end"   "bins.n"    "prop"      "x.tr.lo"   "x.tr.hi"  
  # [7] "y.tr.lo"   "y.tr.hi"   "mue"       "sig"       "r2"        "rank.r2"  
  #[13] "y.RL1"     "y.RL2"

  #  Change some names to TMC standard
  tab.qqw.names <- names(tab.qqw)
  tab.qqw.names[which(tab.qqw.names=="mue")]      <- "mue.qqw"
  tab.qqw.names[which(tab.qqw.names=="sig")]      <- "sigma.qqw"
  tab.qqw.names[which(tab.qqw.names=="prop")]     <- "prop.qqw"
  tab.qqw.names[which(tab.qqw.names=="y.RL1")]    <- "y.RL1.qqw"
  tab.qqw.names[which(tab.qqw.names=="y.RL2")]    <- "y.RL2.qqw"

  names(tab.qqw) <- tab.qqw.names

  #  Is there any information in tab.qqw? (May be not if too few bins)
  if (is.na(tab.qqw[1, "bin.start"]))
  { # No information, no further action, return a single NA
    #cat("\n [TIQQ] No result from TIQQ for proportion ", 
    #        x.tr.prop.min, " - ", x.tr.prop.max, "\n")
    #print(tab.qqw)
    tab.qqw <- NA
  } else
  {
    #  At least 1 line with information
    tab.qqw   <- RowToMatrix(tab.qqw)
    tab.qqw.n <- nrow(tab.qqw)

    #  Add the lambda used to tab.qqw
    tab.qqw <- data.frame(lambda.qqw=rep(lambda, times=nrow(tab.qqw)), tab.qqw,
                        stringsAsFactors=FALSE)

    #  Transform breaks to linear scale
    tab.qqw[ , "x.tr.lo"] <- BoxCoxInv(tab.qqw[ ,"y.tr.lo"], lambda)
    tab.qqw[ , "x.tr.hi"] <- BoxCoxInv(tab.qqw[ ,"y.tr.hi"], lambda)

    #  Transform RLs to linear scale
    tab.qqw <- data.frame(tab.qqw, 
                        x.RL1.qqw=BoxCoxInv(tab.qqw[ ,"y.RL1.qqw"], lambda),
                        stringsAsFactors=FALSE)
    tab.qqw <- data.frame(tab.qqw, 
                        x.RL2.qqw=BoxCoxInv(tab.qqw[ ,"y.RL2.qqw"], lambda),
                        stringsAsFactors=FALSE)

    #  Calculate xl.n xc.n, xr.n and corresponding prevalence: moved to QQW(),
    #  because here potentially a reduced dataset is used

    #  Show first 3 QQ plot solutions, if there are that many
    if (!is.na(figE))
    { dev.set(figE)
      if (tab.qqw.n >= 1) { abline(c(tab.qqw[1, "mue.qqw"], 
                                     tab.qqw[1, "sigma.qqw"]), col="blue") }
      if (tab.qqw.n >= 2) { abline(c(tab.qqw[2, "mue.qqw"], 
                                     tab.qqw[2, "sigma.qqw"]), col="orange") }
      if (tab.qqw.n >= 3) { abline(c(tab.qqw[3, "mue.qqw"], 
                                     tab.qqw[3, "sigma.qqw"]), col="red") }
      legend("topleft", 
           c("QQ plot solution 1", "QQ plot solution 3", "QQ plot solution 2"),
           col=c("red", "orange", "yellow"), lwd=3, lty=1, cex=0.8) 
    }

    #  If there is information on the distribution parameters from other 
    #  sources (data generation, other methods): calculate the error in the 
    #  false positive rate (in percentage points)
    #  lambda.gen, mue.gen, sig.gen: external information!
    if (!is.na(lambda.gen))
    {
      x.RL1.FP.err.qqw <- cdf.PN(tab.qqw[ , "x.RL1.qqw"], lambda.gen, mue.gen, 
                               sig.gen) 
      x.RL1.FP.err.qqw <- 100*(x.RL1.FP.err.qqw - RL1.p)/RL1.p
      x.RL2.FP.err.qqw <- cdf.PN(tab.qqw[ , "x.RL2.qqw"], lambda.gen, mue.gen,
                               sig.gen)
      x.RL2.FP.err.qqw <- 100*(x.RL2.FP.err.qqw - RL2.p)/RL2.p
      x.RL.FP.err.qqw  <- abs(x.RL1.FP.err.qqw) + abs(x.RL2.FP.err.qqw)

      tab.qqw <- data.frame(tab.qqw, x.RL1.FP.err.qqw=x.RL1.FP.err.qqw,
                          stringsAsFactors=FALSE)
      tab.qqw <- data.frame(tab.qqw, x.RL2.FP.err.qqw=x.RL2.FP.err.qqw,
                          stringsAsFactors=FALSE)
      tab.qqw <- data.frame(tab.qqw, x.RL.FP.err.qqw=x.RL.FP.err.qqw,
                          stringsAsFactors=FALSE)
    }

    #  Show QQ plot based parameter estimates
    if (!is.na(figF))
    { dev.set(figF)
      #bringToTop(figF)

      plot(tab.qqw[ , "mue.qqw"], tab.qqw[ , "sigma.qqw"], type="p", col=qqcol,
         xlab="mue", ylab="sig", 
         main=paste("QQ plot parameter estimates, lambda=",lambda, sep=""),
         sub=subtitle, cex.sub=0.8)
    }

    #  Show QQ plot based RL estimates
    if (!is.na(figG))
    { dev.set(figG)
      #bringToTop(figG)

      plot(tab.qqw[ , "x.RL1.qqw"], tab.qqw[ , "x.RL2.qqw"], type="p", col=qqcol,
        xlab="x.RL1.qqw", ylab="x.RL2.qqw", 
        main=paste("QQ plot RL estimates, lambda=",lambda, sep=""),
        sub=subtitle, cex.sub=0.8)

      if (!is.na(lambda.gen))
      {
        points(x.RL1.gen, x.RL2.gen, pch=16, col=gencol, cex=1.5)
      }
    }

    #  Show QQ plot based relative errors in RL estimates
    if (!is.na(figH) & !is.na(lambda.gen))
    { dev.set(figH)
      #bringToTop(figH)

      plot(tab.qqw[ , "x.RL1.FP.err.qqw"], tab.qqw[ , "x.RL2.FP.err.qqw"], 
         type="p", col=qqcol,
         xlab="x.RL1 FP error in %", ylab="x.RL2 FP error in %", 
         main=paste("FP errors from QQ plot estimates (qqw), lambda=",lambda,
                  sep=""),
       sub=subtitle, cex.sub=0.8)
    }

  }

  return(tab.qqw)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


# ****************************************************************************
#  F_tmc.master.R
#
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Master function for TMC estimation. Function version, replaces source
#  version TMC_seg102_tmc_master5.R. This function does the TMC estimation 
#  for a single dataset. QQW results are also produced (and reported) as they
#  serve as initial values for TMC.
#  
#  Input is the data vector plus parameters controlling the estimation.
#  Action here:
#
#  - find initial values via a QQ plot approach in a loop over lambda values
#  - update the QQ plot estimate via TMC, thereby estimating a possibly new
#    lambda, but keeping the truncation interval fixed  
#  - Do bias reduction by tmu, switched on/ off by do.tmu
#  - Reported results: from QQ plot (tab.qqw) and TMC (tab.tmc). Note that
#    qqw here does not mean the same as in earlier versions. 
#  
#  CHANGE HISTORY
#
#  12.05.2021 Use of fixed parameter components installed
#  13.01.2012 Call changed, parameters for RL@.pen added, q.fact now 
#             alphabetically
#  06.12.2020 Call changed, x.kde.mode removed, x.lt.mode, x.ge.mode changed
#  04.12.2020 Changed components for CalcPrev() call
#  04.12.2020 Components for CalcPrev() added to call 
#  19.11.2020 x removed from tmc1 call
#  14.11.2020 xsupp added to call
#  09.10.2020 do.tmu removed
#             Only tab.tmc is returned
#  08.12.2020 RB3 removed from function call
#  26.09.2020 Plot devices re-organized, colours re-organized
#  09.09.2020 QQ plot estimation + tmc + tmu 
#  31.08.2020 Developed from TMC_seg102_tmc_master5.R

# ============================================================================

tmc.master <- function(x, x.hist, tab.qqw, 
                       theta.fix, theta.est, idx.fix, idx.est,
                       lambda.min, lambda.max, 
                       df.est, df.con, p.fit.min,
                       l.fact, p.fact, r.fact, w.fact,
                       x.Q1, x.Q2, RL1.p, RL2.p, fastnull, fastnull.chi,                     
                       print.log.message,
                       x.RL1.npa, x.RL2.npa, 
                       xsupp, gencol, kdecol,tmccol, figJ)
{  
  #  INPUT

  #  tab.qqw     matrix containing the QQW results, initial values for TMC
  #  ...


  #  OUTPUT
  #  tab.tmc  Vector containing the TMC result. For components see
  #           tab.tmc.names
 

  # ==========================================================================

  # --------------------------------------------------------------------------
  if (print.log.message) { cat("%%%   tmc.master   Start\n") }
  # --------------------------------------------------------------------------

  # ==========================================================================
  #  Dimensions
  
  x.n   <- length(x)  

  # ==========================================================================
  #  Table for TMC results

  #  Indices ilo, ihi refer to components of counts 
  tab.tmc.names <- c("ilo", "ihi", "x.tr.lo", "x.tr.hi", "x.tr.n", "x.tr.prop",
                     "rc", "iter", 
                     "lambda.tmc", "mue.tmc", "sigma.tmc", 
                     "x.RL1.tmc", "x.RL2.tmc", "x.RL1.cilo", "x.RL1.cihi", 
                      "x.RL2.cilo", "x.RL2.cihi", "reldist.tmc", "p.fit",
                     "opt.crit", "prev.l.tmc", "prev.c.tmc", "prev.r.tmc",
                     "prev.l.tmc.pen", "prev.r.tmc.pen", "xc.n.tmc",
                     "chi2.total", "chi2.total.df", "chi2.trun",
                     "chi2.trun.df", "chi2.path")
  tab.tmc.names.n <- length(tab.tmc.names)

  #  Matrix containing candidates for the final TMC results 
  tab.tmc.can           <- matrix(NA, nrow=1, ncol=length(tab.tmc.names))
  colnames(tab.tmc.can) <- tab.tmc.names

  # ==========================================================================
  #  Error code
  error.tmc.master <- FALSE

  # =========================================================================
  #  Loop over QQW  results

  tab.qqw.n <- nrow(tab.qqw)

  for (iqq in 1:tab.qqw.n)
  {                
    # -----------------------------------------------------------------------
    #  Set constraints for TMC estimation. Take initial values from previous 
    #  QQ analysis. Transform all components, including those that are fixed. 
    theta.ini <- c(lambda=tab.qqw[iqq, "lambda.qqw"],
                   mue=tab.qqw[iqq, "mue.qqw"],
                   sigma=tab.qqw[iqq, "sigma.qqw"])

    y.ini     <- BoxCox(x, theta.ini[1])

    RB3           <- matrix(NA,nrow=3,ncol=2)
    dimnames(RB3) <- list(c("lambda","mue","sigma"),c("min","max"))

    RB3[1, ]  <- c( lambda.min-fastnull,lambda.max+fastnull)
    if (y.ini[1] <  0) { RB3[2, 1] <- 2*y.ini[1] } 
    if (y.ini[1] == 0) { RB3[2, 1] <- 2*fastnull }   #  25.05.2020
    if (y.ini[1] >  0) { RB3[2, 1] <- 0.5*y.ini[1] } 
    RB3[2, 2]  <- 2.0*tail(y.ini,n=1)                # (y is sorted ascending)

    RB3[3, ]  <- c( 0.001*theta.ini[3], 2*theta.ini[3])

    #  Transform initial parameter to real scale
    ttheta.ini <- o2R(theta.ini, RB3, fastnull=fastnull)

    # -----------------------------------------------------------------------
    #  Do TMC analysis. Returns vector TMC. Results go preliminary
    #  in tab.tmc.can.can. The best of these candidates goes into
    #  tab.tmc.rep.

    ilo <- tab.qqw[iqq, "bin.start"] 
    ihi <- tab.qqw[iqq, "bin.end"] 

    x.tr.lo <- x.hist$breaks[ilo]
    x.tr.hi <- x.hist$breaks[ihi+1]
    x.tr.n  <- sum(x.hist$counts[ilo:ihi])   #  cannot be taken from qqw!
    x.lt.tr.n <- sum(x <  x.tr.lo)
    x.ge.tr.n <- sum(x >= x.tr.hi)

    tab.tmc.can.can   <- tmc1(ttheta.ini, RB3, idx.fix, idx.est, 
                              x.hist, x.tr.n, x.tr.lo, x.tr.hi,
                              x.lt.tr.n, x.ge.tr.n,
                              l.fact, p.fact, r.fact, w.fact, 
                              x.Q1, x.Q2, RL1.p, RL2.p,
                              df.est, df.con, fastnull, fastnull.chi,
                              ilo, ihi, x.RL1.npa, x.RL2.npa, 
                              tab.tmc.names)

    if (iqq == 1)
    { tab.tmc.can <- tab.tmc.can.can} else
    { tab.tmc.can <- rbind(tab.tmc.can, tab.tmc.can.can) }

  }  #  iqq ...

  tab.tmc.can <- RowToMatrix(tab.tmc.can)
  row.names(tab.tmc.can) <- NULL

  # =========================================================================

  #cat("\n [tmc.master] All results from TMC, not ordered\n")
  #print(tab.tmc.can)
  
  if (nrow(tab.tmc.can) > 1)
  {
    tab.tmc.can <- tab.tmc.can[order(-tab.tmc.can[ ,"p.fit"]), ]
  }

  #cat("\n [tmc.master] All results from TMC, sorted by p.fit\n")
  #print(tab.tmc.can)

  #  Show TMC estimates for each QQ plot based solution
  if (!is.na(figJ))
  { dev.set(figJ)
    #bringToTop(figJ)

    plot(x.hist,
         border=bordercol2, col=histcol2, freq=FALSE, lty=1,
         main="Collapsed histogram and TMC candidates fits",
         sub=subtitle, cex.sub=0.7)

    for (i in 1:tab.qqw.n) 
    {
      x.pdf.tmc <-  tab.tmc.can[i, "prev.c.tmc"] *
                    pdf.PN(xsupp, tab.tmc.can[i, "lambda.tmc"],  
                          tab.tmc.can[i, "mue.tmc"], tab.tmc.can[i, "sigma.tmc"])
      lines(xsupp, x.pdf.tmc, col=tmccol, lty=2)
      abline(v=c(tab.tmc.can[i, "x.RL1.tmc"], tab.tmc.can[i, "x.RL2.tmc"]),
             col=tmccol, lty=2)      
      #locator(1)
    }
  }

  # -------------------------------------------------------------------------
  #  Select optimal TMC solution from these candidates

  tab.tmc <- tab.tmc.can[1, ]

  # -------------------------------------------------------------------------
  #  Add asymptotic confidence intervals
  xc.RL1.tmc.CI <- CIQuant.PNV(tab.tmc["x.tr.n"], RL1.p, 
                               tab.tmc["x.RL1.tmc"],
                               tab.tmc["lambda.tmc"], 
                               tab.tmc["mue.tmc"], 
                               tab.tmc["sigma.tmc"],
                               alpha=alpha)
  tab.tmc["x.RL1.cilo"] <- xc.RL1.tmc.CI[1]
  tab.tmc["x.RL1.cihi"] <- xc.RL1.tmc.CI[2]

  xc.RL2.tmc.CI <- CIQuant.PNV(tab.tmc["x.tr.n"], RL1.p, 
                               tab.tmc["x.RL2.tmc"],
                               tab.tmc["lambda.tmc"], 
                               tab.tmc["mue.tmc"], 
                               tab.tmc["sigma.tmc"],
                               alpha=alpha)
  tab.tmc["x.RL2.cilo"] <- xc.RL2.tmc.CI[1]
  tab.tmc["x.RL2.cihi"] <- xc.RL2.tmc.CI[2]

  # -------------------------------------------------------------------------

  #  Plot result
  x.pdf.tmc <-  tab.tmc["prev.c.tmc"] *
                pdf.PN(xsupp, tab.tmc["lambda.tmc"],  
                       tab.tmc["mue.tmc"], tab.tmc["sigma.tmc"])

  if (!is.na(figJ))
  { dev.set(figJ)
    lines(xsupp, x.pdf.tmc, col=tmccol, lty=1, lwd=2)
    abline(v=c(tab.tmc["x.RL1.tmc"], tab.tmc["x.RL2.tmc"]),
           col=tmccol, lty=1, lwd=2)
  }      
  
  # --------------------------------------------------------------------------
  if (print.log.message) { cat("%%%   tmc.master   End\n") }
  # --------------------------------------------------------------------------

  return(tab.tmc=tab.tmc)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_tmc0_V3.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Function executing tmc estimation (the inner loop of the complete etimation)

#  CHANGE HISTORY
#  Changes in functions are documented in the functions!

#  13.01.2012 Call changed, parameters for RL@.pen added, q.fact now 
#             alphabetically
#  06.12.2020 Call changed, x.kde.mode removed, x.lt.mode, x.ge.mode changed 
#  04.12.2020 Components for CalcPrev() added to call
#  19.11.1010 x removed from call
#  09.10.2020 do.tmu removed
#  28.09.2020 ok range changed to { , )
#  08.09.2020 tmu installed, execution via do.tmu.  
#  30.03.2020 eps removed from parameter list
#  17.09.2019 p.fact as additional parameter
#  15.03.2019 Installation from directory ProgWW17_RH

# ============================================================================

tmc0 <- function(ttheta.ini, RB3, idx.fix, idx.est,
                 x.hist, x.tr.n, x.tr.lo, x.tr.hi,
                 x.lt.tr.n, x.ge.tr.n,
                 l.fact, p.fact, r.fact, w.fact, 
                 x.Q1, x.Q2, RL1.p, RL2.p,
                 df.est, df.con, 
                 fastnull, fastnull.chi)
{ 
  # ttheta.ini    Initial estimates, transformed to whole real line
  # RB3           Constraints for parameter estimation
  # x             Data, complete, original scale
  # x.hist        Histogram of x, original scale
  # x.tr.lo       Lower limit of truncation interval, original scale
  # x.tr.hi       Upper limit of truncation interval, original scale
  # l.fact        Factor "weight for truncation interval size"
  # w.fact        Factor "too high prediction" in penalty term 
  # p.fact        Factor "negative prevalence" in penalty term 
  # fastnull      Absicherung gegen kleine Rundungsfehler
  # fastnull.chi  Absicherung gegen kleine Erwartungswerte beim chi^2-Test

  # --------------------------------------------------------------------------

  # ==========================================================================
  #if (print.log.message) { cat("%%%   F_tmc0_V3.R  Start\n") }
  # ==========================================================================

  # ttheta contains all 3 components of theta (lambda, mue, sigma),
  # but possibly some of them are fixed by the user. Pass only components
  # to nlm that are to be estimated. 

  ttheta.ini.est <- Extract(ttheta.ini, idx.est)

  if (is.null(idx.fix[1]))
  { # No fixed component
    ttheta.ini.fix <- NULL
  }  else
  { # >= 1 components in idx.fix
    ttheta.ini.fix <- ttheta.ini[idx.fix]
  }

  tmc <- nlm(chi2.PNV.tr.lms, ttheta.ini.est, iterlim=500, print.level=0,
             gradtol=1.e-8, steptol=1.e-8, 
             RB=RB3, 
             ttheta.ini.fix=ttheta.ini.fix, idx.fix=idx.fix, idx.est=idx.est,
             x.hist=x.hist, x.tr.lo=x.tr.lo, x.tr.hi=x.tr.hi,
             x.lt.tr.n=x.lt.tr.n, 
             x.ge.tr.n=x.ge.tr.n,
             l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, w.fact=w.fact,
             x.Q1=x.Q1, x.Q2=x.Q2, RL1.p=RL1.p, RL2.p=RL2.p,
             df.est=df.est, df.con=df.con, 
             opt.crit.only=TRUE, fastnull=fastnull.chi)

  #  Transform estimated and fixed parameters to original scale
  #  First build complete parameter vector
  ttheta.hut <- Compose(ttheta.ini.fix, tmc$estimate, idx.fix, idx.est)

  theta.hut <- R2o(ttheta.hut,RB3)
   
  result <- c(theta.hut[1], theta.hut[2], theta.hut[3], 
              tmc$minimum, tmc$code, tmc$iterations)
  names(result) <- c("lambda.tmc", "mue.tmc", "sigma.tmc",
                     "chi2.total.min", "rc", "iter")
  #  Meaning of "code" ("rc"):
  #  1: relative gradient is close to zero, current iterate is 
  #     probably solution.
  #  2: successive iterates within tolerance, current iterate is 
  #     probably solution.
  #  3: last global step failed to locate a point lower than estimate. 
  #     Either estimate is an approximate local minimum of the function 
  #     or steptol is too small.
  #  4: iteration limit exceeded.
  #  5: maximum step size stepmax exceeded five consecutive times. 
  #     Either the function is unbounded below, becomes asymptotic 
  #     to a finite value from above in some direction or stepmax is too small.

  # ==========================================================================
  #if (print.log.message) { cat("%%%   F_tmc0_V3.R  End\n") }
  # ==========================================================================

  return(result)
}   
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_tmc1.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Function for TMC parameter estimation plus the calculation of estimate
#  properties

#  CHANGE HISTORY

#  13.01.2012 Call changed, parameters for RL@.pen added, q.fact now 
#             alphabetically
#  06.12.2020 Call changed, x.kde.mode removed, x.lt.mode, x.ge.mode changed 
#  04.12.2020 Components for CalcPrev() added to call
#  19.11.2020 x removed from tmc1 call
#             x removed from tmc0( ...
#             length(x) replaced by sum(x.hist$counts)
#  09.10.2020 do.tmu removed
#  29.03.2020 chi2.per.df replaced by opt.crit
#  12.03.2020 RL1.p, RL2.p added to parameter list and quantile calculation
#  18.02.2020 Start

# ============================================================================


tmc1 <- function(ttheta.ini, RB3, idx.fix, idx.est,
                 x.hist, x.tr.n, x.tr.lo, x.tr.hi, 
                 x.lt.tr.n, x.ge.tr.n,
                 l.fact, p.fact, r.fact, w.fact, 
                 x.Q1, x.Q2, RL1.p, RL2.p,
                 df.est, df.con, fastnull, fastnull.chi,
                 ilo, ihi, x.RL1.npa, x.RL2.npa, result.names)
{ 
  # ==========================================================================
  #if (print.log.message) { cat("%%%   tmc1   Start\n") }
  # ==========================================================================

  #  Prepare result vector
  result <- rep(NA, times=length(result.names))
  names(result) <- result.names

  x.n <- sum(x.hist$counts)

  #  Do parameter estimation

  TMC <- tmc0(ttheta.ini, RB3, idx.fix, idx.est,
              x.hist, x.tr.n, x.tr.lo, x.tr.hi,
              x.lt.tr.n, x.ge.tr.n,
              l.fact, p.fact, r.fact, w.fact, 
              x.Q1, x.Q2, RL1.p, RL2.p,
              df.est, df.con, fastnull, fastnull.chi)

  #  RLs and other derived quantities
  x.RL1.tmc <- q.PN(RL1.p, TMC["lambda.tmc"],  TMC["mue.tmc"],  TMC["sigma.tmc"])
  x.RL2.tmc <- q.PN(RL2.p, TMC["lambda.tmc"],  TMC["mue.tmc"],  TMC["sigma.tmc"])

  #  Relative distance between (P025.c, P975.c) and (x.RL1.tmc, RL2.tmc)
  #  This is available only for test data!

  if (is.na(x.RL1.npa))
  { # No test data
    reldist.tmc <- NA
  } else
  { #  Test data  
    reldist.tmc <- reldist(c(x.RL1.tmc, x.RL2.tmc), c(x.RL1.npa, x.RL2.npa))
  }

  #  Fill result vector
  result["ilo"]         <- ilo
  result["ihi"]         <- ihi
  result["x.tr.lo"]     <- x.tr.lo
  result["x.tr.hi"]     <- x.tr.hi
  result["x.tr.n"]      <- x.tr.n
  result["x.tr.prop"]   <- x.tr.n / x.n
  result["rc"]          <- TMC["rc"]
  result["iter"]        <- TMC["iter"]
  result["lambda.tmc"]  <- TMC["lambda.tmc"]
  result["mue.tmc"]     <- TMC["mue.tmc"]
  result["sigma.tmc"]   <- TMC["sigma.tmc"]
  result["x.RL1.tmc"]   <- x.RL1.tmc
  result["x.RL2.tmc"]   <- x.RL2.tmc
  result["reldist.tmc"] <- reldist.tmc

  # Calculate properties of the actual estimate
  TMC.prop <-  chi2trunc(x.hist,x.tr.lo,x.tr.hi,
                         x.lt.tr.n, x.ge.tr.n, 
                         x.Q1, x.Q2, RL1.p, RL2.p,
                         TMC["lambda.tmc"],TMC["mue.tmc"],TMC["sigma.tmc"],
                         df.est,df.con,
                         l.fact, p.fact, r.fact, w.fact, 
                         opt.crit.only=FALSE, fastnull=fastnull)

  result["p.fit"]          <- TMC.prop$res["chi2.trun.p"]
  result["opt.crit"]       <- TMC.prop$res["opt.crit"]
  result["prev.l.tmc"]     <- TMC.prop$res["prev.l.tmc"]
  result["prev.c.tmc"]     <- TMC.prop$res["prev.c.tmc"]
  result["prev.r.tmc"]     <- TMC.prop$res["prev.r.tmc"]
  result["prev.l.tmc.pen"] <- TMC.prop$res["prev.l.tmc.pen"]
  result["prev.r.tmc.pen"] <- TMC.prop$res["prev.r.tmc.pen"]
  result["xc.n.tmc"]       <- TMC.prop$res["xc.n.tmc"]

  result["chi2.total"]     <- TMC.prop$res["chi2.total"]
  result["chi2.total.df"]  <- TMC.prop$res["chi2.total.df"]
  result["chi2.trun"]      <- TMC.prop$res["chi2.trun"]
  result["chi2.trun.df"]   <- TMC.prop$res["chi2.trun.df"]
  result["chi2.path"]      <- TMC.prop$res["chi2.path"]

  # ==========================================================================
  # if (print.log.message) { cat("%%%   tmc1   End\n") }
  # ==========================================================================

  return(result)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_TransformKDE.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Transform the kde of x to the kde of y <- f(x)
 
#  #  CHANGE HISTORY
#
# 18.11.2020 Start
# ============================================================================

TransformKDE <- function(x.kde, y)
{
  #  INPUT 
  #  x.kde          kernel density estimate of x, has components $x, $y
  #  y              transformed x.kde$x

  #  OUTPUT 
  #  y.kde          kernel density estimate of y, has components $x, $y
  # =========================================================================

  kernel.n <- length(x.kde$x)

  #  Calculate area of each trapezoid in x.kde
  x.area <- diff(x.kde$x) * (x.kde$y[1:(kernel.n-1)] + x.kde$y[2:kernel.n])/2
  # cat("\n sum(x.area): ", sum(x.area), "\n")

  y.diff <- diff(y)
  y.kde.y <- rep(NA, times=kernel.n)
  y.kde.y[1] <- x.kde$y[1]     #  should be zero in both cases

  for (i in 2:kernel.n)
  { 
    y.kde.y[i] <- 2*x.area[i-1]/y.diff[i-1] -  y.kde.y[i-1]
  }
  y.area <- y.diff * (y.kde.y[1:(kernel.n-1)] + y.kde.y[2:kernel.n])/2
  # cat("\n sum(y.area): ", sum(y.area), "\n")
  
  return(y.kde=list(x=y, y=y.kde.y))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  Test
#x.kde <- list(x=seq(6, 14, by=0.1), y=dnorm(seq(6, 14, by=0.1), mean=10, sd=1))
#kern.n <- length(x.kde$x)

#y <- log(x.kde$x)
##y <- x.kde$x

#y.kde <- TransformKDE(x.kde, y) 

#cbind(x.kde$x, y.kde$x)
#cbind(x.kde$y, y.kde$y)

#z <- exp(y.kde$x)
#z.kde <- TransformKDE(y.kde, z)

#cbind(x.kde$y, z.kde$y)
#  Hin- und Rücktransformation ist ok
# ............................................................................

#x <- rnorm(1000, mean=10, sd=2)
#x.kde <- density(x)

#y <- log(x)
#y.kde1 <- density(y)

#y.kde2 <- TransformKDE(x.kde, log(x.kde$x))

#plot(y.kde1, type="l", col="red")
#lines(y.kde2$x, y.kde2$y, type="l", col="blue")

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_TukeyWW.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  'Tukey' method of removing 'outliers'
#  Find a subset of the input data which has data only in the range
#  Q1 - 1.5*IQR, Q3 + 1.5*IQR, calculate mean and SD, adjust SD estimate
#  because of using a truncated sample

#  Source for 'Tukey' method: Horn, Pesce: Reference intervals: an update.
#  Clinica Chimica Acta 334(2003) 5-23, p.13
#  Expansion by WW: iterated solution
#  Source for truncation adjustment: WW

#  'Fences' in Tukey terminology correspond to truncation interval limits  

#  NEEDS
#  Quantile()   function to calculate quantiles
 
#  #  CHANGE HISTORY
#  22.05.2021 Iteration stop criterion x.tr.prop.min added
#  16.05.2021 Iteration made optional, default 1 step
#  25.04.2021 Iteration added
#  22.04.2021 Start  
# ============================================================================

TukeyWW <- function(x, xlabel, IQR.fact=1.5, iter.max=1, x.tr.prop.min=0.50, 
                    RL1.p=0.025, RL2.p=0.975, 
                    print.details=FALSE, figA=NA)
{
  #  INPUT 
  #  x              data vector. NAs are not checked, not allowed
  #  IQR.fact       IQR range factor 
  #  iter.max       controls iteration: 1 : basic (original) Tukey approach
  #                                    >1 : extended version
  #  x.tr.prop.min  iteration emergency break, prevents generation of too small
  #                 Tukey data sets  
  #  RL1.p, RL2.p   probabiilities for RLs       
  #  print.details  print details of the calculation and iteration process

  #  OUTPUT 
  #  mue.hat        estimate for mean
  #  sigma.hat      estimate for standard deviation 
  #  --> see more terms in the return statement
  # ---------------------------------------------------------------------------
  
  #  For a standard normal distribution, the fences below correspond to a 
  #  probability of ...
  Q1.z    <- qnorm(0.25)
  Q3.z    <- -Q1.z
  IQR.z   <- Q3.z - Q1.z
  z.lo    <- Q1.z - IQR.fact*IQR.z     #  lower fence, theoretical
  z.hi    <- Q3.z + IQR.fact*IQR.z     #  upper fence, theoretical
  p.trunc <- pnorm(z.hi) - pnorm(z.lo) #  probability between fences
  range.z <- z.hi - z.lo

  #  Calculate fences for input data
  x.n <- length(x)
  ok  <- rep(TRUE, times=x.n)

  #  Make the loop start at all
  out.lo    <- 1
  out.hi    <- 1
  x.tr.prop <- 1
  iter      <- 0

  #  Iterate until no more cases lie outside the fences or the data 
  #  set is reduced to less than x.tr.prop.min

  #  Problem: with rounded data, the decrease in IQR can be smaller than
  #  the rounding unit, causing that no cases are removed and IQR does no more
  #  change. 
  while ((out.lo > 0 | out.hi > 0) & (iter < iter.max) & 
          (x.tr.prop > x.tr.prop.min))
  #while (                            (iter < iter.max) & 
  #        (x.tr.prop > x.tr.prop.min))
  {
    iter <- iter + 1
    Q1 <- Quantile(x[ok], probs=0.25)
    Q3 <- Quantile(x[ok], probs=0.75)
    IQR <- Q3 - Q1

    #  The fences x.lo and x.hi are part of the accepted interval
    #  Note that these limits may coincide with  existing values if the data
    #  is rounded.
    x.lo <- Q1 - IQR.fact*IQR
    x.hi <- Q3 + IQR.fact*IQR

    #  Shift the theoretical limits to really existing values
    # x.lo <- min(x[x >= x.lo])
    # x.hi <- max(x[x <= x.hi])
    #  No good idea. Changes the proportion of data between fences (p.trunc)
    #  to an unknown level

    out.lo <- sum(x[ok] < x.lo)
    out.hi <- sum(x[ok] > x.hi)
    ok <- (x.lo <= x) & (x <= x.hi)
    x.tr.prop <- sum(ok) / x.n

    if (print.details)
    {
      cat("\n Iteration ", iter, " x.lo", x.lo, " x.hi", x.hi,
                                 " out.lo", out.lo, " out.hi", out.hi, 
                                 " x.tr.prop", x.tr.prop)
    }
  }
  cat("\n")

  #  No more values outside the fences
  mean.hat  <- mean(x[ok])
  sigma.raw <- sd(x[ok])

  # Adjusted estimate of sigma: leave the mean as it is and look for the sd
  # that assigns a probability p.trunc to (x.lo, x.hi)
  
  range.x   <- x.hi - x.lo
  sigma.adj <- range.x/range.z  

  # Calculate xc.n. Needs total numbers of removed.
  out.lo <- sum(x < x.lo)
  out.hi <- sum(x > x.hi)
  xc.n   <- (x.n - out.lo - out.hi)/p.trunc

  #  Distribution based  estimate for RLx
  x.RL1 <- qnorm(RL1.p, mean=mean.hat, sd=sigma.adj)
  x.RL2 <- qnorm(RL2.p, mean=mean.hat, sd=sigma.adj)

  if (print.details)
  {
    cat("\n Total # of values        :", x.n,
        "\n Range of input data      :", min(x), max(x),  
        "\n Range of truncated data  :", min(x[ok]), max(x[ok]),  
        "\n Tukey fences             :", x.lo, x.hi,  
        "\n # of values removed below:", out.lo, 
        "\n # of values removed above:", out.hi,
        "\n Tukey estimate for mean  :", mean.hat,
        "\n Raw Tukey estimate for sd:", sigma.raw,
        "\n Q1.z, Q3.z               :", Q1.z, Q3.z,
        "\n p.trunc                  :", p.trunc,
        "\n range.z                  :", range.z,
        "\n range.x                  :", range.x,
        "\n Adjusted estimate for sd :", sigma.adj,
        "\n Estimated xc.n           :", xc.n,
        "\n Allowed # of iterations  :", iter.max,
        "\n")
  }

  return(c(mean=mean.hat, 
           sd.raw=sigma.raw, sd.adj=unname(sigma.adj), 
           x.tr.lo=unname(x.lo), x.tr.hi=unname(x.hi), 
           z.tr.lo=unname(z.lo), z.tr.hi=unname(z.hi), 
           out.lo=out.lo, out.hi=out.hi,            
           x.RL1=x.RL1, x.RL2=x.RL2, 
           xc.n=xc.n)) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test: see ../dev/Test_Tukey.R 


# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#  F_Untie.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Replace tied data (due to rounding) by replacing identical values by an 
#  equidistant sequence in the range of the rounding unit. Suitable for
#  a common rounding over the whole data set.
 
#  #  CHANGE HISTORY
#  27.08.2020 Start
#
# ==============================================================================

Untie <- function(x, round.unit)
{
  #  INPUT 
  #  x          data vector, possibly with ties
  #  round.unit unit by which the data was rounded

  #  OUTPUT 
  #  x.untie    data vector with ties replaced
  # ---------------------------------------------------------------------------

  #  Replace ties by an increasing sequence
  x.n     <- length(x)
  x.table <- table(x)
  x.val   <- as.numeric(names(x.table)) 
  x.val.n <- length(x.val)
  x.untie <- rep(NA, times=x.n)

  #  i  Position in x.val and x.table
  #  j  Position in x.untie
 
  #  Go to through all unique values and resolve ties
  #  Special case: x = 0, because replacing values must all be >= 0

  i <- 1
  if (x.val[i] == 0)
  { j <- x.table[i]
    if (x.table[i] == 1)
    { x.untie[i] <- x.val[i] } else
    { delta <- round.unit/(2 * x.table[i])
      x.untie[1:j] <- delta/2 + (0:(x.table[i]-1))*delta 
    }
  } else
  { # x.val[i] > 0
    j <- x.table[i]
    if (x.table[i] == 1)
    { x.untie[i] <- x.val[i] } else
    { delta <- round.unit/x.table[i]
      x.untie[1:j] <- x.val[i] - round.unit/2 + delta/2 + (0:(x.table[i]-1))*delta 
    }
  }

  for (i in 2:(x.val.n)) 
  { # x.val[i] > 0
    if (x.table[i] == 1)
    { j1 <- j + 1
      j2 <- j1
      x.untie[j1] <- x.val[i] } else
    { delta <- round.unit/x.table[i]
      j1 <- j + 1
      j2 <- j1 + x.table[i] - 1
      x.untie[j1:j2] <- x.val[i] - round.unit/2 + delta/2 + (0:(x.table[i]-1))*delta 
    }
    j <- j2
  }
  return(x.untie) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test

#round.unit <- 1
#x <- c(0,0,1,1,2,2,2,3,3,3,3,5,5,5,5,5)
#x <- c(1,1,2,2,2,3,3,3,3,5,5,5,5,5)
#x <- rep(0:10, each=3)

#x.untie <- Untie(x, round.unit)

#print(cbind(x, Round(x.untie, round.unit), x.untie))
#mean(x)
#mean(x.untie)
# *****************************************************************************
# F_x.mean.PN.R

# ----------------------------------------------------------------------       

x.mean.PN <- function(lambda,mue,sigma,fastnull=1.e-10)
{
  #  Mean of the power normal distribution, numerically
  #  Parameters of the PND: lambda, mue, sigma
  #  siehe Freeman & Modarres, SPL, 2006, p. 766
  #  Dichte ist für x>0 positiv, Null für x = 0, sonst Null

  #  #  CHANGE HISTORY
  #  03.03.2021 q.X replaced by q.PN
  #             q.PN replaced by 0, Inf
  #  31.07.     Integration limits fixed from quantiles
  #  27.07.2019 Installation
  # ===========================================================================

  #  Calculate integral numerically

  #  Using 0, Inf as integration limits fails sometimes
  #  Use quantiles instead
  #x.lo <- q.PN(0.0001,lambda,mue,sigma)
  #x.hi <- q.PN(1-0.0001,lambda,mue,sigma)

  EX <- integrate(x.pdf.PN, 0, Inf, lambda=lambda, mue=mue, sigma=sigma)

  return(EX[[1]]) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
#  F_x.mode.PN.R

# --------------------------------------------------------------------------
x.mode.PN <- function(lambda,mue,sigma,fastnull=1.e-10)
{
  #  Mode of x with x ~ PND
  #  Parameters of the PND: lambda, mue, sigma
  #  siehe Freeman & Modarres, SPL, 2006, p. 766
  #  Numerical calculation (no case distinction necessary

  #  #  CHANGE HISTORY
  #  27.07.2019 Installation
  # ===========================================================================

  #  Determine integration interval

  x.min <- q.X(  0.001, lambda, mue, sigma)
  x.max <- q.X(1-0.001, lambda, mue, sigma)

  x <- seq(x.min, x.max, length.out=1000)
  f <- pdf.PN(x, lambda, mue, sigma)
  x.mode <- x[which.max(f)]
  if (length(x.mode) > 1)
  { 
    cat("\n x has multiple modes - first is used \n")
    print(x.mode)
    x.mode <- x.mode[1] 
  } 
  return(x.mode) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_x.pdf.PN.R

# -----------------------------------------------------------------------------

x.pdf.PN <- function(x, lambda,mue,sigma,fastnull=1.e-10)
{
  #  Auxiliary function for mean.PN
  #  Parameters of the PND: lambda, mue, sigma
  #  siehe Freeman & Modarra, SPL, 2006, p. 766

  #  #  CHANGE HISTORY
  #  27.07.2019 Installation
  # ===========================================================================

  return(x * pdf.PN(x,lambda,mue,sigma))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
#  F_x.Var.PN.R

# ----------------------------------------------------------------------       

x.Var.PN <- function(EX, lambda,mue,sigma,fastnull=1.e-10)
{
  #  Variance of the power normal distribution, numerically
  #  Parameters of the PND: lambda, mue, sigma
  #  siehe Freeman & Modarra, SPL, 2006, p. 766
  #  Dichte ist für x>0 positiv, Null für x = 0, sonst Null

  #  #  CHANGE HISTORY
  #  03.03.2021 q.X replaced by q.PN
  #             Quantiles replaced by 0, Inf
  #             EX required as argument 
  #  31.07.     Integration limits fixed from quantiles
  #  27.07.2019 Installation
  # ===========================================================================

  #  Calculate integral numerically

  #  Using 0, Inf as integration limits fails sometimes
  #  Use quantiles instead
  # x.lo <- q.PN(0.0001,lambda,mue,sigma)
  # x.hi <- q.PN(1-0.0001,lambda,mue,sigma)

  VarX <- integrate(xx.pdf.PN, 0, Inf, EX=EX, lambda=lambda, mue=mue, 
                    sigma=sigma)
  return(VarX[[1]]) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  xx.pdf.PN.R

xx.pdf.PN <- function(x, EX, lambda,mue,sigma,fastnull=1.e-10)
{
  #  Auxiliary function for Var.PN
  #  Parameters of the PND: lambda, mue, sigma
  #  siehe Freeman & Modarra, SPL, 2006, p. 766

  #  #  CHANGE HISTORY
  #  27.07.2019 Installation
  # ===========================================================================
  
  return( (x-EX)^2 * pdf.PN(x,lambda,mue,sigma))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

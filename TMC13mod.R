# =============================================================================
#  TMC.R - Truncated Minimum Chi-square estimation of reference limits
#
#  Algorithm: TMC13, revision 2022-11-15a, W. Wosniok, University of Bremen
#             https://user.math.uni-bremen.de/c05c/TMC13_public_2022-11-15a/
#             Wosniok W, Haeckel R. Clin Chem Lab Med 2019;57:1933-47,
#             https://doi.org/10.1515/cclm-2018-1341
#             modified by G. Brandhorst for https://github.com/gubruol/ReferenceRangeR/
#  Usage:     source("TMC.R")
#             res <- tmc(x)            # x numeric or character, "<DL" allowed
#             res$RL1; res$RL2         # estimated reference limits
#             print(res$myplot)        # or grDevices::replayPlot(res$myplot)
#             tmc(x, RL1.p = 0.025, RL2.p = 0.975, s.fact = 0.67)  # defaults
# =============================================================================

tmc <- local({

  EPS <- 1e-10
  LSEQ <- c(0, 0.25, 0.5, 0.75, 1)

  COL <- list(hist = "lightgray", border = "black", tmc = "darkgreen", kde = "black",
              dif = adjustcolor("red", alpha.f = 0.10), difborder = "red",
              exp = rgb(0, 139, 0, 16, maxColorValue = 255),
              expborder = rgb(0, 0, 255, 64, maxColorValue = 255))

  bc <- function(x, l) if (l > EPS) (x^l - 1)/l else log(x)

  cdf.pn <- function(x, l, m, s) {
    if (l <= EPS) return(plnorm(x, m, s))
    T <- 1/(l*s) + m/s; K <- pnorm(T); cdf <- rep(0, length(x)); g <- x > EPS
    cdf[g] <- if (K > EPS) (1/K) * (pnorm((bc(x[g], l) - m)/s) - pnorm(-T)) else 1
    cdf
  }

  pdf.pn <- function(x, l, m, s) {
    pdf <- rep(0, length(x)); g <- x > EPS
    if (l <= EPS) pdf[g] <- dlnorm(x[g], m, s) else {
      pdf[g] <- (1/pnorm(1/(l*s) + m/s)) * (1/(s*sqrt(2*pi))) * (x[g])^(l - 1) *
                exp(-0.5*((bc(x[g], l) - m)/s)^2)
      pdf[g & !is.finite(pdf)] <- EPS
    }
    pdf
  }

  q.pn <- function(p, l, m, s) {
    if (l < EPS) return(qlnorm(p, m, s))
    qV <- qnorm(1 - (1 - p)*pnorm(1/(l*s) + m/s)); qV[qV == -Inf] <- -10; qV[qV == Inf] <- 10
    W <- l*(s*qV + m) + 1
    ifelse(W > 0, W^(1/l), 0)
  }

  Round <- function(x, u) u * trunc((abs(x) + 0.5*u)/u + 1e-15) * sign(x)

  quant <- function(x, p) {
    tb <- table(x); v <- as.numeric(names(tb)); y <- cumsum(tb)/length(x)
    if (length(v) == 1) return(rep(v, length(p)))
    ru <- min(diff(v)); xx <- v + ru/2
    if (any(p < y[1])) { xx <- c(xx[1] - ru/2, xx); y <- c(0, y) }
    k <- findInterval(p, y); n <- length(y); q <- rep(NA_real_, length(p)); i <- k >= 1 & k < n
    q[i] <- xx[k[i]] + (xx[k[i] + 1] - xx[k[i]]) * (p[i] - y[k[i]]) / (y[k[i] + 1] - y[k[i]])
    q[k == n & p == y[n]] <- xx[n]
    q
  }

  proc.dl <- function(v, s.fact) {
    num <- suppressWarnings(as.numeric(v)); DL <- NA
    lt <- !is.na(v) & is.na(num) & grepl("<", v, fixed = TRUE)
    if (any(lt)) {
      DL <- suppressWarnings(max(as.numeric(sub("<", " ", unique(v[lt]), fixed = TRUE)),
                                 na.rm = TRUE))
      if (is.finite(DL)) {
        sur <- as.numeric(as.character(s.fact * DL))
        num[!is.na(num) & num < DL] <- sur; num[lt] <- sur
      } else DL <- NA
    }
    list(x = num, DL = DL, lt = lt)
  }

  find.ru <- function(v) {
    v <- v[is.finite(v) & v > 0]
    if (!length(v)) return(1)
    for (k in 6:-6) { r <- v/10^k; if (all(abs(r - round(r)) < 1e-6)) return(10^k) }
    1e-6
  }

  smooth.hist <- function(h) {
    d <- h$density; nb <- length(d)
    if (nb < 5) return(h)
    k <- if (nb < 20) 2 else 3; w <- rep(1, 2*k + 1)/(2*k + 1); d2 <- d
    for (i in (k + 1):(nb - k)) {
      s <- (i - k):(i + k); xs <- h$mids[s]; xs2 <- xs^2; ys <- d[s]
      d2[i] <- max(0, predict(lm(ys ~ xs + xs2, weights = w))[k + 1])
    }
    br <- h$breaks; bw <- diff(br); d2 <- d2/sum(d2*bw)
    structure(list(breaks = br, counts = sum(h$counts) * d2 * bw, density = d2,
                   mids = (br[-1] + br[-length(br)])/2), class = "histogram")
  }

  collapse.bins <- function(br, cnt, lo, hi, nmin) {
    n <- sum(cnt)
    while (cnt[1] == 0) { br <- br[-1]; cnt <- cnt[-1]; lo <- lo - 1; hi <- hi - 1 }
    while (cnt[length(cnt)] == 0) { br <- br[-length(br)]; cnt <- cnt[-length(cnt)] }
    L <- br[-length(br)]; U <- br[-1]; N <- length(cnt)
    oL <- L[lo:hi]; oC <- cnt[lo:hi]; oU <- U[hi]; i <- lo - 1
    while (i >= 1) {
      s <- cnt[i]
      while (s < nmin && i > 1) { i <- i - 1; s <- s + cnt[i] }
      if (s >= nmin) { oL <- c(L[i], oL); oC <- c(s, oC) }
      else { oL[1] <- L[i]; oC[1] <- oC[1] + s }
      i <- i - 1
    }
    i <- hi + 1
    while (i <= N) {
      s <- cnt[i]; j <- i
      while (s < nmin && j < N) { j <- j + 1; s <- s + cnt[j] }
      if (s >= nmin) { oL <- c(oL, L[i]); oC <- c(oC, s) } else oC[length(oC)] <- oC[length(oC)] + s
      oU <- U[j]; i <- j + 1
    }
    br <- c(oL, oU)
    structure(list(breaks = br, counts = oC, density = oC/(diff(br)*n),
                   mids = (br[-length(br)] + br[-1])/2, xname = "x",
                   equidist = length(unique(diff(br))) == 1), class = "histogram")
  }

  find.histo <- function(x, x.val, ru, DL, nmin) {
    n <- length(x); x.max <- max(x)
    b0 <- if (is.na(DL)) max(x.val[1] - ru/2, 0.05*ru) else 2*EPS
    eqd <- function(h, b2) {
      hh <- hist(x, breaks = c(b0, seq(b2, x.max + 1.01*h*ru, by = h*ru)),
                 right = FALSE, plot = FALSE)
      cdf <- cumsum(hh$counts)/n; lo <- min(which(cdf > 0.2)); hi <- max(which(cdf < 0.8))
      list(h = hh, lo = lo, hi = hi,
           go = (8 < hi - lo + 1 & !all(hh$counts[lo:hi] >= nmin)) | (hi - lo + 1 > 50))
    }
    h <- 1; e <- eqd(1, x.val[1] + ru/2); go <- e$go
    while (go) {
      h <- h + 1; e1 <- eqd(h, b0 + h*ru/2); go <- e1$go
      if (go) e <- e1
    }
    hh <- e$h; lo <- e$lo; hi <- e$hi
    if (!is.na(DL)) {
      hh <- hist(x, breaks = c(hh$breaks[1], hh$breaks[hh$breaks >= DL]),
                 right = FALSE, plot = FALSE)
      cdf <- cumsum(hh$counts)/n; lo <- min(which(cdf > 0.2)); hi <- max(lo, which(cdf < 0.8))
    }
    hh <- smooth.hist(hh)
    collapse.bins(hh$breaks, hh$counts, lo, hi, nmin/2)
  }

  proportions <- function(cnt, mi, pmin) {
    N <- length(cnt); n <- sum(cnt); cdf <- cumsum(cnt)/n
    idx <- vapply(c(0, seq(0.04, 0.96, by = 0.03), 1), function(p) which.min(abs(p - cdf))[1], 1L)
    D <- unique(data.frame(idx = idx, cdf = cdf[idx])); D <- D[order(D$cdf), ]; P <- nrow(D)
    out <- NULL
    for (a in seq_len(P - 1)) for (b in (a + 1):P) {
      lo <- D$idx[a]; hi <- D$idx[b]; ntr <- sum(cnt[lo:hi]); pr <- ntr/n
      if (hi - lo + 1 >= 6 && pmin <= pr) {
        nlt <- if (lo > 1) sum(cnt[1:(lo - 1)]) else 0
        nge <- if (hi < N) sum(cnt[(hi + 1):N]) else 0
        out <- rbind(out, c(ilo = lo, ihi = hi, prop.lo = nlt/n, x.tr.prop = pr, prop.hi = nge/n,
                            x.lt.tr.n = nlt, x.ge.tr.n = nge))
      }
    }
    if (is.null(out)) return(NULL)
    mp <- cdf[mi]; i3 <- tail(which(cdf <= mp - 0.15), 1); j3 <- which(cdf >= mp + 0.15)[1]
    q <- ifelse(out[, "ilo"] <= mi & mi <= out[, "ihi"], 2, 1)
    if (length(i3) && !is.na(j3)) q[out[, "ilo"] <= i3 & out[, "ihi"] >= j3] <- 3
    out <- data.frame(out, TI.qual = q)
    out[order(-out$TI.qual, -out$x.tr.prop), ]
  }

  wawo <- function(s) {
    n <- length(s); n1 <- sum(s == 1); n0 <- n - n1
    if (n0 == 0 || n1 == 0) return(list(z = qnorm(EPS), p = EPS))
    r <- 1 + sum(s[-1] != s[-n]); mu <- 2*n0*n1/n + 1
    z <- (r - mu)/sqrt((mu - 1)*(mu - 2)/(n - 1))
    list(z = z, p = pnorm(z))
  }

  oc2 <- function(chi2, df, node.p = 0.1) {
    node <- qchisq(1 - node.p, df)
    if (chi2 <= node) return(pchisq(chi2, df))
    sl <- dchisq(node, df); 1 - node.p - sl*node + sl*chi2
  }

  chi2trunc <- function(h, ilo, ihi, nlt, nge, l, m, s, full = FALSE) {
    cnt <- h$counts; N <- length(cnt); n <- sum(cnt); ntr <- sum(cnt[ilo:ihi])
    F.lo <- cdf.pn(h$breaks[ilo], l, m, s); F.hi <- cdf.pn(h$breaks[ihi + 1], l, m, s)
    c.n <- min(ntr/max(1e-6, F.hi - F.lo), n)
    p <- diff(cdf.pn(h$breaks, l, m, s)); p[p < 1e-10] <- 1e-10; p <- p/sum(p)
    e <- c.n * p; d <- cnt - e; TI <- ilo:ihi
    ww <- wawo(1 * (d[TI] > 0))
    chi.TI <- sum(d[TI]^2/e[TI]); df <- ihi - ilo - 3
    out <- function(i, nobs) {
      if (!length(i)) return(NA); eo <- sum(e[i]); dif <- nobs - eo; c2 <- dif^2/eo
      if (dif >= 0) 0 else c2 * pchisq(c2, 1)
    }
    chi.path <- sum(c(out(seq_len(ilo - 1), nlt), out(if (ihi < N) (ihi + 1):N, nge)), na.rm = TRUE)
    oc <- oc2(chi.TI + chi.path + ww$z^2 * (ww$z < 0), df + 1)
    if (!full) return(oc)
    list(oc = oc, p.fit = 1 - pchisq(chi.TI, df), p.rt = ww$p,
         prev = c(nlt - c.n*F.lo, c.n, nge - c.n*(1 - F.hi))/n,
         tab = data.frame(x.lo = h$breaks[-(N + 1)], x.hi = h$breaks[-1], count = cnt,
                          e = e, d = d))
  }

  kde.mode <- function(x, ru) {
    k <- density(x, n = 4096); from <- min(k$x[k$x > 0]); to <- max(k$x); adj <- 1
    k <- density(x, n = 4096, from = from, to = to)
    if (k$bw < 0.75*ru) adj <- 0.75*ru/k$bw
    repeat {
      k <- density(x, n = 4096, adjust = adj, from = from, to = to)
      cdf <- cumsum(k$y); cdf <- cdf/tail(cdf, 1); inner <- 0.025 <= cdf & cdf < 0.975
      xx <- k$x[inner]; yy <- k$y[inner]; m <- length(xx)
      mid <- yy[2:(m - 1)]; xm <- xx[2:(m - 1)]; is.max <- yy[1:(m - 2)] <= mid & mid > yy[3:m]
      if (!any(is.max)) return(NA)
      if (sum(is.max) == 1) return(xm[is.max])
      dx <- abs(diff(xm[is.max]))
      if (!any(0.55*ru <= dx & dx <= 1.45*ru)) return(xm[which.max(mid)])
      adj <- adj*1.05
    }
  }

  draw <- function(h, x.val, r, tab, mode, xlim) {
    op <- par(las = 1, tcl = 0.5, mar = c(3.1, 5.6, 1.1, 2.1)); on.exit(par(op))
    ymax <- 1.1 * max(h$density); ymin <- -0.05 * ymax
    plot(h, freq = FALSE, xlim = xlim, ylim = c(ymin, ymax), col = COL$hist, border = COL$border,
         xlab = "", ylab = "", main = "")
    title(ylab = "Density", line = 4)
    if (length(x.val) > 1) points(x.val, rep(0.5*ymin, length(x.val)), type = "h", col = "gray10")
    xs <- seq(xlim[1], xlim[2], length.out = 101)
    lines(xs, r$prev[2] * pdf.pn(xs, r$lambda, r$mue, r$sigma), col = COL$tmc, lwd = 2)
    bw <- tab$x.hi - tab$x.lo; nn <- sum(tab$count)
    hst <- function(cnt, dens) structure(list(breaks = c(tab$x.lo, tail(tab$x.hi, 1)),
      counts = cnt, density = dens, mids = (tab$x.lo + tab$x.hi)/2, xname = "", equidist = FALSE),
      class = "histogram")
    plot(hst(tab$e, tab$e/(nn*bw)), border = COL$expborder, col = COL$exp, add = TRUE)
    plot(hst(tab$d, pmax(tab$d/(nn*bw), ymin)), border = COL$difborder, col = COL$dif, add = TRUE)
    if (!is.na(mode)) abline(v = mode, col = COL$kde, lty = 2)
    for (k in 1:2) {
      RL <- c(r$RL1, r$RL2)[k]; lines(c(RL, RL), c(0, 0.85*ymax), col = COL$tmc)
      text(RL, 0.85*ymax, pos = c(2, 4)[k], labels = format(RL, digits = 3), cex = 0.8,
           col = COL$tmc)
    }
    lines(r$TI[1:2], rep(0.9*ymin, 2), col = COL$tmc, lwd = 2)
  }

  function(x, RL1.p = 0.025, RL2.p = 0.975, s.fact = 0.67) {
    if (is.factor(x)) x <- as.character(x)
    dl <- proc.dl(x, s.fact); DL <- dl$DL
    v <- dl$x[!is.na(dl$x) & dl$x >= 0]
    if (length(v) < 300) stop(sprintf("TMC needs at least 300 values (got %d).", length(v)),
                              call. = FALSE)
    num <- suppressWarnings(as.numeric(x))
    ru <- max(find.ru(num[!dl$lt & !is.na(num) & num >= 0]), 10^ceiling(log10(max(v)/1e6)))
    v <- Round(v, ru); v[v == 0] <- 0.1*ru; x <- sort(v); n <- length(x)
    pct.DL <- if (is.na(DL)) 0 else 100*sum(x <= DL)/n
    if (pct.DL > 65)
      stop(sprintf("TMC: %.1f%% of values are below the detection limit (max. 65%%).", pct.DL),
           call. = FALSE)
    Q2 <- quant(x, RL2.p)
    x.val <- as.numeric(names(table(x)))
    h <- find.histo(x, x.val, ru, DL, max(Round(n/150, 1), 50))
    TI <- proportions(h$counts, which.max(h$density), if (is.na(DL)) 0.55 else 0.40)
    if (is.null(TI)) stop("TMC: no truncation interval candidates found.", call. = FALSE)
    x.red <- if (n < 1000) x else quant(x, seq(0.005, 0.995, by = 0.005))
    m <- length(x.red); red.cdf <- (1:m)/m; qx <- qnorm(ppoints(m))
    Y <- vapply(LSEQ, function(l) bc(x.red, l), numeric(m))

    fit.ti <- function(t) {
      crit <- function(th, full = FALSE)
        chi2trunc(h, t$ilo, t$ihi, t$x.lt.tr.n, t$x.ge.tr.n, th[1], th[2], th[3], full)
      sub <- t$prop.lo < red.cdf & red.cdf < (1 - t$prop.hi)
      qq <- vapply(seq_along(LSEQ), function(j) {
        f <- lm.fit(cbind(1, qx[sub]), Y[sub, j]); fv <- f$fitted.values
        mss <- sum((fv - mean(fv))^2); r2 <- mss/(mss + sum(f$residuals^2))
        c(1 - (1 - r2) * ((sum(sub) - 1)/f$df.residual), f$coefficients)
      }, numeric(3))
      j <- which.max(qq[1, ]); th <- c(LSEQ[j], qq[2, j], qq[3, j])
      y1 <- bc(x[1], th[1]); yn <- bc(x[n], th[1])
      RB <- rbind(c(-EPS, 1 + EPS),
                  c(if (y1 < 0) 2*y1 else if (y1 == 0) 2*EPS else 0.5*y1, 2*yn),
                  c(0.001, 2)*th[3])
      oc0 <- crit(th)
      if (oc0 > 2) {
        f <- c(0.950, 0.970, 0.980, 0.990, 0.995, 1.005, 1.010, 1.020, 1.030, 1.050)
        g <- cbind(rep(f*th[2], each = 10), rep(f*th[3], times = 10))
        ocg <- apply(g, 1, function(p) crit(c(th[1], p))); k <- which.min(ocg)
        if (ocg[k] < oc0) th[2] <- g[k, 1]
      }
      R2o <- function(p) RB[, 1] + (RB[, 2] - RB[, 1])/(1 + exp(-pmax(-50, pmin(p, 50))))
      nl <- nlm(function(p) crit(R2o(p)), log((th - RB[, 1])/(RB[, 2] - th)),
                iterlim = 500, gradtol = 1e-8, steptol = 1e-8)
      th <- R2o(nl$estimate); RL <- q.pn(c(RL1.p, RL2.p), th[1], th[2], th[3])
      pr <- crit(th, full = TRUE)
      c(ilo = t$ilo, ihi = t$ihi, nlt = t$x.lt.tr.n, nge = t$x.ge.tr.n, prop = t$x.tr.prop,
        lambda = th[1], mue = th[2], sigma = th[3], RL1 = RL[1], RL2 = RL[2], oc = nl$minimum,
        p.fit = pr$p.fit, p.rt = pr$p.rt,
        score = sum(c(TRUE, Q2 + ru >= RL[2], -0.05 <= pr$prev & pr$prev <= 1.05)))
    }

    R <- NULL
    for (i in seq_len(nrow(TI))) {
      R <- rbind(R, fit.ti(TI[i, ]))
      acc <- isTRUE(R[i, "score"] == 5 && R[i, "p.fit"] >= 0.025 && R[i, "p.rt"] >= 0.025)
      if (acc) break
    }
    b <- as.list(R[if (acc) i else order(-R[, "score"], R[, "oc"])[1], ])
    pr <- chi2trunc(h, b$ilo, b$ihi, b$nlt, b$nge, b$lambda, b$mue, b$sigma, full = TRUE)
    z <- qnorm(0.975)
    ci <- function(p, Q)
      Q + c(-1, 1) * z * sqrt(p*(1 - p)/(pdf.pn(Q, b$lambda, b$mue, b$sigma)^2 * n))
    res <- list(RL1 = b$RL1, RL2 = b$RL2, lambda = b$lambda, mue = b$mue, sigma = b$sigma,
                RL1.ci = ci(RL1.p, b$RL1), RL2.ci = ci(RL2.p, b$RL2),
                prev = setNames(pr$prev, c("l", "c", "r")),
                TI = c(lo = h$breaks[b$ilo], hi = h$breaks[b$ihi + 1], prop = b$prop),
                p.fit = b$p.fit, p.rt = b$p.rt, opt.crit = b$oc, score = b$score, accepted = acc,
                n = n, DL = DL, pct.lt.DL = pct.DL, round.unit = ru, myplot = NULL)
    if (abs(res$lambda) < EPS) res$lambda <- 0

    q <- stats::quantile(x, c(0.1, 0.9), names = FALSE)
    xlim <- c(max(0, q[1] - (q[2] - q[1])/1.3), q[2] + (q[2] - q[1])/1.3)
    mode <- tryCatch(kde.mode(x, ru), error = function(e) NA)
    old <- grDevices::dev.cur()
    grDevices::pdf(NULL, width = 9, height = 6); grDevices::dev.control("enable")
    dev <- grDevices::dev.cur()
    res$myplot <- tryCatch({ draw(h, x.val, res, pr$tab, mode, xlim)
                             grDevices::recordPlot() },
                           error = function(e) { warning("TMC plot failed: ", conditionMessage(e))
                                                 NULL })
    grDevices::dev.off(dev); if (old > 1) grDevices::dev.set(old)
    res
  }
})

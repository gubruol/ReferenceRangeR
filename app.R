library(shiny)
library(bslib)
library(ggplot2)

BASE_PATH <- switch(Sys.info()[["sysname"]],
                    Windows = "Z:/R/kc.uol.de/referenceranger/",
                    Linux   = "/srv/shiny-server/referenceranger/",
                    "./referenceranger/")
options(shiny.maxRequestSize = 200 * 1024^2)
ASYNC <- requireNamespace("mirai", quietly = TRUE)
if (ASYNC) { mirai::daemons(2); onStop(function() mirai::daemons(0)) }


### CORE Function Definition Start ########################################################

CORE <- quote({
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  MALE_LABELS    <- c("male", "männlich", "Mann", "M", "m")
  FEMALE_LABELS  <- c("female", "weiblich", "Frau", "F", "f", "W", "w")
  DIVERSE_LABELS <- c("D", "d", "diverse", "Diverse")
  SEX_LEVELS     <- c("F", "M", "D", "X")
  MIN_N <- 500L; WARN_N_SMALL <- 2000L; WARN_N_LARGE <- 100000L; MIN_GROUP_N <- 100L
  RI_METHODS <- c(refineR = "refineR", tmc = "TMC", tml = "TML",
                  kosmic = "kosmic", reflimr = "reflimR")
  METHOD_REFERENCES <- c(
    refineR = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8346497/",
    tmc     = "https://doi.org/10.1515/cclm-2018-1341",
    tml     = "https://doi.org/10.1515/CCLM.2007.249",
    kosmic  = "https://www.nature.com/articles/s41598-020-58749-2",
    reflimr = "https://doi.org/10.1515/labmed-2023-0042")
  RRR_GITHUB     <- "https://github.com/gubruol/ReferenceRangeR"
  RRR_HELP       <- "help.pdf"   # served from www/
  RRR_DISCLAIMER <- "https://kc.uol.de/disclaimer/"
  RRR_REMOTE     <- "https://raw.githubusercontent.com/gubruol/ReferenceRangeR/refs/heads/main/"
  TML_PACKAGES <- c("MASS", "nlme", "msm", "geoR", "modeest", "admisc", "date")
  TMC_PACKAGES <- c(TML_PACKAGES, "snpar", "tseries", "mgcv", "stringr")
  
  need_package <- function(pkg, purpose) if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("The '%s' package is needed to %s. Please install it.", pkg, purpose), call. = FALSE)
  
  as_number <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))
  
  age_limits_valid <- function(from, to) isTRUE(!is.na(from) && !is.na(to) && to > 0 && to > from)
  
  too_few <- function(n, what) list(ok = FALSE, n = n, message = sprintf(
    "The selection contains %s %s. At least %s are needed.",
    format(n, big.mark = ","), what, format(MIN_N, big.mark = ",")))
  
  read_table_file <- function(path, filename = path) {
    out <- if (tolower(tools::file_ext(filename)) %in% c("xlsx", "xls")) {
      need_package("readxl", "read Excel files")
      as.data.frame(readxl::read_excel(path, col_types = "text", guess_max = 100000))
    } else {
      need_package("data.table", "read text files")
      data.table::fread(path, data.table = FALSE, showProgress = FALSE,
                        colClasses = "character", na.strings = c("", "NA"))
    }
    if (nrow(out) == 0L) stop("The file contains no data rows.", call. = FALSE)
    out
  }
  
  read_pasted_text <- function(txt) {
    lines <- strsplit(trimws(txt %||% ""), "\r?\n")[[1]]
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines) < 2L) stop("Please paste a header line and at least one data row.", call. = FALSE)
    seps <- c("\t", ";", ",", " +")
    n <- vapply(seps, function(s) { m <- gregexpr(s, lines[1])[[1]]
    if (m[1] == -1L) 0L else length(m) }, integer(1))
    if (max(n) == 0L) stop("No column separator found in the header line.", call. = FALSE)
    parts <- strsplit(lines, seps[which.max(n)])
    width <- max(lengths(parts))
    m <- do.call(rbind, lapply(parts, function(p) c(p, rep(NA_character_, width - length(p)))))
    out <- as.data.frame(m[-1, , drop = FALSE], stringsAsFactors = FALSE)
    nms <- trimws(m[1, ]); bad <- is.na(nms) | !nzchar(nms)
    nms[bad] <- paste0("V", seq_along(nms))[bad]
    names(out) <- nms
    out
  }
  
  guess_columns <- function(df) {
    nms <- names(df); low <- tolower(nms)
    pick <- function(patterns) {
      for (p in patterns) { hit <- which(grepl(p, low)); if (length(hit)) return(nms[hit[1]]) }
      NULL
    }
    result <- pick(c("^result$", "^value$", "^wert$", "result", "value", "wert",
                     "conc", "analyte", "measure"))
    if (is.null(result)) {
      share <- vapply(df, function(x) mean(!is.na(as_number(x))), numeric(1))
      if (any(share > 0.8)) result <- nms[which.max(share)]
    }
    list(result = result,
         age = pick(c("^age$", "^alter$", "age", "alter", "years", "jahre")),
         sex = pick(c("^sex$", "^gender$", "^geschlecht$", "sex", "gender", "geschlecht")),
         trimester = pick(c("trimester", "trimenon", "gestation", "ssw")))
  }
  
  normalise_sex <- function(x, female = NULL, male = NULL, diverse = NULL) {
    x <- trimws(as.character(x)); x[!nzchar(x)] <- NA_character_
    out <- ifelse(is.na(x), NA_character_, "X")
    out[x %in% c(FEMALE_LABELS,  female )] <- "F"
    out[x %in% c(MALE_LABELS,    male   )] <- "M"
    out[x %in% c(DIVERSE_LABELS, diverse)] <- "D"
    factor(out, levels = SEX_LEVELS)
  }
  
  prepare_data <- function(raw, cols, sex_labels = list()) {
    if (is.null(cols$result) || !cols$result %in% names(raw))
      stop("Please choose which column holds the measured result.", call. = FALSE)
    n_input <- nrow(raw)
    res_chr <- trimws(as.character(raw[[cols$result]]))
    res_chr[!nzchar(res_chr)] <- NA_character_
    res_chr <- gsub(",", ".", res_chr, fixed = TRUE)
    censored <- !is.na(res_chr) & grepl("^<", res_chr)
    res_num <- suppressWarnings(as.numeric(res_chr))
    discarded <- c("Empty"                 = sum(is.na(res_chr)),
                   "Not a number"          = sum(!is.na(res_chr) & is.na(res_num) & !censored),
                   "Zero or negative"      = sum(!is.na(res_num) & res_num <= 0),
                   "Below limit of quantification" = sum(censored))
    res_num[!is.na(res_num) & res_num <= 0] <- NA_real_
    age_chr <- if (!is.null(cols$age) && cols$age %in% names(raw))
      trimws(as.character(raw[[cols$age]])) else rep(NA_character_, n_input)
    age_given <- !is.na(age_chr) & nzchar(age_chr)
    age <- if (!is.null(cols$age) && cols$age %in% names(raw))
      as_number(raw[[cols$age]]) else rep(NA_real_, n_input)
    sex_raw <- if (!is.null(cols$sex) && cols$sex %in% names(raw))
      as.character(raw[[cols$sex]]) else rep(NA_character_, n_input)
    out <- data.frame(result = res_num, result_chr = res_chr, age = age,
                      sex = normalise_sex(sex_raw, sex_labels$female, sex_labels$male,
                                          sex_labels$diverse),
                      stringsAsFactors = FALSE)
    has_tri <- !is.null(cols$trimester) && cols$trimester %in% names(raw)
    tri_given <- rep(FALSE, n_input)
    if (has_tri) {
      tri_chr <- trimws(as.character(raw[[cols$trimester]]))
      tri_given <- !is.na(tri_chr) & nzchar(tri_chr)
      tri <- suppressWarnings(as.integer(tri_chr))
      tri[!tri %in% 1:3] <- NA_integer_
      out$trimester <- tri
    }
    keep <- !is.na(out$result) | censored
    out <- out[keep, , drop = FALSE]
    rownames(out) <- NULL
    age_given <- age_given[keep]
    tri_given <- tri_given[keep]
    seen <- trimws(sex_raw[!is.na(sex_raw) & nzchar(trimws(sex_raw))])
    custom <- setdiff(seen, c(FEMALE_LABELS, MALE_LABELS, DIVERSE_LABELS))
    unknown <- setdiff(custom, unlist(sex_labels))
    list(data = out, report = list(
      n_input = n_input, n_usable = sum(!is.na(out$result)),
      n_with_age = sum(!is.na(out$age)),
      n_bad_age = sum(age_given & is.na(out$age)),
      n_with_sex = sum(!is.na(out$sex) & out$sex != "X"),
      n_unknown_sex = sum(!is.na(out$sex) & out$sex == "X"),
      n_with_trimester = if (has_tri) sum(!is.na(out$trimester)) else 0L,
      n_bad_trimester = if (has_tri) sum(tri_given & is.na(out$trimester)) else 0L,
      discarded = discarded[discarded > 0], has_trimester = has_tri,
      custom_labels = sort(unique(custom)),
      unknown_labels = sort(unique(unknown))))
  }
  
  subset_data <- function(df, sex = "A", age_from = NA, age_to = NA, trimester = NULL) {
    if (sex %in% c("M", "F", "D")) df <- df[!is.na(df$sex) & df$sex == sex, , drop = FALSE]
    if (!is.null(trimester) && "trimester" %in% names(df) && isTRUE(trimester %in% 1:3))
      df <- df[!is.na(df$trimester) & df$trimester == trimester, , drop = FALSE]
    if (age_limits_valid(age_from, age_to))
      df <- df[!is.na(df$age) & df$age >= age_from & df$age <= age_to, , drop = FALSE]
    df
  }
  
  quiet_reflim <- function(values, ...) tryCatch(
    reflimR::reflim(stats::na.omit(values), plot.it = FALSE, print.n = FALSE, ...),
    error = function(e) NULL, warning = function(w) NULL)
  
  estimate_pu_percent <- function(values, reflim_result = NULL) {
    rl <- reflim_result %||% quiet_reflim(values)
    if (is.null(rl)) return(NA_real_)
    lo <- rl$limits[1]; hi <- rl$limits[2]
    if (!is.finite(lo) || !is.finite(hi) || lo <= 0 || hi <= 0) return(NA_real_)
    2.39 * (-0.25 + 100 * (-1 + exp(((log(hi) - log(lo)) / 3.92)^2))^0.5)^0.5
  }
  
  plot_ylim <- function(values) {
    q <- stats::quantile(stats::na.omit(values), probs = c(0.1, 0.9))
    unname(c(q[1] - (q[2] - q[1]) / 1.3, q[2] + (q[2] - q[1]) / 1.3))
  }
  
  plot_comparison_limits <- function(estimated_low, estimated_high, reference_low, reference_high) {
    pU <- reflimR::permissible_uncertainty(reference_low, reference_high)
    inside <- function(x, lo, hi) if (x > lo && x < hi) "#10D010" else "red"
    col1 <- inside(estimated_low, pU[1], pU[2]); col2 <- inside(estimated_high, pU[3], pU[4])
    graphics::rect(pU[1], 0, pU[2], 99999, col = grDevices::adjustcolor(col1, alpha.f = 0.05), border = FALSE)
    graphics::rect(pU[3], 0, pU[4], 99999, col = grDevices::adjustcolor(col2, alpha.f = 0.05), border = FALSE)
    graphics::abline(v = reference_low,  col = grDevices::adjustcolor(col1, alpha.f = 0.40), lwd = 2, lty = 3)
    graphics::abline(v = reference_high, col = grDevices::adjustcolor(col2, alpha.f = 0.40), lwd = 2, lty = 3)
    graphics::abline(v = estimated_low,  col = "#10D010", lwd = 3)
    graphics::abline(v = estimated_high, col = "#10D010", lwd = 3)
  }
  
  record_plot <- function(expr) {
    grDevices::pdf(NULL, width = 9, height = 6)
    on.exit(grDevices::dev.off(), add = TRUE)
    grDevices::dev.control("enable")
    res <- force(expr)
    if (inherits(res, c("gg", "ggplot", "trellis"))) print(res)
    grDevices::recordPlot()
  }
  
  rrr_plot_theme <- function(base_size = 14) ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "#E9ECEF"),
                   axis.text  = ggplot2::element_text(colour = "#495057"),
                   axis.title = ggplot2::element_text(colour = "#495057"),
                   plot.margin = ggplot2::margin(8, 8, 8, 8))
  
  method_script <- function(file, base_path) {
    local <- file.path(base_path, file)
    if (file.exists(local)) return(local)
    cache <- tryCatch(tools::R_user_dir("ReferenceRangeR", "cache"),
                      error = function(e) file.path(tempdir(), "ReferenceRangeR"))
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    cached <- file.path(cache, file)
    if (!file.exists(cached) || file.size(cached) == 0) {
      ok <- tryCatch({ utils::download.file(paste0(RRR_REMOTE, file), cached, mode = "wb", quiet = TRUE)
        file.exists(cached) && file.size(cached) > 0 }, error = function(e) FALSE)
      if (!ok) { unlink(cached)
        stop(sprintf("%s is not in %s and could not be downloaded from the project repository.",
                     file, base_path), call. = FALSE) }
    }
    cached
  }
  
  load_method_scripts <- function(files, base_path, packages) {
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop(sprintf("This method needs the following packages: %s.",
                                      paste(missing, collapse = ", ")), call. = FALSE)
    for (p in packages) suppressPackageStartupMessages(
      library(p, character.only = TRUE, quietly = TRUE))
    for (f in files) source(method_script(f, base_path))
    invisible(TRUE)
  }
  
  with_method_device <- function(expr) {
    owd <- setwd(tempdir()); devices <- grDevices::dev.list()
    grDevices::pdf(NULL, width = 9, height = 6)
    grDevices::dev.control("enable")
    on.exit({ setwd(owd)
      for (d in setdiff(grDevices::dev.list(), devices)) try(grDevices::dev.off(d), silent = TRUE) },
      add = TRUE)
    force(expr)
  }
  
  compute_reference_interval <- function(values, values_chr = NULL, method = "refineR",
                                         opts = list(), base_path = ".") {
    values <- as.numeric(stats::na.omit(values))
    chr <- values_chr[!is.na(values_chr) & nzchar(values_chr)]
    n <- if (method == "tmc") length(chr) else length(values)
    if (n < MIN_N) return(c(too_few(n, "usable results"), list(method = method)))
    out <- list(ok = TRUE, method = method, n = n, limits = c(NA_real_, NA_real_),
                ci = NULL, lambda = NA_real_, lambda_note = "", note = "",
                plot = NULL)
    
    if (method == "refineR") {
      need_package("refineR", "use the refineR method")
      fit <- refineR::findRI(Data = values,
                             model = if (isTRUE(opts$modboxcox)) "modBoxCox" else "BoxCox",
                             NBootstrap = opts$nbootstrap %||% 0)
      ri <- refineR::getRI(fit)
      out$limits <- c(ri[1, 2], ri[2, 2]); out$lambda <- fit$Lambda
      if ((opts$nbootstrap %||% 0) > 0) out$ci <- c(ri[1, 3], ri[1, 4], ri[2, 3], ri[2, 4])
      out$plot <- record_plot(plot(fit))
      
    } else if (method == "tmc") {
      res <- with_method_device({
        load_method_scripts(c("TMC.settings.R", "TMC.functions.R"), base_path, TMC_PACKAGES)
        tmc(chr)
      })
      out$limits <- unname(c(res$RL1, res$RL2)); out$lambda <- unname(res$lambda)
      out$plot <- res$myplot
      
    } else if (method == "tml") {
      q90 <- stats::quantile(values, probs = 0.90)
      decimals <- if (q90 >= 1000) 0 else if (q90 >= 100) 1 else if (q90 >= 10) 2 else 3
      if (isTRUE(opts$fast_tml) && decimals > 0) decimals <- decimals - 1
      v <- round(values, decimals)
      rl <- quiet_reflim(v)
      path_right <- is.null(rl) || sum(v > rl$limits[2]) > sum(v < rl$limits[1])
      res <- with_method_device({
        load_method_scripts("TML.R", base_path, TML_PACKAGES)
        tml(v, path_right)
      })
      out$limits <- unname(c(res$DL25, res$DL975)); out$lambda <- unname(res$lambda)
      out$plot <- res$myplot
      
    } else if (method == "kosmic") {
      need_package("tidykosmic", "use the kosmic method")
      if (!any((values %% 1) > 0)) return(list(ok = FALSE, n = n, method = method, message = paste(
        "kosmic needs results with at least one decimal place.",
        "The selected data are all whole numbers.")))
      fit <- tidykosmic::kosmic(values, decimals = 1)
      s <- summary(fit)
      out$limits <- unname(c(s[1], s[3])); out$plot <- record_plot(plot(fit))
      
    } else if (method == "reflimr") {
      caught <- new.env(parent = emptyenv())
      out$plot <- record_plot(caught$fit <- reflimR::reflim(values, targets = opts$targets))
      out$limits <- unname(c(caught$fit$limits[1], caught$fit$limits[2]))
      out$lambda <- 1 - as.numeric(isTRUE(caught$fit$lognormal))
      out$lambda_note <- "(binary result)"
      
    } else stop("Unknown method: ", method, call. = FALSE)
    out
  }
  
  compute_sex_difference <- function(df) {
    df <- df[!is.na(df$result) & !is.na(df$sex), , drop = FALSE]
    df$sex <- droplevels(df$sex)
    n_before <- nrow(df)
    df <- df[df$sex %in% names(which(table(df$sex) >= MIN_GROUP_N)), , drop = FALSE]
    df$sex <- droplevels(df$sex)
    n <- nrow(df)
    if (n < MIN_N) {
      res <- too_few(n, "results with a recorded sex")
      if (n < n_before) res$message <- paste(res$message, sprintf(
        "(%s were dropped because their sex group held fewer than %s members.)",
        format(n_before - n, big.mark = ","), MIN_GROUP_N))
      return(res)
    }
    n_groups <- nlevels(df$sex); p_value <- NA_real_; test_used <- ""
    if (n_groups == 2L) {
      p_value <- stats::wilcox.test(result ~ sex, data = df)$p.value
      test_used <- "Wilcoxon rank-sum test"
    } else if (n_groups > 2L) {
      p_value <- stats::kruskal.test(result ~ sex, data = df)$p.value
      test_used <- "Kruskal-Wallis test"
    }
    medians <- stats::aggregate(result ~ sex, data = df, FUN = stats::median)
    median_table <- data.frame(sex = as.character(medians$sex),
                               n = as.integer(table(df$sex)[as.character(medians$sex)]),
                               median = medians$result, stringsAsFactors = FALSE)
    list(ok = TRUE, n = n, data_removed = n < n_before, n_groups = n_groups,
         p_value = p_value, test_used = test_used, median_table = median_table,
         median_diff = round((max(median_table$median) - min(median_table$median)) /
                               min(median_table$median) * 100, 2),
         pu_percent = estimate_pu_percent(df$result),
         plot = ggplot2::ggplot(df, ggplot2::aes(.data$sex, .data$result)) +
           ggplot2::geom_violin(fill = "#40668d", alpha = 0.18, colour = NA) +
           ggplot2::geom_boxplot(fill = "#40668d", colour = "#1d3557",
                                 outlier.shape = NA, width = 0.12, alpha = 0.9) +
           ggplot2::coord_cartesian(ylim = plot_ylim(df$result)) +
           ggplot2::labs(x = NULL, y = NULL) + rrr_plot_theme())
  }
  
  compute_age_drift <- function(df, max_groups = 3L, max_points = 10000L) {
    df <- df[!is.na(df$result) & !is.na(df$age) & df$age > 0, , drop = FALSE]
    n_cases <- nrow(df)
    if (n_cases < MIN_N) return(too_few(n_cases, "results with a recorded age"))
    need_package("qgam", "check the age drift")
    if (n_cases > max_points) df <- df[sample.int(n_cases, max_points), , drop = FALSE]
    
    pu_absolute <- estimate_pu_percent(df$result) * stats::median(df$result) / 100
    if (!is.finite(pu_absolute)) return(list(ok = FALSE, n = n_cases,
                                             message = "The permissible uncertainty could not be estimated for these data."))
    
    age_digits <- 3 - floor(log10(max(df$age)))
    step <- 10^-age_digits
    from <- seq(round(min(df$age), age_digits), round(max(df$age), age_digits), by = step)
    groups <- data.frame(from = from, to = from + step)
    
    fit <- qgam::qgam(result ~ s(age), data = df, qu = 0.5)
    pred <- stats::predict(fit, newdata = data.frame(age = (groups$from + groups$to) / 2), se = TRUE)
    groups$median <- pred$fit
    ci_low <- pred$fit - 1.96 * pred$se.fit; ci_high <- pred$fit + 1.96 * pred$se.fit
    bin <- findInterval(df$age, groups$from)
    groups$count <- tabulate(bin[bin >= 1L], nbins = nrow(groups))
    keep <- groups$count > 0L
    groups <- groups[keep, , drop = FALSE]
    
    drift_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = df, ggplot2::aes(.data$age, .data$result),
                          colour = "grey60", alpha = 0.10, size = 0.8) +
      ggplot2::geom_ribbon(data = data.frame(x = from[keep], lo = ci_low[keep], hi = ci_high[keep]),
                           ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
                           fill = "#00406b", alpha = 0.20) +
      ggplot2::geom_line(data = groups, ggplot2::aes(.data$from, .data$median),
                         colour = "#00406b", linewidth = 1) +
      ggplot2::coord_cartesian(ylim = plot_ylim(df$result)) +
      ggplot2::labs(x = "age", y = NULL) + rrr_plot_theme()
    
    ord <- order(df$age); ages <- df$age[ord]; vals <- df$result[ord]
    window_median <- function(lo, hi) {
      i <- findInterval(lo, ages, left.open = TRUE) + 1L
      j <- findInterval(hi, ages)
      if (j < i) NA_real_ else stats::median(vals[i:j])
    }
    repeat {
      if (nrow(groups) <= 1L) break
      gaps <- diff(groups$median)
      idx <- which.min(gaps)
      if (gaps[idx] >= (if (nrow(groups) > max_groups) 0.50 else 0.25) * pu_absolute) break
      groups$to[idx]     <- groups$to[idx + 1L]
      groups$count[idx]  <- groups$count[idx] + groups$count[idx + 1L]
      groups$median[idx] <- window_median(groups$from[idx], groups$to[idx])
      groups <- groups[-(idx + 1L), , drop = FALSE]
    }
    groups <- groups[groups$count >= 50L, , drop = FALSE]
    if (nrow(groups) == 0L) return(list(ok = FALSE, n = n_cases,
                                        message = "No age group large enough to report was found."))
    
    groups$from <- round(groups$from, age_digits - 1)
    groups$to   <- round(groups$to,   age_digits - 1)
    groups$median <- round(groups$median, if (max(groups$median) > 100) 0
                           else if (max(groups$median) > 10) 1 else 2)
    groups$size <- round(100 * groups$count / sum(groups$count), 0)
    if (nrow(groups) >= 2L) {
      shade <- scales::col_numeric(palette = c("red", "red", rep("#10D010", 8)),
                                   domain = c(0, 100))(groups$size)
      drift_plot <- drift_plot +
        ggplot2::geom_rect(data = transform(groups, shade = shade),
                           ggplot2::aes(xmin = .data$from, xmax = .data$to,
                                        ymin = -Inf, ymax = Inf, fill = .data$shade),
                           alpha = 0.10, inherit.aes = FALSE) +
        ggplot2::scale_fill_identity() + ggplot2::theme(legend.position = "none")
    }
    list(ok = TRUE, n = n_cases, n_used = nrow(df), subsampled = n_cases > max_points,
         groups = groups[, c("from", "to", "median", "size")],
         stratification_needed = nrow(groups) >= 2L, plot = drift_plot)
  }
  
  format_report_limit <- function(x) if (!isTRUE(is.finite(x))) "-" else
    formatC(x, format = "f", big.mark = ",",
            digits = if (abs(x) > 100) 0 else if (abs(x) > 10) 1 else 2)
  
  method_settings_text <- function(method, opts = list()) {
    nb <- opts$nbootstrap %||% 0
    switch(method,
           refineR = sprintf("%s, %s", if (isTRUE(opts$modboxcox)) "modified Box-Cox" else "Box-Cox",
                             if (nb > 0) paste(nb, "bootstrap iterations") else "no bootstrap"),
           tml = if (isTRUE(opts$fast_tml)) "fast mode (3 significant figures)" else "full precision",
           kosmic = "1 decimal place", "default settings")
  }
  
  wrap_text <- function(txt, width, cex) {
    words <- strsplit(txt, " +")[[1]]
    lines <- character(); current <- ""
    for (w in words) {
      candidate <- if (nzchar(current)) paste(current, w) else w
      if (nzchar(current) && graphics::strwidth(candidate, cex = cex) > width) {
        lines <- c(lines, current); current <- w
      } else current <- candidate
    }
    c(lines, current)
  }
  
  plot_raster <- function(p, width_in, height_in, res = 200, overlay = NULL) {
    if (is.null(p) || !requireNamespace("png", quietly = TRUE)) return(NULL)
    tmp <- tempfile(fileext = ".png")
    on.exit(unlink(tmp), add = TRUE)
    previous <- grDevices::dev.cur()
    opened <- tryCatch({ grDevices::png(tmp, width = width_in, height = height_in,
                                        units = "in", res = res, bg = "white"); TRUE },
                       error = function(e) FALSE)
    if (!opened) return(NULL)
    ours <- grDevices::dev.cur()
    drawn <- tryCatch({ grDevices::replayPlot(p); TRUE }, error = function(e) FALSE)
    if (drawn && is.function(overlay)) try(overlay(), silent = TRUE)
    grDevices::dev.off(ours)
    if (previous > 1L && previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
    if (!drawn) return(NULL)
    tryCatch(png::readPNG(tmp), error = function(e) NULL)
  }
  
  export_ri_pdf <- function(path, result, settings = list(), source_label = NULL,
                            generated = Sys.time()) {
    navy <- "#00406b"; grey <- "#5A6672"; line <- "#D6DBE1"
    dash <- function(x) if (is.null(x) || !nzchar(as.character(x))) "-" else as.character(x)
    W <- 8.27; H <- 11.69; L <- 0.055; R <- 0.945; COL2 <- 0.52; VAL <- 0.125
    
    cmp <- settings$comparison
    overlay <- if (!is.null(cmp) && result$method %in% c("refineR", "tmc", "tml"))
      function() plot_comparison_limits(result$limits[1], result$limits[2], cmp[1], cmp[2])
    
    grDevices::pdf(path, width = W, height = H, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    
    y <- 0.965
    put <- function(txt, x = L, cex = .72, font = 1, col = "black", adj = 0)
      graphics::text(x, y, txt, adj = c(adj, .5), cex = cex, font = font, col = col)
    down <- function(dy) y <<- y - dy
    rule <- function() { down(.014); graphics::segments(L, y, R, y, col = line); down(.026) }
    section <- function(title, x = L) { put(toupper(title), x = x, cex = .66, font = 2, col = navy)
      down(.036) }
    fit <- function(txt, x, cex, limit) {
      txt <- as.character(txt)
      if (x + graphics::strwidth(txt, cex = cex) <= limit) return(txt)
      while (nchar(txt) > 4L && x + graphics::strwidth(paste0(txt, "..."), cex = cex) > limit)
        txt <- substr(txt, 1L, nchar(txt) - 1L)
      paste0(txt, "...")
    }
    field <- function(label, value, x, limit) {
      put(label, x = x, cex = .70, col = grey)
      put(fit(dash(value), x + VAL, .70, limit), x = x + VAL, cex = .70)
      down(.026)
    }
    
    put("ReferenceRangeR", cex = 1.3, font = 2, col = navy)
    put(format(generated, "%d %B %Y"), x = R, cex = .8, col = grey, adj = 1)
    down(.018); rule()
    
    section("Estimated reference interval"); down(.004)
    put(sprintf("%s  –  %s", format_report_limit(result$limits[1]),
                format_report_limit(result$limits[2])), cex = 1.7, font = 2, col = navy)
    down(.048)
    if (!is.null(result$ci)) {
      put(sprintf("95%% confidence interval:  lower %s – %s,  upper %s – %s",
                  format_report_limit(result$ci[1]), format_report_limit(result$ci[2]),
                  format_report_limit(result$ci[3]), format_report_limit(result$ci[4])),
          cex = .72, col = grey)
      down(.028)
    }
    if (!is.null(result$note) && nzchar(result$note)) {
      put(result$note, cex = .72, col = grey); down(.028)
    }
    if (isTRUE(is.finite(result$lambda))) {
      lam <- sprintf("%.2f", result$lambda)
      if (nzchar(result$lambda_note %||% "")) lam <- paste(lam, result$lambda_note)
      put(bquote("Estimated" ~ lambda * ":" ~ .(lam)), cex = .72, col = grey)
      down(.028)
    }
    down(.014)
    y_cols <- y
    
    section("Selection")
    field("Data source", source_label, L, COL2 - 0.02)
    field("Results used", format(result$n, big.mark = ","), L, COL2 - 0.02)
    field("Sex", switch(settings$sex %||% "A", A = "no selection", M = "male",
                        F = "female", D = "non-binary"), L, COL2 - 0.02)
    field("Age", if (age_limits_valid(settings$age_from, settings$age_to))
      sprintf("%s to %s", settings$age_from, settings$age_to) else "no selection",
      L, COL2 - 0.02)
    if (isTRUE(settings$trimester %in% 1:3))
      field("Trimester", paste0(settings$trimester, "."), L, COL2 - 0.02)
    y_left_end <- y
    
    y <- y_cols
    section("Method", COL2)
    field("Method", RI_METHODS[[result$method]], COL2, 0.955)
    field("Settings", method_settings_text(result$method, settings$opts), COL2, 0.955)
    put("Reference", x = COL2, cex = .70, col = grey); down(.022)
    put(fit(METHOD_REFERENCES[[result$method]], COL2, .56, 0.955), x = COL2, cex = .56, col = grey)
    down(.026)
    if (!is.null(cmp)) {
      down(.012)
      pU <- reflimR::permissible_uncertainty(cmp[1], cmp[2])
      verdict <- function(x, lo, hi) if (isTRUE(x > lo && x < hi))
        "within permissible uncertainty" else "outside permissible uncertainty"
      section("Comparison with known limits", COL2)
      field("Known limits", sprintf("%s – %s", format_report_limit(cmp[1]),
                                    format_report_limit(cmp[2])), COL2, 0.955)
      field("Lower limit", verdict(result$limits[1], pU[1], pU[2]), COL2, 0.955)
      field("Upper limit", verdict(result$limits[2], pU[3], pU[4]), COL2, 0.955)
    }
    y_right_end <- y
    
    y <- 0.128
    rule()
    for (ln in wrap_text(paste(
      "For research purpose only. This application does not constitute a",
      "(in-vitro diagnostic) medical device."), R - L, .62)) {
      put(ln, cex = .62, col = grey); down(.021)
    }
    down(.004)
    for (ln in wrap_text(paste0("Source code : ", RRR_GITHUB,
                                " - Disclaimer: ", RRR_DISCLAIMER), R - L, .62)) {
      put(ln, cex = .62, col = grey); down(.021)
    }
    
    PLOT <- c(x0 = L - 0.012, x1 = R + 0.032, y0 = 0.148,
              y1 = min(y_left_end, y_right_end) - 0.020)
    img <- plot_raster(result$plot, (PLOT[["x1"]] - PLOT[["x0"]]) * W,
                       (PLOT[["y1"]] - PLOT[["y0"]]) * H, overlay = overlay)
    if (!is.null(img))
      graphics::rasterImage(img, PLOT[["x0"]], PLOT[["y0"]], PLOT[["x1"]], PLOT[["y1"]])
    else if (!is.null(result$plot)) {
      graphics::par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
      drawn <- tryCatch({ grDevices::replayPlot(result$plot); TRUE }, error = function(e) FALSE)
      if (drawn && is.function(overlay)) try(overlay(), silent = TRUE)
    }
    invisible(path)
  }
  
  make_demo_data <- function(n = 50000L, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    age <- round(stats::runif(n, min = 1, max = 99), 2)
    base <- stats::SSlogis(age, 5, 60, 6)
    result <- base + stats::rnorm(n, mean = 37.5, sd = 3.75)
    k <- round(n / 100)
    for (b in list(c(1, k, 20, 5), c(k + 1, 3 * k, 28, 3),
                   c(3 * k + 1, 8 * k, 50, 4.5), c(8 * k + 1, 10 * k, 60, 10))) {
      i <- seq.int(b[1], b[2])
      result[i] <- base[i] + stats::rnorm(length(i), b[3], b[4])
    }
    data.frame(result = as.character(round(result, 2)), age = age,
               sex = sample(c("M", "F"), n, replace = TRUE), stringsAsFactors = FALSE)
  }
  
})

eval(CORE)

### Layout Definitions Start ########################################################

rrr_theme <- bs_theme(version = 5, preset = "shiny", primary = "#40668d",
                      "border-radius" = "0.75rem", "card-cap-padding-y" = "0.5rem")

UMO_SVG <- local({
  f <- "www/umo.svg"
  if (!file.exists(f)) return(NULL)
  s <- sub("<\\?xml[^>]*\\?>", "", paste(readLines(f, warn = FALSE), collapse = "\n"))
  HTML(gsub("#ffffff", "#ffffff", s, ignore.case = TRUE))
})

rrr_head <- tagList(
  useBusyIndicators(),
  tags$style(HTML("
    :root {
      --rrr-page: #ced7e0; --rrr-surface: #FFFFFF; --rrr-line: rgba(16,24,40,.09); --rrr-header: #00406b;
      --rrr-brand: #ffffff; --bs-body-bg: #F4F6FA; --rrr-tile: #e8ecf1;
      --rrr-shadow: 0 1px 2px rgba(16,24,40,.04), 0 6px 16px -6px rgba(16,24,40,.10);
      --rrr-nav-row1-bg: #40668d;  --rrr-nav-row1-h: 64px;
      --rrr-nav-row2-bg: #40668d;  --rrr-nav-row2-h: 40px;
      --rrr-text-light: #838891; --rrr-text-dark: #40444c;
      --rrr-details-light: #ced7e0;
      --rrr-green-text: #3c862e;
    }
    
    body { background: var(--rrr-page); }
    .card { --bs-card-bg: var(--rrr-surface); --bs-card-border-color: var(--rrr-line);
            box-shadow: var(--rrr-shadow); }
    .card-header { background: transparent; border-bottom-color: var(--rrr-line);
                   font-weight: 600; letter-spacing: -0.005em; }
    .navbar { background: var(--rrr-header) !important;
              border-bottom: 1px solid var(--rrr-line) !important; }
    .navbar-brand { color: var(--rrr-brand) !important; font-weight: 660;
                    letter-spacing: -0.02em; display: inline-flex; align-items: center; }
    .navbar-brand img { width: 50px; height: 50px;
                        object-fit: cover; margin-right: .55rem;
                        padding: 3px; background: var(--rrr-nav-row1-bg); box-sizing: border-box; }
                        
    .standard-text-light { color: var(--rrr-text-light) !important; font-weight: 375; font-size: .9rem;
                            letter-spacing: -0.02em; display: inline-flex; align-items: center; }
    .standard-text-green-bold { color: var(--rrr-green-text) !important; font-weight: 500; font-size: .9rem;
                            letter-spacing: -0.02em; display: inline-flex; align-items: center; }
    .standard-text-dark-bold { color: var(--rrr-text-dark) !important; font-weight: 500; font-size: .9rem;
                            letter-spacing: -0.02em; display: inline-flex; align-items: center; }
    .standard-text-dark { color: var(--rrr-text-dark) !important; font-weight: 400; font-size: .9rem;
                            letter-spacing: -0.02em; display: inline-flex; align-items: center; }
    .standard-text-white-bold { color: var(--rrr-brand) !important; font-weight: 550; font-size: .8rem;
                            letter-spacing: -0.02em; display: inline-flex; align-items: center; }

    
    .navbar { background: transparent !important; padding: 0 !important;
              border-bottom: 1px solid var(--rrr-line) !important; }
    .navbar > .container-fluid, .navbar > .container {
      flex-wrap: wrap !important; padding: 0; row-gap: 0; }
    
    /* Row 1: brand + UMO */
    .rrr-nav-row1 {
      flex: 0 0 100%; display: flex; align-items: center;
      background: var(--rrr-nav-row1-bg); min-height: var(--rrr-nav-row1-h);
      padding: 0 1rem; }
    .rrr-nav-row1 .rrr-nav-spacer { flex: 1 1 auto; }
    
    /* Row 2: tab links + gear */
    .rrr-nav-row2 {
      flex: 0 0 100%; display: flex; align-items: center;
      background: var(--rrr-nav-row2-bg); min-height: var(--rrr-nav-row2-h);
      padding: 0 1rem; border-top: 1px solid var(--rrr-line); }
    .rrr-nav-row2 .navbar-collapse { display: flex !important; padding: 0;
                                     flex: 0 0 auto !important; }
    .rrr-nav-row2 .rrr-nav-spacer { flex: 1 1 auto; }
    .rrr-nav-row2 > .nav-item.dropdown { color: var(--rrr-details-light); flex: 0 0 auto; list-style: none; position: relative; }
    .rrr-nav-row2 > .nav-item.dropdown .dropdown-menu { right: 0; left: auto; }

    .navbar > .container-fluid,
    .navbar > .container {
      flex-wrap: wrap !important; align-items: center; row-gap: .0rem; }
      
    /* brand: row 1, left (flex:0 prevents stretching) */
    
    .navbar-brand { flex: 0 0 auto; }
    
    /* UMO + gear: after JS move they are direct children of .container-fluid */
    
    .navbar > .container-fluid > .bslib-nav-spacer,
    .navbar > .container-fluid > .bslib-nav-item,
    .navbar > .container-fluid > .nav-item.dropdown { flex: 0 0 auto; }
    
    /* spacer between brand and UMO/gear */
    
    .navbar > .container-fluid > .bslib-nav-spacer { flex: 1 1 auto; }
    
    /* Row 2: collapse fills full width → tabs inside stay left-aligned */
    
    .navbar .navbar-collapse {
      order: 1; flex: 0 0 100%;
      display: flex !important; align-items: center; padding: 0; }
    .navbar .navbar-nav {
      flex-direction: row; align-items: center; flex-wrap: nowrap; gap: .25rem; }
    .navbar .navbar-nav .dropdown-menu { position: absolute; }
    .navbar .nav-underline { --bs-nav-underline-gap: 1rem; }
    .nav-underline .nav-link { border-bottom-width: 2px; margin-bottom: 0; }
    .navbar .nav-underline .nav-link,
    .navbar .nav-underline :where(ul.nav.navbar-nav > li) > .nav-link {
      padding-top: .25rem; padding-bottom: .25rem; margin-bottom: 0; color: var(--rrr-details-light); }
    .nav-underline .nav-link.active, .nav-underline .show > .nav-link,
    .nav-underline .nav-link:hover, .nav-underline .nav-link:focus {
      color: var(--rrr-brand); border-bottom-color: var(--rrr-brand); }
      
    .card-header .nav-underline .nav-link { padding-bottom: .25rem; }
    
    .navbar .nav-item.dropdown > .nav-link { display: flex; align-items: center;
                                             padding-top: .25rem; padding-bottom: .25rem;
                                             margin-bottom: 0; }
                                             
    .rrr-field > .form-label { margin-bottom: .1rem !important; }
    .rrr-field .bslib-grid, .rrr-field .form-group { margin-bottom: 0 !important; }
    .rrr-field > .form-text { margin-top: .15rem !important; }
    
    .navbar .bslib-nav-item { display: flex; align-items: center; }
    .navbar .bslib-nav-item .shiny-input-container,
    .navbar .bslib-input-switch { margin-bottom: 0; }
    
    .rrr-umo { display: flex; align-items: center; padding: 0 .6rem 0 .2rem; }
    .navbar { padding-top: 7px; padding-bottom: 7px; }
    .rrr-umo svg { height: 50px; width: auto; display: block; }
    
    .navbar .dropdown-menu { padding-top: .35rem; }
    .navbar .dropdown-menu .bslib-nav-item { padding: .1rem 1rem .35rem; }
    .navbar .dropdown-menu .form-switch { padding-left: 1rem; margin: 0; display: flex;
                                          align-items: center; gap: .6rem; }
    .navbar .dropdown-menu .form-switch .form-check-input { float: none; margin: 0; order: 2; }
    .navbar .dropdown-menu .form-switch .form-check-label { order: 1; }
    .rrr-canvas { background: #FFFFFF; border: 1px solid var(--rrr-line);
                  border-radius: .5rem; padding: .25rem; }
                  
    .rrr-tile { background: var(--rrr-tile) !important; color: var(--bs-body-color) !important;
            border: 1px solid var(--rrr-line) !important;
            height: 100px !important;          /* adjust here */
            max-height: none !important;
            min-height: 0 !important;
            flex: 0 0 auto !important; }
    .rrr-tile .value-box-area { padding: 0 .2rem !important;
                            justify-content: center; gap: 0; }
    .rrr-tile .value-box-title { font-weight: 500; opacity: .75;font-size: .9rem; }
    .rrr-tile .value-box-value { font-weight: 640; letter-spacing: -0.02em; font-size: 1.1rem; }
    
    #preview table th, #preview table td { text-align: center; }
    .bslib-sidebar-layout > .sidebar { font-size: .9rem; }
    .bslib-sidebar-layout > .sidebar .form-label { margin-bottom: .2rem; }
    .limit-value { font-variant-numeric: tabular-nums; }
    .report-list { font-size: .875rem; margin-bottom: 0; }
    .report-list dt { font-weight: 500; }
    table.compact-table { font-size: .875rem; margin-bottom: .5rem; }
    .invalid-range input { border-color: var(--bs-danger) !important; }

    /* One scrollbar only: every zone grows to fit its content and the page
       is the single scroller. bslib gives cards, sidebars and tab panes
       'overflow: auto' plus a height borrowed from the viewport, which is
       what produced a second scrollbar inside the preview table and the
       column selector. */
    .bslib-page-navbar > .container-fluid,
    .bslib-page-navbar > .container,
    .bslib-page-navbar .tab-content,
    .bslib-page-navbar .tab-pane { height: auto !important; flex: 0 0 auto; }
    
    .card:not(.rrr-tile) > .card-body { justify-content: flex-start; align-items: stretch; }
    .card:not(.rrr-tile), .card:not(.rrr-tile) > .card-body, .bslib-sidebar-layout > .main {
      height: auto !important; max-height: none !important; overflow: visible !important; }
    /* Above the mobile breakpoint the sidebar is a plain column and may grow
       with the page; below it bslib turns it into a collapsible overlay that
       needs its own height rules, so leave that alone. */
    @media (min-width: 576px) {
      .bslib-sidebar-layout > .sidebar,
      .bslib-sidebar-layout > .sidebar > .sidebar-content {
        height: auto !important; max-height: none !important; overflow: visible !important; }
    }
    /* ...except a card blown up with the expand button, which really is
       supposed to fill the screen. */
    .bslib-card[data-full-screen=\"true\"],
    .bslib-card[data-full-screen=\"true\"] > .card-body {
      height: 100% !important; overflow: auto !important; }
    /* htmlwidget outputs (the preview table) otherwise keep the 400px
       flex-basis the fill system gave them, leaving a gap or a clip. */
    .html-widget-output { height: auto !important; flex: 0 0 auto !important; }
    .bslib-sidebar-layout { --_vert-border: none; }
    /* the placeholder shown before a plot exists needs its own inset: the
       canvas itself is padding-free so the plot can use the whole box */
    .rrr-canvas > .shiny-output-error { padding: 1.1rem; }
    /* 420px is the canvas's natural height, but when a taller sidebar stretches
       the card, bslib would centre the canvas with auto margins and leave a gap
       above and below its border. Grow into the card instead, so the canvas
       border sits as close to the card as it does at the sides. */
    .rrr-canvas { height: 420px !important; flex: 1 0 auto !important;
                  margin-top: 0 !important; margin-bottom: 0 !important; }
    .bslib-card[data-full-screen=\"true\"] .rrr-canvas { height: 100% !important; }
    
    .card-header .nav-underline .nav-link[data-value=Upload].active { color: var(--rrr-header); border-bottom-color: var(--rrr-header);}
    .card-header .nav-underline .nav-link[data-value=Upload]:hover { color: var(--rrr-header); border-bottom-color: var(--rrr-header);}
    .card-header .nav-underline .nav-link[data-value=Upload]:focus { color: var(--rrr-header); border-bottom-color: var(--rrr-header);}
    .card-header .nav-underline .nav-link[data-value=Paste].active { color: var(--rrr-header); border-bottom-color: var(--rrr-header);}
    .card-header .nav-underline .nav-link[data-value=Paste]:hover { color: var(--rrr-header); border-bottom-color: var(--rrr-header);}
    .card-header .nav-underline .nav-link[data-value=Paste]:focus { color: var(--rrr-header); border-bottom-color: var(--rrr-header);}
    .card-header .nav-underline .nav-link[data-value=Demo].active { color: var(--rrr-header); border-bottom-color: var(--rrr-header);}
    .card-header .nav-underline .nav-link[data-value=Demo]:hover { color: var(--rrr-header); border-bottom-color: var(--rrr-header);}
    .card-header .nav-underline .nav-link[data-value=Demo]:focus { color: var(--rrr-header); border-bottom-color: var(--rrr-header);}
    
    #focus
    .card-header .nav-underline .nav-link[data-value=Upload] { color: var(--rrr-text-light); }
    .card-header .nav-underline .nav-link[data-value=Paste] { color: var(--rrr-text-light)}
    .card-header .nav-underline .nav-link[data-value=Demo] { color: var(--rrr-text-light); }
    
  ")),
  tags$script(HTML("
    (function() {
      function rrr_two_row_nav() {
        var c = document.querySelector('.navbar > .container-fluid, .navbar > .container');
        if (!c || c.querySelector('.rrr-nav-row1')) return;
        // with collapsible = FALSE there is no .navbar-collapse, only the ul
        var col = c.querySelector('.navbar-collapse');
        var nav = c.querySelector('.navbar-nav');
        if (!nav) return;

        var row1 = document.createElement('div'); row1.className = 'rrr-nav-row1';
        var row2 = document.createElement('div'); row2.className = 'rrr-nav-row2';

        // row 1: brand, spacer, UMO logo
        var brand = c.querySelector('.navbar-brand');
        if (brand) row1.appendChild(brand);
        var spacer = document.createElement('div'); spacer.className = 'rrr-nav-spacer';
        row1.appendChild(spacer);
        nav.querySelectorAll(':scope > .bslib-nav-spacer').forEach(function(el) { el.remove(); });
        nav.querySelectorAll(':scope > .bslib-nav-item').forEach(function(el) { row1.appendChild(el); });

        // gear dropdown: goes to row 2, right-aligned
        var gear = nav.querySelector(':scope > .nav-item.dropdown');

        // row 2: tab links, then spacer + gear
        row2.appendChild(col || nav);
        var spacer2 = document.createElement('div'); spacer2.className = 'rrr-nav-spacer';
        row2.appendChild(spacer2);
        if (gear) row2.appendChild(gear);
        c.prepend(row2); c.prepend(row1);
      }
      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', rrr_two_row_nav);
      } else { rrr_two_row_nav(); }
    })();
    Shiny.addCustomMessageHandler('rrr-mark-range', function(msg) {
      var el = document.getElementById(msg.id);
      if (el) el.classList.toggle('invalid-range', msg.invalid);
    });"))
)

### Buttons/Function Definitions Start ########################################################

age_range_input <- function(id_from, id_to, label = "Age range") div(
  class = "rrr-field",
  tags$label(label, class = "form-label rrr-tight-label"),
  div(id = paste0(id_from, "_wrap"), class = "rrr-tight-range",
      layout_columns(col_widths = c(6, 6), gap = "0.4rem", fill = FALSE,
                     numericInput(id_from, NULL, value = NA, min = 0),
                     numericInput(id_to, NULL, value = NA, min = 0))),
  div(class = "form-text mt-0", "from / to (leave empty for all ages)"))

task_button <- function(id, label, ...) input_task_button(
  id, label, ...,
  icon_busy = tags$span(class = "spinner-border spinner-border-sm",
                        role = "status", `aria-hidden` = "true"))

sex_radio <- function(id) radioButtons(id, "Sex", inline = FALSE,
                                       choiceNames = c("all", "male", "female"), choiceValues = c("A", "M", "F"))

METHOD_LABELS <- unname(Map(function(name, url)
  HTML(sprintf("<b>%s</b> <a href='%s' target='_blank'>(lit.)</a>", name, url)),
  RI_METHODS, METHOD_REFERENCES[names(RI_METHODS)]))

tile <- function(title, value, note = NULL, note_class = "small")
  value_box(title, value, class = "rrr-tile",
            if (!is.null(note)) span(class = note_class, note))

plot_canvas <- function(id) card_body(class = "rrr-canvas", padding = 0, fill = FALSE,
                                      plotOutput(id, height = "100%"))

html_table <- function(df, digits = 2) {
  cell <- function(x) if (is.numeric(x)) format(round(x, digits), trim = TRUE) else as.character(x)
  tags$table(class = "table table-sm compact-table",
             tags$thead(tags$tr(lapply(names(df), tags$th))),
             tags$tbody(lapply(seq_len(nrow(df)), function(i)
               tags$tr(lapply(seq_along(df), function(j) tags$td(cell(df[i, j])))))))
}

# The age groups double as a shortcut: every row can hand its from/to to the
# age-range inputs, so a group can be carried straight into a calculation.
drift_group_table <- function(groups) {
  cols <- c("from", "to", "median")
  tags$table(class = "table table-sm compact-table",
             tags$thead(tags$tr(lapply(cols, tags$th), tags$th())),
             tags$tbody(lapply(seq_len(nrow(groups)), function(i)
               tags$tr(lapply(cols, function(cn) tags$td(format(groups[[cn]][i], trim = TRUE))),
                       tags$td(class = "text-end py-1",
                               tags$button(type = "button", class = "btn btn-outline-primary btn-sm py-0 px-2",
                                           onclick = sprintf(
                                             "Shiny.setInputValue('use_age_group', %d, {priority: 'event'})", i),
                                           "use"))))))
}

html_dl <- function(rows, dt_class = "col-5", dd_class = "col-7")
  tags$dl(class = "report-list row",
          lapply(rows, function(r) tagList(tags$dt(class = dt_class, r[1]),
                                           tags$dd(class = dd_class, r[2]))))

hint <- function(...) div(class = "text-secondary", ...)

format_limit <- function(x) if (!isTRUE(is.finite(x))) "-" else
  formatC(x, format = "f", big.mark = ",",
          digits = if (abs(x) > 100) 0 else if (abs(x) > 10) 1 else 2)

format_p <- function(p) if (!isTRUE(is.finite(p))) "n/a" else
  if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

big <- function(x) format(x, big.mark = ",", scientific = FALSE)

### UI Definition Start ########################################################

ui <- page_navbar(
  title = tagList(img(src = "www/rrr.svg", alt = ""), "ReferenceRangeR"),
  id = "nav", theme = rrr_theme, header = rrr_head,
  navbar_options = navbar_options(underline = TRUE, collapsible = FALSE),
  fillable = FALSE,
  
  nav_panel("Input Data", value = "data",
            layout_columns(col_widths = c(5, 7),
                             navset_card_underline(
                               nav_panel("Upload", div(class = "pt-2",
                                                      fileInput("file", NULL, width = "100%",
                                                                 accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
                                                                 buttonLabel = "Browse", placeholder = "CSV, TSV or Excel"),
                                                       div(class = "standard-text-light",
                                                           "The separator and the decimal sign are detected automatically.", br(),
                                                           "Values below the limit of quantification may be written as \"<5\"."))),
                              nav_panel(value="Paste", title="Paste Data from Clipboard", div(class = "pt-2",
                                                                          textAreaInput("pasted", NULL, rows = 5, width = "100%",
                                                                                        placeholder = "e.g.\nresult  age sex\n5.4  42  F\n6.1  37  M"),
                                                                          div(class = "standard-text-light", "Paste with a header line, e.g. copied from Excel.",br()),
                                                                          actionButton("use_paste", "Use pasted data", class = "btn-primary w-100"))),
                              nav_panel("Demo", div(class = "pt-2",
                                                     p(class = "standard-text-light", "50,000 simulated results with an age-dependent median ",
                                                       "and about 10% pathological values."),
                                                     actionButton("use_demo", "Generate demo data", class = "btn-primary w-100")))),
                           card(card_header(
                             class = "d-flex align-items-center",
                             span("Options"),
                             span(class = "ms-auto",
                                  input_switch("advanced", div(class = "standard-text-light", "Advanced"), value = FALSE)
                             ),
                             tags$style(HTML("
    .card-header .form-check { margin-bottom: 0; }
    .card-header .shiny-input-container { width: auto !important; }
  "))
                           ),
                                uiOutput("welcome_ui"),
                                uiOutput("column_ui"),
                                div(hr(),uiOutput("next_step")))),
            uiOutput("data_summary"),
            uiOutput("preview_controls")),
  
  nav_panel("Check Data", value = "explore",
            layout_sidebar(
              # bslib only blends the layout into the page while the tab pane is a fill
              # container; the panes are plain document flow now, so switch its frame
              # off explicitly instead of letting it draw a border and rounded corners.
              border = FALSE, border_radius = FALSE,
              sidebar = sidebar(width = 250, open = "always", resizable = FALSE,
                                sex_radio("ex_sex"), age_range_input("ex_age_from", "ex_age_to"),
                                conditionalPanel("input.advanced",
                                                 sliderInput("max_groups",div(class="standard-text-dark", "Maximum number of age groups"),
                                                             min = 2, max = 12, value = 3, step = 1, ticks = FALSE)),
                                task_button("run_sex", div(class="standard-text-white-bold", "Check sex differences"), class = "btn-primary w-100 mt-1"),
                                task_button("run_drift",div(class="standard-text-white-bold", "Check age drift"), class = "btn-primary w-100 mt-1")),
              layout_columns(col_widths = c(7, 5),
                             card(full_screen = TRUE, min_height = "380px", plot_canvas("explore_plot")),
                             card(card_header("Result"), min_height = "380px", uiOutput("explore_summary"),
                                  uiOutput("explore_next"))))),
  
  nav_panel("Calculate Reference Interval", value = "ri",
            layout_sidebar(
              border = FALSE, border_radius = FALSE,
              sidebar = sidebar(width = 250, open = "always", resizable = FALSE,
                                sex_radio("ri_sex"),
                                conditionalPanel("output.has_trimester && input.advanced",
                                                 radioButtons("ri_trimester", "Trimester", inline = FALSE,
                                                              choiceNames = c("any", "1.", "2.", "3."),
                                                              choiceValues = c("0", "1", "2", "3"))),
                                age_range_input("ri_age_from", "ri_age_to"),
                                radioButtons("method", "Method", choiceNames = METHOD_LABELS,
                                             choiceValues = names(RI_METHODS)),
                                conditionalPanel("input.method == 'refineR' && input.advanced",
                                                 checkboxInput("modboxcox", "Use modified Box-Cox", value = FALSE),
                                                 sliderInput("nbootstrap", "Bootstrap iterations",
                                                             min = 0, max = 50, value = 0, step = 5, ticks = FALSE),
                                                 div(class = "form-text mb-2", "0 = no confidence intervals (fastest)")),
                                conditionalPanel("input.method == 'tml' && input.advanced",
                                                 checkboxInput("fast_tml", "Fast mode (3 significant figures)", value = TRUE)),
                                checkboxInput("compare", div(class = "standard-text-dark","Compare with known limits"), value = FALSE),
                                conditionalPanel("input.compare",
                                                 div(class = "rrr-field",
                                                     layout_columns(col_widths = c(6, 6), gap = "0.4rem", fill = FALSE,
                                                                    numericInput("cmp_low", NULL, value = NA, step = 0.1),
                                                                    numericInput("cmp_high", NULL, value = NA, step = 0.1)),
                                                     div(class = "form-text mt-0", "lower / upper limit"))),
                                task_button("run_ri", "Calculate", class = "btn-primary w-100 mt-2"),
                                uiOutput("ri_export")),
              uiOutput("ri_boxes"),
              layout_columns(col_widths = c(7, 5),
                             card(full_screen = TRUE, min_height = "380px", plot_canvas("ri_plot")),
                             card(card_header("Details"), min_height = "380px", uiOutput("ri_summary"))))),
  
  nav_spacer(),
  nav_item(div(class = "rrr-umo", UMO_SVG)),
  nav_menu(title = icon("circle-question"), align = "right",
           nav_item(a("Help", href = RRR_HELP, target = "_blank")),
           nav_item(a("GitHub", href = RRR_GITHUB, target = "_blank")),
           nav_item(a("Disclaimer", href = RRR_DISCLAIMER, target = "_blank")))
)

### Server Definition Start ########################################################

server <- function(input, output, session) {
  
  raw <- reactiveVal(NULL); guessed <- reactiveVal(NULL); source_name <- reactiveVal(NULL)
  
  load_raw <- function(df, label) {
    raw(df); guessed(guess_columns(df)); source_name(label)
    showNotification(sprintf("%s loaded: %s rows, %s columns.",
                             label, big(nrow(df)), ncol(df)), duration = 4)
  }
  notify_error <- function(e, prefix = NULL)
    showNotification(paste(c(prefix, conditionMessage(e)), collapse = " "),
                     type = "error", duration = 10)
  
  observeEvent(input$file, tryCatch(
    load_raw(read_table_file(input$file$datapath, input$file$name), input$file$name),
    error = function(e) notify_error(e, "Could not read the file:")))
  
  observeEvent(input$use_paste, tryCatch(
    load_raw(read_pasted_text(input$pasted), "Pasted data"), error = notify_error))
  
  observeEvent(input$use_demo, {
    load_raw(make_demo_data(), "Demo data")
    for (id in c("ex_age_from", "ri_age_from")) updateNumericInput(session, id, value = 18)
    for (id in c("ex_age_to", "ri_age_to")) updateNumericInput(session, id, value = 45)
  })
  

  output$column_ui <- renderUI({
    df <- raw()
    if (is.null(df)) return(NULL)
    g <- guessed(); cols <- names(df)
    optional <- function(id, label, sel)
      selectInput(id, label, choices = c("(none)" = "", cols), selected = sel %||% "",
                  width = "100%")
    div(
      card_body(padding = c(5, 10), gap = "0.25rem", fillable = FALSE, class = "px-2",
        div(class = "standard-text-light","Adjust column selection if necessary"),
         layout_columns(col_widths = c(4, 4, 4), fill = FALSE,
                        selectInput("col_result", div(class = "standard-text-dark","Result"), choices = cols, selected = g$result,
                                    width = "100%"),
                        optional("col_age",div(class = "standard-text-dark", "Age"), g$age), optional("col_sex",div(class = "standard-text-dark", "Sex"), g$sex),
                        conditionalPanel("input.advanced",
                                         optional("col_trimester", div(class = "standard-text-dark","Trimester"), g$trimester))),
         uiOutput("sex_label_ui")))
  })
  
  output$welcome_ui <- renderUI({
    if (!is.null(raw())) return(NULL)
    div(class = "standard-text-light",
                "Load a file, paste data or generate demo data to continue.")
  })
    
  col_map <- reactive({
    blank <- function(x) if (is.null(x) || !nzchar(x)) NULL else x
    list(result = blank(input$col_result), age = blank(input$col_age),
         sex = blank(input$col_sex), trimester = blank(input$col_trimester))
  })
  
  NO_SEX_LABELS <- list(female = NULL, male = NULL)
  sex_labels <- reactiveVal(NO_SEX_LABELS)
  
  observe(sex_labels(list(female = input$lab_female, male = input$lab_male)))
  
  observeEvent(input$reset_sex_labels, {
    sex_labels(NO_SEX_LABELS)
    for (id in c("lab_female", "lab_male"))
      updateSelectInput(session, id, selected = character(0))
  })
  
  prepared <- reactive({
    df <- raw(); cm <- col_map()
    if (is.null(df) || is.null(cm$result) || !all(unlist(cm) %in% names(df))) return(NULL)
    tryCatch(prepare_data(df, cm, sex_labels()),
             error = function(e) { notify_error(e); NULL })
  })
  
  dat <- reactive(prepared()$data)
  
  output$has_trimester <- reactive(isTRUE(prepared()$report$has_trimester))
  outputOptions(output, "has_trimester", suspendWhenHidden = FALSE)
  
  output$sex_label_ui <- renderUI({
    r <- prepared()$report
    if (!length(r$custom_labels)) return(NULL)
    choices <- c("(none)" = "", r$custom_labels)
    # The panel is redrawn whenever an assignment changes, so carry the current
    # picks over; isolate() keeps that read from re-triggering the redraw.
    picked <- function(id) isolate(input[[id]])
    tagList(hr(),
            div(class = "standard-text-light",
                if (length(r$unknown_labels))
                  sprintf("Unrecognised sex labels: %s. Assign them below.",
                          paste(r$unknown_labels, collapse = ", "))
                else "All unrecognised sex labels are assigned."),
            layout_columns(col_widths = c(6, 6), fill = FALSE,
                           selectInput("lab_female",div(class = "standard-text-dark", "female"), choices = choices, multiple = TRUE,
                                       selected = picked("lab_female")),
                           selectInput("lab_male", div(class = "standard-text-dark","male"), choices = choices, multiple = TRUE,
                                       selected = picked("lab_male"))),
            actionButton("reset_sex_labels", "Reset labels",
                         class = "btn-outline-primary btn-sm"))
  })
  
  output$data_summary <- renderUI({
    r <- prepared()$report
    if (is.null(r)) return(NULL)
    warn <- function(cond) if (cond) "small text-warning" else "small text-secondary"
    n_discarded <- r$n_input - r$n_usable
    tiles <- list(
      tile("Total Data", big(r$n_input)),
      tile("Usable results", big(r$n_usable),
           if (n_discarded > 0) sprintf("%s discarded / <LOQ", big(n_discarded)),
           note_class = warn(n_discarded > 0 || r$n_usable < MIN_N)),
      tile("With age", big(r$n_with_age),
           if (r$n_bad_age > 0) sprintf("%s unrecognised", big(r$n_bad_age)),
           note_class = warn(r$n_bad_age > 0)),
      tile("With sex", big(r$n_with_sex),
           if (r$n_unknown_sex > 0) sprintf("%s unrecognised", big(r$n_unknown_sex)),
           note_class = warn(r$n_unknown_sex > 0)),
      if (isTRUE(input$advanced) && isTRUE(r$has_trimester))
        tile("With trimester", big(r$n_with_trimester),
             if (r$n_bad_trimester > 0) sprintf("%s unrecognised", big(r$n_bad_trimester)),
             note_class = warn(r$n_bad_trimester > 0)))
    # A dropped tile must not be handed over as NULL: layout_column_wrap would
    # still emit an empty grid cell for it, which keeps the unused column alive
    # and stops the remaining tiles from filling the row.
    tagList(
      div(class = "d-flex flex-nowrap gap-3 mb-3",
          lapply(Filter(Negate(is.null), tiles),
                 function(t) div(class = "flex-fill", style = "min-width: 0;", t))),
      if (length(r$discarded)) card(
        card_header("What happened to the discarded rows"),
        html_dl(Map(function(label, n) c(label, paste0(big(n),
                                                       if (label == "Below limit of quantification")
                                                         " <LOQ - only used by TMC")),
                    names(r$discarded), r$discarded), dt_class = "col-sm-3", dd_class = "col-sm-9")))
  })
  
  output$next_step <- renderUI({
    if (is.null(prepared()$report)) return(NULL)
    div(class = "d-flex justify-content-end gap-2 mb-3",
        actionButton("go_check", "Proceed to data check →",
                     class = "btn-outline-primary"),
        actionButton("go_calc", "Proceed to calculation →", class = "btn-primary"))
  })
  
  output$explore_next <- renderUI({
    if (is.null(explore_res())) return(NULL)
    div(class = "d-flex justify-content-end mt-3",
        actionButton("go_calc2", "Proceed to calculation →", class = "btn-primary"))
  })
  
  observeEvent(input$use_age_group, {
    r <- explore_res()
    if (is.null(r) || !identical(r$kind, "drift") || !isTRUE(r$value$ok)) return()
    groups <- r$value$groups
    i <- as.integer(input$use_age_group)
    if (is.na(i) || i < 1L || i > nrow(groups)) return()
    from <- groups$from[i]; to <- groups$to[i]
    updateNumericInput(session, "ex_age_from", value = from)
    updateNumericInput(session, "ex_age_to", value = to)
    updateNumericInput(session, "ri_age_from", value = from)
    updateNumericInput(session, "ri_age_to", value = to)
    showNotification(sprintf("Age range set to %s - %s.",
                             format(from, trim = TRUE), format(to, trim = TRUE)),
                     duration = 3)
  })
  
  observeEvent(input$go_check, nav_select("nav", "explore"))
  observeEvent(input$go_calc, nav_select("nav", "ri"))
  observeEvent(input$go_calc2, nav_select("nav", "ri"))
  
  # The preview shows either the cleaned data or the mapped columns exactly as
  # they were read, so unparseable results and unrecognised sex labels stay
  # visible. Loading a new source returns to the cleaned view.
  preview_original <- reactiveVal(FALSE)
  observeEvent(input$toggle_preview, preview_original(!preview_original()))
  observeEvent(raw(), preview_original(FALSE))
  
  output$preview_controls <- renderUI({
    if (is.null(prepared()$report)) return(NULL)
    original <- preview_original()
    card(full_screen = TRUE, min_height = "320px",
         card_header(
           class = "d-flex align-items-center justify-content-between",
           "Preview",
           actionButton("toggle_preview",
                        if (original) "View cleaned data" else "View original data",
                        class = "btn-outline-secondary btn-sm flex-shrink-0")),
         card_body(
           padding = c(6, 16, 12, 16), gap = "0.25rem",
           div(class = "form-text mt-0 mb-2", if (original)
             "Showing the first 100 rows as they were read, including discarded ones"
             else "Showing the first 100 usable entries",
           DT::DTOutput("preview"))))
  })
  
  preview_table <- reactive({
    if (!preview_original()) {
      d <- utils::head(req(dat()), 100)
      show <- data.frame(result = d$result_chr, age = d$age, sex = as.character(d$sex),
                         stringsAsFactors = FALSE)
      if ("trimester" %in% names(d)) show$trimester <- d$trimester
      return(show)
    }
    df <- req(raw())
    cols <- Filter(function(x) !is.null(x) && x %in% names(df),
                   col_map()[c("result", "age", "sex", "trimester")])
    if (!length(cols)) return(NULL)
    d <- utils::head(df[, unlist(cols), drop = FALSE], 100)
    stats::setNames(as.data.frame(lapply(d, as.character), stringsAsFactors = FALSE),
                    names(cols))
  })
  
  output$preview <- DT::renderDT({
    DT::datatable(req(preview_table()), rownames = FALSE, style = "bootstrap5",
                  class = "table-sm table-striped",
                  options = list(dom = "t", pageLength = 100, ordering = TRUE)
                  )
  }, server = TRUE)
  
  run_job <- function(what, args) {
    job <- quote(switch(what,
                        sex   = compute_sex_difference(args$df),
                        drift = compute_age_drift(args$df, args$max_groups),
                        ri    = compute_reference_interval(args$values, args$values_chr,
                                                           args$method, args$opts, args$base_path)))
    if (ASYNC) mirai::mirai({ eval(core); eval(job) },
                            core = CORE, job = job, what = what, args = args)
    else promises::promise_resolve(eval(job))
  }
  
  new_task <- function(what, button)
    bind_task_button(ExtendedTask$new(function(args) run_job(what, args)), button)
  sex_task <- new_task("sex", "run_sex")
  drift_task <- new_task("drift", "run_drift")
  ri_task <- new_task("ri", "run_ri")
  
  collect <- function(task, target, kind, settings = function() NULL)
    observeEvent(task$status(), {
      if (task$status() == "success")
        target(list(kind = kind, value = task$result(), settings = isolate(settings())))
      else if (task$status() == "error")
        notify_error(tryCatch(task$result(), error = function(e) e), "Calculation failed:")
    }, ignoreInit = TRUE)
  
  result_or <- function(r, placeholder) {
    if (is.null(r)) return(if (is.null(placeholder)) NULL else hint(placeholder))
    if (!isTRUE(r$value$ok)) return(div(class = "text-danger", r$value$message))
    NULL
  }
  
  explore_res <- reactiveVal(NULL)
  collect(sex_task, explore_res, "sex")
  collect(drift_task, explore_res, "drift")
  
  explore_subset <- reactive(subset_data(req(dat()), sex = input$ex_sex,
                                         age_from = input$ex_age_from, age_to = input$ex_age_to))
  
  sync_to_ri <- function() {
    updateRadioButtons(session, "ri_sex", selected = input$ex_sex)
    updateNumericInput(session, "ri_age_from", value = input$ex_age_from)
    updateNumericInput(session, "ri_age_to", value = input$ex_age_to)
  }
  
  observeEvent(input$run_sex, {
    d <- explore_subset(); sync_to_ri()
    sex_task$invoke(list(df = d[, c("result", "sex")]))
  })
  
  observeEvent(input$run_drift, {
    d <- explore_subset(); sync_to_ri()
    drift_task$invoke(list(df = d[, c("result", "age")], max_groups = input$max_groups %||% 3))
  })
  
  output$explore_plot <- renderPlot({
    r <- explore_res()
    validate(need(!is.null(r), "Choose a selection and start an analysis."))
    validate(need(isTRUE(r$value$ok), r$value$message))
    r$value$plot
  })
  
  output$explore_summary <- renderUI({
    r <- explore_res()
    early <- result_or(r, div(class = "standard-text-light", paste(
      "Two checks are available:",
      "whether the results differ between the sexes,",
      "and whether the median drifts with age. Both suggest whether the",
      "reference interval should be stratified.")))
    if (!is.null(early)) return(early)
    if (r$kind == "sex") sex_summary_ui(r$value) else drift_summary_ui(r$value)
  })

  sex_summary_ui <- function(v) {
    pu <- v$pu_percent
    colour <- if (!is.finite(pu)) "text-secondary" else if (v$median_diff > pu) "text-danger"
    else if (v$median_diff > 0.5 * pu) "text-warning" else "text-success"
    verdict <- if (v$n_groups < 2 || !isTRUE(is.finite(v$p_value)))
      div(class = "standard-text-dark mb-1","More than one group is needed for a statistical comparison.")
    else if (v$p_value < 0.05) tagList(
      div(class = "text-danger fw-semibold",
          sprintf("Significant difference (%s, %s)", format_p(v$p_value), v$test_used)),
      div(class = "standard-text-light mb-0",
          "With large samples this test can reach significance for differences ",
          "of no practical relevance. Compare the median deviation below."))
    else div(class = paste("standard-text-green-bold mb-1"),
             sprintf("No significant difference (%s, %s)", format_p(v$p_value), v$test_used))
    tagList(
      html_table(v$median_table),
      p(class = "standard-text-dark mb-1", sprintf("n = %s selected", big(v$n)),
        if (isTRUE(v$data_removed))
          span(br(), sprintf("Groups with fewer than %s members were removed.", MIN_GROUP_N))),
      verdict, hr(class = "my-2"),
      p(class = "standard-text-dark mb-1", "Largest deviation between medians: ",
        span(class = paste("fw-bold limit-value", colour), sprintf("%.2f%%", v$median_diff))),
      p(class = "standard-text-light mb-0",
        sprintf("Estimated permissible uncertainty: %s",
                if (is.finite(pu)) sprintf("%.2f%%", pu) else "not available")))
  }
  
  drift_summary_ui <- function(v) tagList(
    p(class = "standard-text-dark", sprintf("n = %s selected", big(v$n))),
      if (isTRUE(v$subsampled))
        p(class = "standard-text-light",  span(sprintf("The curve was fitted on a random subsample of %s.", big(v$n_used)))),
    if (!v$stratification_needed)
      div(class = "standard-text-green-bold mb-1", "No stratification by age necessary.")
    else tagList(
      div(class = "standard-text-dark-bold mb-1", "Age groups with comparable medians"),
      drift_group_table(v$groups)))
  
  ri_res <- reactiveVal(NULL); ri_pending <- reactiveVal(NULL)
  collect(ri_task, ri_res, "ri", settings = ri_pending)
  
  comparison_limits <- reactive({
    lo <- input$cmp_low; hi <- input$cmp_high
    if (!isTRUE(input$compare) || is.na(lo) || is.na(hi) || hi <= 0 || hi <= lo) NULL
    else c(lo, hi)
  })
  
  observeEvent(input$run_ri, {
    d <- req(dat())
    tri <- if (isTRUE(prepared()$report$has_trimester)) as.integer(input$ri_trimester)
    d <- subset_data(d, sex = input$ri_sex, age_from = input$ri_age_from,
                     age_to = input$ri_age_to, trimester = tri)
    method <- input$method
    n <- sum(!is.na(if (method == "tmc") d$result_chr else d$result))
    if (n >= MIN_N && n < WARN_N_SMALL)
      showNotification(sprintf("Only %s results selected - Consider increasing the sample size.",
                               big(n)), type = "warning", duration = 8)
    else if (n > WARN_N_LARGE)
      showNotification(sprintf("%s results selected - this may take a while.", big(n)),
                       type = "warning", duration = 8)
    cmp <- comparison_limits()
    opts <- list(nbootstrap = input$nbootstrap %||% 0, modboxcox = isTRUE(input$modboxcox),
                 fast_tml = isTRUE(input$fast_tml), targets = cmp)
    ri_pending(list(sex = input$ri_sex, trimester = tri, age_from = input$ri_age_from,
                    age_to = input$ri_age_to, comparison = cmp, opts = opts,
                    source = source_name()))
    ri_task$invoke(list(values = d$result, values_chr = d$result_chr, method = method,
                        base_path = BASE_PATH, opts = opts))
  })
  
  output$ri_export <- renderUI({
    r <- ri_res()
    if (is.null(r) || !isTRUE(r$value$ok)) return(NULL)
    downloadButton("export_pdf", "Export as PDF", class = "btn-outline-primary w-100 mt-2")
  })
  
  output$export_pdf <- downloadHandler(
    filename = function() sprintf("ReferenceRangeR_%s_%s.pdf",
                                  RI_METHODS[[isolate(ri_res())$value$method]],
                                  format(Sys.Date(), "%Y-%m-%d")),
    contentType = "application/pdf",
    content = function(file) {
      r <- isolate(ri_res())
      export_ri_pdf(file, r$value, r$settings, r$settings$source)
    })
  
  output$ri_boxes <- renderUI({
    r <- ri_res()
    if (is.null(r) || !isTRUE(r$value$ok)) return(NULL)
    v <- r$value
    layout_column_wrap(width = "220px", fill = FALSE, class = "mb-3",
                       tile("Lower limit", span(class = "limit-value", format_limit(v$limits[1]))),
                       tile("Upper limit", span(class = "limit-value", format_limit(v$limits[2]))),
                       tile("Method", RI_METHODS[[v$method]], sprintf("n = %s", big(v$n))))
  })
  
  output$ri_plot <- renderPlot({
    r <- ri_res()
    validate(need(!is.null(r), "Choose the settings on the left and press Calculate."))
    validate(need(isTRUE(r$value$ok), r$value$message))
    v <- r$value
    ok <- tryCatch({ grDevices::replayPlot(v$plot); TRUE }, error = function(e) FALSE)
    validate(need(ok, "The plot for this method could not be displayed."))
    cmp <- r$settings$comparison
    if (!is.null(cmp) && v$method %in% c("refineR", "tmc", "tml"))
      plot_comparison_limits(v$limits[1], v$limits[2], cmp[1], cmp[2])
  })
  
  output$ri_summary <- renderUI({
    r <- ri_res()
    early <- result_or(r, NULL)
    if (!is.null(early) || is.null(r)) return(early)
    v <- r$value; s <- r$settings
    rows <- list(c("n", big(v$n)),
                 c("Sex", switch(s$sex, A = "no selection", M = "male", F = "female", D = "non-binary")),
                 c("Age", if (age_limits_valid(s$age_from, s$age_to))
                   sprintf("%s to %s", s$age_from, s$age_to) else "no selection"))
    if (isTRUE(s$trimester %in% 1:3))
      rows <- c(rows, list(c("Trimester", paste0(s$trimester, "."))))
    rows <- c(rows, list(c("Method", RI_METHODS[[v$method]])))
    tagList(html_dl(rows), hr(class = "my-2"),
            p(class = "mb-1 standard-text-dark-bold",
              sprintf("Estimated reference interval: %s - %s",
                      format_limit(v$limits[1]), format_limit(v$limits[2]))),
            if (!is.null(v$ci)) p(class = "small mb-1",
                                  sprintf("95%% confidence intervals: lower %s - %s, upper %s - %s",
                                          format_limit(v$ci[1]), format_limit(v$ci[2]),
                                          format_limit(v$ci[3]), format_limit(v$ci[4]))),
            if (!is.null(s$comparison)) p(class = "standard-text-light mb-1",
                                          sprintf("Comparison limits: %s - %s",
                                                  format_limit(s$comparison[1]), format_limit(s$comparison[2]))),
            if (isTRUE(is.finite(v$lambda))) p(class = "standard-text-light mb-0",
                                               paste0(sprintf("Estimated \u03bb: %.2f", v$lambda),
                                                      if (nzchar(v$lambda_note %||% "")) paste0(" ", v$lambda_note))),
            if (nzchar(v$note)) p(class = "standard-text-light mb-0", v$note))
  })
  
  for (ids in list(c("ex_age_from", "ex_age_to"), c("ri_age_from", "ri_age_to"))) local({
    id <- ids
    observe(session$sendCustomMessage("rrr-mark-range",
                                      list(id = paste0(id[1], "_wrap"), invalid = isTRUE(input[[id[1]]] > input[[id[2]]]))))
  })
}

shinyApp(ui, server)

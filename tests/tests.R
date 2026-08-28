suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
})

HERE <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  dirname(if (length(f)) normalizePath(f[1], winslash = "/") else "tests/tests.R")
})
ROOT <- normalizePath(file.path(HERE, ".."), winslash = "/")
APP <- file.path(ROOT, "app.R")
FIX <- file.path(HERE, "fixtures")
BASELINE <- file.path(FIX, "baseline.csv")
SCRIPTS <- Sys.getenv("RRR_BASE_PATH", "Z:/R/kc.uol.de/referenceranger/")

PACKAGES <- c(
  shiny = "app", bslib = "app", ggplot2 = "app", htmltools = "app",
  data.table = "CSV/TXT upload", readxl = "XLSX/XLS upload", DT = "data preview",
  refineR = "refineR method", reflimR = "reflimR method + permissible uncertainty",
  tidykosmic = "kosmic method [github divinenephron/tidykosmic]",
  qgam = "age drift", mgcv = "age drift + TMC", scales = "age drift shading",
  promises = "synchronous calculation", MASS = "TML", nlme = "TML", msm = "TML",
  geoR = "TML", modeest = "TML", admisc = "TML", date = "TML",
  snpar = "TMC [github debinqiu/snpar]", tseries = "TMC", stringr = "TMC")
OPTIONAL <- c("mirai", "png")

CASES <- list(
  list(id = "demo_refiner",  file = "demo.csv", method = "refineR"),
  list(id = "demo_reflimr",  file = "demo.csv", method = "reflimr"),
  list(id = "demo_kosmic",   file = "demo.csv", method = "kosmic"),
  list(id = "demo_tmc",      file = "demo.csv", method = "tmc"),
  list(id = "demo_tml",      file = "demo.csv", method = "tml",
       opts = list(fast_tml = TRUE)),
  list(id = "demo_tml_full", file = "demo.csv", method = "tml",
       opts = list(fast_tml = FALSE)))

CASE_DEFAULTS <- list(sex = "A", age_from = NA_real_, age_to = NA_real_,
                      col_result = NULL, col_age = NULL, col_sex = NULL,
                      opts = list(), tol = 2)

field <- function(case, name) {
  v <- case[[name]]
  if (is.null(v)) CASE_DEFAULTS[[name]] else v
}

S <- new.env()
S$tier <- ""
S$tally <- list()
S$fail <- character(0)
S$notes <- character(0)
S$measured <- list()

tier <- function(name, expr) {
  S$tier <- name
  if (is.null(S$tally[[name]])) S$tally[[name]] <- c(0L, 0L, 0L)
  tryCatch(force(expr), error = function(e)
    note(2L, sprintf("%s | ABORTED: %s", name, conditionMessage(e))))
}

note <- function(slot, msg = NULL) {
  t <- S$tally[[S$tier]]
  t[slot] <- t[slot] + 1L
  S$tally[[S$tier]] <- t
  if (!is.null(msg)) S$fail <- c(S$fail, msg)
}

brief <- function(x) {
  s <- tryCatch(paste(utils::capture.output(str(x, vec.len = 4, give.attr = FALSE)),
                      collapse = " "), error = function(e) "?")
  s <- trimws(gsub("\\s+", " ", s))
  if (nchar(s) > 110) paste0(substr(s, 1, 107), "...") else s
}

ok <- function(label, cond) {
  if (isTRUE(cond)) note(1L)
  else note(2L, sprintf("%s | %s", S$tier, label))
}

eq <- function(label, actual, expected, tol = NULL) {
  good <- tryCatch(
    if (is.null(tol)) identical(actual, expected)
    else isTRUE(all.equal(actual, expected, tolerance = tol)),
    error = function(e) FALSE)
  if (good) note(1L)
  else note(2L, sprintf("%s | %s\n      got      %s\n      expected %s",
                        S$tier, label, brief(actual), brief(expected)))
}

fails <- function(label, expr, pattern) {
  m <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  if (!is.na(m) && grepl(pattern, m)) note(1L)
  else note(2L, sprintf("%s | %s\n      expected an error matching '%s', got: %s",
                        S$tier, label, pattern, if (is.na(m)) "no error" else m))
}

skip <- function(reason) note(3L, NULL)

load_core <- function(app = APP) {
  env <- new.env(parent = globalenv())
  for (e in parse(app))
    if (is.call(e) && length(e) == 3L && is.name(e[[2]]) &&
        as.character(e[[2]]) == "CORE") { eval(e, env); break }
  if (is.null(env$CORE)) stop("no CORE block in ", app)
  eval(env$CORE, env)
  for (e in parse(app))
    if (is.call(e) && length(e) == 3L && is.name(e[[2]]) &&
        is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))
      try(eval(e, env), silent = TRUE)
  env
}

E <- load_core()
SRC <- readLines(APP, warn = FALSE)
fixture <- function(...) file.path(FIX, ...)
read_semi <- function(p) utils::read.table(p, sep = ";", header = TRUE,
                                           colClasses = "character",
                                           stringsAsFactors = FALSE)
installed <- function(p) nzchar(system.file(package = p))
have <- function(p) suppressMessages(suppressPackageStartupMessages(
  requireNamespace(p, quietly = TRUE)))
remark <- function(...) S$notes <- c(S$notes, paste0(...))
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))

quiet <- function(expr) {
  depth <- sink.number()
  spare <- file(tempfile(), open = "wt")
  keep <- file(tempfile(), open = "wt")
  sink(keep)
  sink(spare)
  on.exit({
    while (sink.number() > depth) sink()
    close(spare)
    close(keep)
  }, add = TRUE)
  suppressWarnings(suppressPackageStartupMessages(suppressMessages(force(expr))))
}

run_method <- function(values, values_chr, method, opts = list()) {
  op <- options()
  on.exit(options(op), add = TRUE)
  quiet(E$compute_reference_interval(values, values_chr, method, opts, SCRIPTS))
}

record <- function(id, method, n, lower, upper)
  S$measured[[length(S$measured) + 1L]] <-
    data.frame(id = id, method = method, n = as.integer(n),
               lower = as.numeric(lower), upper = as.numeric(upper),
               stringsAsFactors = FALSE)

upload <- function(name)
  list(name = name, size = file.size(fixture(name)), type = "text/csv",
       datapath = fixture(name))

demo_prepared <- function() {
  raw <- E$read_table_file(fixture("demo.csv"))
  E$prepare_data(raw, E$guess_columns(raw))
}


tier("environment", {
  missing <- names(PACKAGES)[!vapply(names(PACKAGES), installed, logical(1))]
  eq("all required packages installed", missing, character(0))
  for (p in OPTIONAL) if (!installed(p))
    remark("optional package '", p, "' is not installed")

  named <- unlist(regmatches(SRC, gregexpr(
    '(?:requireNamespace|library|loadNamespace)\\(\\s*"[A-Za-z][A-Za-z0-9.]*"', SRC, perl = TRUE)))
  named <- gsub('.*"([A-Za-z][A-Za-z0-9.]*)"', "\\1", named)
  qualified <- unlist(regmatches(SRC, gregexpr(
    "\\b[A-Za-z][A-Za-z0-9.]*(?=:::?[A-Za-z.])", SRC, perl = TRUE)))
  bundled <- rownames(installed.packages(priority = c("base", "recommended")))
  found <- setdiff(unique(c(named, qualified, E$TML_PACKAGES, E$TMC_PACKAGES)),
                   c(bundled, "p"))
  eq("no undeclared dependency", sort(setdiff(found, c(names(PACKAGES), OPTIONAL))),
     character(0))

  ok("CORE block present", is.call(E$CORE))
  for (fn in c("compute_reference_interval", "compute_sex_difference",
               "compute_age_drift", "prepare_data", "subset_data", "read_table_file",
               "read_pasted_text", "guess_columns", "normalise_sex", "export_ri_pdf",
               "make_demo_data"))
    ok(paste(fn, "defined inside CORE"), is.function(E[[fn]]))
  body_txt <- paste(deparse(E$CORE), collapse = "\n")
  for (bad in c("input$", "output$", "session$"))
    ok(paste("CORE free of", bad), !grepl(bad, body_txt, fixed = TRUE))

  eq("MIN_N", E$MIN_N, 500L)
  eq("MIN_GROUP_N", E$MIN_GROUP_N, 100L)
  eq("SEX_LEVELS", unname(E$SEX_LEVELS), c("F", "M", "D", "X"))
  eq("methods", sort(names(E$RI_METHODS)),
     sort(c("refineR", "tmc", "tml", "kosmic", "reflimr")))
  eq("method references", sort(names(E$METHOD_REFERENCES)), sort(names(E$RI_METHODS)))

  need <- unique(c("demo.csv", "demo_dot.xlsx", "demo_comma.xlsx", "bad_data.csv",
                   vapply(CASES, function(c) c$file, character(1))))
  eq("fixtures present", need[!file.exists(fixture(need))], character(0))
})


tier("ingest", {
  eq("as_number decimal comma", E$as_number(c("1,5", "2.5", "x")), c(1.5, 2.5, NA_real_))
  eq("as_number blanks", E$as_number(c("", " ")), c(NA_real_, NA_real_))
  eq("as_number negative comma", E$as_number("-0,25"), -0.25)

  ok("age window valid", E$age_limits_valid(18, 45))
  ok("age window inverted", !E$age_limits_valid(45, 18))
  ok("age window empty", !E$age_limits_valid(18, 18))
  ok("age window NA", !E$age_limits_valid(NA, 45))
  ok("age window negative", !E$age_limits_valid(-5, -1))

  for (sep in c("\t", ";", ",")) {
    d <- E$read_pasted_text(paste0("result", sep, "age\n1.5", sep, "20\n2.5", sep, "30"))
    eq(paste("paste separator", encodeString(sep)), names(d), c("result", "age"))
    eq(paste("paste rows", encodeString(sep)), nrow(d), 2L)
  }
  d <- E$read_pasted_text("result   age\n1.5   20\n2.5   30")
  eq("paste multi-space", names(d), c("result", "age"))

  d <- E$read_pasted_text("result;age;sex\n1.5;20;F\n2.5;30")
  eq("ragged row padded", ncol(d), 3L)
  ok("ragged cell is NA", is.na(d[2, 3]))
  eq("blank header named", names(E$read_pasted_text("result;;sex\n1.5;20;F")),
     c("result", "V2", "sex"))

  fails("paste needs two lines", E$read_pasted_text("result;age"),
        "header line and at least one data row")
  fails("paste needs a separator", E$read_pasted_text("resultage\n1.5"),
        "No column separator found")

  p <- tempfile(fileext = ".csv")
  writeLines(c("result,age,sex", "1.5,20,F", "2.5,30,M"), p)
  d <- E$read_table_file(p)
  eq("csv columns", names(d), c("result", "age", "sex"))
  ok("csv read as character", all(vapply(d, is.character, logical(1))))

  p2 <- tempfile()
  writeLines(c("result;age", "1,5;20", "2,5;30"), p2)
  eq("extension taken from filename", nrow(E$read_table_file(p2, "uploaded.csv")), 2L)

  p3 <- tempfile(fileext = ".csv")
  writeLines("result,age", p3)
  fails("empty file rejected", E$read_table_file(p3), "no data rows")

  g <- E$guess_columns(data.frame(result = "1", age = "2", sex = "F"))
  eq("guess english", c(g$result, g$age, g$sex), c("result", "age", "sex"))
  g <- E$guess_columns(data.frame(Wert = "1", Alter = "2", Geschlecht = "F"))
  eq("guess german", c(g$result, g$age, g$sex), c("Wert", "Alter", "Geschlecht"))
  eq("guess trimester", E$guess_columns(data.frame(Trimenon = "1", Wert = "2"))$trimester,
     "Trimenon")
  eq("exact header beats substring",
     E$guess_columns(data.frame(result_flag = "x", result = "1"))$result, "result")
  eq("numeric fallback",
     E$guess_columns(data.frame(label = c("a", "b", "c"),
                                measurement = c("1.5", "2.5", "3.5")))$result, "measurement")
  ok("no guess when nothing numeric",
     is.null(E$guess_columns(data.frame(label = c("a", "b"), note = c("x", "y")))$result))

  eq("sex levels", levels(E$normalise_sex("F")), c("F", "M", "D", "X"))
  ok("female labels", all(as.character(E$normalise_sex(E$FEMALE_LABELS)) == "F"))
  ok("male labels", all(as.character(E$normalise_sex(E$MALE_LABELS)) == "M"))
  ok("diverse labels", all(as.character(E$normalise_sex(E$DIVERSE_LABELS)) == "D"))
  eq("unknown sex", as.character(E$normalise_sex("Nonsense")), "X")
  ok("blank sex is NA", is.na(E$normalise_sex("")))
  eq("custom female label", as.character(E$normalise_sex("Dame", female = "Dame")), "F")
  eq("sex trimmed", as.character(E$normalise_sex(" F ")), "F")

  raw <- data.frame(result = c("1,5", "", "abc", "0", "-3", "<0.10", "2.5"),
                    age = as.character(20:26),
                    sex = c("F", "M", "F", "M", "F", "M", "F"), stringsAsFactors = FALSE)
  p <- E$prepare_data(raw, list(result = "result", age = "age", sex = "sex"))
  eq("discard empty", unname(p$report$discarded[["Empty"]]), 1L)
  eq("discard junk", unname(p$report$discarded[["Not a number"]]), 1L)
  eq("discard non-positive", unname(p$report$discarded[["Zero or negative"]]), 2L)
  eq("discard censored", unname(p$report$discarded[["Below limit of quantification"]]), 1L)
  eq("usable count", p$report$n_usable, 2L)

  p <- E$prepare_data(data.frame(result = c("1.5", "<0.10", "2.5"),
                                 stringsAsFactors = FALSE), list(result = "result"))
  eq("censored kept as text", p$data$result_chr, c("1.5", "<0.10", "2.5"))
  ok("censored has no number", is.na(p$data$result[2]))

  p <- E$prepare_data(data.frame(result = c("1,5", "2,5"), sex = c("weiblich", "Zwitter"),
                                 stringsAsFactors = FALSE),
                      list(result = "result", sex = "sex"))
  eq("comma normalised", p$data$result, c(1.5, 2.5))
  eq("unknown label reported", p$report$unknown_labels, "Zwitter")

  p <- E$prepare_data(data.frame(result = c("1,5", "2,5"), sex = c("weiblich", "Zwitter"),
                                 stringsAsFactors = FALSE),
                      list(result = "result", sex = "sex"), list(diverse = "Zwitter"))
  eq("assigned label applied", as.character(p$data$sex), c("F", "D"))
  eq("assigned label no longer unknown", p$report$unknown_labels, character(0))
  eq("assigned label stays selectable", p$report$custom_labels, "Zwitter")

  fails("result column required", E$prepare_data(data.frame(a = "1"), list()),
        "which column holds the measured result")
  fails("result column must exist",
        E$prepare_data(data.frame(a = "1"), list(result = "nope")),
        "which column holds the measured result")

  raw <- data.frame(result = as.character(1:4), trimester = c("1", "2", "3", "9"),
                    stringsAsFactors = FALSE)
  p <- E$prepare_data(raw, list(result = "result", trimester = "trimester"))
  ok("trimester detected", p$report$has_trimester)
  eq("trimester counted", p$report$n_with_trimester, 3L)
  eq("unrecognised trimester counted", p$report$n_bad_trimester, 1L)

  pt <- E$prepare_data(data.frame(result = as.character(1:4), trimester = c("1", "", "", ""),
                                  stringsAsFactors = FALSE),
                       list(result = "result", trimester = "trimester"))
  eq("blank trimester is not unrecognised", pt$report$n_bad_trimester, 0L)
  eq("no trimester column means no count",
     unlist(E$prepare_data(data.frame(result = "1", stringsAsFactors = FALSE),
                           list(result = "result"))$report[c("n_with_trimester",
                                                             "n_bad_trimester")]),
     c(n_with_trimester = 0L, n_bad_trimester = 0L))
  eq("trimester out of range dropped", p$data$trimester, c(1L, 2L, 3L, NA_integer_))
  ok("no trimester column when unmapped",
     !E$prepare_data(raw, list(result = "result"))$report$has_trimester)

  df <- data.frame(result = 1:6, age = c(10, 20, 30, 40, 50, 60),
                   sex = E$normalise_sex(c("F", "M", "F", "M", "F", "M")),
                   trimester = c(1L, 1L, 2L, 2L, 3L, 3L))
  eq("subset all", nrow(E$subset_data(df)), 6L)
  eq("subset female", nrow(E$subset_data(df, sex = "F")), 3L)
  eq("subset absent group", nrow(E$subset_data(df, sex = "D")), 0L)
  eq("subset age window", nrow(E$subset_data(df, age_from = 20, age_to = 40)), 3L)
  eq("subset trimester", nrow(E$subset_data(df, trimester = 2)), 2L)
  eq("subset unusable window ignored",
     nrow(E$subset_data(df, age_from = 40, age_to = 20)), 6L)
})


tier("fixtures", {
  on_disk <- E$read_table_file(fixture("demo.csv"))
  fresh <- E$make_demo_data(3000L, seed = 42L)
  eq("demo.csv rows", nrow(on_disk), 3000L)
  eq("demo.csv is make_demo_data(3000, seed = 42)", on_disk$result,
     as.character(fresh$result))
  eq("demo.csv ages", as.numeric(on_disk$age), fresh$age, tol = 1e-9)
  eq("demo.csv sexes", on_disk$sex, fresh$sex)

  dot <- E$read_table_file(fixture("demo_dot.xlsx"), "demo_dot.xlsx")
  comma <- E$read_table_file(fixture("demo_comma.xlsx"), "demo_comma.xlsx")
  eq("xlsx columns", names(dot), c("result", "age", "sex"))
  eq("xlsx rows", nrow(dot), 9L)
  eq("german and english xlsx are identical", dot, comma)
  ok("xlsx holds no commas", !any(grepl(",", unlist(dot), fixed = TRUE)))

  pd <- E$prepare_data(dot, E$guess_columns(dot))
  pc <- E$prepare_data(comma, E$guess_columns(comma))
  eq("xlsx prepared identically", pd$data, pc$data)
  eq("xlsx usable", pd$report$n_usable, 9L)
  eq("xlsx first result", pd$data$result[1], 23.96)
  eq("xlsx matches demo.csv numerically", as.numeric(dot$result),
     as.numeric(head(on_disk$result, 9L)))
  eq("excel float expansion is exact as a number", pd$data$result[7], 18.65)
  eq("excel float expansion survives in text", pd$data$result_chr[7], "18.649999999999999")

  as_text <- data.frame(result = c("23,96", "18,91"), age = c("90,65", "92,83"),
                        sex = c("F", "M"), stringsAsFactors = FALSE)
  pt <- E$prepare_data(as_text, list(result = "result", age = "age", sex = "sex"))
  eq("comma in text cells", pt$data$result, c(23.96, 18.91))

  raw <- E$read_table_file(fixture("bad_data.csv"), "bad_data.csv")
  eq("bad_data rows", nrow(raw), 16L)
  eq("bad_data header, BOM stripped", names(raw), c("result", "age", "sex"))
  p <- E$prepare_data(raw, E$guess_columns(raw))
  d <- p$report$discarded
  eq("bad_data below limit of quantification", unname(d[["Below limit of quantification"]]), 3L)
  eq("bad_data empty", unname(d[["Empty"]]), 1L)
  eq("bad_data not a number", unname(d[["Not a number"]]), 2L)
  ok("bad_data has no non-positive", !("Zero or negative" %in% names(d)))
  eq("bad_data usable", p$report$n_usable, 10L)
  eq("bad_data kept", nrow(p$data), 13L)
  eq("bad_data censored kept", sum(startsWith(p$data$result_chr, "<")), 3L)
  ok("greater-than is not censored", !any(startsWith(p$data$result_chr, ">")))
  eq("bad_data unparseable ages", sum(is.na(p$data$age)), 2L)
  eq("bad_data unknown sex", p$report$n_unknown_sex, 1L)
  eq("bad_data unknown label", p$report$unknown_labels, "P")
  eq("bad_data recognised sexes", p$report$n_with_sex, 12L)
  ok("umlaut label recognised", "M" %in% as.character(p$data$sex))

  r <- E$compute_reference_interval(p$data$result, p$data$result_chr, "refineR",
                                    base_path = SCRIPTS)
  ok("bad_data too small for a method", !isTRUE(r$ok))
  eq("bad_data method n", r$n, 10L)
})


tier("guards", {
  small <- seq(1, 2, length.out = 10)
  for (m in names(E$RI_METHODS)) {
    r <- E$compute_reference_interval(small, as.character(small), m, base_path = SCRIPTS)
    ok(paste(m, "refuses small n"), !isTRUE(r$ok))
    eq(paste(m, "small n message"), r$message,
       "The selection contains 10 usable results. At least 500 are needed.")
  }
  fails("unknown method errors",
        E$compute_reference_interval(rnorm(600, 100, 10), method = "nonsense"),
        "Unknown method")

  if (have("tidykosmic")) {
    r <- E$compute_reference_interval(rep(seq(80L, 120L), length.out = 600), method = "kosmic",
                                      base_path = SCRIPTS)
    ok("kosmic refuses whole numbers", !isTRUE(r$ok) &&
         grepl("at least one decimal place", r$message))
  } else skip("tidykosmic")

  df <- data.frame(result = c(rnorm(600, 100, 10), rnorm(50, 130, 10)),
                   sex = E$normalise_sex(c(rep("F", 600), rep("D", 50))))
  r <- E$compute_sex_difference(df)
  eq("tiny sex group dropped", r$n, 600L)
  ok("drop is reported", r$data_removed)
  eq("one group left", r$n_groups, 1L)

  df <- data.frame(result = c(rnorm(450, 100, 10), rnorm(60, 130, 10)),
                   sex = E$normalise_sex(c(rep("F", 450), rep("D", 60))))
  r <- E$compute_sex_difference(df)
  ok("below MIN_N after dropping", !isTRUE(r$ok))
  ok("drop explained", grepl("dropped because their sex group held fewer than 100", r$message))

  set.seed(11)
  r <- E$compute_sex_difference(data.frame(
    result = rnorm(600, 100, 10), sex = E$normalise_sex(rep(c("F", "M"), each = 300))))
  eq("two groups use wilcoxon", r$test_used, "Wilcoxon rank-sum test")
  r <- E$compute_sex_difference(data.frame(
    result = rnorm(900, 100, 10), sex = E$normalise_sex(rep(c("F", "M", "D"), each = 300))))
  eq("three groups use kruskal", r$test_used, "Kruskal-Wallis test")

  r <- E$compute_age_drift(data.frame(result = rnorm(100, 100, 10), age = runif(100, 20, 60)))
  ok("drift refuses small n", !isTRUE(r$ok))

  tbl <- as.character(E$drift_group_table(
    data.frame(from = c(1, 10), to = c(10, 20), median = c(5, 15), size = c(50, 50))))
  ok("group table drops the size column", !grepl("size", tbl, fixed = TRUE))
  for (h in c("from", "to", "median"))
    ok(paste("group table keeps", h), grepl(paste0("<th>", h, "</th>"), tbl, fixed = TRUE))
  eq("group table offers one use button per row",
     length(gregexpr(">use</button>", tbl, fixed = TRUE)[[1]]), 2L)
  ok("use button reports its row",
     grepl("use_age_group&#39;, 2", tbl, fixed = TRUE))

  lo <- 30; hi <- 50
  eq("permissible uncertainty formula",
     E$estimate_pu_percent(NULL, list(limits = c(lo, hi))),
     2.39 * (-0.25 + 100 * (-1 + exp(((log(hi) - log(lo)) / 3.92)^2))^0.5)^0.5)
  ok("pu rejects zero limit", is.na(E$estimate_pu_percent(NULL, list(limits = c(0, 50)))))
  ok("pu rejects NA limit", is.na(E$estimate_pu_percent(NULL, list(limits = c(30, NA)))))

  q <- quantile(1:100, probs = c(0.1, 0.9))
  eq("plot ylim", E$plot_ylim(1:100),
     unname(c(q[1] - (q[2] - q[1]) / 1.3, q[2] + (q[2] - q[1]) / 1.3)))

  eq("format large", E$format_report_limit(1234.6), "1,235")
  eq("format banker rounding", E$format_report_limit(1234.5), "1,234")
  eq("format medium", E$format_report_limit(45.67), "45.7")
  eq("format small", E$format_report_limit(1.234), "1.23")
  eq("format missing", E$format_report_limit(NA), "-")
  eq("format infinite", E$format_report_limit(Inf), "-")
  for (x in c(0, 1.234, 9.99, 45.67, 100.4, 1234.5, -55.5))
    eq(paste("ui formatter agrees at", x), E$format_limit(x), E$format_report_limit(x))
  eq("thousands separator", E$big(1234567), "1,234,567")
  eq("p very small", E$format_p(0.0001), "p < 0.001")
  eq("p ordinary", E$format_p(0.0432), "p = 0.043")
  eq("p missing", E$format_p(NA), "n/a")
  for (m in names(E$RI_METHODS))
    ok(paste("settings text for", m), is.character(E$method_settings_text(m, list())))

  result <- list(ok = TRUE, method = "reflimr", n = 600L, limits = c(80.2, 119.8),
                 ci = NULL, lambda = NA_real_, note = "test", plot = NULL)
  stamp <- as.POSIXct("2026-01-01 12:00:00", tz = "UTC")
  a <- tempfile(fileext = ".pdf"); b <- tempfile(fileext = ".pdf")
  E$export_ri_pdf(a, result, generated = stamp)
  E$export_ri_pdf(b, result, generated = stamp)
  ok("pdf written", file.exists(a) && file.size(a) > 1000)
  eq("pdf magic bytes", readBin(a, "raw", 4L), charToRaw("%PDF"))
  ra <- readBin(a, "raw", file.size(a))
  rb <- readBin(b, "raw", file.size(b))
  eq("pdf size stable for fixed inputs", length(ra), length(rb))
  ok("pdf differs only in its embedded creation date",
     length(ra) == length(rb) && sum(ra != rb) <= 64L)
  unlink(c(a, b))
})


tier("app", suppressPackageStartupMessages({
  testServer(ROOT, {
    session$setInputs(file = upload("demo.csv"))
    eq("upload rows", nrow(raw()), 3000L)
    eq("upload name", source_name(), "demo.csv")
    g <- guessed()
    eq("upload guessed columns", c(g$result, g$age, g$sex), c("result", "age", "sex"))
    ok("prepared is NULL before a column is chosen", is.null(prepared()))
    session$setInputs(col_result = g$result, col_age = g$age, col_sex = g$sex,
                      col_trimester = "")
    eq("column map", col_map()$result, "result")
    eq("prepared input rows", prepared()$report$n_input, 3000L)
    eq("prepared usable", prepared()$report$n_usable, 3000L)
    eq("prepared with sex", prepared()$report$n_with_sex, 3000L)
    ok("no trimester", !prepared()$report$has_trimester)
    eq("data rows", nrow(dat()), 3000L)
    ok("summary renders", !is.null(output$data_summary))
    ok("preview renders", !is.null(output$preview))
  })

  testServer(ROOT, {
    session$setInputs(pasted = "result;age;sex\n1,5;20;F\n2,5;30;M\n3,5;40;F")
    session$setInputs(use_paste = 1)
    eq("paste rows", nrow(raw()), 3L)
    eq("paste label", source_name(), "Pasted data")
    session$setInputs(col_result = "result", col_age = "age", col_sex = "sex",
                      col_trimester = "")
    eq("paste values", prepared()$data$result, c(1.5, 2.5, 3.5))
    eq("paste sexes", as.character(prepared()$data$sex), c("F", "M", "F"))
  })

  testServer(ROOT, {
    session$setInputs(pasted = "no-separator-here")
    session$setInputs(use_paste = 1)
    ok("bad paste leaves data empty", is.null(raw()))
  })

  testServer(ROOT, {
    session$setInputs(use_demo = 1)
    eq("demo rows", nrow(raw()), 50000L)
    eq("demo label", source_name(), "Demo data")
    eq("demo columns", names(raw()), c("result", "age", "sex"))
    session$setInputs(col_result = "result", col_age = "age", col_sex = "sex",
                      col_trimester = "")
    ok("demo mostly usable", prepared()$report$n_usable > 40000L)
  })

  testServer(ROOT, {
    session$setInputs(file = upload("bad_data.csv"))
    session$setInputs(col_result = "result", col_age = "age", col_sex = "sex",
                      col_trimester = "")
    d <- prepared()$report$discarded
    eq("upload bad_data censored", unname(d[["Below limit of quantification"]]), 3L)
    eq("upload bad_data empty", unname(d[["Empty"]]), 1L)
    eq("upload bad_data junk", unname(d[["Not a number"]]), 2L)
    eq("upload bad_data usable", prepared()$report$n_usable, 10L)
    eq("upload bad_data unknown label", prepared()$report$unknown_labels, "P")
    ok("messy summary renders", !is.null(output$data_summary))
    ok("messy preview renders", !is.null(output$preview))
  })

  for (f in c("demo_dot.xlsx", "demo_comma.xlsx"))
    testServer(ROOT, {
      session$setInputs(file = upload(f))
      eq(paste("upload", f, "rows"), nrow(raw()), 9L)
      session$setInputs(col_result = "result", col_age = "age", col_sex = "sex",
                        col_trimester = "")
      eq(paste("upload", f, "usable"), prepared()$report$n_usable, 9L)
      eq(paste("upload", f, "first result"), prepared()$data$result[1], 23.96)
      eq(paste("upload", f, "first age"), prepared()$data$age[1], 90.65)
    })

  testServer(ROOT, {
    session$setInputs(file = upload("demo.csv"))
    session$setInputs(col_result = "age", col_age = "", col_sex = "sex",
                      col_trimester = "")
    eq("manual override wins", col_map()$result, "age")
    ok("unmapped age is empty", all(is.na(prepared()$data$age)))
    eq("override reads the chosen column", prepared()$data$result[1], 90.65)
  })

  testServer(ROOT, {
    session$setInputs(pasted = "result;sex\n1.5;Dame\n2.5;Herr\n3.5;Dame")
    session$setInputs(use_paste = 1)
    session$setInputs(col_result = "result", col_age = "", col_sex = "sex",
                      col_trimester = "")
    eq("labels unknown at first", prepared()$report$unknown_labels, c("Dame", "Herr"))
    session$setInputs(lab_female = "Dame", lab_male = "Herr")
    eq("custom labels applied", as.character(prepared()$data$sex), c("F", "M", "F"))
    eq("assigned labels stay listed", prepared()$report$custom_labels, c("Dame", "Herr"))
    ok("selector stays visible once assigned", !is.null(output$sex_label_ui))
    session$setInputs(reset_sex_labels = 1)
    eq("reset clears the assignments", sex_labels(), list(female = NULL, male = NULL))
    eq("reset makes the labels unknown again", prepared()$report$unknown_labels,
       c("Dame", "Herr"))
    eq("reset undoes the coding", as.character(prepared()$data$sex), c("X", "X", "X"))
  })

  testServer(ROOT, {
    session$setInputs(pasted = "result;age;sex;trimester
5.4;42;F;1
6.1;37;M;2
7.2;51;F;3")
    session$setInputs(use_paste = 1)
    session$setInputs(col_result = "result", col_age = "age", col_sex = "sex",
                      col_trimester = "trimester")
    # a dropped tile must leave no empty grid cell behind, otherwise the unused
    # column stays alive and the remaining tiles stop short of the full width
	count_matches <- function(pattern, x) {
	  m <- gregexpr(pattern, x, fixed = TRUE)[[1]]
	  if (identical(m[1], -1L)) 0L else length(m)
	}
	cells <- function() count_matches("rrr-tile", output$data_summary$html)
	session$setInputs(advanced = FALSE)
	eq("four tiles when the trimester one is dropped", cells(), 4L)
	session$setInputs(advanced = TRUE)
    eq("five tiles in advanced mode", cells(), 5L)
  })

  testServer(ROOT, {
    session$setInputs(pasted = "result;age;sex
5.4;42;F
abc;37;M
;51;weiblich")
    session$setInputs(use_paste = 1)
    session$setInputs(col_result = "result", col_age = "age", col_sex = "sex",
                      col_trimester = "")
    eq("preview starts on the cleaned data", preview_original(), FALSE)
    eq("cleaned preview keeps only usable rows", nrow(preview_table()), 1L)
    eq("cleaned preview columns", names(preview_table()), c("result", "age", "sex"))
    session$setInputs(toggle_preview = 1)
    ok("toggle switches to the original data", preview_original())
    eq("original preview keeps every row", nrow(preview_table()), 3L)
    ok("original preview keeps unparsed text", "abc" %in% preview_table()$result)
    ok("original preview keeps the raw sex label", "weiblich" %in% preview_table()$sex)
    session$setInputs(toggle_preview = 2)
    ok("toggle switches back to the cleaned data", !preview_original())
    session$setInputs(use_demo = 1)
    ok("a new source returns to the cleaned view", !preview_original())
  })

  testServer(ROOT, {
    session$setInputs(file = upload("demo.csv"))
    session$setInputs(col_result = "result", col_age = "age", col_sex = "sex",
                      col_trimester = "")
    session$setInputs(ex_sex = "F", ex_age_from = 18, ex_age_to = 45)
    s <- explore_subset()
    ok("subset is female only", all(as.character(s$sex) == "F"))
    ok("subset within age window", all(s$age >= 18 & s$age <= 45))
    ok("subset is smaller", nrow(s) > 0L && nrow(s) < 3000L)
    session$setInputs(ex_sex = "A", ex_age_from = NA, ex_age_to = NA)
    eq("no filter selects everything", nrow(explore_subset()), 3000L)
  })

  testServer(ROOT, {
    session$setInputs(compare = FALSE, cmp_low = 30, cmp_high = 50)
    ok("comparison off", is.null(comparison_limits()))
    session$setInputs(compare = TRUE)
    eq("comparison on", comparison_limits(), c(30, 50))
    session$setInputs(cmp_low = 50, cmp_high = 30)
    ok("comparison inverted", is.null(comparison_limits()))
    session$setInputs(cmp_low = 30, cmp_high = NA)
    ok("comparison incomplete", is.null(comparison_limits()))
  })

  testServer(ROOT, {
    session$setInputs(file = upload("demo.csv"))
    session$setInputs(col_result = "result", col_age = "age", col_sex = "sex",
                      col_trimester = "")
    session$setInputs(ri_sex = "A", ri_age_from = NA, ri_age_to = NA, method = "reflimr",
                      compare = FALSE, nbootstrap = 0, modboxcox = FALSE, fast_tml = FALSE)
    session$setInputs(run_ri = 1)
    for (i in 1:200) {
      later::run_now(0.05)
      session$flushReact()
      if (identical(ri_task$status(), "error")) break
      if (!identical(ri_task$status(), "running") && !is.null(ri_res())) break
    }
    eq("task succeeded", ri_task$status(), "success")
    r <- ri_res()
    eq("task result kind", r$kind, "ri")
    ok("task result ok", isTRUE(r$value$ok))
    eq("task method", r$value$method, "reflimr")
    eq("task n", r$value$n, 3000L)
    ok("task limits ordered", all(is.finite(r$value$limits)) &&
         r$value$limits[1] < r$value$limits[2])
    eq("task settings sex", r$settings$sex, "A")
    eq("task settings source", r$settings$source, "demo.csv")
  })
}))


tier("stratification", {
  p <- demo_prepared()
  r <- E$compute_sex_difference(p$data)
  ok("sex difference ok", isTRUE(r$ok))
  eq("sex n", r$n, 3000L)
  eq("sex groups", r$n_groups, 2L)
  eq("sex test", r$test_used, "Wilcoxon rank-sum test")
  eq("sex group sizes sum", sum(r$median_table$n), 3000L)
  ok("sex plot built", inherits(r$plot, "ggplot"))
  record("demo_sex", "sex", r$n, r$p_value, r$median_diff)

  if (have("qgam")) {
    set.seed(1)
    r <- quiet(E$compute_age_drift(p$data))
    ok("drift ok", isTRUE(r$ok))
    eq("drift n", r$n, 3000L)
    ok("drift not subsampled", !isTRUE(r$subsampled))
    ok("drift suggests stratification", r$stratification_needed)
    eq("drift columns", names(r$groups), c("from", "to", "median", "size"))
    ok("drift groups ascend", all(diff(r$groups$from) > 0))
    ok("drift groups are intervals", all(r$groups$to > r$groups$from))
    ok("drift covers the data", abs(sum(r$groups$size) - 100) <= 2)
    ok("drift plot built", inherits(r$plot, "ggplot"))
    record("demo_drift", "drift", nrow(r$groups), r$groups$to[1],
           r$groups$from[nrow(r$groups)])

    set.seed(99); a <- quiet(E$compute_age_drift(p$data, max_points = 1000L))
    set.seed(99); b <- quiet(E$compute_age_drift(p$data, max_points = 1000L))
    ok("drift subsamples above max_points", isTRUE(a$subsampled) && a$n_used == 1000L)
    eq("seeded subsample is repeatable", a$groups, b$groups)
  } else skip("qgam")
})


tier("methods", {
  ids <- vapply(CASES, function(c) c$id, character(1))
  methods <- vapply(CASES, function(c) c$method, character(1))
  ok("case ids unique", !anyDuplicated(ids))
  eq("case methods known", sort(setdiff(methods, names(E$RI_METHODS))), character(0))

  for (case in CASES[order(methods %in% c("tmc", "tml"))]) {
    id <- case$id
    method <- case$method

    if (method == "kosmic" && !have("tidykosmic")) { skip(id); next }
    files <- switch(method, tmc = c("TMC.settings.R", "TMC.functions.R"), tml = "TML.R", NULL)
    if (!is.null(files) && !all(file.exists(file.path(SCRIPTS, files)))) {
      remark("skipped ", id, ": ", paste(files, collapse = " + "), " not in ", SCRIPTS)
      skip(id)
      next
    }

    path <- fixture(case$file)
    if (!file.exists(path)) { skip(id); next }
    raw <- E$read_table_file(path, path)
    cols <- E$guess_columns(raw)
    for (nm in c("result", "age", "sex")) {
      o <- field(case, paste0("col_", nm))
      if (!is.null(o)) cols[[nm]] <- o
    }
    prep <- E$prepare_data(raw, cols)
    d <- E$subset_data(prep$data, sex = field(case, "sex"),
                       age_from = field(case, "age_from"),
                       age_to = field(case, "age_to"))

    res <- run_method(d$result, d$result_chr, method, field(case, "opts"))
    if (!isTRUE(res$ok)) {
      note(2L, sprintf("methods | %s (%s) produced no result: %s", id, method,
                       if (is.null(res$message)) "unknown" else res$message))
      next
    }
    note(1L)
    ok(paste(id, "limits finite and ordered"),
       all(is.finite(res$limits)) && res$limits[1] < res$limits[2])
    record(id, method, res$n, res$limits[1], res$limits[2])
  }
})


measured <- if (length(S$measured)) do.call(rbind, S$measured) else NULL
recorded_now <- FALSE

if (!is.null(measured)) {
  if (!file.exists(BASELINE)) {
    utils::write.table(measured, BASELINE, sep = ";", row.names = FALSE, quote = FALSE)
    recorded_now <- TRUE
  } else {
    base <- read_semi(BASELINE)
    tol <- function(id) {
      hit <- Filter(function(c) c$id == id, CASES)
      if (!length(hit)) 1e-6 else field(hit[[1]], "tol")
    }
    S$tier <- "baseline"
    S$tally[["baseline"]] <- c(0L, 0L, 0L)
    for (i in seq_len(nrow(measured))) {
      r <- measured[i, ]
      b <- base[base$id == r$id & base$method == r$method, , drop = FALSE]
      if (nrow(b) != 1L) {
        note(2L, sprintf("baseline | %s (%s) is not in baseline.csv", r$id, r$method))
        next
      }
      eq(paste(r$id, "n"), r$n, as.integer(b$n[1]))
      for (side in c("lower", "upper")) {
        old <- num(b[[side]][1]); new <- r[[side]]
        delta <- if (isTRUE(is.finite(old)) && old != 0) 100 * (new - old) / abs(old) else NA_real_
        if (isTRUE(abs(delta) <= tol(r$id))) note(1L)
        else note(2L, sprintf("baseline | %s %s moved %+.3f%% (%.6g to %.6g), tolerance %g%%",
                              r$id, side, delta, old, new, tol(r$id)))
      }
    }
  }
}

RULE <- strrep("-", 76)
cat("\n", RULE, "\n", sep = "")
cat(sprintf("ReferenceRangeR tests    R %s.%s    %s\n",
            R.version$major, R.version$minor, APP))
cat(RULE, "\n\n", sep = "")

total_pass <- 0L; total_fail <- 0L; total_skip <- 0L
for (nm in names(S$tally)) {
  t <- S$tally[[nm]]
  total_pass <- total_pass + t[1]; total_fail <- total_fail + t[2]
  total_skip <- total_skip + t[3]
  cat(sprintf("  %-16s %s  %4d ok%s\n", nm, if (t[2] == 0L) "PASS" else "FAIL", t[1],
              if (t[3] > 0L) sprintf(", %d skipped", t[3]) else ""))
}

if (!is.null(measured)) {
  base <- if (file.exists(BASELINE) && !recorded_now) read_semi(BASELINE) else NULL
  cat("\n")
  cat(sprintf("  %-14s %-8s %6s  %-21s %-21s\n", "id", "method", "n", "lower", "upper"))
  for (i in seq_len(nrow(measured))) {
    r <- measured[i, ]
    b <- if (is.null(base)) NULL else
      base[base$id == r$id & base$method == r$method, , drop = FALSE]
    show <- function(side) {
      v <- formatC(r[[side]], format = "g", digits = 6)
      if (is.null(b) || nrow(b) != 1L) v
      else sprintf("%s (%s)", v, formatC(num(b[[side]][1]), format = "g", digits = 6))
    }
    cat(sprintf("  %-14s %-8s %6d  %-21s %-21s\n", r$id, r$method, r$n,
                show("lower"), show("upper")))
  }
}

if (length(S$notes)) {
  cat("\n")
  for (n in S$notes) cat("  note: ", n, "\n", sep = "")
}

cat("\n", RULE, "\n", sep = "")
if (recorded_now)
  cat(sprintf("  baseline recorded: %s\n", BASELINE))
if (total_fail == 0L) {
  cat(sprintf("  ALL PASS    %d checks%s\n", total_pass,
              if (total_skip > 0L) sprintf(", %d skipped", total_skip) else ""))
} else {
  cat(sprintf("  %d FAILED of %d checks\n\n", total_fail, total_pass + total_fail))
  for (f in S$fail) cat("   ", f, "\n", sep = "")
  cat("\n  to accept new reference limits, delete", BASELINE, "and run again\n")
}
cat(RULE, "\n\n", sep = "")
quit(status = if (total_fail == 0L) 0L else 1L, save = "no")

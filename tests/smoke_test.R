#!/usr/bin/env Rscript

APP <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  here <- if (length(f)) dirname(normalizePath(f[1], winslash = "/")) else getwd()
  for (p in c(file.path(here, "app.R"), file.path(getwd(), "app.R"),
              "/srv/shiny-server/referenceranger/app.R"))
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  stop("app.R not found next to this script, in the working directory, or in /srv/shiny-server/referenceranger", call. = FALSE)
})
APPDIR <- dirname(APP)

PASS <- 0L
FAIL <- character(0)
NOTE <- character(0)

check <- function(label, cond, detail = "") {
  good <- isTRUE(tryCatch(cond, error = function(e) { detail <<- conditionMessage(e); FALSE }))
  cat(sprintf("  [%s] %s%s\n", if (good) "ok  " else "FAIL", label,
              if (!good && nzchar(detail)) paste0("  --  ", detail) else ""))
  if (good) PASS <<- PASS + 1L else FAIL <<- c(FAIL, label)
  invisible(good)
}

remark <- function(...) {
  msg <- paste0(...)
  cat(sprintf("  [note] %s\n", msg))
  NOTE <<- c(NOTE, msg)
}

section <- function(name) cat(sprintf("\n%s\n", name))

top_value <- function(name) {
  env <- new.env(parent = globalenv())
  for (e in parse(APP))
    if (is.call(e) && length(e) == 3L && is.name(e[[2]]) &&
        as.character(e[[2]]) == name) { eval(e, env); break }
  env[[name]]
}

cat(strrep("=", 72), "\n", sep = "")
cat("ReferenceRangeR smoke test\n")
cat(sprintf("  %s\n", R.version.string))
cat(sprintf("  %s on %s\n", APP, Sys.info()[["sysname"]]))
cat(strrep("=", 72), "\n", sep = "")


section("R packages")

REQUIRED <- c("shiny", "bslib", "ggplot2", "htmltools", "data.table", "readxl", "DT",
              "refineR", "reflimR", "tidykosmic", "qgam", "mgcv", "scales", "promises",
              "MASS", "nlme", "msm", "geoR", "modeest", "admisc", "date",
              "snpar", "tseries", "stringr")
OPTIONAL <- c("mirai", "png")

broken <- character(0)
for (p in REQUIRED) {
  loaded <- suppressWarnings(suppressMessages(suppressPackageStartupMessages(
    tryCatch({ library(p, character.only = TRUE); TRUE }, error = function(e) FALSE))))
  if (!loaded) broken <- c(broken, p)
}
check(sprintf("all %d required packages load", length(REQUIRED)), length(broken) == 0L,
      paste("cannot load:", paste(broken, collapse = ", ")))
for (p in OPTIONAL)
  if (!nzchar(system.file(package = p)))
    remark(sprintf("optional package '%s' missing%s", p,
                   if (p == "mirai") " -- calculations run synchronously"
                   else " -- PDF export loses its page layout"))


section("app.R")

parsed <- tryCatch({ parse(APP); TRUE }, error = function(e) FALSE)
check("app.R parses", parsed)

E <- new.env(parent = globalenv())
core_ok <- tryCatch({
  for (e in parse(APP))
    if (is.call(e) && length(e) == 3L && is.name(e[[2]]) &&
        as.character(e[[2]]) == "CORE") { eval(e, E); break }
  eval(E$CORE, E)
  TRUE
}, error = function(e) FALSE)
check("CORE block evaluates", core_ok)

check("compute functions defined",
      all(vapply(c("compute_reference_interval", "compute_sex_difference",
                   "compute_age_drift", "prepare_data", "read_table_file",
                   "make_demo_data"), function(f) is.function(E[[f]]), logical(1))))

check("upload limit is set", grepl("shiny.maxRequestSize",
                                   paste(readLines(APP, warn = FALSE), collapse = " ")))


section("method scripts")

BASE <- tryCatch(top_value("BASE_PATH"), error = function(e) NA_character_)
cat(sprintf("  BASE_PATH = %s\n", BASE))
check("BASE_PATH directory exists", !is.na(BASE) && dir.exists(BASE),
      "the tmc and tml methods will try to download their scripts from GitHub")

for (f in c("TMC.settings.R", "TMC.functions.R", "TML.R"))
  check(sprintf("%s present", f), !is.na(BASE) && file.exists(file.path(BASE, f)))

if (!is.na(BASE) && all(file.exists(file.path(BASE, c("TMC.settings.R", "TMC.functions.R"))))) {
  loaded <- tryCatch({
    suppressWarnings(suppressMessages(E$load_method_scripts(
      c("TMC.settings.R", "TMC.functions.R"), BASE, E$TMC_PACKAGES)))
    TRUE
  }, error = function(e) FALSE)
  check("TMC sources and defines tmc()", loaded && is.function(get0("tmc")))
} else check("TMC sources and defines tmc()", FALSE, "scripts missing")

if (!is.na(BASE) && file.exists(file.path(BASE, "TML.R"))) {
  loaded <- tryCatch({
    suppressWarnings(suppressMessages(E$load_method_scripts("TML.R", BASE, E$TML_PACKAGES)))
    TRUE
  }, error = function(e) FALSE)
  check("TML sources and defines tml()", loaded && is.function(get0("tml")))
} else check("TML sources and defines tml()", FALSE, "TML.R missing")


section("static assets")

assets <- c("www/rrr.webp", "www/umo.svg")
gone <- assets[!file.exists(file.path(APPDIR, assets))]
if (length(gone)) {
  remark(sprintf("%s missing -- the app runs, but the navbar logo will not render",
                 paste(gone, collapse = " and ")))
} else {
  cat("  [ok  ] static assets present\n")
}


section("user interface")

ui_ok <- tryCatch({
  u <- new.env(parent = globalenv())
  old <- setwd(APPDIR)
  on.exit(setwd(old), add = TRUE)
  suppressWarnings(suppressMessages(suppressPackageStartupMessages(
    sys.source(APP, envir = u, toplevel.env = u))))
  nchar(as.character(htmltools::renderTags(u$ui)$html)) > 10000
}, error = function(e) FALSE)
check("UI builds", ui_ok)


section("end to end")

e2e <- tryCatch({
  set.seed(1)
  d <- E$make_demo_data(1000L)
  p <- E$prepare_data(d, E$guess_columns(d))
  list(n = p$report$n_usable,
       r = E$compute_reference_interval(p$data$result, method = "reflimr", base_path = BASE))
}, error = function(e) list(n = 0L, r = list(ok = FALSE, message = conditionMessage(e))))

check("demo data prepares", isTRUE(e2e$n >= 1000L))
check("reflimR returns a usable interval",
      isTRUE(e2e$r$ok) && all(is.finite(e2e$r$limits)) && e2e$r$limits[1] < e2e$r$limits[2],
      if (isTRUE(e2e$r$ok)) "" else paste(e2e$r$message))
if (isTRUE(e2e$r$ok))
  cat(sprintf("         %d values -> %.4g to %.4g\n", e2e$r$n, e2e$r$limits[1], e2e$r$limits[2]))

check("small selections are refused",
      !isTRUE(E$compute_reference_interval(rnorm(10), method = "reflimr")$ok))


section("writable locations")

tmp <- tempfile(fileext = ".pdf")
check("tempdir is writable",
      tryCatch({ grDevices::pdf(tmp); grDevices::dev.off(); file.exists(tmp) },
               error = function(e) FALSE))
unlink(tmp)

cache <- tryCatch(tools::R_user_dir("ReferenceRangeR", "cache"), error = function(e) NA)
if (!is.na(cache) && !dir.exists(cache))
  remark(sprintf("method-script cache %s does not exist yet (created on demand)", cache))


cat("\n", strrep("=", 72), "\n", sep = "")
if (length(FAIL) == 0L) {
  cat(sprintf("  READY    %d checks passed%s\n", PASS,
              if (length(NOTE)) sprintf(", %d note(s)", length(NOTE)) else ""))
} else {
  cat(sprintf("  NOT READY    %d of %d checks failed\n\n", length(FAIL), PASS + length(FAIL)))
  for (f in FAIL) cat("    ", f, "\n", sep = "")
}
cat(strrep("=", 72), "\n\n", sep = "")
quit(status = if (length(FAIL) == 0L) 0L else 1L, save = "no")

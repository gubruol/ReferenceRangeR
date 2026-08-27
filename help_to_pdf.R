# Renders help.md to www/help.pdf.
# Rscript help_to_pdf.R [input.md] [output.pdf]

args   <- commandArgs(trailingOnly = TRUE)
HERE   <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "."
INPUT  <- if (length(args) >= 1) args[1] else file.path(HERE, "help.md")
OUTPUT <- if (length(args) >= 2) args[2] else file.path(HERE, "www", "help.pdf")

NAVY <- "#003B73"; GREY <- "#5A6672"; LINE <- "#D6DBE1"
W <- 8.27; H <- 11.69                      # A4 portrait, inches
ML <- 0.85; MR <- 0.85; MT <- 0.80; MB <- 0.75
TW <- W - ML - MR                          # text width

PT <- 1 / 72.27                            # cex is relative to 12pt on this device
size <- function(pt) pt / 12

state <- new.env(parent = emptyenv())
state$y <- MT
state$page <- 0L
state$maxx <- 0            # widest point reached, checked against the margin at the end

# ---------------------------------------------------------------- inline runs
# Splits a line into styled runs. face: 1 plain, 2 bold, 3 italic; mono for `code`.
parse_inline <- function(txt) {
  runs <- list(); buf <- ""; face <- 1L; mono <- FALSE
  # `code` is rendered in italic rather than a monospace family: only the sans
  # faces are embedded in the PDF, so a mono request has no font to fall back on
  # and its glyphs come out broken.
  push <- function() {
    if (nzchar(buf)) runs[[length(runs) + 1L]] <<-
      list(text = buf, face = if (mono) 3L else face, mono = FALSE)
    buf <<- ""
  }
  i <- 1L; n <- nchar(txt)
  while (i <= n) {
    two <- substr(txt, i, i + 1L); one <- substr(txt, i, i)
    if (two == "**") {
      push(); face <- if (face == 2L) 1L else 2L; i <- i + 2L
    } else if (one == "*") {
      push(); face <- if (face == 3L) 1L else 3L; i <- i + 1L
    } else if (one == "`") {
      push(); mono <- !mono; i <- i + 1L
    } else {
      buf <- paste0(buf, one); i <- i + 1L
    }
  }
  push()
  if (!length(runs)) runs <- list(list(text = "", face = 1L, mono = FALSE))
  runs
}

run_family <- function(mono) "sans"

measure <- function(txt, cex, face, mono)
  graphics::strwidth(txt, units = "inches", cex = cex, font = face,
                     family = run_family(mono))

# ------------------------------------------------------------------- paging
new_page <- function() {
  if (state$page > 0L) footer()
  graphics::plot.new()          # the first call fills the page pdf() already opened
  state$page <- state$page + 1L
  graphics::par(mai = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  graphics::plot.window(xlim = c(0, W), ylim = c(H, 0))
  state$y <- MT
}

footer <- function() {
  graphics::segments(ML, H - MB + 0.22, W - MR, H - MB + 0.22, col = LINE)
  graphics::text(W - MR, H - MB + 0.40, sprintf("%d", state$page),
                 adj = c(1, 0.5), cex = size(8), col = GREY, family = "sans")
  graphics::text(ML, H - MB + 0.40, "ReferenceRangeR — user documentation",
                 adj = c(0, 0.5), cex = size(8), col = GREY, family = "sans")
}

need <- function(h) if (state$y + h > H - MB) new_page()

# --------------------------------------------------------------- text output
# Draws styled runs with word wrapping. Returns invisibly; advances state$y.
draw_runs <- function(runs, cex, x, width, leading, col = "black", first_x = x) {
  # A word separator may sit at the end of one run and the start of the next
  # ("in the **Columns** panel"), so carry the trailing space across runs.
  words <- list(); prev_space <- FALSE
  for (r in runs) {
    lead <- startsWith(r$text, " "); trail <- endsWith(r$text, " ")
    parts <- strsplit(r$text, " ", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts)]
    for (k in seq_along(parts))
      words[[length(words) + 1L]] <- list(
        text = parts[k], face = r$face, mono = r$mono,
        space_before = if (k == 1L) (prev_space || lead) else TRUE)
    prev_space <- if (length(parts)) trail else (prev_space || lead || trail)
  }
  line <- list(); lx <- 0; cur_x <- first_x
  flush <- function(last) {
    if (!length(line)) return(invisible())
    need(leading)
    px <- cur_x
    for (w in line) {
      s <- if (w$space_before && px > cur_x) paste0(" ", w$text) else w$text
      graphics::text(px, state$y, s, adj = c(0, 0.5), cex = cex, font = w$face,
                     col = col, family = run_family(w$mono))
      px <- px + measure(s, cex, w$face, w$mono)
    }
    state$maxx <- max(state$maxx, px)
    state$y <- state$y + leading
    line <<- list(); lx <<- 0; cur_x <<- x
  }
  avail <- function() width - (cur_x - x)
  for (w in words) {
    piece <- function() if (w$space_before && length(line)) paste0(" ", w$text) else w$text
    if (length(line) && lx + measure(piece(), cex, w$face, w$mono) > avail()) flush(FALSE)
    lx <- lx + measure(piece(), cex, w$face, w$mono)
    line[[length(line) + 1L]] <- w
  }
  flush(TRUE)
  invisible()
}

para <- function(txt, cex = size(9.5), leading = 0.165, gap = 0.075, col = "black") {
  draw_runs(parse_inline(txt), cex, ML, TW, leading, col)
  state$y <- state$y + gap
}

heading <- function(txt, level) {
  spec <- switch(as.character(level),
    "1" = list(cex = size(21), col = NAVY, face = 2L, before = 0.00, after = 0.16, lead = 0.34),
    "2" = list(cex = size(13.5), col = NAVY, face = 2L, before = 0.30, after = 0.11, lead = 0.24),
    list(cex = size(10.5), col = "black", face = 2L, before = 0.20, after = 0.06, lead = 0.20))
  state$y <- state$y + spec$before
  need(spec$lead + 0.2)
  runs <- lapply(parse_inline(txt), function(r) { r$face <- spec$face; r })
  draw_runs(runs, spec$cex, ML, TW, spec$lead, spec$col)
  state$y <- state$y + spec$after
}

bullet <- function(txt, marker = NULL) {
  cex <- size(9.5); indent <- 0.22
  need(0.165)
  if (is.null(marker))
    graphics::points(ML + 0.07, state$y, pch = 16, cex = 0.30, col = NAVY)
  else
    graphics::text(ML, state$y, marker, adj = c(0, 0.5), cex = cex, font = 2, col = NAVY,
                   family = "sans")
  draw_runs(parse_inline(txt), cex, ML + indent, TW - indent, 0.165, first_x = ML + indent)
  state$y <- state$y + 0.030
}

rule <- function() {
  state$y <- state$y + 0.10
  need(0.10)
  graphics::segments(ML, state$y, W - MR, state$y, col = LINE)
  state$y <- state$y + 0.16
}

# --------------------------------------------------------------------- tables
strip_md <- function(s) gsub("[*`]", "", s)

split_row <- function(ln) {
  ln <- sub("^\\s*\\|", "", ln); ln <- sub("\\|\\s*$", "", ln)
  trimws(strsplit(ln, "|", fixed = TRUE)[[1]])
}

draw_table <- function(rows) {
  header <- rows[[1]]; body <- rows[-1]
  ncol <- length(header)
  cex <- size(8.5); pad <- 0.10; leading <- 0.155
  natural <- vapply(seq_len(ncol), function(j) {
    cells <- c(header[j], vapply(body, function(r) if (length(r) >= j) r[j] else "", character(1)))
    max(vapply(cells, function(s) measure(strip_md(s), cex, 1L, FALSE), numeric(1))) + 2 * pad
  }, numeric(1))
  widths <- if (sum(natural) <= TW) natural * (TW / sum(natural)) else natural * (TW / sum(natural))
  xs <- ML + c(0, cumsum(widths)[-ncol])

  wrap_cell <- function(s, w) {
    words <- strsplit(strip_md(s), " +")[[1]]
    if (!length(words)) return("")
    out <- character(); cur <- ""
    for (wd in words) {
      cand <- if (nzchar(cur)) paste(cur, wd) else wd
      if (nzchar(cur) && measure(cand, cex, 1L, FALSE) > w - 2 * pad) {
        out <- c(out, cur); cur <- wd
      } else cur <- cand
    }
    c(out, cur)
  }

  row_out <- function(cells, face) {
    lines <- lapply(seq_len(ncol), function(j)
      wrap_cell(if (length(cells) >= j) cells[j] else "", widths[j]))
    hgt <- max(vapply(lines, length, integer(1))) * leading
    need(hgt + 0.06)
    y0 <- state$y
    for (j in seq_len(ncol)) {
      yy <- y0
      for (ln in lines[[j]]) {
        graphics::text(xs[j] + pad, yy, ln, adj = c(0, 0.5), cex = cex, font = face,
                       col = if (face == 2L) NAVY else "black", family = "sans")
        state$maxx <- max(state$maxx, xs[j] + pad + measure(ln, cex, face, FALSE))
        yy <- yy + leading
      }
    }
    state$y <- y0 + hgt
  }

  state$y <- state$y + 0.06
  need(0.5)
  row_out(header, 2L)
  state$y <- state$y + 0.03
  graphics::segments(ML, state$y, ML + sum(widths), state$y, col = LINE)
  state$y <- state$y + 0.09
  for (r in body) {
    row_out(r, 1L)
    state$y <- state$y + 0.03
    graphics::segments(ML, state$y, ML + sum(widths), state$y, col = LINE)
    state$y <- state$y + 0.09
  }
  state$y <- state$y + 0.06
}

# ---------------------------------------------------------------------- main
md <- readLines(INPUT, encoding = "UTF-8", warn = FALSE)
Encoding(md) <- "UTF-8"

dir.create(dirname(OUTPUT), showWarnings = FALSE, recursive = TRUE)
# cairo_pdf is preferred over pdf(): R's own pdf device remaps character 45 to
# the minus glyph, which would leave every hyphen in a URL or DOI uncopyable.
if (capabilities("cairo")) {
  grDevices::cairo_pdf(OUTPUT, width = W, height = H, onefile = TRUE,
                       family = "sans", fallback_resolution = 300)
} else {
  grDevices::pdf(OUTPUT, width = W, height = H, onefile = TRUE, encoding = "WinAnsi.enc",
                 title = "ReferenceRangeR")
}
new_page()

i <- 1L
while (i <= length(md)) {
  ln <- md[i]; trimmed <- trimws(ln)

  if (!nzchar(trimmed)) {
    i <- i + 1L; next
  }
  if (grepl("^---+$", trimmed)) {
    rule(); i <- i + 1L; next
  }
  if (grepl("^#{1,6} ", trimmed)) {
    lvl <- nchar(sub("^(#+) .*$", "\\1", trimmed))
    heading(sub("^#+ ", "", trimmed), lvl); i <- i + 1L; next
  }
  if (startsWith(trimmed, "|")) {
    rows <- list(); j <- i
    while (j <= length(md) && startsWith(trimws(md[j]), "|")) {
      cells <- split_row(md[j])
      if (!all(grepl("^:?-{2,}:?$", cells))) rows[[length(rows) + 1L]] <- cells
      j <- j + 1L
    }
    draw_table(rows); i <- j; next
  }
  if (grepl("^[-*] ", trimmed) || grepl("^[0-9]+\\. ", trimmed)) {
    numbered <- grepl("^[0-9]+\\. ", trimmed)
    marker <- if (numbered) sub("^([0-9]+\\.) .*$", "\\1", trimmed) else NULL
    txt <- sub("^([-*]|[0-9]+\\.) ", "", trimmed)
    j <- i + 1L                       # absorb continuation lines
    while (j <= length(md) && grepl("^\\s{2,}\\S", md[j]) &&
           !grepl("^\\s*([-*]|[0-9]+\\.) ", md[j])) {
      txt <- paste(txt, trimws(md[j])); j <- j + 1L
    }
    bullet(txt, marker); i <- j; next
  }
  # paragraph: join until a blank line or a block marker
  txt <- trimmed; j <- i + 1L
  while (j <= length(md) && nzchar(trimws(md[j])) &&
         !grepl("^(#{1,6} |[-*] |[0-9]+\\. |\\||---+$)", trimws(md[j]))) {
    txt <- paste(txt, trimws(md[j])); j <- j + 1L
  }
  para(txt); i <- j
}

footer()
grDevices::dev.off()
cat(sprintf("wrote %s (%d pages, %.0f kB)\n", OUTPUT, state$page,
            file.size(OUTPUT) / 1024))
if (state$maxx > W - MR + 0.01)
  warning(sprintf("content reaches %.2f in, past the right margin at %.2f in",
                  state$maxx, W - MR), call. = FALSE)

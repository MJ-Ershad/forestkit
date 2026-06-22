#' forestkit: publication-quality meta-analysis figures
#'
#' Two figure helpers that share one philosophy: compute and save the numeric
#' result to disk *first*, then draw the figure, so a plotting failure can never
#' cost you the analysis. Every figure is written twice — a 300-DPI PNG for
#' submission and a vector PDF for infinite sharpness.

# internal: device dimensions that grow with the number of rows
.fig_height <- function(k, per_row = 0.34, base = 1.0, floor = 3.0)
  max(floor, base + per_row * k)

#' Publication-quality forest plot for a continuous between-group pool
#'
#' Fits a three-level random-effects model and draws a forest plot, saving a
#' result CSV, a 300-DPI PNG, and a vector PDF under \code{dir}.
#'
#' @param data    Data frame with group means/SDs/Ns and id columns.
#' @param label   Outcome label used in titles and the result row.
#' @param stem    File stem for the three output files.
#' @param dir     Output directory (created if needed). Default ".".
#' @param cols    Named list mapping required roles to column names. Defaults to
#'   \code{m_e, sd_e, n_e, m_c, sd_c, n_c, study, cluster, slab}.
#' @param measure Effect-size measure passed to \code{metafor::escalc}. "MD" or "SMD".
#' @return Invisibly, the one-row result data frame.
#' @export
pub_forest <- function(data, label, stem, dir = ".",
                       cols = list(m_e = "m_e", sd_e = "sd_e", n_e = "n_e",
                                   m_c = "m_c", sd_c = "sd_c", n_c = "n_c",
                                   study = "study_id", cluster = "cluster",
                                   slab = "slab"),
                       measure = "MD") {
  if (!requireNamespace("metafor", quietly = TRUE))
    stop("forestkit needs the 'metafor' package.", call. = FALSE)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  g <- function(role) data[[cols[[role]]]]
  d <- data.frame(m_e = g("m_e"), sd_e = g("sd_e"), n_e = g("n_e"),
                  m_c = g("m_c"), sd_c = g("sd_c"), n_c = g("n_c"),
                  study = g("study"), cluster = g("cluster"),
                  slab = if (!is.null(cols$slab) && cols$slab %in% names(data))
                           g("slab") else as.character(g("study")))
  d <- d[stats::complete.cases(d[, c("m_e","sd_e","n_e","m_c","sd_c","n_c")]), ]
  d$cluster <- ifelse(is.na(d$cluster) | d$cluster == "",
                      paste0("SOLO_", d$study), as.character(d$cluster))
  if (nrow(d) < 2) stop("Need at least 2 rows with complete data.", call. = FALSE)

  es <- metafor::escalc(measure = measure, m1i = m_e, sd1i = sd_e, n1i = n_e,
                        m2i = m_c, sd2i = sd_c, n2i = n_c, data = d)
  es$effect <- seq_len(nrow(es))
  m <- metafor::rma.mv(yi, vi, random = ~ 1 | cluster / study / effect,
                       data = es, method = "REML", test = "t", dfs = "contain")

  # ---- result FIRST ----
  res <- data.frame(label = label, measure = measure, k = nrow(es),
                    estimate = as.numeric(m$b), ci_lo = m$ci.lb, ci_hi = m$ci.ub,
                    p = m$pval, tau2 = sum(m$sigma2), stringsAsFactors = FALSE)
  utils::write.csv(res, file.path(dir, paste0(stem, "_result.csv")), row.names = FALSE)

  # ---- figure (PNG @300dpi + vector PDF) ----
  win <- 9.2; hin <- .fig_height(nrow(es))
  xl  <- paste0(label, " (", measure, ")")
  draw <- function() metafor::forest(
    m, slab = d$slab, header = c("Study", paste0(measure, " [95% CI]")),
    xlab = xl, cex = 0.85, mgp = c(2, 0.6, 0))

  for (dev in c("png", "pdf")) {
    fn <- file.path(dir, paste0(stem, "_forest.", dev))
    tryCatch({
      if (dev == "png") grDevices::png(fn, width = win, height = hin, units = "in", res = 300)
      else grDevices::pdf(fn, width = win, height = hin)
      draw(); grDevices::dev.off()
      message("saved ", fn)
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      message("skipped ", fn, " (", conditionMessage(e), ") — result.csv already written")
    })
  }
  invisible(res)
}

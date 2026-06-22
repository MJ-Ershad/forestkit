#' Restricted natural-spline dose-response figure with a robust CI band
#'
#' Models a between-group mean difference as a smooth function of a continuous
#' dose moderator using a natural cubic spline basis, with cluster-robust
#' inference, and draws the predicted curve plus a confidence band, the study
#' points (sized by precision), and the knot positions. As with
#' \code{pub_forest}, the numeric summary is written before the figure.
#'
#' @param data    Data frame with effect-size inputs and a continuous \code{dose}.
#' @param label   Outcome label for the title and result row.
#' @param stem    Output file stem.
#' @param dir     Output directory. Default ".".
#' @param df      Spline degrees of freedom (number of basis columns). Default 3.
#' @param cols    Named column map; see \code{pub_forest}. Adds \code{dose}.
#' @return Invisibly, the one-row result data frame.
#' @export
dose_spline <- function(data, label, stem, dir = ".", df = 3,
                        cols = list(m_e = "m_e", sd_e = "sd_e", n_e = "n_e",
                                    m_c = "m_c", sd_c = "sd_c", n_c = "n_c",
                                    study = "study_id", cluster = "cluster",
                                    dose = "dose")) {
  for (pkg in c("metafor", "clubSandwich"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("forestkit needs the '", pkg, "' package.", call. = FALSE)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  g <- function(role) data[[cols[[role]]]]
  d <- data.frame(m_e = g("m_e"), sd_e = g("sd_e"), n_e = g("n_e"),
                  m_c = g("m_c"), sd_c = g("sd_c"), n_c = g("n_c"),
                  study = g("study"), cluster = g("cluster"),
                  dose = suppressWarnings(as.numeric(g("dose"))))
  d <- d[stats::complete.cases(d), ]
  if (nrow(d) < df + 2) stop("Too few complete rows for df = ", df, ".", call. = FALSE)
  d$cluster <- factor(ifelse(is.na(d$cluster) | d$cluster == "",
                             paste0("SOLO_", d$study), as.character(d$cluster)))
  es <- metafor::escalc(measure = "MD", m1i = m_e, sd1i = sd_e, n1i = n_e,
                        m2i = m_c, sd2i = sd_c, n2i = n_c, data = d)
  es$effect <- seq_len(nrow(es))

  # natural cubic spline basis on dose (base 'splines', no rms dependency)
  B <- splines::ns(es$dose, df = df)
  knots <- attr(B, "knots"); bnd <- attr(B, "Boundary.knots")
  colnames(B) <- paste0("s", seq_len(ncol(B)))

  m_sp  <- metafor::rma.mv(yi, vi, mods = B, random = ~ 1 | cluster / study / effect,
                           data = es, method = "REML", test = "t")
  m_lin <- metafor::rma.mv(yi, vi, mods = ~ dose, random = ~ 1 | cluster / study / effect,
                           data = es, method = "ML", test = "t")
  m_spml <- metafor::rma.mv(yi, vi, mods = B, random = ~ 1 | cluster / study / effect,
                            data = es, method = "ML", test = "t")
  lrt <- metafor::anova.rma(m_lin, m_spml)
  rob <- clubSandwich::coef_test(m_sp, vcov = "CR2", cluster = es$cluster)
  # overall test of the nonlinear part = joint test of the spline columns beyond linear
  nonlinear_p <- tryCatch(lrt$pval, error = function(e) NA_real_)

  # ---- result FIRST ----
  res <- data.frame(label = label, k = nrow(es),
                    clusters = nlevels(es$cluster), df = df,
                    dose_min = min(es$dose), dose_max = max(es$dose),
                    nonlinear_LRT_p = round(as.numeric(nonlinear_p), 4),
                    AIC_linear = round(stats::AIC(m_lin), 2),
                    AIC_spline = round(stats::AIC(m_spml), 2),
                    stringsAsFactors = FALSE)
  utils::write.csv(res, file.path(dir, paste0(stem, "_spline_result.csv")), row.names = FALSE)

  # ---- predicted curve with robust band ----
  V <- clubSandwich::vcovCR(m_sp, type = "CR2", cluster = es$cluster)
  grid <- seq(min(es$dose), max(es$dose), length.out = 200)
  Bg <- predict(B, grid)
  Xg <- cbind(1, Bg)
  pred <- as.numeric(Xg %*% as.numeric(m_sp$beta))
  se <- sqrt(pmax(0, rowSums((Xg %*% V) * Xg)))
  tcrit <- stats::qt(0.975, df = max(2, nlevels(es$cluster) - ncol(Xg)))
  lo <- pred - tcrit * se; hi <- pred + tcrit * se

  for (dev in c("png", "pdf")) {
    fn <- file.path(dir, paste0(stem, "_spline.", dev))
    if (dev == "png") grDevices::png(fn, width = 9.2, height = 5.4, units = "in", res = 300)
    else grDevices::pdf(fn, width = 9.2, height = 5.4)
    graphics::par(mar = c(4.6, 4.8, 3.0, 1.2))
    graphics::plot(NA, xlim = range(grid), ylim = range(c(lo, hi, es$yi)),
                   xlab = "Dose (continuous moderator)",
                   ylab = "Treatment - control mean difference",
                   main = sprintf("%s\nnatural spline (df = %d); nonlinearity LRT p = %s",
                                  label, df, format(round(as.numeric(nonlinear_p), 3))),
                   cex.main = 0.95)
    graphics::abline(h = 0, col = "grey55", lty = 2)
    graphics::polygon(c(grid, rev(grid)), c(lo, rev(hi)),
                      col = grDevices::rgb(0.20, 0.40, 0.70, 0.18), border = NA)
    graphics::lines(grid, pred, col = grDevices::rgb(0.12, 0.30, 0.62), lwd = 2.4)
    wt <- 1 / es$vi
    cexv <- 0.8 + 2.0 * (wt - min(wt)) / (max(wt) - min(wt) + 1e-9)
    graphics::points(es$dose, es$yi, pch = 21, cex = cexv, col = "grey20",
                     bg = grDevices::rgb(0.85, 0.45, 0.10, 0.6))
    graphics::rug(c(bnd, knots), lwd = 2, col = "grey30")
    grDevices::dev.off()
    message("saved ", fn)
  }
  invisible(res)
}

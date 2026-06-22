#' Dose-response restricted cubic spline plot
#'
#' Fits a restricted cubic spline of \code{y} on \code{dose} (via \pkg{rms} if
#' available, otherwise a natural spline through \code{stats::lm}) and plots the
#' fitted curve with a pointwise confidence band over the observed points.
#' Writes a 300-dpi PNG plus a vector PDF.
#'
#' As with \code{pub_forest()}, the fitted prediction grid is written to
#' \code{<file>_fit.csv} \emph{before} the plotting device is opened, so the
#' numbers survive even if rendering fails.
#'
#' @param dose Numeric predictor (e.g. cumulative drug or exposure dose).
#' @param y Numeric outcome.
#' @param knots Number of spline knots (default 3) or a numeric vector of knot
#'   locations.
#' @param file Output path stem (without extension).
#' @param xlab,ylab Axis labels.
#' @param level Confidence level for the band (default 0.95).
#' @return Invisibly, the prediction grid data frame.
#' @export
dose_spline <- function(dose, y, knots = 3, file = "spline",
                        xlab = "Dose", ylab = "Outcome", level = 0.95) {
  ok <- stats::complete.cases(dose, y)
  dose <- dose[ok]; y <- y[ok]
  stopifnot(length(dose) >= 4)

  grid <- seq(min(dose), max(dose), length.out = 200)

  use_rms <- requireNamespace("rms", quietly = TRUE)
  if (use_rms) {
    dd <- data.frame(dose = dose, y = y)
    kn <- if (length(knots) == 1) knots else knots
    fit <- rms::ols(y ~ rms::rcs(dose, kn), data = dd)
    pr  <- rms::Predict(fit, dose = grid, conf.int = level)
    pred <- data.frame(dose = pr$dose, fit = pr$yhat,
                       lower = pr$lower, upper = pr$upper)
    nl_p <- tryCatch(stats::anova(fit)["dose", "P"][[1]],
                     error = function(e) NA_real_)
  } else {
    if (!requireNamespace("splines", quietly = TRUE)) {
      stop("Need either 'rms' or the base 'splines' package.", call. = FALSE)
    }
    df_spline <- if (length(knots) == 1) knots else length(knots)
    fit <- stats::lm(y ~ splines::ns(dose, df = df_spline))
    p   <- stats::predict(fit, newdata = data.frame(dose = grid),
                          interval = "confidence", level = level)
    pred <- data.frame(dose = grid, fit = p[, "fit"],
                       lower = p[, "lwr"], upper = p[, "upr"])
    nl_p <- NA_real_
  }

  # --- write the fitted numbers FIRST -------------------------------------- #
  tryCatch(
    utils::write.csv(pred, paste0(file, "_fit.csv"), row.names = FALSE),
    error = function(e)
      warning("could not write fit csv: ", conditionMessage(e), call. = FALSE)
  )

  width <- 9.2; height <- 5.5
  draw <- function() {
    op <- graphics::par(mar = c(4.5, 4.5, 1, 1)); on.exit(graphics::par(op))
    plot(dose, y, pch = 16, col = grDevices::adjustcolor("grey40", 0.6),
         xlab = xlab, ylab = ylab, bty = "l")
    graphics::polygon(c(pred$dose, rev(pred$dose)),
                      c(pred$lower, rev(pred$upper)),
                      col = grDevices::adjustcolor("steelblue", 0.20),
                      border = NA)
    graphics::lines(pred$dose, pred$fit, col = "steelblue", lwd = 2.2)
    if (!is.na(nl_p)) {
      graphics::legend("topright", bty = "n",
                       legend = sprintf("nonlinearity p = %.3f", nl_p))
    }
  }

  grDevices::png(paste0(file, ".png"), width = width, height = height,
                 units = "in", res = 300)
  draw(); grDevices::dev.off()
  grDevices::pdf(paste0(file, ".pdf"), width = width, height = height)
  draw(); grDevices::dev.off()

  message(sprintf("forestkit: wrote %s.png, %s.pdf, %s_fit.csv",
                  file, file, file))
  invisible(pred)
}

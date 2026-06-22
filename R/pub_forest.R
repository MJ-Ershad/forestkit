#' Publication-quality forest plot
#'
#' Draws a forest plot of point estimates with confidence intervals and writes
#' both a 300-dpi PNG and a vector PDF. A pooled diamond is drawn at the bottom
#' if \code{pooled} is supplied.
#'
#' The plot height scales with the number of rows so labels never collide:
#' \code{height = max(3.0, 1.0 + 0.34 * k)} inches, width ~9.2 inches.
#'
#' Safety rule: the numeric table behind the figure is written to
#' \code{<file>_data.csv} \emph{before} the device is opened, inside a
#' \code{tryCatch}, so a plotting crash can never lose the underlying numbers.
#'
#' @param estimate Numeric vector of point estimates.
#' @param lower,upper Numeric vectors of CI bounds (same length).
#' @param label Character vector of row labels.
#' @param file Output path stem (without extension).
#' @param pooled Optional named list \code{list(estimate=, lower=, upper=,
#'   label=)} drawn as a summary diamond.
#' @param xlab Axis label.
#' @param ref Reference line (default 0).
#' @return Invisibly, the data frame that was plotted.
#' @export
pub_forest <- function(estimate, lower, upper, label,
                       file = "forest",
                       pooled = NULL, xlab = "Effect size", ref = 0) {
  stopifnot(length(estimate) == length(lower),
            length(lower) == length(upper),
            length(upper) == length(label))

  df <- data.frame(label = label, estimate = estimate,
                   lower = lower, upper = upper,
                   stringsAsFactors = FALSE)

  # --- write numbers FIRST, before any device is opened -------------------- #
  tryCatch(
    utils::write.csv(df, paste0(file, "_data.csv"), row.names = FALSE),
    error = function(e)
      warning("could not write data csv: ", conditionMessage(e), call. = FALSE)
  )

  k <- nrow(df)
  width  <- 9.2
  height <- max(3.0, 1.0 + 0.34 * k)

  draw <- function() {
    op <- graphics::par(mar = c(4, 14, 1, 2)); on.exit(graphics::par(op))
    ys <- rev(seq_len(k))
    yp <- if (is.null(pooled)) NULL else 0
    xr <- range(c(lower, upper, pooled$lower, pooled$upper, ref), na.rm = TRUE)
    pad <- diff(xr) * 0.08
    plot(NA, xlim = c(xr[1] - pad, xr[2] + pad),
         ylim = c(if (is.null(pooled)) 0.5 else -0.5, k + 0.5),
         yaxt = "n", ylab = "", xlab = xlab, bty = "n")
    graphics::abline(v = ref, col = "grey60", lty = 2)
    graphics::axis(2, at = c(ys, yp), labels = c(df$label, pooled$label),
                   las = 1, tick = FALSE, cex.axis = 0.9)
    graphics::segments(df$lower, ys, df$upper, ys, lwd = 1.6)
    graphics::points(df$estimate, ys, pch = 15, cex = 1.3)
    if (!is.null(pooled)) {
      d <- pooled
      graphics::polygon(c(d$lower, d$estimate, d$upper, d$estimate),
                        c(0, 0.22, 0, -0.22), col = "black", border = NA)
    }
  }

  grDevices::png(paste0(file, ".png"), width = width, height = height,
                 units = "in", res = 300)
  draw(); grDevices::dev.off()

  grDevices::pdf(paste0(file, ".pdf"), width = width, height = height)
  draw(); grDevices::dev.off()

  message(sprintf("forestkit: wrote %s.png, %s.pdf, %s_data.csv (k = %d)",
                  file, file, file, k))
  invisible(df)
}

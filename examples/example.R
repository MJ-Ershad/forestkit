# forestkit examples ---------------------------------------------------------
# install.packages("remotes"); remotes::install_github("MJ-Ershad/forestkit")
library(forestkit)

## 1) Forest plot with a pooled diamond -------------------------------------
set.seed(42)
k    <- 8
est  <- rnorm(k, 0.3, 0.25)
se   <- runif(k, 0.08, 0.18)
pub_forest(
  estimate = est,
  lower    = est - 1.96 * se,
  upper    = est + 1.96 * se,
  label    = paste("Study", LETTERS[1:k]),
  pooled   = list(estimate = mean(est),
                  lower = mean(est) - 0.10,
                  upper = mean(est) + 0.10,
                  label = "Pooled (RE)"),
  xlab = "Standardized mean difference",
  file = "demo_forest"
)
# -> demo_forest.png, demo_forest.pdf, demo_forest_data.csv

## 2) Dose-response spline ---------------------------------------------------
dose <- runif(120, 0, 50)
y    <- 0.2 + 0.04 * dose - 0.0007 * dose^2 + rnorm(120, 0, 0.4)
dose_spline(
  dose = dose, y = y, knots = 4,
  xlab = "Cumulative dose", ylab = "Outcome",
  file = "demo_spline"
)
# -> demo_spline.png, demo_spline.pdf, demo_spline_fit.csv

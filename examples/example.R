# examples/example.R — forestkit on a small synthetic dataset
# Run from the package root:  Rscript examples/example.R
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

set.seed(42)
k <- 24
dat <- data.frame(
  study_id = paste0("s", seq_len(k)),
  cluster  = rep(paste0("lab", 1:8), each = 3),
  slab     = paste0("Study ", seq_len(k)),
  n_e = sample(20:60, k, TRUE), n_c = sample(20:60, k, TRUE),
  m_e = rnorm(k, 52, 4), m_c = rnorm(k, 50, 4),
  sd_e = runif(k, 5, 9),  sd_c = runif(k, 5, 9),
  dose = round(runif(k, 1, 25), 1)
)

dir.create("outputs", showWarnings = FALSE)
pub_forest(dat, label = "Example outcome", stem = "demo", dir = "outputs")
dose_spline(dat, label = "Example outcome", stem = "demo", dir = "outputs", df = 3)
cat("\nWrote forest + spline figures and result CSVs to outputs/\n")

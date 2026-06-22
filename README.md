# forestkit

Two small R helpers that turn meta-analysis output into **publication-quality
figures** — and, importantly, write the numbers to disk *before* drawing so a
rendering crash never costs you the results.

- `pub_forest()` — forest plot with CIs and an optional pooled diamond
- `dose_spline()` — restricted cubic spline dose-response curve with a CI band

Every call produces a **300-dpi PNG**, a **vector PDF**, and a **CSV** of the
underlying values. Figure height scales with the number of rows
(`max(3.0, 1.0 + 0.34 * k)` inches) so labels never overlap.

## Install

```r
# install.packages("remotes")
remotes::install_github("MJ-Ershad/forestkit")
```

`rms` is optional: if present, `dose_spline()` uses `rms::rcs()` and reports a
nonlinearity p-value; otherwise it falls back to a base `splines::ns()` fit.

## Forest plot

```r
library(forestkit)
pub_forest(
  estimate = est, lower = lo, upper = hi, label = study_names,
  pooled   = list(estimate = 0.31, lower = 0.21, upper = 0.41,
                  label = "Pooled (RE)"),
  xlab = "Standardized mean difference",
  file = "fig2_forest"
)
#> forestkit: wrote fig2_forest.png, fig2_forest.pdf, fig2_forest_data.csv (k = 8)
```

## Dose-response spline

```r
dose_spline(
  dose = cumulative_dose, y = outcome, knots = 4,
  xlab = "Cumulative dose", ylab = "Outcome",
  file = "fig3_spline"
)
#> forestkit: wrote fig3_spline.png, fig3_spline.pdf, fig3_spline_fit.csv
```

## The "numbers before plot" rule

Both functions write their data/fit CSV inside a `tryCatch` **before** any
graphics device is opened. If the device or a downstream draw call errors, the
CSV is already on disk — you keep the estimates and just re-render.

See [`examples/example.R`](examples/example.R) for runnable demos.

## License

MIT (c) 2026 Mohamadjavad Ershadmanesh

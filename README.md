# forestkit

Publication-quality forest and dose-response spline figures for meta-analysis in R. The guiding rule: write the numeric result to disk first, then draw the picture. Every figure is saved as a 300-DPI PNG plus a vector PDF.

Two functions: pub_forest builds three-level forest plots, and dose_spline fits a natural cubic spline dose-response curve with cluster-robust confidence bands. Both depend on metafor and clubSandwich.

MIT licensed, 2026 Mohamadjavad Ershadmanesh.

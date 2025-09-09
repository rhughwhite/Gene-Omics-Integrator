# Gene-Omics-Integrator

Multimodal statistical utilities for cancer genomics, centered on fast mutual‑information (MI) estimation and publication‑ready visualizations across heterogeneous data types (continuous and categorical).

These functions were used to support analyses appearing in:
- **Nature Genetics (2025)**: s41588-025-02128-y
- **Nature (2021)**: s41586-021-03850-3

> Note: This repository surfaces the reusable R helpers extracted from internal analyses so others can reproduce similar MI‑based summaries and figures.

---

## Highlights

- Compute MI between two vectors with flexible discretization and optional normalization.
- Pairwise MI across *lists* of matrices/data frames (e.g., gene × sample blocks across data types).
- Permutation‑based significance estimation (wrapper).
- Compact summary plots: violin plots, summary heatmaps, and multi‑panel layouts.
- Simple rank‑based inverse normal transform (RNOmni) helper.
- Colour‑legend helper compatible with `BoutrosLab.plotting.general` styles.

## Installation

This was developed as a standard R package. From the repo root:

```r
# install dependencies you may need:
# install.packages(c('doParallel', 'foreach'))
# remotes::install_github('BoutrosLab/plotting.general') # or your local fork

# install this package from a local path
remotes::install_local('.')
```

> Requires R ≥ 3.0.2. Plots use `BoutrosLab.plotting.general` (≥ 6.0.3). Parallelization uses `foreach`/`doParallel`.

---

## Function quick reference & snippets

Below are excerpts adapted from the package documentation (`man/*.Rd`). Usage blocks are shown verbatim where possible and minimal examples provided for orientation.

### `mutual.information` — MI between two vectors

**Usage**

```r
mutual.information(
  data.1,
  data.2,
  data.1.categorical = FALSE,
  data.2.categorical = FALSE,
  n.bins,
  FD.n.bins = FALSE,
  equal.frequency.bins = FALSE,
  normalize.mutual.info = TRUE,
  use.shrinkage = FALSE,
  simple.sample.size.correction = FALSE
)
```

**Value**: Mutual information (optionally normalized).  
Used internally by higher‑level wrappers.

---

### `get.mutual.information` — wrapper with permutations

**Usage**

```r
get.mutual.information(
  data.1,
  data.2,
  data.1.categorical = FALSE,
  data.2.categorical = FALSE,
  n.bins,
  FD.n.bins = FALSE,
  equal.frequency.bins = FALSE,
  calculate.pvalue,
  n.permutations = 10^5,
  normalize.mutual.info = TRUE,
  use.shrinkage = FALSE,
  simple.sample.size.correction = FALSE
)
```

**Notes**: Wraps `mutual.information`; handles missing data, categorical/continuous inputs, and optional permutation‑based p‑values.

**Minimal example**

```r
set.seed(1)
x <- rnorm(200)
y <- x + rnorm(200, sd = 0.5)
mi <- get.mutual.information(x, y, n.bins = 10, calculate.pvalue = FALSE)
mi
```

---

### `calculate.mutual.information.pairwise.data.types` — MI across data‑type matrices

Compute MI for every pair of rows across aligned matrices (e.g., gene‑wise MI between expression and copy‑number).

**Usage**

```r
calculate.mutual.information.pairwise.data.types(
  omics.data.list,
  data.types.categorical = setNames(rep(FALSE, length(omics.data.list)),
    names(omics.data.list)),
  comparisons = NULL,
  match.rows = TRUE,
  match.columns = TRUE,
  filter.rows = FALSE,
  filter.threshold = NULL,
  normalize.input = FALSE,
  n.permutations,
  calculate.pvalue = FALSE,
  FDR.threshold = 0.05,
  n.cores = 1,
  n.bins = NULL,
  FD.n.bins = FALSE,
  sample.size.bins = TRUE,
  equal.frequency.bins = FALSE,
  normalize.mutual.info = TRUE,
  use.shrinkage = FALSE,
  simple.sample.size.correction = FALSE
)
```

**Workflow sketch**

```r
data('simple.data')  # toy inputs included
omic.list <- list(
  expression = simple.data$expression,
  cnv        = simple.data$cnv
)

MI <- calculate.mutual.information.pairwise.data.types(
  omic.list,
  filter.rows = TRUE,
  filter.threshold = 0.0,   # keep all rows
  normalize.mutual.info = TRUE
)
str(MI)  # list of matrices: observed MI, expected (perm), p.value, etc (see docs)
```

---

### Visualization helpers

#### `create.mutual.information.violinplot`

**Usage**

```r
create.mutual.information.violinplot(
  MI,
  log.scale,
  ylabel = "Mutual Information",
  yaxis.cex = 0.55,
  y.min.log = 0.01,
  style = "Nature",
  yaxis.lab = NULL,
  yat = NULL,
  ylimits = NULL
)
```

```r
vp <- create.mutual.information.violinplot(MI, log.scale = TRUE, ylabel = "Mutual Information")
# draw.plot(vp, filename = 'mi_violinplot.tiff', resolution = 300)  # using BoutrosLab.plotting.general
```

#### `create.mutual.information.summary.statistic.heatmap`

**Usage**

```r
create.mutual.information.summary.statistic.heatmap(
  summary.statistic,
  statistic.name,
  colour.scheme
)
```

```r
sum_stat <- matrix(runif(100), nrow = 10)
hm <- create.mutual.information.summary.statistic.heatmap(sum_stat, 'mean(MI)', colour.scheme = c('white','steelblue'))
```

#### `create.mutual.information.multipanelplot`

**Usage**

```r
create.mutual.information.multipanelplot(
  MI.violinplot,
  data.types.heatmap,
  num.genes.heatmap,
  num.sig.genes.heatmap,
  perc.sig.genes.heatmap,
  plot.objects.heights = c(0.56, 0.14, 0.1, 0.1, 0.1),
  data.types.legend,
  num.genes.legend,
  num.sig.genes.legend,
  perc.sig.genes.legend,
  return.plot = FALSE,
  filename = NULL,
  plot.width = NULL,
  plot.height = NULL,
  y.spacing = -0.5,
  style = "Nature"
)
```

```r
# Combine violin, summary heatmap, and legends into a single multi‑panel figure
mp <- create.mutual.information.multipanelplot(
  MI.violinplot = vp$plot,
  summary.statistic.heatmap = hm$summary.statistic.heatmap,
  summary.statistic.legend  = hm$summary.statistic.legend,
  # ... optional legends for gene sets or significance
)
```

#### `create.data.type.colourlegend.heatmap`

**Usage**

```r
create.data.type.colourlegend.heatmap(
  data.type.pairs,
  data.type.colours,
  data.type.legend.labels
)
```

---

### Utilities

#### `normalize.data`

```r
normalize.data(input.matrix)
```

Applies a rank‑based inverse normal transform (see `RNOmni::RankNorm`); returns a matrix of normalized rows.

#### `bins.fd`

```r
bins.fd(values)
```

Scott/Freedman–Diaconis‑style helper to pick sensible `n.bins` for discretization.

---

## Example end‑to‑end

```r
data('simple.data')

# 1) Pairwise MI between data types
MI <- calculate.mutual.information.pairwise.data.types(
  list(
    expression = simple.data$expression,
    cnv        = simple.data$cnv,
    splicing   = simple.data$splicing
  ),
  normalize.mutual.info = TRUE
)

# 2) Visualize
vp <- create.mutual.information.violinplot(MI, log.scale = TRUE)
sum_stat <- sapply(MI$mutual.information, rowMeans)  # toy summary
hm <- create.mutual.information.summary.statistic.heatmap(sum_stat, 'rowMeans(MI)', c('white','steelblue'))

# 3) Compose into a figure
mp <- create.mutual.information.multipanelplot(
  MI.violinplot = vp$plot,
  summary.statistic.heatmap = hm$summary.statistic.heatmap,
  summary.statistic.legend  = hm$summary.statistic.legend
)
```

---

## Citing

If this code contributes to your work, please cite the relevant papers:

- Nature Genetics (2025): s41588-025-02128-y  
- Nature (2021): s41586-021-03850-3

And consider citing `BoutrosLab.plotting.general` for plotting.
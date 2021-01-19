
# scPOP

scPOP is a lightweight, low dependency R package which brings together
multiple metrics to evaluate batch correction for single cell RNA-seq.
The package includes the Local Simpson Index (LISI) and Average
Silhouette Width (ASW) metrics from Harmony and kBET, respectively, as
well as the Adjusted Rand Index (ARI) and Normalized Mutual Information
(NMI) algorithms.

# Usage

We provide a toy dataset in .h5ad format

## Bioconductor framework example

For some metric, a distance matrix mist be calculated, which requires
significant memory usage. we provide the function `stratified_sample` to
downsample data based on a grouping variable

``` r
library(zellkonverter,quietly = T)
library(SingleCellExperiment)
library(scPOP)
sce <- zellkonverter::readH5AD('scEiaD_all_anndata_mini_ref.h5ad')
idx <- stratified_sample(colnames(sce), sce$batch)
sce_ds <- sce[,idx]
sil_score <- silhouette_width(reducedDim(sce_ds, 'X_scvi'), sce_ds$batch)
```

Both the NMI and ARI metrics compare sets of labels, and thus can scale
well to large data sets

``` r
ari_score <- ari(sce$batch, sce$CellType_predict)
nmi_score <- nmi(sce$batch, sce$CellType_predict)
```

Lisi requires about 10gb of memory for 100K cells and scales linearly
with number of cells

``` r
lisi_score <- lisi(X = reducedDim(sce, 'X_scvi'), meta_data = as.data.frame(colData(sce)), label_colnames = 'CellType_predict' )
```

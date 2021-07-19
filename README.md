
# scPOP

`scPOP` is a lightweight, low dependency R package which brings together
multiple metrics to evaluate batch correction for single cell RNA-seq.
The package includes the Local Simpson Index (LISI) and Average
Silhouette Width (ASW) metrics from Harmony and kBET, respectively, as
well as the Adjusted Rand Index (ARI) and Normalized Mutual Information
(NMI) algorithms.

# Installation

Install with the following:

    library(devtools)
    devtools::install_github('vinay-swamy/scPOP')

Note that to install this package, you may require additional software
to compile Rcpp code:

-   [xcode CLI](https://www.youtube.com/watch?v=Z01lzHNrSdU) for MacOS
-   [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows

# Metrics

The metrics we include are :

-   Adjusted Rand Index(ARI): The amount of overlap between two sets of
    labels.
-   Normalized Mutual Information: The amount of overlap between two
    sets of labels.
-   Silhouette Width: Average difference between within label distance
    and nearest label distance. This requires calculation of distance
    matrix, which does not scale well with larger samplesizes. We
    provide methods for easily subsampling data within `scPOP`
-   Local Inverse Simpson Index (LISI): (avg) Number of labels present
    in the N nearest cells

It’s important to note that based on the type of label being evaluated,
the “optimal” score a given metric may change. For example, when
calculating LISI based on batch, a highscore is better(multiple batches
close together), but when calculating based on Cell Type, a low score is
better( the same celltypes are close together.)

# Usage

The ideal use case for `scPOP` to generate metrics on dataset for which
multiple rounds of batch correction have been calculated. These metrics
can be used to rank different

We provide a toy dataset in .h5ad format

``` r
download.file('https://hpc.nih.gov/~mcgaugheyd/scEiaD/colab/scEiaD_all_anndata_mini_ref.h5ad', 'scEiaD_all_anndata_mini_ref.h5ad')
```

We recommend calculating all metrics at once using `run_all_metrics`.
This function requires a matrix of reduced dimensions, a data.frame
containing metadata, and the names of 3 columns

-   `batch_key`: column corresponds to batch for each cell
-   `label1_key`: primary label for each cell, ie Cell Type
-   `secondary_label`: for each cell, ie Cluster

We recommend the `zellkonverter` for reading .h5ad formatted into R. The
example we provide uses data in the
[SingleCellExperiment](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)
format, but as scPOP only requires vectors of labels and matrices of
reduced dimensions, data from other frameworks like
[Seurat](https://satijalab.org/seurat/) can be easily used.

``` r
library(zellkonverter,quietly = T)
library(SingleCellExperiment,quietly = T)
library(scPOP)
sce <- zellkonverter::readH5AD('scEiaD_all_anndata_mini_ref.h5ad')
sce
```

    ## class: SingleCellExperiment 
    ## dim: 15114 27350 
    ## metadata(0):
    ## assays(1): X
    ## rownames(15114): ENSG00000000005 ENSG00000000419 ... ENSG00000288602
    ##   ENSG00000288642
    ## rowData names(5): vst.mean vst.variance vst.variance.expected
    ##   vst.variance.standardized vst.variable
    ## colnames(27350): CTTTGCGAGATGTGGC_ERS2852885
    ##   GATCGTATCGAGAGCA_ERS2852885 ... ACGATACCAAGCTGTT_SRS6424737
    ##   AACTCAGAGCCCAGCT_SRS6424747
    ## colData names(32): nCount_RNA nFeature_RNA ... doublet_score_scran
    ##   CellType_predict
    ## reducedDimNames(2): X_scvi X_scviumap
    ## altExpNames(0):

Running scPOP

``` r
metrics <-run_all_metrics(reduction = reducedDim(sce, 'X_scvi'), 
                          metadata = colData(sce),
                          batch_key = 'batch',
                          label1_key = 'CellType_predict',
                          label2_key = 'cluster', 
                          run_name = 'example')
```

    ## Calculating LISI...

    ## Done
    ## Calculating Silhoette width...

    ## Done
    ## Calculating ARI...

    ## Done
    ## Calculating NMI...

    ## Done

``` r
metrics
```

    ##       run ari_label nmi_label lisi_batch lisi_CellType_predict lisi_cluster
    ## 1 example 0.4459971 0.6276339   2.400331              1.137047     1.205916
    ##   silWidth_batch silWidth_CellType_predict silWidth_cluster
    ## 1     -0.1615347                 0.2745293         0.138743

These metrics can be applied to multiple integration runs to determine
the optimal integration method/parameters. To illustrate this, we’ll
generate some fake data.

``` r
multi_run_example <-  lapply(c(23232, 23423423, 66774, 2341345, 56733), function(i){
    set.seed(i)
    sce_shuffled <- sce
    sce_shuffled$batch <- sample(sce_shuffled$batch, ncol(sce))
    sce_shuffled$CellType_predict <- sample(sce_shuffled$CellType_predict, ncol(sce))
    sce_shuffled$cluster <- sample(sce_shuffled$cluster, ncol(sce))
    run_all_metrics(reduction = reducedDim(sce_shuffled, 'X_scvi'), 
                  metadata = colData(sce_shuffled),
                  batch_key = 'batch',
                  label1_key = 'CellType_predict',
                  label2_key = 'cluster', 
                  run_name = as.character(i), 
                  sil_width_prop = .25, 
                  sil_width_group_key = 'CellType_predict', 
                  quietly=T)
    
})

run_metrics <-  do.call(rbind,  multi_run_example)
run_metrics
```

    ##        run     ari_label   nmi_label lisi_batch lisi_CellType_predict
    ## 1    23232 -5.892982e-04 0.007564143   8.174771              4.797056
    ## 2 23423423 -3.130038e-04 0.007987631   8.174047              4.806773
    ## 3    66774  9.051521e-06 0.007990235   8.210295              4.819635
    ## 4  2341345 -1.652940e-04 0.008120943   8.187121              4.801828
    ## 5    56733 -6.915233e-05 0.007969419   8.178482              4.804363
    ##   lisi_cluster silWidth_batch silWidth_CellType_predict silWidth_cluster
    ## 1     7.933616    -0.18476700                -0.2956764       -0.3619559
    ## 2     7.925464    -0.16575566                -0.3013349       -0.2507972
    ## 3     7.972977    -0.06813073                -0.2794385       -0.4588756
    ## 4     7.945181    -0.13626901                -0.2239636       -0.3859973
    ## 5     8.006362    -0.17093947                -0.2862585       -0.4382246

We provide a method, `calc_sumZscore`, to aggregate these metrics
together across multiple runs to generate a single score for each run

``` r
run_metrics$sumZscore <-  calc_sumZscore(run_metrics, 'batch')
run_metrics[,c('run', 'sumZscore')]
```

    ##        run sumZscore
    ## 1    23232 -1.711125
    ## 2 23423423  1.283290
    ## 3    66774 -1.835191
    ## 4  2341345  3.599619
    ## 5    56733 -1.336593

We also provide functions for running each function individually

``` r
ari_score <- ari(sce$batch, sce$CellType_predict)
ari_score
```

    ## [1] 0.1371779

``` r
nmi_score <- nmi(sce$batch, sce$CellType_predict)
nmi_score
```

    ## [1] 0.3246491

Lisi requires about 10gb of memory for 100K cells and scales linearly
with number of cells

``` r
lisi_score <- lisi(X = reducedDim(sce, 'X_scvi'), meta_data = as.data.frame(colData(sce)), label_colnames = 'CellType_predict' )
head(lisi_score)
```

    ##                             CellType_predict
    ## CTTTGCGAGATGTGGC_ERS2852885         1.000000
    ## GATCGTATCGAGAGCA_ERS2852885         1.000000
    ## CCTCTGACACCCAGTG_ERS2852885         1.003521
    ## GCGCAACAGAAGCCCA_ERS2852885         1.000000
    ## GACGTTATCGGTCCGA_ERS2852885         1.000000
    ## TGCGCAGGTACCAGTT_ERS2852885         1.000000

For some the `silhouette_width`, a distance matrix mist be calculated,
which requires significant memory usage. We provide the function
`stratified_sample` to downsample data based on a grouping variable

``` r
idx <- stratified_sample(colnames(sce), sce$batch)
sce_ds <- sce[,idx]
sil_score <- silhouette_width(reduction = reducedDim(sce_ds, 'X_scvi'), 
                              meta.data = colData(sce_ds),  
                              keys ='CellType_predict')
```

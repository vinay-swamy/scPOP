#' Compute Local Inverse Simpson's Index (LISI)
#'
#' Use this function to compute LISI scores of one or more labels.
#'
#' @param X A matrix with cells (rows) and features (columns).
#' @param meta_data A data frame with one row per cell.
#' @param label_colnames Which variables to compute LISI for.
#' @param perplexity The effective number of each cell's neighbors.
#' @param nn_eps Error bound for nearest neighbor search with \code{RANN:nn2()}.
#' Default of 0.0 implies exact nearest neighbor search.
#'
#' @return A data frame of LISI values. Each row is a cell and each
#' column is a different label variable.
#'
#' @importFrom RANN nn2
#' @export
#'
#' @examples
#' data(sceiad_subset_data)
#' features <- sceiad_subset_data[, paste0('scviDim_', 1:8)]
#' metadata <- sceiad_subset_data[,c('Barcode', 'cluster',  'subcluster',
#'                                      'CellType', 'CellType_predict')]
#' lisi_scores <- lisi(features, metadata, c('CellType_predict'))
#' head(lisi_scores)
#'

lisi <- function(
    X, meta_data, label_colnames, perplexity = 30, nn_eps = 0
) {
    N <- nrow(meta_data)
    dknn <- nn2(X, k = perplexity * 3, eps = nn_eps)
    lisi_df <- data.frame(matrix(NA, N, length(label_colnames)))
    lisi_df <- Reduce(cbind, lapply(label_colnames, function(label_colname) {
        labels <- data.frame(meta_data)[, label_colname, drop = TRUE]
        if (any(is.na(labels))) {
            message('Cannot compute LISI on missing values')
            return(rep(NA, N))
        } else {
            ## don't count yourself in your neighborhood
            dknn$nn.idx <- dknn$nn.idx[, 2:ncol(dknn$nn.idx)]
            dknn$nn.dists <- dknn$nn.dists[, 2:ncol(dknn$nn.dists)]
            labels <- as.integer(factor(labels)) - 1
            n_batches <- length(unique(labels))
            simpson <- compute_simpson_index(
                t(dknn$nn.dists),
                t(dknn$nn.idx) - 1,
                labels,
                n_batches,
                perplexity
            )
            return(1 / simpson)
        }
    }))
    lisi_df <- as.data.frame(lisi_df)
    colnames(lisi_df) <- label_colnames
    row.names(lisi_df) <- row.names(meta_data)
    return(lisi_df)
}


## Borrowed from https://github.com/jchiquet/aricode

#' @importFrom Matrix sparseMatrix
#'
sortPairs <- function(c1, c2, spMat=FALSE){
    if (anyNA(c1) | anyNA(c2))
        stop("NA are not supported.")

    if (((!is.vector(c1) & !is.factor(c1)) | is.list(c1)) | ((!is.vector(c2) & !is.factor(c2)) | is.list(c2)))
        stop("c1 and c2 must be vectors or factors but not lists.")

    if (length(c1) != length(c2))
        stop("the two vectors must have the same length.")

    n <- length(c1)

    ## if c1 and c2 are integer
    if (is.integer(c1) & is.integer(c2)) {
        ## getRank is O(n) if max(c1)-min(c1) and max(c2)-min(c2) is of order length(c1)=length(c2)
        ## NOTE: getRank does not assume c1 and c2 are between 0 and n
        res1 <- getRank(c1)
        res2 <- getRank(c2)
        mylevels <- list(c1=res1$index, c2=res2$index)
        c1 <- res1$translated  # here ranks are in [0, n)
        c2 <- res2$translated  # here ranks are in [0, n)
    } else if (is.factor(c1) & is.factor(c2)) {
        mylevels <- list(c1 = levels(c1), c2 = levels(c2))
        c1 <- as.integer(c1) - 1L
        c2 <- as.integer(c2) - 1L
    } else {
        ## if neither a factor nor an integer or different of types force to factor then integer
        mylevels <- list(c1 = unique(c1), c2 = unique(c2))
        c1 <- as.integer(factor(c1, levels = mylevels$c1)) - 1L
        c2 <- as.integer(factor(c2, levels = mylevels$c2)) - 1L
    }


    i_order <- order(c1, c2, method="radix") - 1L
    out <- countPairs(c1, c2, i_order)

    if (spMat) {
        spOut <- sparseMatrix(i=out$pair_c1,
                              j=out$pair_c2,
                              x=out$pair_nb,
                              dims=sapply(mylevels,length),
                              dimnames = mylevels, index1=FALSE)
    } else {
        spOut <- NULL
    }

    res <- list(spMat = spOut,
                levels = mylevels,
                nij = out$pair_nb,
                ni. = out$c1_nb,
                n.j = out$c2_nb,
                pair_c1 = out$pair_c1,
                pair_c2 = out$pair_c2
    )
    res
}

#' Adjusted Rand Index
#'
#' A function to compute the adjusted rand index between two classifications
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a scalar with the adjusted rand index.

#' @export
#'
#' @examples
#' ## calculate Adjusted Rand Index on two  sets of labels
#' data(sceiad_subset_data)
#' ari(sceiad_subset_data$CellType_predict, sceiad_subset_data$cluster)
ari <- function(c1, c2){

    ## get pairs using C
    ## ensure that values of c1 and c2 are between 0 and n1
    res <- sortPairs(c1, c2)

    ## get ARI using pairs
    N <- length(c1)

    stot <- sum(choose(res$nij, 2), na.rm=TRUE)
    srow <- sum(choose(res$ni., 2), na.rm=TRUE)
    scol <- sum(choose(res$n.j, 2), na.rm=TRUE)

    expectedIndex <-(srow*scol)/(choose(N,2))
    maximumIndex <- (srow+scol)/2

    if (expectedIndex == maximumIndex & stot != 0) {
        res <- 1
    } else {
        res <- (stot-expectedIndex)/(maximumIndex-expectedIndex)
    }
    res
}


entropy <- function(c1, c2){
    res <- sortPairs(c1, c2)

    N <- length(c1)

    H.UV <- - sum(res$nij * log(res$nij))/N + log(N)
    H.U  <- - sum(res$ni. * log(res$ni.))/N + log(N)
    H.V  <- - sum(res$n.j * log(res$n.j))/N + log(N)

    res <- list(UV = H.UV, U = H.U, V = H.V, sortPairs = res)
    res
}


#' Normalized mutual information (NMI)
#'
#' A function to compute the NMI between two classifications
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @param variant a string in ("max", "min", "sqrt", "sum", "joint"): different variants of NMI. Default use "max".
#' @return a scalar with the normalized mutual information .
#'
#' @export
#' @examples
#'  ## calculate Normalized Mutal Information score for two sets of labels
#'  data(sceiad_subset_data)
#'  nmi(sceiad_subset_data$CellType_predict, sceiad_subset_data$cluster)
nmi <- function(c1, c2, variant = c("max", "min", "sqrt", "sum", "joint")) {

    variant <- match.arg(variant)

    H   <- entropy(c1,c2)

    MI  <- - H$UV + H$U + H$V

    D <- switch(variant,
                "max" = max(H$U, H$V),
                "sqrt" = sqrt(H$U * H$V),
                "min" = min(H$U, H$V),
                "sum" = .5*(H$U + H$V),
                "joint" = H$UV)
    res <- MI / D
    res
}

#' batch_sil
#'
#' @description
#'     Determine batch/bio effect using the silhouette
#'     coefficient (adopted from scone):
#' @param reduction a matrix of reduced dimensions
#' @param meta.data dataframe with meta.data associated with reduction
#' @param keys columns in meta.data to calculate silhoette for
#'     to use (default: all)
#' @return
#'     The average silhouette width for all clusters.
#'     For batch effect, the smaller the better.
#'     For biological effect, the larger the better.
#'
#' @importFrom cluster silhouette
#' @importFrom stats dist
#' @export
#'
#' @examples
#'
#' ## calculate the the silhoeuette width score on two sets of labels
#' ## NOTE: this requires computation of a distance matrix, so does not
#' ##       scale well to large datasets
#' data(sceiad_subset_data)
#' features <- sceiad_subset_data[, paste0('scviDim_', 1:8)]
#' metadata <- sceiad_subset_data[,c('Barcode', 'cluster',
#'               'subcluster', 'CellType', 'CellType_predict')]
#' silhouette_width(features, metadata, 'CellType_predict')
#'
silhouette_width <- function(reduction, meta.data, keys){
    # in scone, they use svd to compute principal components.
    # For now, we'll keep the PCA object created in
    # prcomp to be consistent
    # proj <- svd(scale(t(expr),center = TRUE,scale = TRUE),
    #             nu = 3, nv =0)$u

    dd <- dist(reduction)
    sil_width <- sapply(keys, function(x){
        summary(cluster::silhouette(as.numeric(as.factor(meta.data[[x]])),
                                    dd))$avg.width
    }   )
    names(sil_width) <- keys
    return(sil_width)
}

#' Calc_sumZscore
#' @description
#'     Aggregate multiple integration metrics across multiple integration runs, ie from different
#'     batch correction algorithms, or different parameters for the same algorithms
#'
#' @param metric_df_list a list of data.frames generated by applying \code{\link{run_all_metrics}} to multiple
#'     sets of integrations
#' @param batch_key name of batch column in metadata used when generating run_all_metrics
#'
#' @return a vector of aggregated, z-scored metrics
#' @export
#'
#' @examples
#' \donttest{
#' library(scPOP)
#' data(sceiad_subset_data)
#'
#' features <- sceiad_subset_data[, paste0('scviDim_', 1:8)]
#' metadata_1 <- sceiad_subset_data[,c('Barcode', 'cluster',  'subcluster',
#'                                   'batch', 'CellType', 'CellType_predict')]
#'
#' ## scramble example dataset to generate multiple integration runs
#' metadata_2 <- metadata_1
#' metadata_2$batch <- sample(metadata_2$batch, length(metadata_2$batch))
#' metadata_2$CellType_predict <- sample(metadata_2$CellType_predict,
#'                                    length(metadata_2$CellType_predict))
#' metadata_2$cluster <- sample(metadata_2$cluster, length(metadata_2$cluster))
#'
#' metadata_3 <- metadata_1
#' metadata_3$batch <- sample(metadata_3$batch, length(metadata_3$batch))
#' metadata_3$CellType_predict <- sample(metadata_3$CellType_predict,
#'                                   length(metadata_3$CellType_predict))
#' metadata_3$cluster <- sample(metadata_3$cluster, length(metadata_3$cluster))
#' integration_data_list <- list( metadata_1, metadata_2, metadata_3)
#' metric_df_list <- lapply(integration_data_list, function(x)
#'                           run_all_metrics(reduction = features,
#'                           metadata = x,
#'                           batch_key = 'batch',
#'                           label1_key = 'CellType_predict',
#'                           label2_key = 'cluster',
#'                           run_name = 'example',
#'                           quietly =TRUE
#'                           )
#'                         )
#'
#'  calc_sumZscore(metric_df_list,'batch' )
#'  }
calc_sumZscore <- function(metric_df_list, batch_key){
    stopifnot( 'Must have at least two sets of metrics to compare' = length(metric_df_list)>1)
    df <- do.call(rbind, metric_df_list)
    id <- df[,1]
    df <- df[,-1]
    batch_cols <- colnames(df)[grepl(batch_key, colnames(df))]
    label_cols <- colnames(df)[!grepl(batch_key, colnames(df))]
    negative_cases <- c(  batch_cols[grepl('silWidth', batch_cols)],
                          label_cols[grepl('lisi', label_cols)]
        )
    positive_cases <- c(  batch_cols[!grepl('silWidth', batch_cols)],
                          label_cols[!grepl('lisi', label_cols)]
    )
    nc_scaled <- apply(df[,negative_cases], 2, scale) * -1
    pc_scaled <- apply(df[,positive_cases], 2, scale)
    pc_scaled[is.nan(nc_scaled)] <- 0
    nc_scaled[is.nan(nc_scaled)] <- 0
    sumZ <- rowSums(cbind(pc_scaled, nc_scaled))
    return(sumZ)
}

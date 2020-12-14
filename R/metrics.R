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
#'
#' ## Example with 400 cells.
#' library(lisi)
#' library(dplyr)
#' library(tidyr)
#' library(ggplot2)
#' library(magrittr)
#'
#' head(lisi::meta_data)
#'
#' ## Let's color cells by labels. For label 1, there are mixed and non-mixed
#' ## groups. For label 2, all cells are well mixed.
#' lisi::X %>%
#'   cbind(lisi::meta_data) %>%
#'   sample_frac(1L, FALSE) %>%
#'   gather(key, val, label1, label2) %>%
#'   ggplot(aes(X1, X2, color = val)) +
#'     geom_point(shape = 21) +
#'     facet_wrap(~key)
#'
#' ## Now to compute and plot the LISI values for each label.
#' lisi_res <- compute_lisi(lisi::X, lisi::meta_data, c('label1', 'label2'))
#' head(lisi_res)
#'
#' lisi::X %>%
#'   cbind(lisi_res) %>%
#'   sample_frac(1L, FALSE) %>%
#'   gather(key, lisi_value, label1, label2) %>%
#'   ggplot(aes(X1, X2, color = lisi_value)) +
#'     geom_point(shape = 21) +
#'     facet_wrap(~key)
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
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' ARI(cl,iris$Species)
#' @export
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
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' NMI(cl,iris$Species)
#' @export
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
#' @param pca.data a list as created by \code{prcomp},
#'     \code{batch_sil.R} needs \code{$x}:
#'     the principal components (PCs,
#'     correctly: the rotated data)
#' @param batch vector with the batch
#'     covariate (for each cell)
#' @param nPCs the number of principal components
#'     to use (default: 3)
#' @return
#'     The average silhouette width for all clusters.
#'     For batch effect, the smaller the better.
#'     For biological effect, the larger the better.
#' @examples
#'     testdata <- create_testset_multibatch(n.genes=1000,
#'                   n.batch=3, plattform='any')
#'     pca.data <- prcomp(testdata$data, center=TRUE)
#'     batch.silhouette <- batch_sil(pca.data, testdata$batch)
#'
#' @importFrom cluster silhouette
#' @importFrom stats dist
#' @export
silhouette_width <- function(pca.data, batch, nPCs = 3){
    # in scone, they use svd to compute principal components.
    # For now, we'll keep the PCA object created in
    # prcomp to be consistent
    # proj <- svd(scale(t(expr),center = TRUE,scale = TRUE),
    #             nu = 3, nv =0)$u

    dd <- as.matrix(dist(pca.data$x[, seq_len(nPCs)]))
    summary(silhouette(as.numeric(batch), dd))$avg.width
}


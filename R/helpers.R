#' Generate a stratified subsample for a vector given a grouping
#'
#' Use this function to compute LISI scores of one or more labels.
#'
#' @param indexer A vector containing cell barcodes/labels to subsample
#' @param grouping A vector containg a groups to stratify by ( same size as indexer)
#' @param sample_proportion proportion to sample data (default: .1)
#' @param min_count Minimum number of samples in a group to keep
#' @param seed seed value for set.seed
#'
#' @return A subsampled vector generated from indexer
#' @export
stratified_sample <- function(indexer, grouping,sample_proportion=.1, min_count=0, seed=424242){
    df <- data.frame(indexer=indexer, label=grouping)
    dfl <- split(df, df$label)
    lengths <- sapply(dfl, nrow)
    keep <- lengths > min_count
    dfl <- dfl[keep]
    set.seed(seed)
    indices <- sapply(dfl, function(x) {
        size <- as.integer(nrow(x)*sample_proportion)
        sample(x$indexer,size)})
    names(indices) <- NULL
    idx <- unlist(indices)
    return(idx)
}

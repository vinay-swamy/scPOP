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
#'
#' @examples
#' data(sceiad_subset_data)
#' rownames(sceiad_subset_data) <- sceiad_subset_data$Barcode
#' res  = stratified_sample(sceiad_subset_data$Barcode, sceiad_subset_data$cluster)
#' dim(sceiad_subset_data[res, ])
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

qmessage <- function(msg, quiet){
    if(!quiet) message(msg)
}

#' Running All Metrics
#' @param reduction A matrix of reduced dimensions
#' @param metadata A data.frame containing information like batch, cell type, etc
#' @param batch_key Name of column in metadata corresponding to batch
#' @param label1_key Name of column in metadata corresponding to primary cell label, eg Cell type
#' @param label2_key Name of column in metadata corresponding to secondary cell label, eg cluster identity
#' @param sil_width_prop (optional) proportion of data to use for \code{\link{silhouette_width}}
#' @param sil_width_group_key (optional) which column in metadata to use for stratified sampling of data
#' @param run_name (optional) name to refer to dataset
#' @param quietly (optional) if TRUE dont print anything
#'
#' @return A one row \code{data.frame} of calculated metrics
#' @export
run_all_metrics <- function(reduction, metadata, batch_key, label1_key, label2_key, run_name=NULL,
                            sil_width_prop=1, sil_width_group_key=NULL, quietly=F){
    if(is.null(run_name))run_name <- sample(letters, 12, replace = T)
    metadata <- as.data.frame( metadata)
    keys <- c(batch_key, label1_key, label2_key)
    qmessage('Calculating LISI...',quietly)
    lisi <- lisi(reduction, metadata, keys )
    lisi <- lapply(lisi, mean)
    names(lisi) <- paste0('lisi_', keys)
    lisi <- as.data.frame(lisi)
    qmessage('Done\nCalculating Silhoette width...',quietly)
    if(sil_width_prop < 1){
        if(is.null(sil_width_group_key) ) sil_width_group_key <- label1_key
        idx <- stratified_sample(rownames(reduction), metadata[[sil_width_group_key]])
        rd_ds <- reduction[idx, ]
        md_ds <- metadata[idx, ]
        sw <- silhouette_width(rd_ds, md_ds, keys)
    } else{
        sw <- silhouette_width(reduction, metadata, keys)
    }
    names(sw) <- paste0('silWidth_', keys)
    sw <- as.data.frame(t(sw))
    qmessage('Done\nCalculating ARI...',quietly)
    ari_batch <- ari(metadata[[batch_key]], metadata[[label1_key]])
    ari_label <- ari(metadata[[label1_key]], metadata[[label2_key]])
    qmessage('Done\nCalculating NMI...',quietly)
    nmi_batch <- nmi(metadata[[batch_key]], metadata[[label1_key]])
    nmi_label <- nmi(metadata[[label1_key]], metadata[[label2_key]])
    qmessage('Done',quietly)
    scores <- data.frame(run=run_name,
               #ari_batch= ari_batch,
               ari_label= ari_label,
               #nmi_batch= nmi_batch,
               nmi_label= nmi_label)
    scores <- do.call(cbind, list( scores, lisi, sw) )
    return(scores)


}


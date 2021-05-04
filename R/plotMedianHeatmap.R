#' @title Visualizes the median pixel intensities per spot and channel
#'
#' @description Helper function for estimating the spillover matrix. This
#' function visualizes the median pixel intensities per spot (columns) and per
#' channel (rows).
#'
#' @param object a \code{SingleCellExperiment} object containing pixel
#' intensities per channel. Individual pixels are stored as columns and
#' channels are stored as rows.
#' @param spot_id character string indicating which \code{colData(object)} entry
#' stores the isotope names of the spotted metal. Entries should be of the 
#' form (mt)(mass) (e.g. Sm152 for Samarium isotope with the atomic mass 152).
#' @param channel_id character string indicating which \code{rowData(object)} 
#' entry contains the isotope names of the acquired channels. 
#' @param assay_type character string indicating which assay to use.
#' @param statistics the statistic to use when aggregating channels per spot
#' @param log should the aggregated pixel intensities be \code{log10(x + 1)} 
#' transformed?
#' @param theshold single numeric indicating the aggregated pixel threshold.
#' This facilitates the identification of spots with low aggregated intensities. 
#' @param order_metals should the metals be ordered based on spotted mass?
#' @param ... arguments passed to \code{pheatmap}.
#'
#' @return returns a SCE object
#'
#' @examples
#'
#' @importFrom pheatmap pheatmap
#' @importFrom scuttle aggregateAcrossCells
#' @importFrom viridis viridis
#' @importFrom stringr str_extract
#' @export
plotSpotHeatmap <- function(object, 
                            spot_id = "sample_id",
                            channel_id = "channel_name",
                            assay_type = "counts",
                            statistic = "median",
                            log = TRUE,
                            threshold = NULL,
                            order_metals = TRUE,
                            ...){
    
    .validSpotHeatmapInput(object, spot_id, channel_id, assay_type, 
                           log, threshold, order_metals)
    
    cur_out <- aggregateAcrossCells(object, object[[spot_id]], 
                                    statistics = statistic,
                                    use.assay.type = assay_type)
    
    if (log) {
        cur_mat <- log10(assay(cur_out, assay_type) + 1)
    } else {
        cur_mat <- assay(cur_out, assay_type)
    }
    
    if (!is.null(threshold)) {
        cur_mat <- (cur_mat > threshold) * 1
    }
    
    colnames(cur_mat) <- cur_out[[spot_id]] 
    rownames(cur_mat) <- rowData(cur_out)[[channel_id]] 
    
    # Overwrite default colour
    args <- list(...)
    
    if ("color" %in% names(args)) {
        cur_col <- args$color 
    } else {
        cur_col <- viridis(100)
    }
    
    # Order rows and cols based on spot metal
    if (order_metals) {
        cur_spots <- colnames(cur_mat)
        cur_mass <- as.numeric(str_extract(cur_spots, "[0-9]{2,3}$"))
        cur_spots <- cur_spots[order(cur_mass)]
        
        cur_channels <- rownames(cur_mat)
        cur_mass <- as.numeric(str_extract(cur_channels, "[0-9]{2,3}"))
        cur_channels <- cur_channels[order(cur_mass)]
        cur_isotope <- str_extract(cur_channels, "[A-Za-z]{1,2}[0-9]{2,3}")
        cur_rownames <- c(cur_channels[cur_isotope %in% cur_spots],
                          cur_channels[!cur_isotope %in% cur_spots])
        
        cur_mat <- cur_mat[cur_rownames,cur_spots]
        
        cur_cluster_cols <- FALSE
        cur_cluster_rows <- FALSE
    } else {
        cur_cluster_cols <- ifelse("cluster_cols" %in% names(args), args$cluster_cols, TRUE)
        cur_cluster_rows <- ifelse("cluster_rows" %in% names(args), args$cluster_rows, TRUE)
    }
    
    pheatmap(cur_mat, color = cur_col, 
             cluster_cols = cur_cluster_cols,
             cluster_rows = cur_cluster_rows)
}

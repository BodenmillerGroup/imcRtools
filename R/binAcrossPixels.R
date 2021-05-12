#' @title Aggregate consecutive pixels per spot
#'
#' @description Helper function for estimating the spillover matrix. Per metal 
#' spot, consecutive pixels a aggregated (default: summed).
#'
#' @param object a \code{SingleCellExperiment} object containing pixel
#' intensities per channel. Individual pixels are stored as columns and
#' channels are stored as rows.
#' @param bin_size single numeric indicating how many consecutive pixels should
#' be aggregated.
#' @param spot_id character string indicating which \code{colData(object)} entry
#' stores the isotope names of the spotted metal. 
#' @param assay_type character string indicating which assay to use.
#' @param statistic character string indicating the statistic to use for
#' aggregating consecutive pixels.
#' @param ... additional arguments passed to \code{aggregateAcrossCells}
#'
#' @return returns a SCE object
#'
#' @examples
#'
#' @importFrom scuttle aggregateAcrossCells
#' @export
binAcrossPixels <- function(object, 
                            bin_size,
                            spot_id = "sample_id",
                            assay_type = "counts",
                            statistic = "sum",
                            ...){
    
    .valid.binAcrossPixels.input(object, bin_size, spot_id, assay_type)
    
    cur_split <- split(object[[spot_id]], f = object[[spot_id]])
    cur_split <- lapply(cur_split, function(x){ceiling(seq_along(x)/bin_size)})
    
    cur_df <- DataFrame(spot_id = object[[spot_id]],
                         bin = unlist(cur_split))
    
    cur_out <- aggregateAcrossCells(object, cur_df, 
                                    statistics = statistic,
                                    use.assay.type = assay_type,
                                    ...)
    
    return(cur_out)
    
}

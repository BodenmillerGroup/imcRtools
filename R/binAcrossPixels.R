#' @title Aggregate consecutive pixels per single-metal spot
#'
#' @description Helper function for estimating the spillover matrix. Per metal 
#' spot, consecutive pixels a aggregated (default: summed).
#'
#' @param object a \code{SingleCellExperiment} object containing pixel
#' intensities for all channels. Individual pixels are stored as columns and
#' channels are stored as rows.
#' @param bin_size single numeric indicating how many consecutive pixels per
#' spot should be aggregated.
#' @param spot_id character string indicating which \code{colData(object)} entry
#' stores the isotope names of the spotted metal. 
#' @param assay_type character string indicating which assay to use.
#' @param statistic character string indicating the statistic to use for
#' aggregating consecutive pixels.
#' @param ... additional arguments passed to \code{aggregateAcrossCells}
#'
#' @return returns the binned pixel intensities in form of a 
#' \code{SingleCellExperiment} object
#'
#' @examples
#' path <- system.file("extdata/spillover", package = "imcRtools")
#' # Read in .txt files
#' sce <- readSCEfromTXT(path)
#' dim(sce)
#' 
#' # Visualizes heatmap before aggregation
#' plotSpotHeatmap(sce)
#' 
#' # Sum consecutive pixels
#' sce <- binAcrossPixels(sce, bin_size = 10)
#' dim(sce)
#' 
#' # Visualizes heatmap after aggregation
#' plotSpotHeatmap(sce)
#' 
#' @seealso \code{\link[scuttle]{aggregateAcrossCells}} for the aggregation
#' function
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
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
    
    if (!statistic %in% c("sum", "mean", "median")) {
        stop("'statistic' must be 'sum', 'mean' or 'median'.")
    }
    
    cur_split_tmp <- split(object[[spot_id]], f = object[[spot_id]])
    cur_split <- lapply(cur_split_tmp, function(x){ceiling(seq_along(x)/bin_size)})
    
    if (!isTRUE(all.equal(as.vector(unlist(cur_split_tmp)), object[[spot_id]]))) {
        stop("Spot IDs of pixels within 'object' are not ordered alphabetically.")
    }
    
    cur_df <- DataFrame(spot_id = object[[spot_id]],
                        bin = unlist(cur_split))
    
    cur_out <- aggregateAcrossCells(object, cur_df, 
                                    statistics = statistic,
                                    use.assay.type = assay_type,
                                    ...)
    
    return(cur_out)
    
}

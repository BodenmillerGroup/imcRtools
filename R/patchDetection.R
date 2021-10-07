#' @title Function to detect patches containing defined cell types
#'
#' @description Function to detect spatial clusters of defined types of cells.
#' By defining a certain distance threshold, all cells within the vicinity
#' of these clusters are detected as well.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param colPairName single character indicating the \code{colPair(object)}
#' entry containing the neighbor information.

#' @return
#' 
#' @examples
#'   
#' @author Tobias Hoch 
#' @author adapted by Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#'
#' @export
patchDetection <- function(object, 
                           pattern,
                           colPairName,
                           expand_by = 0,
                           min_patch_size = 1,
                           name = "patch_id"){
    
    # .valid.patchDetection.input(object, colPairName, expand_by)
    
    cur_graph <- graph_from_edgelist(as.matrix(colPair(object[,pattern], colPairName)))
    cur_clusters <- components(cur_graph)$membership
    cur_out <- vector(mode = "character", length = ncol(object))
    cur_out[!pattern] <- NA
    cur_out[pattern] <- cur_clusters
    colData(object)[[name]] <- cur_out
    
    return(object)
}

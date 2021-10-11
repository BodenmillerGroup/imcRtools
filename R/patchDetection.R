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
                           patch_cells,
                           colPairName,
                           min_patch_size = 1,
                           name = "patch_id",
                           expand_by = 0,
                           coords = c("Pos_X", "Pos_Y"),
                           convex = FALSE,
                           img_id = NULL,
                           BPPARAM = SerialParam()){
    
    .valid.patchDetection.input(object, patch_cells, colPairName, 
                                min_patch_size, name, expand_by, coords,
                                convex, img_id)
    
    cur_graph <- graph_from_edgelist(as.matrix(colPair(object[,patch_cells], 
                                                       colPairName)))
    cur_components <- components(cur_graph)
    cur_clusters <- cur_components$membership
    
    if (min_patch_size > 1) {
        cur_clusters[!(cur_clusters %in% 
                           which(cur_components$csize >= min_patch_size))] <- NA
    } 
    
    cur_out <- vector(mode = "character", length = ncol(object))
    cur_out[!patch_cells] <- NA
    cur_out[patch_cells] <- cur_clusters
    colData(object)[[name]] <- cur_out
    
    if (expand_by > 0) {
        object <- .expand_patch(object, name = name, 
                                expand_by = expand_by, 
                                coords = coords,
                                convex = convex,
                                img_id = img_id,
                                BPPARAM = BPPARAM)
    }
    
    return(object)
}

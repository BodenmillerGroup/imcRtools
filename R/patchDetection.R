#' @title Function to detect patches containing defined cell types
#'
#' @description Function to detect spatial clusters of defined types of cells.
#' By defining a certain distance threshold, all cells within the vicinity
#' of these clusters are detected as well.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param patch_cells logical vector of length equal to the number of cells
#' contained in \code{object}. \code{TRUE} entries define the cells to consider
#' for patch detection (see Details).
#' @param colPairName single character indicating the \code{colPair(object)}
#' entry containing the neighbor information.
#' @param min_patch_size single integer indicating the minimum number of 
#' connected cells that make up a patch before expansion.
#' @param name single character specifying the \code{colData} entry storing
#' the patch IDs in the returned object.
#' @param expand_by single numeric indicating in which vicinity range cells
#' should be considered as belonging to the patch (see Details).
#' @param coords character vector of length 2 specifying the names of the
#' \code{colData} (for a \code{SingleCellExperiment} object) or the
#' \code{spatialCoords} entries of the cells' x and y locations.
#' @param convex should the convex hull be computed before expansion? Default:
#' the concave hull is computed.
#' @param img_id single character indicating the \code{colData(object)} entry
#' containing the unique image identifiers.
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam-class}} object
#' defining how to parallelize computations.
#' 
#' @section Detecting patches of defined cell types:
#' This function works as follows:
#' 
#' 1. Only cells defined by \code{patch_cells} are considered for patch 
#' detection.
#' 
#' 2. Patches of connected cells are detected. Here, cell-to-cell connections
#' are defined by the interaction graph stored in 
#' \code{colPair(object, colPairName)}. At this point, patches that contain 
#' fewer than \code{min_patch_size} cells are removed.
#' 
#' 3. If \code{expand_by > 0}, a concave (default) or convex hull is constructed
#' around each patch. This is is then expanded by \code{expand_by} and cells
#' within the expanded hull are detected and assigned to the patch. This 
#' expansion only works if a patch contains at least 3 cells.
#' 
#' The returned object contains an additional entry 
#' \code{colData(object)[[name]]}, which stores the patch ID per cell. \code{NA}
#' indicate cells that are not part of a patch. 
#' 
#' @return An object of \code{class(object)} containing a patch ID for each 
#' cell in \code{colData(object)[[name]]}. 
#' 
#' @examples
#' library(cytomapper)
#' data(pancreasSCE)
#' 
#' # Visualize cell types
#' plotSpatial(pancreasSCE, img_id = "ImageNb", node_color_by = "CellType")
#'
#' # Build interaction graph
#' pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
#'                                  type = "expansion", threshold = 20)
#' 
#' # Detect patches of "celltype_B" cells
#' pancreasSCE <- patchDetection(pancreasSCE, 
#'                               patch_cells = pancreasSCE$CellType == "celltype_B",
#'                               colPairName = "expansion_interaction_graph")
#'                               
#' plotSpatial(pancreasSCE, img_id = "ImageNb", node_color_by = "patch_id")
#' 
#' # Include cells in vicinity
#' pancreasSCE <- patchDetection(pancreasSCE, 
#'                               patch_cells = pancreasSCE$CellType == "celltype_B",
#'                               colPairName = "expansion_interaction_graph",
#'                               expand_by = 20, 
#'                               img_id = "ImageNb")
#' 
#' plotSpatial(pancreasSCE, img_id = "ImageNb", node_color_by = "patch_id")
#'   
#' @author Tobias Hoch 
#' @author adapted by Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#' 
#' @references
#' \href{https://www.biorxiv.org/content/10.1101/2021.07.29.454093v1}{
#' Hoch, T. et al., Multiplexed Imaging Mass Cytometry of Chemokine Milieus in 
#' Metastatic Melanoma Characterizes Features of Response to Immunotherapy., 
#' bioRxiv 2021}
#' 
#' @importFrom igraph graph_from_edgelist components 
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

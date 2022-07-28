#' @title Function to compute the area of each patch
#'
#' @description This function constructs a polygons around patch cells and
#' computes their area.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param pach_name single character indicating the \code{colData(object)} entry
#' containing the patch cell identifiers.
#' @param coords character vector of length 2 specifying the names of the
#' \code{colData} (for a \code{SingleCellExperiment} object) or the
#' \code{spatialCoords} entries of the cells' x and y locations.
#' @param convex should the convex hull be computed before expansion? Default:
#' the concave hull is computed.
#' 
#' @return A DataFrame object containing the patch identifier, the constructed
#' polygon and the polygon size.
#' 
#' @examples
#' library(cytomapper)
#' data(pancreasSCE)
#'
#' # Build interaction graph
#' pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
#'                                  type = "expansion", threshold = 20)
#' 
#' # Detect patches of "celltype_B" cells
#' pancreasSCE <- patchDetection(pancreasSCE, 
#'                               patch_cells = pancreasSCE$CellType == "celltype_B",
#'                               expand_by = 5, img_id = "ImageNb",
#'                               colPairName = "expansion_interaction_graph")
#'                               
#' # Compute the patch area
#' patchSize(pancreasSCE)
#'   
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#' 
#' @importFrom sf st_area
#' @export
patchSize <- function(object, 
                      patch_name = "patch_id",
                      coords = c("Pos_X", "Pos_Y"),
                      convex = FALSE){
    
    data <- polygon <- NULL
    
    if (is(object, "SpatialExperiment")) {
        out <- cbind(colData(object), spatialCoords(object)) %>% 
            as_tibble %>%
            filter(!is.na(!!sym(patch_name))) %>%
            nest_by(!!sym(patch_name)) %>%
            summarize(
                polygon = list(.polygon_function(x = data,
                                                 coords = coords,
                                                 convex = convex)),
                size = ifelse(is.na(polygon), NA, st_area(polygon)))
    } else {
        out <- colData(object) %>% as_tibble %>%
            filter(!is.na(!!sym(patch_name))) %>%
            nest_by(!!sym(patch_name)) %>%
            summarize(polygon = list(.polygon_function(x = data,
                                                       coords = coords,
                                                       convex = convex)),
                      size = ifelse(is.na(polygon), NA, st_area(polygon)))
    }
    
    return(as(out, "DataFrame"))
}

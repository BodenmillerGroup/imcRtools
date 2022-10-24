#' @title Function to calculate minimal distance to cells of interest 
#' 
#' @description Function to return the distance of the closest cell of interest
#' for each cell in the data. In the case of patched/clustered cells negative
#' distances are returned by default which indicate the distance of the cells
#' of interest to the closest cell that is not of the type of cells of
#' interest.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param x_cells logical vector of length equal to the number of cells
#' contained in \code{object}. \code{TRUE} entries define the cells to which
#' distances will be calculated.
#' @param name character specifying the name of the \code{colData} entry to safe
#' the distances in.
#' @param coords character vector of length 2 specifying the names of the
#' \code{colData} (for a \code{SingleCellExperiment} object) or the
#' \code{spatialCoords} entries of the cells' x and y locations.
#' @param img_id single character indicating the \code{colData(object)} entry
#' containing the unique image identifiers.
#' @param return_neg logical indicating whether negative distances are to be
#' returned for the distances of patched/spatially clustered cells.
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam-class}} object
#' defining how to parallelize computations.
#' 
#' @return returns an object of \code{class(object)} containing a new column 
#' entry to \code{colData(object)[[name]]}.
#' 
#' @examples
#' library(cytomapper)
#' data(pancreasSCE)
#' 
#' # Build interaction graph
#' pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb",
#' type = "expansion",threshold = 20)
#' 
#' # Detect patches of "celltype_B" cells
#' pancreasSCE <- patchDetection(pancreasSCE,
#'                              img_id = "ImageNb",
#'                              patch_cells = pancreasSCE$CellType == "celltype_B",
#'                              colPairName = "expansion_interaction_graph",
#'                              min_patch_size = 20,
#'                              expand_by = 1)
#'
#' plotSpatial(pancreasSCE, 
#'             img_id = "ImageNb", 
#'             node_color_by = "patch_id",
#'             scales = "free")
#'
#' # Distance to celltype_B patches
#' pancreasSCE <- minDistToCells(pancreasSCE,
#'                              x_cells = !is.na(pancreasSCE$patch_id),
#'                              coords = c("Pos_X","Pos_Y"),
#'                              img_id = "ImageNb")
#'
#' plotSpatial(pancreasSCE,
#'             img_id = "ImageNb",
#'             node_color_by = "distToCells",
#'             scales = "free")
#'
#' @author Daniel Schulz (\email{daniel.schulz@@uzh.ch})
#' @importFrom distances distances distance_columns
#' @importFrom MatrixGenerics rowMins
#' @export
minDistToCells <- function(object,
                           x_cells,
                           img_id,
                           name = "distToCells",
                           coords = c("Pos_X","Pos_Y"),
                           return_neg = TRUE,
                           BPPARAM = SerialParam()){
  
  .valid.minDistToCells.input(object,x_cells,name,coords,img_id,return_neg)
  
  cur_meta <- metadata(object)
  metadata(object) <- list()
  
  cur_intmeta <- int_metadata(object)
  
  object$x_cells <- x_cells

  cur_out <- bplapply(
    unique(colData(object)[[img_id]]),
    function(x){
      
      cur_obj <- object[,as.character(colData(object)[[img_id]]) == x]
      
      cur_obj[[name]] <- NA
      if (sum(cur_obj$x_cells) == 0) {
        return(cur_obj)
      }
      
      patch_cells <- which(cur_obj$x_cells)
      non_patch_cells <- which(!cur_obj$x_cells)
      
      if (is(object, "SpatialExperiment")) {
        dist_mat <- distances(spatialCoords(cur_obj))
      } else {
        dist_mat <- distances(as.matrix(colData(cur_obj)[,coords]))
      }
      
      pos_dist <- distance_columns(dist_mat,column_indices = patch_cells)
      dist_to_patch <- rowMins(pos_dist)
      neg_dist <- distance_columns(dist_mat,column_indices = non_patch_cells)
      dist_from_patch <- rowMins(neg_dist)
      
      # cells that had a 0 distance to the cells of interest can be substitutes
      # with the negative distances from the cells of interest
      if(return_neg == TRUE) {
        dist_to_patch[dist_to_patch == 0] <- -dist_from_patch[dist_to_patch == 0]
      }
      cur_obj[[name]] <- dist_to_patch
      
      return(cur_obj)
    }, BPPARAM = BPPARAM)
  
  cur_out <- do.call("cbind", cur_out)
  
  cur_out$x_cells <- NULL
  
  metadata(cur_out) <- cur_meta
  int_metadata(cur_out) <- cur_intmeta
  
  return(cur_out)
}
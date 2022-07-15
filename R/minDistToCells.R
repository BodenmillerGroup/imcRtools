#' @title Function to calculate minimal distance of cells of interest 
#' 
#' @description Function to return the distance of the closest cell of interest
#'   for each cell in the data. In the case of patched/clustered cells negative
#'   distances can be returned which indicate the closest cell that is not of
#'   the type of cells of interest.
#' 
#' function to calculate the distance to the border of patches of cells.
# positive values are indicative of cells that are outside of an object and negative values are indicative of cells inside of an object.

#' @param object single cell object
#' @param x_cells logical vector of length equal to the number of cells
#' contained in \code{object}. \code{TRUE} entries define the cells to consider
#' for patch detection (see Details).
#' @param name character specifying the name of the colData entry to safe the distances in.
#' @param coords character vector of length 2 specifying the names of the \code{colData} (for a \code{SingleCellExperiment} object) or the \code{spatialCoords} entries of the cells' x and y locations.
#' @param img_id single character indicating the \code{colData(object)} entry containing the unique image identifiers.
#' @param return_neg logical indicating whether negative distances are to be returned for the distances of patched/spatially clustered cells.
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam-class}} object defining how to parallelize computations.
minDistToCells <- function(object,
                           x_cells,
                           name,
                           coords,
                           img_id,
                           return_neg = TRUE,
                           BPPARAM = SerialParam()){
  
  cur_meta <- metadata(object)
  metadata(object) <- list()
  
  cur_intmeta <- int_metadata(object)
  
  object$x_cells <- x_cells
  
  #cells <- colnames(object[,x_cells])
  #other_cells <- colnames(object[,!x_cells])
  
  cur_out <- bplapply(
    unique(colData(object)[[img_id]]),
    function(x){
      
      
      # get one image and all cells within the mask
      cur_obj <- object[,as.character(colData(object)[[img_id]]) == x]
      
      cur_obj[[name]] <- NA
      if (sum(cur_obj$x_cells) == 0) {
        return(cur_obj)
      }
      
      # get cells of interest and other cells
      patch_cells <- which(cur_obj$x_cells)
      non_patch_cells <- which(!cur_obj$x_cells)
      
      # calculate the distances for all against all cells
      dist_mat <- distances(as.matrix(colData(cur_obj)[,coords]))
      
      # select only those columns (cells) that are part of the patch
      pos_dist <- distances::distance_columns(dist_mat,column_indices = patch_cells)
      # for each row (cell) get the minimal distance to a cell of the patch (columns)
      dist_to_patch <- rowMins(pos_dist)
      
      # select only those columns (cells) that are NOT part of the patch
      neg_dist <- distances::distance_columns(dist_mat,column_indices = non_patch_cells)
      # # for each row (cell) get the minimal distance to a cell NOT part of the patch (columns)
      dist_from_patch <- rowMins(neg_dist)
      
      # cells that had a 0 distance to the patch can be substitutes with the negative distances from the patch
      if(return_neg == TRUE) {
        dist_to_patch[dist_to_patch == 0] <- -dist_from_patch[dist_from_patch != 0]
      }
      cur_obj[[name]] <- dist_to_patch
      
      return(cur_obj)
    }, BPPARAM = BPPARAM)
  
  cur_out <- do.call("cbind", cur_out)
  
  metadata(cur_out) <- cur_meta
  int_metadata(cur_out) <- cur_intmeta
  
  return(cur_out)
}

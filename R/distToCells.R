#' @title Function to calculate distance to cells of interest
#'
#' @description Function to return the min, max, mean or median distance to the cells of interest
#' for each cell in the data. In the case of patched/clustered cells negative
#' distances are returned by default which indicate the distance of the cells
#' of interest to the cells that are not of the type of cells of
#' interest.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param x_cells logical vector of length equal to the number of cells
#' contained in \code{object}. \code{TRUE} entries define the cells to which
#' distances will be calculated.
#' @param name character specifying the name of the \code{colData} entry to safe
#' the distances in.
#' @param statistics one of "min", "max", "mean" or "meadian" specifying the distance statistics to use when computing
#' the distances.
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
#' @section Ordering of the output object:
#' The \code{minDistToCells} function operates on individual images.
#' Therefore the returned object is grouped by entries in \code{img_id}.
#' This means all cells of a given image are grouped together in the object.
#' The ordering of cells within each individual image is the same as the ordering
#' of these cells in the input object.
#'
#' @return returns an object of \code{class(object)} containing a new column
#' entry to \code{colData(object)[[name]]}. Cells in the object are grouped
#' by entries in \code{img_id}.
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
#' pancreasSCE <- distToCells(pancreasSCE,
#'                              x_cells = !is.na(pancreasSCE$patch_id),
#'                              coords = c("Pos_X","Pos_Y"),
#'                              statistics = "min",
#'                              img_id = "ImageNb")
#'
#' plotSpatial(pancreasSCE,
#'             img_id = "ImageNb",
#'             node_color_by = "distToCells",
#'             scales = "free")
#'
#' @author Daniel Schulz & Bruno Palau (\email{daniel.schulz@@uzh.ch})
#' @importFrom distances distances distance_columns
#' @importFrom MatrixGenerics rowMins rowMaxs rowMeans rowMedians
#' @export
distToCells <- function (object,
                        x_cells,
                        img_id,
                        name = "distToCells",
                        coords = c("Pos_X","Pos_Y"),
                        statistics = "min",
                        return_neg = TRUE,
                        BPPARAM = SerialParam()){
  .valid.distToCells.input(object, x_cells, name, coords, statistics,
                              img_id, return_neg)
  cur_meta <- metadata(object)
  metadata(object) <- list()
  cur_intmeta <- int_metadata(object)
  object$x_cells <- x_cells
  cur_out <- bplapply(
    unique(colData(object)[[img_id]]),
    function(x) {
      cur_obj <- object[, as.character(colData(object)[[img_id]]) == x]
      cur_obj[[name]] <- NA
      if (sum(cur_obj$x_cells) == 0 | sum(cur_obj$x_cells) ==
          ncol(cur_obj)) {
        return(cur_obj)
      }
      patch_cells <- which(cur_obj$x_cells)
      non_patch_cells <- which(!cur_obj$x_cells)
      if (is(object, "SpatialExperiment")) {
        dist_mat <- distances(spatialCoords(cur_obj))
      }
      else {
        dist_mat <- distances::distances(as.matrix(colData(cur_obj)[,
                                                                    coords]))
      }
      pos_dist <- distance_columns(dist_mat, column_indices = patch_cells)
      neg_dist <- distance_columns(dist_mat, column_indices = non_patch_cells)
      dist_to_patch <- switch(statistics,
                              mean = rowMeans(pos_dist),
                              median = rowMedians(pos_dist),
                              min = rowMins(pos_dist),
                              max = rowMaxs(pos_dist))


      if (return_neg == TRUE) {
        dist_from_patch <- switch(statistics,
                                  mean = rowMeans(neg_dist),
                                  median = rowMedians(neg_dist),
                                  min = rowMins(neg_dist),
                                  max = rowMaxs(neg_dist))

        dist_to_patch[cur_obj$x_cells] <- -dist_from_patch[cur_obj$x_cells]
      }
      cur_obj[[name]] <- dist_to_patch
      return(cur_obj)
    }, BPPARAM = BPPARAM)
  cur_out <- do.call("cbind", cur_out)
  cur_out$x_cells <- NULL
  metadata(cur_out) <- cur_meta
  int_metadata(cur_out) <- cur_intmeta
  message("The returned object is ordered by the '", img_id,
          "' entry.")
  return(cur_out)
}

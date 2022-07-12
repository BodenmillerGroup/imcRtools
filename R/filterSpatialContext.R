#' @title Filter spatial contexts
#'
#' @description Function to filter detected spatial contexts (SCs) based on a
#' user-defined threshold for number of group entries and/or cells.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param entry a single character specifying the \code{colData(object)} entry 
#' containing the \code{detectSpatialContext} output. Defaults to 
#' "spatial_context".
#' @param group_by a single character indicating the \code{colData(object)}
#' entry by which SCs are grouped. This is usually the image or patient ID. 
#' Defaults to "sample_id".
#' @param group_threshold a single numeric specifying the minimum number of
#' group entries in which a SC is detected.
#' @param cells_threshold a single numeric specifying the minimum total number
#' of cells in a SC.
#' @param name a single character specifying the name of the output saved in 
#'  \code{colData(object)}. Defaults to "spatial_context_filtered".
#' 
#' @return returns an object of \code{class(object)} containing a new column
#' entry to \code{colData(object)[[name]]} and a new \code{dataframe} entry to
#' \code{metadata(object)[["filterSpatialContext"]]} containing the group and
#' cell counts per SC.
#' 
#' @examples 
#' set.seed(22)
#' library(cytomapper)
#' data(pancreasSCE)
#'
#' ## 1. Cellular neighborhood (CN)
#' sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb",
#'                         type = "knn",
#'                         name = "knn_cn_graph",
#'                         k = 5)
#'
#' sce <- aggregateNeighbors(sce, colPairName = "knn_cn_graph",
#'                          aggregate_by = "metadata",
#'                          count_by = "CellType",
#'                          name = "aggregatedCellTypes")
#'
#' cur_cluster <- kmeans(sce$aggregatedCellTypes, centers = 3)
#' sce$cellular_neighborhood <- factor(cur_cluster$cluster)
#'
#' plotSpatial(sce, img_id = "ImageNb",
#'            colPairName = "knn_cn_graph",
#'            node_color_by = "cellular_neighborhood")
#'
#' ## 2. Spatial context (SC)
#' sce <- buildSpatialGraph(sce, img_id = "ImageNb",
#'                         type = "knn",
#'                         name = "knn_sc_graph",
#'                         k = 15)
#'
#' sce <- aggregateNeighbors(sce, colPairName = "knn_sc_graph",
#'                          aggregate_by = "metadata",
#'                          count_by = "cellular_neighborhood",
#'                          name = "aggregatedNeighborhood")
#'
#' # Detect spatial context
#' sce <- detectSpatialContext(sce, entry = "aggregatedNeighborhood",
#'                            threshold = 0.9)
#'
#' plotSpatial(sce, img_id = "ImageNb",
#'            colPairName = "knn_sc_graph",
#'            node_color_by = "spatial_context")
#'             
#' # Filter spatial context
#' # By group
#' sce <- filterSpatialContext(sce, group_by = "ImageNb", 
#'                             group_threshold = 2)
#' 
#' plotSpatial(sce, img_id = "ImageNb", 
#'             colPairName = "knn_sc_graph", 
#'             node_color_by = "spatial_context_filtered")
#'
#' # By cells
#' sce <- filterSpatialContext(sce, group_by = "ImageNb", 
#'                             cells_threshold = 15)
#'
#' plotSpatial(sce, img_id = "ImageNb", 
#'            colPairName = "knn_sc_graph", 
#'            node_color_by = "spatial_context_filtered")
#'            
#' @seealso 
#' \code{\link[imcRtools]{detectSpatialContext}} for the function to detect
#' spatial contexts
#'
#' \code{\link[imcRtools]{plotSpatialContext}} for the function to plot 
#' spatial context graphs
#' 
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @references
#' \href{https://doi.org/10.1016/j.cels.2021.09.012}{
#' Salil S. et al., Tissue schematics map the specialization of immune tissue 
#' motifs and their appropriation by tumors, Cell Systems, 2022}
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom BiocGenerics table
#' @importFrom dplyr all_of count n filter group_by_at group_by_all pull summarise
#' @importFrom S4Vectors unfactor
#' 
#' @export

filterSpatialContext <- function(object,
                               entry = "spatial_context",
                               group_by = "sample_id",
                               group_threshold = NULL,
                               cells_threshold = NULL,
                               name = "spatial_context_filtered"
                               ){
  
  .valid.filterSpatialContext.input(object, entry, group_by, group_threshold,
                                    cells_threshold, name)
  
  anno <- colData(object) %>% 
    as.data.frame %>% 
    select(all_of(entry), all_of(group_by)) %>% 
    group_by_all() %>% 
    count() %>% 
    group_by_at(entry) %>% 
    summarise(n_cells = sum(n), n_group = n()) %>%
    as.data.frame
  
  if (!is.null(group_threshold)) {
    anno <- anno[anno$n_group >= group_threshold,]
  }
  
  if (!is.null(cells_threshold)) {
      anno <- anno[anno$n_cells >= cells_threshold,]
  }
  
  selected <- anno %>% pull(entry)
  
  out_dat <- colData(object)[[entry]]
  out_dat[!out_dat %in% selected] <- NA 
  
  colData(object)[[name]] <- out_dat
  metadata(object)[["filterSpatialContext"]] <- anno
  
  return(object)
}

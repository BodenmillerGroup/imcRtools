#' @title Filter spatial contexts
#'
#' @description Function to filter detected spatial contexts (SCs) based on a
#' user-defined threshold for number of samples and/or cells.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param entry single character specifying the \code{colData(object)} entry 
#' containing the \code{detectSpatialContext} output. Defaults to 
#' "spatial_context".
#' @param sample_id single character specifying the \code{colData(object)} entry 
#' containing the sample identifiers. Defaults to "sample_id".
#' @param n_cells_threshold single numeric specifying the minimum total number
#' of cells in a SC.
#' @param n_samples_threshold single numeric specifying the minimum number of
#' samples in which a SC is detected.
#' @param name single character specifying the name of the output saved in 
#'  \code{colData(object)}. Defaults to "spatial_context_filtered".
#'
#' @return returns an object of \code{class(object)} containing a new column 
#' entry to \code{colData(object)[[name]]}
#' 
#' @examples 
#' library(cytomapper)
#' data(pancreasSCE)
#' 
#' ## 1. Cellular neighborhood (CN)
#' sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
#'                          type = "knn", 
#'                          name = "knn_cn_graph", 
#'                          k = 5)
#' 
#' sce <- aggregateNeighbors(sce, colPairName = "knn_cn_graph", 
#'                           aggregate_by = "metadata", 
#'                           count_by = "CellType")
#' 
#' cur_cluster <- kmeans(sce$aggregatedNeighbors, centers = 3)
#' sce$cellular_neighborhood <- factor(cur_cluster$cluster)
#' 
#' plotSpatial(sce, img_id = "ImageNb", 
#'             colPairName = "knn_cn_graph", 
#'             node_color_by = "cellular_neighborhood")
#' 
#' ## 2. Spatial context (SC)
#' sce <- buildSpatialGraph(sce, img_id = "ImageNb", 
#'                          type = "knn", 
#'                          name = "knn_sc_graph", 
#'                          k = 15)
#' 
#' sce <- aggregateNeighbors(sce, colPairName = "knn_sc_graph", 
#'                           aggregate_by = "metadata", 
#'                          count_by = "cellular_neighborhood")
#' 
#' # Detect spatial context 
#' sce <- detectSpatialContext(sce, threshold = 0.9)
#' 
#' plotSpatial(sce, img_id = "ImageNb", 
#'             colPairName = "knn_sc_graph", 
#'             node_color_by = "spatial_context")
#'             
#' # Filter spatial context
#' # By samples
#' sce <- filterSpatialContext(sce, sample_id = "ImageNb", 
#'                             n_samples_threshold = 2)
#' 
#' plotSpatial(sce, img_id = "ImageNb", 
#'             colPairName = "knn_sc_graph", 
#'             node_color_by = "spatial_context_filtered")
#'
#' # By cells
#' sce <- filterSpatialContext(sce, sample_id = "ImageNb", 
#'                             n_cells_threshold = 15)
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
#' @importFrom dplyr count filter group_by_at pull summarise
#' @importFrom S4Vectors unfactor
#' 
#' @export

filterSpatialContext <- function(object,
                               entry = "spatial_context",
                               sample_id = "sample_id",
                               n_cells_threshold = NULL,
                               n_samples_threshold = NULL,
                               name = "spatial_context_filtered"
                               ){
  
  .valid.filterSpatialContext.input(object, entry, sample_id, n_cells_threshold,
                                    n_samples_threshold, name)
  
  data <- colData(object)[,colnames(colData(object)) %in% c(entry,sample_id)] %>% 
  table() %>% as.data.frame
  Freq <- as.name("Freq")
  n <- as.name("n")
  
  anno <- data.frame(spatial_context = unfactor(unique(data[,entry])), 
                     n_cells = data %>% group_by_at(entry) %>% 
                       summarise(sum = sum(Freq)) %>% pull(sum), 
                     n_samples = data %>% group_by_at(entry) %>% 
                       filter(Freq != 0) %>% count() %>% pull(n)
                     )
  
  if (!is.null(n_cells_threshold)) {
  anno <- anno[anno$n_cells >= n_cells_threshold,]
  }
  
  if (!is.null(n_samples_threshold)) {
  anno <- anno[anno$n_samples >= n_samples_threshold,]
  }
  
  selected <- anno$spatial_context
  
  out_dat <- ifelse(colData(object)[[entry]] %in% selected, 
                    colData(object)[[entry]], NA)
  
  colData(object)[[name]] <- out_dat
  return(object)
}

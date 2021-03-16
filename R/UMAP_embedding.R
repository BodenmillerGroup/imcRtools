#' @rdname UMAP_embedding
#' @title Computation of a UMAP embedding 
#'
#' @description Computes a UMAP embedding 
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param  path to the Object relationships.csv file produced by the ImcSegmentationPipeline
#'
#' @return Return an updated \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with a list of graph, one per image 
#' @examples 
#' sce = Compute_neighborhood_graph(sce,graph_type = "KNN",graph_parameter = 10)
#'
#'
#' @import uwot 

UMAP_embedding = function(sce,metric="correlation",n_neighbors = 30) {
  
  if (!"Normalized_intensity"%in%names(assays(sce))) {
    stop("The data have not been normalised. Please proceed to normalisation first")
  }
  data_to_project =assays(sce)[["Normalized_intensity"]]
  data_to_project = t(data_to_project)
  
  cat("Computing UMAP embedding...")
  umap_embedding = umap(data_to_project,n_neighbors = n_neighbors,pca = NULL,metric = metric,verbose=F,fast_sgd = TRUE,n_threads = metadata(sce)$N_core )
  cat(" done ! \n")
  
  reducedDim(sce, "UMAP") <- umap_embedding
  return(sce)
}


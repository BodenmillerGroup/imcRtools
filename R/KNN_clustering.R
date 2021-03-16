#' @rdname KNN_clustering
#' @title Perform cell clustering on a computed KNN graph 
#'
#' @description Perform cell clustering on a computed KNN graph 
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param K number of neighbors for the KNN graph computation
#' @param clustering_method method used for graph clustering. Has to be chose among "Louvain","Greedy" and "Infomap"
#' @return Returns an updated \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with an updated colLabels slot
#'
#' @examples 
#' sce = KNN_clustering(sce,K = 30,clustering_method = "Louvain")
#'#'
#' @import SingleCellExperiment 
#' @import igraph 
#' @import N2R 

#' @export

KNN_clustering = function(sce,K=30,clustering_method = "Louvain") {
  if (!"Normalized_intensity"%in%names(assays(sce))) {
    stop("The data have not been normalised. Please proceed to normalisation first")
  }
  
  if (!clustering_method %in%c("Louvain","Greedy","Infomap")) {
    stop("The clustering method required does not exist. Please choose among Louvain,Greedy and Infomap !")
  }
    
  data_to_cluster =assays(sce)[["Normalized_intensity"]]
  data_to_cluster = t(data_to_cluster)
  
  Channel_for_clustering = rowData(sce)$Used_for_clustering
  
  #If no selection of the channels : using all channels
  if (is.null(Channel_for_clustering)) {
    Channel_for_clustering = rep(TRUE,ncol(data_to_cluster))
  }
  
  data_to_cluster = data_to_cluster[,Channel_for_clustering]
  
  cat(paste("Clustering of the data using",as.character(sum(Channel_for_clustering)),"channels \n"))
  
  KNN_graph_matrix =  N2R::Knn(as.matrix(data_to_cluster), K, nThreads=metadata(sce)$N_core, verbose=T, indexType='angular')
  KNN_graph_matrix_symmetric = (KNN_graph_matrix+t(KNN_graph_matrix))/2
  Final_graph <- igraph::graph_from_adjacency_matrix(KNN_graph_matrix_symmetric,mode='undirected',weighted=TRUE)
  
  if (clustering_method == "Louvain") {
    Clustering <- igraph::cluster_louvain(Final_graph)
  }
  
  if (clustering_method == "Greedy") {
    Clustering <- igraph::cluster_fast_greedy(Final_graph,modularity = TRUE)
  }
  
  if (clustering_method == "Infomap") {
    Clustering <- igraph::cluster_infomap(Final_graph,modularity = TRUE)
  }
  
  Clustering_group = membership(Clustering)
  colLabels(sce) = Clustering_group
  cat(paste(as.character(length(unique(Clustering_group))),"clusters have been identified \n"))
  return(sce)
}

#' @rdname Compute_neighborhood_graph
#' @title Computation of a neighborhood graph from an SCE object
#'
#' @description Computes a neighborhood graph using various strategies (KNN, FRNN, Delaunay triangulation)
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param graph_type type of graph to compute. Has to be chosen among the following values : "KNN","FRNN" or "Delaunay"
#' @param graph_parameter parameter for the graph computation. For the KNN it corresponds to the K parameter and for the FRNN to the radius value.
#' @return Return an updated \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with a list of graph, one per image 
#' @examples 
#' sce = Compute_neighborhood_graph(sce,graph_type = "KNN",graph_parameter = 10)
#'
#'
#' @import igraph 
#' @import geometry 
#' @import dbscan 
#' @import N2R 


Compute_neighborhood_graph = function(sce,graph_type,graph_parameter = 10) {
  
  List_graph = list()
  
  if (!graph_type%in%c("KNN","FRNN","Delaunay")) {
    stop("The type of graph required is not available. Please choose among KNN,Epsilon or Delaunay  ! \n")
  }
  
  
  if (graph_type == "KNN") {
    
    for (i in unique(sce$ImageNumber)) {
      Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==i],
                                      Y = sce$Location_Center_Y[sce$ImageNumber==i])
      
      N = nrow(Temp_location_data)
      
      
      KNN_graph_matrix =  N2R::Knn(as.matrix(Temp_location_data),k = graph_parameter, nThreads=1, indexType='L2')
      Temp_graph = igraph::graph_from_adjacency_matrix(KNN_graph_matrix,mode = "undirected")
      Temp_graph = igraph::simplify(Temp_graph,remove.multiple = TRUE,remove.loops = TRUE)
      
      List_graph[[i]] = Temp_graph
    }
  }
  
  
  if (graph_type == "FRNN") {
    
    for (i in unique(sce$ImageNumber)) {
      Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==i],
                                      Y = sce$Location_Center_Y[sce$ImageNumber==i])
      
      FRNN_list_edges = frNN(x = Temp_location_data,eps = graph_parameter)
      FRNN_list_edges = FRNN_list_edges$id
      FRNN_graph = igraph::graph_from_adj_list(FRNN_list_edges,mode = "all",duplicate = TRUE)
      List_graph[[i]] = FRNN_graph
    }
  }
  
  if (graph_type == "Delaunay") {
    
    for (i in unique(sce$ImageNumber)) {
      Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==i],
                                      Y = sce$Location_Center_Y[sce$ImageNumber==i])
      
      
      Delaunay_triangulation = delaunayn(Temp_location_data,output.options = "Fn")
      Delaunay_triangulation_graph = igraph::graph_from_adj_list(Delaunay_triangulation$neighbours,mode = "all",duplicate = TRUE)
      Delaunay_triangulation_graph =  igraph::simplify(Delaunay_triangulation_graph)
      List_graph[[i]] = Delaunay_triangulation_graph
    }
  }
  
  sce@metadata$List_graphs = List_graph
  return(sce)
  
}

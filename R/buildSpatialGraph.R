#' @title Builds an interaction graph based on the cells' locations
#'
#' @description Function to define cell-cell interactions via distance-based
#' expansion, delauney triangulation or k nearest neighbour detection.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param img_id single character indicating the \code{colData(object)} entry
#' containing the unique image identifiers.
#' @param type single character specifying the type of graph to be build.
#' Supported entries are \code{"expansion"} (default) to find interacting
#' cells via distance thresholding; \code{"delauney"} to find interactions via
#' delauney triangulation; \code{"knn"} to find the k nearest neighbouring
#' cells.
#' @param threshold (when \code{type = "expansion"}) single numeric specifying
#' the maximum distance for considering neighbours
#' @param k (when \code{type = "knn"}) single numeric integer defining the
#' number of nearest neighbours to search for.
#' @param coords character vector of length 2 specifying the names of the
#' \code{colData} (for a \code{SingleCellExperiment} object) or the
#' \code{spatialCoords} entries of the cells' x and y locations.
#' @param name single character specifying the name of the graph.
#' @param directed should the returned graph be directed? Only effects the k
#' nearest neighbour graph.
#' @param BNPARAM a \code{\link[BiocNeighbors]{BiocNeighborParam} object
#' defining the algorithm to use.}
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam-class}} object
#' defining how to parallelise computations.
#' @param ... additional parameters passed to the
#' \code{\link[BiocNeighbors]{findNeighbors}} function (\code{type =
#' "expansion"}), the \code{\link[RTriangle]{triangulate}} function
#' (\code{type = "delauney"}) or the \code{\link[BiocNeighbors]{findKNN}}
#' function (\code{type = "knn"})).
#' 
#' @return returns a \code{SpatialExperiment} or \code{SingleCellExperiment}
#' containing the graph in form of a \code{SelfHits} object in
#' \code{colPair(object, name)}.
#'
#' @section Building an interaction graph:
#' This function defines interacting cells in different ways. They are based 
#' on the cells' centroids and do not incorporate cell shape or area.
#' 
#' 1. When \code{type = "expansion"}, all cells within the radius 
#' \code{threshold} are considered interacting cells. 
#' 
#' 2. When \code{type = "delauney"}, interacting cells are found via a delauney
#' triangulation of the cells centroids.
#' 
#' 3. When \code{type = "knn"}, interacting cells are defined as the \code{k}
#' nearest neighbors in the 2D spatial plane. When \code{directed = FALSE},
#' the adjacency matrix of the graph is symmetric and edges between cells are
#' undirected. 
#' 
#' The graph is stored in form of a \code{SelfHits} object in
#' \code{colPair(object, name)}. This object can be regarded as an edgelist
#' and coerced to an \code{igraph} object via 
#' \code{graph_from_edgelist(as.matrix(colPair(object, name)))}.
#'
#' @section Choosing the graph construction method:
#' When finding interactions via \code{expansion} of \code{knn}, the 
#' \code{\link[BiocNeighbors]{findNeighbors}} or 
#' \code{\link[BiocNeighbors]{findKNN}} functions are used. Both functions
#' accept the \code{BNPARAM} via which the graph construction method can be 
#' defined (default \code{\link[BiocNeighbors]{KmknnParam}}). For an overview
#' on the different algorithms, see 
#' \code{\link[BiocNeighbors]{BiocNeighborParam}}. Within the 
#' \code{BiocNeighborParam} object, \code{distance} can be set to
#' \code{"Euclidean"} (default), \code{"Manhattan"} or \code{"Cosine"}.
#'
#' @examples
#' path <- system.file("extdata/mockData/steinbock", package = "imcRtools")
#' spe <- read_steinbock(path)
#' 
#' # Constructing a graph via expansion
#' spe <- buildSpatialGraph(spe, img_id = "sample_id", 
#'                          type = "expansion", threshold = 10)
#' colPair(spe, "expansion_interaction_graph")
#' 
#' # Constructing a graph via delauney triangulation
#' spe <- buildSpatialGraph(spe, img_id = "sample_id", 
#'                          type = "delauney")
#' colPair(spe, "delauney_interaction_graph")
#' 
#' # Constructing a graph via k nearest neighbor search
#' spe <- buildSpatialGraph(spe, img_id = "sample_id", 
#'                          type = "knn", k = 5)
#' colPair(spe, "knn_interaction_graph")
#' 
#' @seealso 
#' \code{\link[BiocNeighbors]{findNeighbors}} for the function finding interactions
#' via expansion
#' 
#' \code{\link[BiocNeighbors]{findKNN}} for the function finding interactions
#' via nearest neighbor search
#' 
#' \code{\link[RTriangle]{triangulate}} for the function finding interactions
#' via delauney triangulation
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#' 
#' @importFrom BiocNeighbors findNeighbors findKNN KmknnParam
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom igraph graph_from_adj_list graph_from_edgelist as.undirected
#' simplify as_edgelist
#' @importFrom RTriangle triangulate pslg
#' @export
buildSpatialGraph <- function(object,
                              img_id,
                              type = c("expansion", "delauney", "knn"),
                              k = NULL,
                              threshold = NULL,
                              coords = c("Pos_X", "Pos_Y"),
                              name = NULL,
                              directed = TRUE,
                              BNPARAM = KmknnParam(),
                              BPPARAM = SerialParam(),
                              ...){
    type <- match.arg(type)
    
    .valid.buildSpatialGraph.input(object, type, img_id, k, threshold, coords,
                                   name, directed)
    
    name <- ifelse(is.null(name), paste0(type, "_interaction_graph"), name)
    
    cur_ind <- unique(as.character(colData(object)[[img_id]]))
    
    cur_out <- bplapply(cur_ind,
                        function(x){
                            
                            cur_obj <- object[,as.character(colData(object)[[img_id]]) == x]
                            
                            # Create coords matrix
                            if (is(cur_obj, "SpatialExperiment")) {
                                cur_coords <- spatialCoords(cur_obj)[,coords]
                            } else {
                                cur_coords <- colData(cur_obj)[,coords]
                            }
                            
                            if (type == "expansion") {
                                cur_graph <- findNeighbors(cur_coords, 
                                                           threshold = threshold, 
                                                           get.distance = FALSE,
                                                           BNPARAM = BNPARAM,
                                                           ...)
                                cur_graph <- graph_from_adj_list(cur_graph$index)
                            } else if (type == "delauney") {
                                cur_graph <- triangulate(pslg(P = cur_coords),
                                                    ...)
                                cur_graph <- graph_from_edgelist(cur_graph$E)
                            } else {
                                cur_graph <- findKNN(cur_coords,
                                                     k = k,
                                                     BNPARAM = BNPARAM,
                                                     get.distance = FALSE,
                                                     ...)
                                cur_graph <- as.list(as.data.frame(t(cur_graph$index)))
                                cur_graph <- graph_from_adj_list(cur_graph)
                            }
                            
                            if (!directed) {
                                cur_graph <- as.undirected(cur_graph, 
                                                           mode = "collapse")
                            } 
                            
                            cur_graph <- simplify(cur_graph)
                            
                            cur_graph <- as_edgelist(cur_graph)
                            cur_graph <- SelfHits(from = cur_graph[,1],
                                                  to = cur_graph[,2], 
                                                  nnode = nrow(cur_coords))
                            colPair(cur_obj, name) <- cur_graph
                            
                            return(cur_obj)
                            
                        }, BPPARAM = BPPARAM)
    
    cur_out <- do.call("cbind", cur_out)
                               
    return(cur_out)
}

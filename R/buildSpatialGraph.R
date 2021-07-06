#' @title Builds an interaction graph based on the cells' locations
#'
#' @description 
#'
#' @param 
#' 
#' @return returns a \code{SpatialExperiment} or \code{SingleCellExperiment}
#' containing the graph in form of a \code{SelfHits} object in colPair(object, name).
#'
#'
#' @examples
#' #TODO
#' 
#' @seealso 
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#' 
#' @importFrom BiocNeighbours findNeighbours
#' @importFrom geometry delauneyn
#' @export
buildSpatialGraph <- function(object,
                              img_id,
                              type = c("expansion", "delauney", "knn"),
                              k = NULL,
                              threshold = NULL,
                              Pos_X = "Pos_X",
                              Pos_Y = "Pos_Y",
                              name = NULL,
                              directed = TRUE,
                              BNPARAM = KmknnParam(),
                              BPPARAM = SerialParam(),
                              ...){
    
    #.valid.buildSpatialGraph.input()
    
    type <- match.arg(type)
    

    name <- ifelse(is.null(name), paste0(type, "_interaction_graph"), name)
    
    cur_ind <- unique(as.character(colData(object)[[img_id]]))
    
    cur_out <- bplapply(cur_ind,
                        function(x){
                            
                            cur_obj <- object[,as.character(colData(object)[[img_id]]) == x]
                            
                            # Create coords matrix
                            if (is(cur_obj, "SpatialExperiment")) {
                                cur_coords <- spatialCoords(cur_obj)[,c(Pos_X, Pos_Y)]
                            } else {
                                cur_coords <- colData(cur_obj)[,c(Pos_X, Pos_Y)]
                            }
                            
                            if (type == "expansion") {
                                cur_graph <- findNeighbors(cur_coords, 
                                                           threshold = threshold, 
                                                           get.distance = FALSE,
                                                           BNPARAM = BNPARAM)
                                cur_graph <- graph_from_adj_list(cur_graph$index)
                            } else if (type == "delauney") {
                                cur_graph <- deldir(cur_coords[,1], cur_coords[,2])
                                cur_graph <- graph_from_edgelist(as.matrix(cur_graph$delsgs[,c("ind1", "ind2")]))
                            } else {
                                cur_graph <- findKNN(cur_coords,
                                                     k = k,
                                                     BNPARAM = BNPARAM,
                                                     get.distance = FALSE)
                                cur_graph <- graph_from_adj_list(as.list(as.data.frame(t(cur_graph$index))))
                            }
                            
                            if (!directed) {
                                cur_graph <- as.undirected(cur_graph, mode = "collapse")
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

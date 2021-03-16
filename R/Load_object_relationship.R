#' @rdname Load_object_relationship
#' @title Loads an object relation csv file produced by the ImcSegmentationPipeline
#'   (https://github.com/BodenmillerGroup/ImcSegmentationPipeline)
#'
#' @description Loads an object relation csv file, constructs the corresponding graph and adds it to a previously defined sce object
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param Path_to_object_relationship path to the Object relationships.csv file produced by the ImcSegmentationPipeline
#'
#' @return Return an updated \code{\link[SingleCellExperiment]{SingleCellExperiment}} object 
#' @examples 
#' sce = Load_object_relationship(sce,Path_to_object_relationship = "Desktop/analysis/cpout/Object relationships.csv")
#'
#'
#' @import readr 
#' @import igraph 



Load_object_relationship = function(sce , Path_to_object_relationship) {
  
  Object_relationship = readr::read_csv(Path_to_object_relationship,progress = TRUE)
  Object_relationship = as.data.frame(Object_relationship)
  
  Object_relationship_reshape = data.frame(ImageNumber = Object_relationship$`First Image Number`,
                                           Cell_number_1 = Object_relationship$`First Object Number`,
                                           Cell_number_2 = Object_relationship$`Second Object Number`)
  
  List_graph = list()
  
  for (i in unique(sce$ImageNumber)) {
    
    Object_relationship_reshape_temp = as.matrix(Object_relationship_reshape[Object_relationship_reshape$ImageNumber==i,-1])
    Temp_graph = igraph::graph_from_edgelist(el = Object_relationship_reshape_temp,directed = FALSE)
    List_graph[[i]] = Temp_graph
  }
  
  sce@metadata$List_graphs = List_graph
  return(sce)
  
}

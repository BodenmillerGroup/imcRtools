#' @title Plot spatial context graph
#'
#' @description Function to splot spatial context graphs based on edge-lists, 
#' which operates on image- and cohort-level. User can specify node, node_label 
#' and edge aesthetics.  
#'
#' @param edges data frame containing a symbolic edge list in the first two 
#' columns. 
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param entry single character specifying the \code{colData(object)} entry 
#' containing the \code{detectSpatialContext} output. If NULL, defaults to 
#' "spatial_context".
#' @param img_id single character specifying the \code{colData(object)} entry 
#' containing the unique image identifiers. If NULL, defaults to "sample_id".
#' @param combined Is the provided edge list created on cohort-level (two 
#' columns) and not image-level (three columns)? If NULL, defaults to TRUE.
#' @param directed should the created graph be directed? Defaults to TRUE.
#' @param node_color_by single character from \code{c("name","freq","n_samples")} 
#' by which the nodes should be colored.
#' @param node_size_by single character from \code{c("freq","n_samples")} 
#' by which the size of the nodes are defined.
#' @param node_color_fix single character specifying the color of all nodes.
#' @param node_size_fix single numeric specifying the size of all nodes.
#' @param node_label_repel should nodes be labelled? Defaults to TRUE.
#' @param node_label_color_by single character from \code{c("name","freq","n_samples")} 
#' by which the node labels should be colored.
#' @param node_label_color_fix single character specifying the color of all node 
#' labels.
#' @param draw_edges should edges be drawn between nodes? Defaults to TRUE.
#' @param edge_color_fix single character specifying the color of all edges.
#' 
#' @return returns a \code{ggplot} object. 
#' 
#' @examples TO DO
#' 
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom stringr str_split
#' @importFrom dplyr count filter group_by_at pull select summarise
#' @importFrom igraph graph_from_data_frame
#' @importFrom cowplot plot_grid
#' @importFrom Biobase listLen
#' @importFrom BiocGenerics table
#' @export

plotSpatialContext <- function(edges,
                               object,
                               entry = NULL,
                               img_id = NULL,
                               combined = NULL,
                               #build graph
                               directed = TRUE,
                               #nodes
                               node_color_by = NULL,#c("name","freq","n_samples"),
                               node_size_by = NULL,#c("freq","n_samples"),
                               node_color_fix = NULL,
                               node_size_fix = NULL,
                               #node labels
                               node_label_repel = TRUE,
                               node_label_color_by = NULL,#c("name","freq","n_samples"),
                               node_label_color_fix = NULL,
                               #plot graph - edges
                               draw_edges = TRUE,
                               edge_color_fix = NULL){
  
  entry <- ifelse(is.null(entry), "spatial_context", entry) #default
  img_id <- ifelse(is.null(img_id), "sample_id", img_id) #default
  combined <- ifelse(is.null(combined), TRUE, combined) #default
  
  .valid.plotSpatialContext.input(edges, object, entry, img_id, combined, directed, node_color_by, node_size_by, node_color_fix, node_size_fix, node_label_repel, node_label_color_by, node_label_color_fix, draw_edges, edge_color_fix)
  
  #data
  data <- colData(object)[,colnames(colData(object)) %in% c(entry,img_id)] %>% table() %>% as.data.frame
  Freq <- as.name("Freq")
  n <- as.name("n")
  
  if(combined == TRUE){ # For combined samples
    anno <- data.frame(spatial_context = unique(data[,entry]), 
                       length = listLen(str_split(unique(data[,entry]),"_")),
                       Freq = data %>% group_by_at(entry) %>% summarise(sum = sum(Freq)) %>% pull(sum), 
                       n_samples = data %>% group_by_at(entry) %>% filter(Freq != 0) %>% count() %>% pull(n) %>% as.character()
    )
    
    g <- graph_from_data_frame(edges, directed = directed,vertices = anno)
    
    #Plot using ggraph
    p <- .generateSpatialContextPlot(graph = g, node_color_by, node_size_by, node_color_fix, node_size_fix, node_label_repel, node_label_color_by, node_label_color_fix, draw_edges, edge_color_fix) #hidden function
    
    return(p)
    
  }else{ # For multiple SC graphs
    anno <- split(data %>% select(-as.name(img_id)), f = data[,img_id])
    
    anno <- lapply(anno, function(x){
      cur_anno <- x %>% filter(Freq !=0)
      cur_anno$length <- listLen(str_split(cur_anno[,entry],"_"))
      return(cur_anno)
    })
    
    #edges
    edges_list <- split(edges %>% select(-as.name(img_id)), f = edges[,img_id])
    
    #generate graph
    g <- mapply(function(x,y){
      g <- graph_from_data_frame(x, directed = TRUE, vertices = y)
      #specify vertical layout using sugiyama
      return(g)}, edges_list, anno)
    
    #generate plots
    all_plots <- lapply(g, function(x){
      p <- .generateSpatialContextPlot(graph = x, node_color_by, node_size_by, node_color_fix, node_size_fix, node_label_repel, node_label_color_by, node_label_color_fix, draw_edges, edge_color_fix) #hidden function
      return(p)
    })
    
    plot_grid(plotlist = all_plots) 
  }
}
#' @title Create edgelist for spatial context graph
#'
#' @description Function to create a symbolic edge list for spatial context 
#' graph construction. Based on single cell spatial context assignments, this 
#' function can operate on image- and cohort-level. 
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param entry single character specifying the \code{colData(object)} entry 
#' containing the \code{detectSpatialContext} output. If NULL, defaults to 
#' "spatial_context".
#' @param img_id single character specifying the \code{colData(object)} entry 
#' containing the unique image identifiers. If NULL, defaults to "sample_id".
#' @param combined should the edge-list be created on a cohort-level and not 
#' image-level? If NULL, defaults to TRUE.
#'
#' @return returns a data frame containing a symbolic edge list in the first
#' two columns ("from" and "to") and, if \code{combined = TRUE}, an additional 
#' column containing the \code{img_id} identifiers.
#' 
#' @examples
#' TO DO
#' 
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom stringr str_split
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider
#' @export

buildEdgeList <- function(object, 
                          entry = NULL,
                          img_id = NULL,
                          combined = NULL){
  
  entry <- ifelse(is.null(entry), "spatial_context", entry) #default
  name <- ifelse(is.null(name), "sample_id", name) #default
  combined <- ifelse(is.null(combined), TRUE, combined) #default
  
  
  .valid.buildEdgeList.input(object, entry, img_id, combined) #validity check
  
  #data  
  data <- colData(object)[,colnames(colData(sce)) %in% c(entry,img_id)] %>% table() %>% as.data.frame
  data_wide <- data %>% pivot_wider(values_from = Freq, names_from = entry) %>% column_to_rownames(img_id)
  
  edges <- if(combined == TRUE){ #Option 1: For all images combined  
    list <- str_split(unique(data$spatial_context), "_")
    list_length <- sapply(list, length)
    edges <- .createEdgeList(list, list_length) #hidden function
    return(edges)
    
  }else{ #Option 2: For each img_id separately
    cur_dat <- apply(data_wide, 1, function(x){
      cur <- x
      names(cur) <- colnames(data_wide)
      dat <- names(cur[cur != 0])
      return(dat)
    })
    
    edges <- lapply(cur_dat, function(z){ 
      list <- str_split(z, "_")
      list_length <- sapply(list, length)
      edges <- .createEdgeList(list, list_length) #hidden function
      return(edges)
    })
    
    edges <- do.call(rbind,edges)
    edges$sample_id <- paste0(str_split(rownames(edges),"\\.", simplify = TRUE)[,1],".",str_split(rownames(edges),"\\.", simplify = TRUE)[,2])
    rownames(edges) <- NULL
    return(edges)
  }
}
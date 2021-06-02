### Helper functions for reading in steinbock data
# Here x is a list of files
#' @importFrom vroom vroom
#' @importFrom BiocParallel bplapply
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SpatialExperiment SpatialExperiment
.read_intensities <- function(x, cell_id, return_as, panel, BPPARAM){
    
    cur_out <-  bplapply(seq_along(x),
        function(y){
            cur_int <- vroom(x[y], progress = FALSE, 
                             col_types = cols())
            cur_counts <- cur_int %>% select(-all_of(cell_id))
                             
            cur_name <- sub("\\.[^.]*$", "", basename(x[y]))
            
            if (return_as == "spe") {
                object <- SpatialExperiment(assays = list(counts = t(as.matrix(cur_counts))),
                                         sample_id = cur_name)
                object$ObjectId <- cur_int[[cell_id]]
            } else {
                object <- SingleCellExperiment(assays = list(counts = t(as.matrix(cur_counts))))     
                object$sample_id <- cur_name
                object$ObjectId <- cur_int[[cell_id]]
            }
            
            return(object)
            
        }, BPPARAM = BPPARAM)
    
    return(cur_out)
}

.read_regionprobs <- function(x, cur_path, cell_id, coords, return_as, BPPARAM){
    
    cur_out <-  bplapply(x,
        function(y){
            cur_sample <- unique(y$sample_id)
            
            cur_file <- list.files(cur_path, pattern = cur_sample, full.names = TRUE)
            
            if (length(test) == 0) {
                return(y)
            }
            
            cur_props <- vroom(cur_file, 
                               progress = FALSE, 
                               col_types = cols()) %>%
                as.data.frame()
            rownames(cur_props) <- cur_props[[cell_id]]
            cur_props <- cur_props[as.character(y$ObjectId),]
                             
            if (return_as == "spe") {
                spatialCoords(y) <- matrix(c(cur_props[[coords[1]]],
                                             cur_props[[coords[2]]]),
                                           ncol = 2, byrow = FALSE,
                                           dimnames = list(as.character(y$ObjectId),
                                                           c("Pos_X", "Pos_Y")))
                colData(y) <- cbind(colData(y),
                                    cur_props[,!colnames(cur_props) %in% c(cell_id, coords)])
            } else {
                colData(y)$Pos_X <- cur_props[[coords[1]]]
                colData(y)$Pos_Y <- cur_props[[coords[2]]]
                colData(y) <- cbind(colData(y),
                                    cur_props[,!colnames(cur_props) %in% c(cell_id, coords)])
            }
            
            return(y)
                             
        }, BPPARAM = BPPARAM)
    return(cur_out)
}

.read_graphs <- function(x, cur_path, return_as, BPPARAM){
    
    cur_out <-  bplapply(x,
        function(y){
            cur_sample <- unique(y$sample_id)
            
            cur_file <- list.files(cur_path, pattern = cur_sample, full.names = TRUE)
            
            if (length(test) == 0) {
                return(y)
            }
            
            cur_graphs <- vroom(cur_file, 
                                progress = FALSE, 
                                col_types = cols()) %>%
                as.data.frame()
            
            cur_hits <- SelfHits(from = cur_graphs[,1],
                                 to = cur_graphs[,2],
                                 nnode = ncol(y))
            
            colPair(y, "neighbourhood") <- cur_hits

            return(y)
    }, BPPARAM = BPPARAM)
    
    return(cur_out)
}

.add_panel <- function(x, path, panel, extract_names_from) {
    
    if (!is.null(panel)) {
        if (file.exists(file.path(path, panel))) {
            cur_panel <- vroom(file.path(path, panel), 
                               progress = FALSE, 
                               col_types = cols())
        } else if (file.exists(panel)) {
            cur_panel <- vroom(panel, 
                               progress = FALSE, 
                               col_types = cols())
        } else {
            warning("'panel' does not exist.")
            return(x)
        }
        
        cur_panel <- as.data.frame(cur_panel)
        rownames(cur_panel) <- cur_panel[,extract_names_from]
        
        cur_panel <- cur_panel[rownames(x),]
        
        rowData(x) <- cur_panel
    }
    
    return(x)

}

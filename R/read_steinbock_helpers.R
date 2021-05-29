### Helper functions for reading in steinbock data
# Here x is a list of files
.read_intensities <- function(x, path, return_as, panel, BPPARAM){
    
    cur_out <-  bplapply(seq_along(int_file_names),
        function(x){
            cur_int <- vroom(int_file_names[x], progress = FALSE, 
                             col_types = cols())
            cur_counts <- cur_int %>% select(-all_of(cell_id))
                             
            cur_name <- sub("\\.[^.]*$", "", basename(int_file_names[x]))
            
            if (return_as == "spe") {
                object <- SpatialExperiment(assays = list(counts = t(as.matrix(cur_counts))),
                                         sample_id = cur_name)
            } else {
                object <- SingleCellExperiment(assays = list(counts = t(as.matrix(cur_counts))))     
                object$sample_id <- cur_name
            }
            
        }, BPPARAM = BPPARAM)
}

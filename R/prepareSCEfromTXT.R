#' @title Generates a SingleCellExperiment from .txt files
#'
#' @description Helper function to process raw .txt files acquired by the
#'   Hyperion system into a \code{\linkS4class{SingleCellExperiment}} container.
#'   This function is mainly used to read-in data generated from a "spillover
#'   slide". Here, each .txt file contains the measurements of multiple pixels
#'   for a single stain across all open channels.
#'
#' @param txt_list a named list containing the read-in txt files for each spot. 
#' Entries to \code{txt_list} need to be coercible to \code{data.frame} objects.
#' @param metadata_cols character vector indicating which column entries of the
#' .txt files should be saved in the \code{colData(sce)} slot.
#' @param verbose logical indicating if additional information regarding the
#' spotted and acquired masses should be shown.

#' @return returns a SCE object
#'
#' @examples
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @export

prepareSCEfromTXT <- function(txt_list, 
                              metadata_cols = c("Start_push", "End_push", 
                                                "Pushes_duration", "X", 
                                                "Y", "Z"),
                              verbose = TRUE){
    
    .validSCEtoTXTinput(txt_list, 
                        metadata_cols,
                        verbose)
   
    # Coerce to data.frame
    txt_list <- lapply(txt_list, as.data.frame)
    
    cur_out <- do.call(rbind, txt_list)
    
    # Construct SCE object
    cell_meta <- DataFrame(cur_out[metadata_cols]) 
    cell_meta$sample_id <- str_extract(rownames(cell_meta), "^[A-Za-z]{1,2}[0-9]{2,3}")
    cell_meta$sample_metal <- str_extract(cell_meta$sample_id, "^[A-Za-z]{1,2}")
    cell_meta$sample_mass <- str_extract(cell_meta$sample_id, "[0-9]{2,3}$")
    
    cur_counts <- cur_out[grepl("[A-Za-z]{1,2}[0-9]{2,3}", colnames(cur_out))]
    cur_counts <- t(cur_counts)
    
    channel_meta <- DataFrame(channel_name = str_extract(rownames(cur_counts), "[A-Za-z]{1,2}[0-9]{2,3}Di"),
                              marker_name = str_extract(rownames(cur_counts), "[A-Za-z]{1,2}[0-9]{2,3}"))
    
    sce <- SingleCellExperiment(assays = list(counts = cur_counts))
    colData(sce) <- cell_meta
    rowData(sce) <- channel_meta
    
    return(sce)
       
}

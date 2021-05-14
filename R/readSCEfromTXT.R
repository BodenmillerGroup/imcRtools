#' @title Generates a SingleCellExperiment from .txt files
#'
#' @description Helper function to process raw .txt files acquired by the
#'   Hyperion system into a \code{\linkS4class{SingleCellExperiment}} container.
#'   This function is mainly used to read-in data generated from a "spillover
#'   slide". Here, each .txt file contains the measurements of multiple pixels
#'   for a single stain across all open channels.
#'
#' @param path path to where the .txt files are located
#' @param pattern pattern to select which files should be read in (default \code{".txt$"})
#' @param metadata_cols character vector indicating which column entries of the
#' .txt files should be saved in the \code{colData(sce)} slot.
#' @param verbose logical indicating if additional information regarding the
#' spotted and acquired masses should be shown.

#' @return returns a SCE object
#' 
#' @section Reading in .txt files:
#' 
#' As part of the mcd folder
#' 
#' This function will try to automatically read in the metal names from the file names
#' For this, the first occurrence of the metal name in form {mt}{123} will be used.
#' The filenames will also be stored in the SingeCellExperiment object
#'
#' @examples
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData<- rowData<-
#' @export
readSCEfromTXT <- function(path, 
                           pattern = ".txt$",
                           metadata_cols = c("Start_push", "End_push", 
                                                "Pushes_duration", "X", 
                                                "Y", "Z"),
                           verbose = TRUE){
    
    txt_list_names <- list.files(path, pattern = pattern, full.names = FALSE)
    cur_names <- 
    
    .valid.prepareSCEfromTXT.input(path
                        metadata_cols,
                        verbose)
   
    # Coerce to data.frame
    txt_list <- list.files(path, pattern = pattern, full.names = TRUE)
    txt_list <- suppressMessages(lapply(txt_list, read_delim, delim = "\t"))
    txt_list <- lapply(txt_list, as.data.frame)
    names(txt_list) <- str_extract(txt_list_names, "[A-Za-z]{1,2}[0-9]{2,3}")
    
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

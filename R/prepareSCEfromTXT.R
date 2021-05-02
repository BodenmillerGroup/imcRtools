#' @title Generates a SingleCellExperiment from .txt files
#'
#' @description Helper function to process raw .txt files acquired by the
#'   Hyperion system into a \code{\linkS4class{SingleCellExperiment}} container.
#'   This function is mainly used to read-in data generated from a "spillover
#'   slide". Here, each .txt file contains the measurements of multiple pixels
#'   for a single stain across all open channels.
#'
#' @param txt_list a named list containing the read-in txt files for each spot
#' @param metadata_cols character vector indicating which column entries of the
#' .txt files should be saved in the \code{colData(sce)} slot.
#' @param verbose logical indicating if additional information regarding the
#' spotted and acquired masses should be shown.

#' @return returns a SCE object
#'
#' @examples
#'
#' @import SingleCellExperiment
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
    
    if (!is.null(names(txt_list))) {
        txt_list <- lapply(names(txt_list), function(x){
            cur_txt <- txt_list[[x]]
            cur_txt$metal_name <- as.character(str_match(x, "[A-Za-z]{2}[0-9]{2,3}"))
            cur_txt$mass <- as.character(str_match(cur_txt$metal_name, "[0-9]{2,3}"))
            return(cur_txt)
        })
    }
   
    cur_out <- do.call(rbind, txt_list)
    
    # Construct SCE object
    cur_col <- 
   
}

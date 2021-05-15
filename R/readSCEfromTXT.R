#' @title Generates a SingleCellExperiment from .txt files
#'
#' @description Helper function to process raw .txt files acquired by the
#' Hyperion system into a \code{\linkS4class{SingleCellExperiment}} object.
#' This function is mainly used to read-in data generated from a "spillover
#' slide". Here, each .txt file contains the measurements of multiple pixels
#' for a single stain across all open channels.
#'
#' @param x input can be of different types:
#' \describe{
#' \item{A path}{Full path to where the single stain .txt files are located.}
#' \item{A list object}{A named list object where each entry is a
#' \code{data.frame} or coercible to one. The names of each entry indicate the
#' spotted metals (see details).}
#' }
#' @param pattern pattern to select which files should be read in (default
#' \code{".txt$"}). Only used when \code{x} is a path.
#' @param metadata_cols character vector indicating which column entries of the
#' .txt files should be saved in the \code{colData(sce)} slot.
#' @param verbose logical indicating if additional information regarding the
#' spotted and acquired masses should be shown.

#' @return returns a SCE object where pixels are stored as columns and acquired
#' channels are stored as rows.
#' 
#' @section Reading in .txt files for spillover correction:
#' 
#' As described in the original publication, single metal spots are acquired
#' using the Hyperion imaging system. Each acquisition corresponds to one spot.
#' All acquisitions are stored in a single .mcd file and individual acquisitions
#' are stored in single .txt files.
#' 
#' This function aggregates these measurements into a single \code{SingleCellExperiment} object.
#' For this, two inputs are possible:
#' 
#' \enumerate{
#' \item \code{x} is a path:
#' By default all .txt files are read in from the specified path. Here, the path
#' should indicate the location of the spillover slide measurement. The file
#' names of the .txt file must contain the spotted metal isotope name in the
#' format \code{(mt)(mass)} (e.g. \code{Sm152} for Samarium isotope with the
#' atomic mass 152). Internally, the first occurrence of such a pattern is read
#' in as the metal isotope name and stored in the \code{colData(sce)$sample_id}
#' slot.
#' 
#' \item \code{x} is a named list:
#' If there are issues with reading in the metal isotope names from the .txt
#' file names, the user can provide a list in which each entry contains the
#' contents of a single .txt file. The names of the list must indicate the
#' spotted metal in the format \code{(mt)(mass)}. These names will be stored in
#' the \code{colData(sce)$sample_id} slot.
#' }
#'
#' @examples
#' # Read files from path
#' path <- system.file("extdata/spillover", package = "imcRtools")
#' 
#' sce <- readSCEfromTXT(path) 
#' sce
#' 
#' # Read files as list
#' cur_file_names <- list.files(path, pattern = ".txt", full.names = TRUE)
#' cur_files <- lapply(cur_file_names, read.delim)
#' names(cur_files) <- sub(".txt", "", basename(cur_file_names))
#' 
#' sce <- readSCEfromTXT(cur_files) 
#' sce
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#' 
#' @references
#' \href{https://www.sciencedirect.com/science/article/pii/S1550413118306910}{Chevrier,
#' S. et al. 2017. “Compensation of Signal Spillover in Suspension and Imaging
#' Mass Cytometry.” Cell Systems 6: 612–20.}
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData<- rowData<-
#' @importFrom stringr str_extract
#' @importFrom readr read_delim
#' @export
readSCEfromTXT <- function(x, 
                           pattern = ".txt$",
                           metadata_cols = c("Start_push", "End_push", 
                                                "Pushes_duration", "X", 
                                                "Y", "Z"),
                           verbose = TRUE){
    
    if (all(is.character(x)) & length(x) == 1) {
        
        if (!dir.exists(x)) {
            stop("Path does not exist.")
        }
        
        cur_names <- list.files(x, pattern = pattern, full.names = FALSE)
        
        if (length(cur_names) == 0) {
            stop("Files could not be read in.")
        }
        
        cur_names <- str_extract(cur_names, "[A-Za-z]{1,2}[0-9]{2,3}")
    
        
        txt_list <- list.files(x, pattern = pattern, full.names = TRUE)
        txt_list <- suppressMessages(lapply(txt_list, read_delim, delim = "\t"))
        txt_list <- lapply(txt_list, as.data.frame)
        names(txt_list) <- cur_names
        
    } else if (is.list(x)) {
        
        if (is.null(names(x))) {
            stop("If 'x' is a list, it needs to be named.")
        }
        cur_names <- names(x)
        txt_list <- lapply(x, as.data.frame)
        
    } else {
        stop("Input 'x' is not of the correct format.")
    }
    
    .valid.readSCEfromTXT.input(txt_list, cur_names,
                                metadata_cols, verbose)
   
    cur_out <- do.call(rbind, txt_list)
    
    # Construct SCE object
    cell_meta <- DataFrame(cur_out[metadata_cols]) 
    cell_meta$sample_id <- str_extract(rownames(cell_meta), 
                                       "^[A-Za-z]{1,2}[0-9]{2,3}")
    cell_meta$sample_metal <- str_extract(cell_meta$sample_id, "^[A-Za-z]{1,2}")
    cell_meta$sample_mass <- str_extract(cell_meta$sample_id, "[0-9]{2,3}$")
    
    cur_counts <- cur_out[grepl("[A-Za-z]{1,2}[0-9]{2,3}", colnames(cur_out))]
    cur_counts <- t(cur_counts)
    
    channel_name <- str_extract(rownames(cur_counts), 
                                "[A-Za-z]{1,2}[0-9]{2,3}Di")
    
    channel_meta <- DataFrame(channel_name = channel_name,
                              marker_name = sub("Di", "", channel_name))
    
    sce <- SingleCellExperiment(assays = list(counts = cur_counts))
    colData(sce) <- cell_meta
    rowData(sce) <- channel_meta
    
    return(sce)
       
}

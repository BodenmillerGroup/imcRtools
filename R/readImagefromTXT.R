#' @title Reads one or multiple .txt files into a CytoImageList object
#'
#' @description Reader function to generate \code{\linkS4class{Image}} objects
#' in form of a \code{\linkS4class{CytoImageList}} container from .txt files.
#'
#' @param path Full path to where the individual .txt files are located. This is
#' usualy the path where the .mcd file is located.
#' @param pattern pattern to select which files should be read in (default
#' \code{".txt$"}). 
#' @param channel_pattern regular expression to select the channel names from
#' the files.
#' 
#' @return returns a \code{\linkS4class{CytoImageList}} object containing one
#' \code{\linkS4class{Image}} object per .txt file.
#' 
#' @section Imaging mass cytometry .txt files:
#' 
#' Explain what they are
#'
#' @examples
#' 
#' @seealso \code{\linkS4class{CytoImageList}} for the container
#' @seealso \code{\linkS4class{Image}} for the multi-channel image object
#' @seealso \code{vignette("cytomapper")} for visualization of multi-channel
#' images
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#' 
#' @importFrom cytomapper CytoImageList channelNames
#' @importFrom abind abind
#' @importFrom EBImage Image
#' @importFrom stringr str_extract
#' @importFrom readr read_delim
#' @export
readSCEfromTXT <- function(path, 
                           pattern = ".txt$",
                           channel_pattern = "[A-Za-z]{1,2}[0-9]{2,3}Di",
                           index_names = c("X", "Y", "Z"),
                           verbose = TRUE){
    
    cur_files <- list.files(path, pattern = pattern, full.names = TRUE)
    
    cur_dat <- suppressMessages(lapply(cur_files, read_delim, delim = "\t"))
    cur_dat <- lapply(cur_dat, function(x){
        cur_mat <- x[,grepl(channel_pattern, colnames(x))]
        cur_ind <- x[,index_names]
        
        n_row <- max(cur_ind[,index_names[1]]) + 1
        n_col <- max(cur_ind[,index_names[2]]) + 1
        
        cur_mat <- lapply(cur_mat, function(y){
            matrix(y, nrow = n_row, ncol = n_col, byrow = FALSE)
        })
        
        cur_mat <- do.call(abind, list(cur_mat, along = 3))
        
        Image(cur_mat)
    })
    
    cur_dat <- CytoImageList(cur_dat)
    
    channelNames(cur_dat) <- str_extract(channelNames(cur_dat), channel_pattern)
    
    return(cur_dat)
}

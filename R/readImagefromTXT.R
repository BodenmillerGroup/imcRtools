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
#' @export
readSCEfromTXT <- function(path, 
                           pattern = ".txt$",
                           channel_pattern = "[A-Za-z]{1,2}[0-9]{2,3}Di",
                           verbose = TRUE){
    
    
    
    
    return(cil)
       
}

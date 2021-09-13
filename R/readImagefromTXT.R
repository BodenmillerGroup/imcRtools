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
#' @param index_names exact names of the columns storing the x and y coordinates 
#' of the image
#' @param BPPARAM parameters for parallelized reading in of images. 
#' This is only recommended for very large images. 
#' 
#' @return returns a \code{\linkS4class{CytoImageList}} object containing one
#' \code{\linkS4class{Image}} object per .txt file.
#' 
#' @section Imaging mass cytometry .txt files:
#' As part of the raw data folder, the Hyperion imaging system writes out
#' one .txt file per acquisition. These files store the ion counts per
#' pixel and channel. 
#' 
#' This function reads these .txt files into a single
#' \code{\linkS4class{CytoImageList}} object for downstream analysis. The
#' \code{pattern} argument allows selection of all .txt files or a specific
#' subset of files. The \code{\link[cytomapper]{channelNames}} of the
#' \code{CytoImageList} object are determined by the \code{channel_pattern}
#' argument.
#'
#' @examples
#' path <- system.file("extdata/mockData/raw", package = "imcRtools")
#' 
#' # Read in all images
#' x <- readImagefromTXT(path)
#' x
#' 
#' # Read in specific files
#' y <- readImagefromTXT(path, pattern = "ROI_002")
#' y
#' 
#' # Read in other channelNames
#' z <- readImagefromTXT(path, channel_pattern = "[A-Za-z]{2}[0-9]{3}")
#' z
#' 
#' @seealso 
#' \code{\linkS4class{CytoImageList}} for the container
#' 
#' \code{\link[BiocParallel]{MulticoreParam}} for parallelized processing
#' 
#' \code{\linkS4class{Image}} for the multi-channel image object
#' 
#' \code{vignette("cytomapper")} for visualization of multi-channel images
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#' 
#' @importFrom cytomapper CytoImageList channelNames<- channelNames
#' @importFrom abind abind
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom EBImage Image
#' @importFrom stringr str_extract
#' @importFrom readr read_delim
#' @export
readImagefromTXT <- function(path, 
                            pattern = ".txt$",
                            channel_pattern = "[A-Za-z]{1,2}[0-9]{2,3}Di",
                            index_names = c("X", "Y"),
                            BPPARAM = SerialParam()){
    
    # Validity checks
    .valid.readImagefromTXT.input(path, pattern)
    
    cur_files <- list.files(path, pattern = pattern, full.names = TRUE)
    cur_names <- list.files(path, pattern = pattern, full.names = FALSE)
    
    cur_dat <- lapply(cur_files, read_delim, delim = "\t", 
                        show_col_types = FALSE)
    cur_dat <- bplapply(cur_dat, function(x){
        
        if (sum(grepl(channel_pattern, colnames(x))) == 0) {
            stop("'channel_pattern' does not match", 
                    " any entries in the .txt files.")
        }
        
        cur_mat <- x[,grepl(channel_pattern, colnames(x))]
        
        if (!all(index_names %in% colnames(x))) {
            stop("'index_names' not in the names of the .txt files.")
        }
        
        cur_ind <- x[,index_names]
        
        n_row <- max(cur_ind[,index_names[1]]) + 1
        n_col <- max(cur_ind[,index_names[2]]) + 1
        
        cur_mat <- lapply(cur_mat, function(y){
            matrix(y, nrow = n_row, ncol = n_col, byrow = FALSE)
        })
        
        cur_mat <- do.call(abind, list(cur_mat, along = 3))
        
        Image(cur_mat)
    }, BPPARAM = BPPARAM)
    
    cur_dat <- CytoImageList(cur_dat)
    names(cur_dat) <- sub("\\.[^.]*$", "", basename(cur_names))
    
    channelNames(cur_dat) <- str_extract(channelNames(cur_dat), channel_pattern)
    
    return(cur_dat)
}

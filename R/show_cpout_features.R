#' @title Helper function to display all features measured by CellProfiler.
#' 
#' @description Searchable datatable object of cell and image features as 
#' extracted by CellProfiler.
#'
#' @param path full path to the CellProfiler output folder
#' @param display single character indicating which features to display.
#' Accepted entries are \code{cell_features} to display extracted single-cell
#' features or \code{image_features} to display extracted image-level features.
#' @param cell_features single character indicating the name of the file storing
#' the extracted cell features.
#' @param image_features single character indicating the name of the file storing
#' the extracted image features.
#' 
#' @return a `DT::datatable` object
#'
#' @examples
#' path <- system.file("extdata/mockData/cpout", package = "imcRtools")
#' 
#' # Display cell features 
#' show_cpout_features(path)
#' 
#' # Display image features 
#' show_cpout_features(path, display = "image_features")
#' 
#' @seealso \code{\link{read_cpout}} for the CellProfiler reader function
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#'
#' @importFrom DT datatable
#' @importFrom vroom vroom
#' @export
show_cpout_features <- function(path,
                                display = c("cell_features", "image_features"),
                                cell_features = "var_cell.csv",
                                image_features = "var_Image.csv"){
    
    display <- match.arg(display)
    
    cur_features <- vroom(paste0(path, eval(parse(text = display))), 
                          col_types = cols(), progress = FALSE)
    
    datatable(cur_features)
}

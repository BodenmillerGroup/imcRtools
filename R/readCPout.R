#' @rdname readCPout
#' @title Read IMC data from CellProfiler output folder
#' 
#' @param cells Path to file containing the features extracted per cell.
#' Rows are cells and columns are features. 
#' This file should contain columns labeled as 'MeanIntensity' which contain the mean ion counts per cell and channel.
#' This should be a .csv file.
#' @param image Path to file containin the image metadata.
#' This should be a .csv file.
#' @param relationships Path to file containing the spatial relationships between cells.
#' This should be a tab-separated file.
#' @param panel Path to file containing the spatial relationships between cells.
#' This should be a .csv file.
#' Please make sure that the rows are ordered based on the channel numbers.
#' Alternatively you can specify the column name the contains the channel numbers by setting \code{channels_name}.
#' @param mean_channel Character that indicates which columns contain the mean ion counts per cell.
#' All other columns will be stored in the \code{rowData} slot of the SingleCellExperiment.
#' @param scale_counts Logical indicating whether the counts should be scaled by the image encoding factor.
#' For example, Cell Profiler scales 16-bit images down by a factor of 2^16.
#' To obtain the original counts, the MeanIntensity values need to be rescaled.   
#' 
#' @return An object of class \code{SingleCellExperiment} containing the raw mean intensities per cell and channel in the \code{counts} slot.
#' Cells are columns and rows are channels
#' 
#' @examples
#' readCPout(system.file("extdata", "cells.csv", package = "imcRtools"),
#'           system.file("extdata", "image.csv", package = "imcRtools"),
#'           system.file("extdata", "Object_relationshipts.txt", package = "imcRtools"),
#'           system.file("extdata", "panel.csv", package = "imcRtools"))
#' 
#' @author Nils Eling \email{nils.eling@@uzh.ch}
#' 
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export

readCPout <- function(cells,
                      image_meta,
                      relationships, 
                      panel,
                      mean_channels = "MeanIntensity_FullStackFiltered",
                      scale_counts = TRUE,
                      channels_name = NULL){
  
  # Initial checks if files exist and are provided in the right format
  errors <- .fileChecks(cells, image_meta, relationships, panel)
  if(length(errors) > 0){stop(errors)}
  
  # Read in the files
  tryCatch(cells.mat <- read.csv(cells, stringsAsFactors = FALSE), 
           error = function(e){"The file containing the cell features could not be read in.\n
             Please check the file format."})
  tryCatch(image.mat <- read.csv(image_meta, stringsAsFactors = FALSE), 
           error = function(e){"The file containing the image metadata could not be read in.\n
             Please check the file format."})
  tryCatch(relationships.mat <- read.delim(relationships, stringsAsFactors = FALSE), 
           error = function(e){"The file containing the cell-cell relationships could not be read in.\n
             Please check the file format."})
  tryCatch(panel.mat <- read.csv(panel, stringsAsFactors = FALSE), 
           error = function(e){"The file containing the cell-cell relationships could not be read in.\n
             Please check the file format."})
  
  # TODO check file formats and if the necessary channels are in the files
  # Also test if the dimensions match
  
  # Obtain the mean expression channels
  image_number <- cells.mat$ImageNumber # check if this exists
  cell_number <- cells.mat$ObjectNumber # check if this exists
  cur_cells <- cells.mat[,grepl(mean_channels, colnames(cells.mat))]
  other_features <- cells.mat[,!grepl(mean_channels, colnames(cells.mat)) & 
                                !colnames(cells.mat) %in% c("ImageNumber", "ObjectNumber")]
  
  # Reorder channels in cells based on panel information
  if(!is.null(channels_name)){
    panel.channels <- paste0("c", panel.mat[,channels_name])
    # check if this exists
    cells.channels <- gsub(".*_", "", colnames(cur_cells))
    
    # Quick check if the channels match
    if(sum(is.na(match(panel.channels, cells.channels))) > 0){
      stop("The cannel information cannot be matched between the panel file and the cells files.")
    }
    
    cur_cells <- cur_cells[,match(panel.channels, cells.channels)]
  }
  else{
    panel.channels <- paste0("c", 1:nrow(panel.mat))
    cells.channels <- gsub(".*_", "", colnames(cur_cells))
    
    # Quick check if the channels match
    if(sum(is.na(match(panel.channels, cells.channels))) > 0){
      stop("The cannel information cannot be matched between the panel file and the cells files.")
    }
    
    cur_cells <- cur_cells[,match(panel.channels, cells.channels)]
  }
  
  # Build the SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=t(cur_cells)))
  
  # Build the row metadata
  row.data <- DataFrame(row.names = colnames(cur_cells),
                        panel.mat)
  
  # Build the column metadata
  col.data <- DataFrame(row.names = rownames(cur_cells),
                        image_number = image_number,
                        cell_number = cell_number,
                        other_features)
  
  # Extend column data by image metadata
  col.data$FileName_FullStack <- image.mat[col.data$image_number,"FileName_FullStack"]
  
  # Scale counts
  
  
}

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
                      scale_counts = TRUE){
  
  # Initial checks if files exist and are provided in the right format
  errors <- .fileChecks(cells, image_meta, relationships, panel)
  if(length(errors) > 0){stop(errors)}
  
  # Read in the files
  tryCatch(cells.mat <- read.csv(cells), 
           error = function(e){"The file containing the cell features could not be read in.\n
             Please check the file format."})
  tryCatch(image.mat <- read.csv(image_meta), 
           error = function(e){"The file containing the image metadata could not be read in.\n
             Please check the file format."})
  tryCatch(relationships.mat <- read.delim(relationships), 
           error = function(e){"The file containing the cell-cell relationships could not be read in.\n
             Please check the file format."})
  tryCatch(panel.mat <- read.csv(panel), 
           error = function(e){"The file containing the cell-cell relationships could not be read in.\n
             Please check the file format."})
  
  errors <- .inputChecks()
}
#' @rdname readCPout
#' @title Read IMC data from CellProfiler output folder
#' 
#' @param cells Path to file containing the features extracted per cell.
#' Rows are cells and columns are features. 
#' This file should contain columns labeled as 'MeanIntensity' which contain the mean ion counts per cell and channel.
#' @param image Path to file containin the image metadata.
#' @param relationships Path to file containing the spatial relationships between cells.
#' @param panel Path to file containing the spatial relationships between cells.
#' @param scale_counts Logical indicating whether the counts should be scaled by the image encoding factor.
#' For example, Cell Profiler scales 16-bit images down by a factor of 2^16.
#' To obtain the original counts, the MeanIntensity values need to be rescaled.   
#' 
#' @return An object of class \code{SingleCellExperiment} containing the raw mean intensities per cell and channel in the \code{counts} slot.
#' Cells are columns and rows are channels
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' @importFrom methods new
#' @importFrom dplyr mutate_all rename
#' @importFrom flowCore colnames exprs exprs<- flowSet 
#'   fsApply identifier isFCSfile keyword read.flowSet
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @export


readCPout <- function(cells,
                      image_meta,
                      relationships, 
                      panel,
                      scale_counts = TRUE)
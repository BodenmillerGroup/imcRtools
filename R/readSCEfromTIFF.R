#' @title Generates a SingleCellExperiment from .tiff files
#'
#' @description Helper function to process .tiff files created with the 
#' steinbock pipeline into a \code{\linkS4class{SingleCellExperiment}}
#' object. This function is mainly used to read-in data generated from a
#' "spillover slide". Here, each .tiff file contains the measurements of
#' multiple pixels for a single stain across all open channels.
#'
#' @param x has to be a path to a folder containing .tiff files.
#' @param image_df_path has to be a path to a images.csv file, 
#' generated with the steinbock pipeline.
#' @param panel_df_path has to be a path to a panel.csv file, 
#' generated with the steinbock pipeline.
#' @param verbose logical indicating if additional information regarding the
#' spotted and acquired masses should be shown.
#' 
#' @return returns a SCE object where pixels are stored as columns and acquired
#' channels are stored as rows.
#' 
#' @section Reading in .tiff files for spillover correction:
#' 
#' As described in the original publication, single metal spots are acquired
#' using the Hyperion imaging system. Each acquisition corresponds to one spot.
#' All acquisitions are stored in a single .mcd file and individual acquisitions
#' are stored in single .tiff files after extraction with the steinbock pipeline.
#' 
#' This function aggregates these measurements into a single
#' \code{SingleCellExperiment} object.
#' 
#' \enumerate{
#' \item \code{x} is a path:
#' All .tiff files are read in from the specified path. Here, the path
#' should indicate the location of the spillover slide measurement. Additionally,
#' the images.csv and panel.csv files generated with the steinbock pipeline 
#' must be passed. The column \code{acquisition_description} in images.csv 
#' as well as the column \code{channel} in panel.csv must contain the 
#' spotted metal isotope name in the format \code{(mt)(mass)} 
#' (e.g. \code{Sm152} for Samarium isotope with the atomic mass 152). 
#'
#' @examples
#' # Read files from path
#' path <- system.file("extdata/spillover_tiff/img", package = "imcRtools")
#' image_df_path <- system.file("extdata/spillover_tiff/images.csv", package = "imcRtools")
#' panel_df_path <- system.file("extdata/spillover_tiff/panel.csv", package = "imcRtools")
#' 
#' sce <- readSCEfromTIFF(path, image_df_path, panel_df_path) 
#' sce
#' 
#' 
#' @author Victor Ibañez (\email{victor.ibanez@@uzh.ch})
#' 
#' @references
#' \href{https://www.sciencedirect.com/science/article/pii/S1550413118306910}{Chevrier,
#' S. et al. 2017. “Compensation of Signal Spillover in Suspension and Imaging
#' Mass Cytometry.” Cell Systems 6: 612–20.}
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom EBImage readImage
#' @importFrom SummarizedExperiment colData<- rowData<-
#' @export
readSCEfromTIFF <- function(x, 
                            image_df_path,
                            panel_df_path,
                            verbose = TRUE) {
  if (!dir.exists(x)) {
    stop("Image folder path does not exist.")
  }
  if (!file.exists(image_df_path)) {
    stop("Image dataframe path does not exist.")
  }
  if (!file.exists(panel_df_path)) {
    stop("Panel dataframe path does not exist.")
  }
  
  tiff_files <- list.files(x, pattern = ".tiff$|.tif$", full.names = TRUE)
  if (length(tiff_files) == 0) {
    stop("Files could not be read in.")
  }
  
  img_df   <- read.csv(image_df_path, stringsAsFactors = FALSE)
  panel_df <- read.csv(panel_df_path, stringsAsFactors = FALSE)
  
  tiff_list <- lapply(seq_along(tiff_files), function(i) {
    img <- readImage(tiff_files[i])
    dims <- dim(img)
    
    mat <- matrix(img, nrow = dims[1]*dims[2], ncol = dims[3])
    
    rownames(mat) <- paste0(basename(img_df$acquisition_description[i]), "_", seq_len(nrow(mat)))
    
    return(mat)
  })
  
  .valid.readSCEfromTIFF.input(img_df, panel_df, verbose)
  
  cur_counts <- do.call(rbind, tiff_list)
  cell_meta <- DataFrame(cur_counts[,1]) 
  cell_meta$sample_id <- str_split(rownames(cell_meta), "\\_", 
                                   simplify = TRUE)[,1]
  cell_meta$sample_metal <- str_extract(cell_meta$sample_id, 
                                        "^[A-Z]{1}[a-z]{0,1}")
  cell_meta$sample_mass <- str_extract(cell_meta$sample_id, "[0-9]{2,3}$")
  
  cur_counts <- t(cur_counts)
  
  channel_meta <- DataFrame(
    channel_name = paste0(panel_df$channel, "Di"),
    marker_name  = panel_df$channel
  )
  
  sce <- SingleCellExperiment(assays = list(counts = cur_counts))
  rowData(sce) <- channel_meta
  colData(sce) <- cell_meta
  sce <- sce[, order(colData(sce)$sample_id, colnames(sce))]
  
  return(sce)
}
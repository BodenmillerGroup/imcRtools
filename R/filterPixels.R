#' @title Filter pixels based on their assigned masses
#'
#' @description Helper function for estimating the spillover matrix. After
#' assigning each pixel to a spotted mass, this function will filter 
#' incorrectly assigned pixels and remove small pixel sets.
#'
#' @param object a \code{SingleCellExperiment} object containing pixel
#' intensities per channel. Individual pixels are stored as columns and
#' channels are stored as rows.
#' @param bc_id character string indicating which \code{colData(object)} entry
#' stores the estimated mass 
#' @param spot_mass character string indicating which \code{colData(object)} 
#' entry stores the true isotope mass of the spotted metal. 
#' @param minevents single numeric indicating the threshold under which pixel
#' sets are excluded from spillover estimation.
#' @param correct_pixels logical indicating if incorrectly assigned pixels 
#' should be excluded from spillover estimation.
#'
#' @return returns a SingleCellExperiment object in which
#' \code{colData(object)$bc_id} has been adjusted based on the filter
#' criteria.
#'
#' @examples
#' path <- system.file("extdata/spillover", package = "imcRtools")
#' sce <- readSCEfromTXT(path)
#' assay(sce, "exprs") <- asinh(counts(sce)/5)
#' 
#' # Pre-process via CATALYST
#' library(CATALYST)
#' 
#' bc_key <- as.numeric(unique(sce$sample_mass))
#' sce <- assignPrelim(sce, bc_key = bc_key)
#' sce <- estCutoffs(sce)
#' sce <- applyCutoffs(sce)
#' 
#' sce <- filterPixels(sce)
#' 
#' table(sce$sample_mass, sce$bc_id)
#' 
#' @author Vito Zanotelli, adapted by 
#' Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#'
#' @importFrom SummarizedExperiment colData<-
#' @export
filterPixels <- function(object, 
                            bc_id = "bc_id",
                            spot_mass = "sample_mass",
                            minevents = 40,
                            correct_pixels = TRUE){
    
    .valid.filterPixels.input(object, bc_id, spot_mass, minevents, 
                                correct_pixels)

    cur_bcs <- as.character(colData(object)[[bc_id]])    
    cur_mass <- as.character(colData(object)[[spot_mass]])   

    if (correct_pixels) {
        cur_bcs[cur_bcs != "0" & cur_bcs != cur_mass] <- "0"
    }
    
    cur_stats <- table(cur_bcs)
    nonfreq <- names(cur_stats)[cur_stats < minevents]
    cur_bcs[cur_bcs %in% nonfreq] <- "0"
    
    colData(object)[[bc_id]] <- cur_bcs
    
    return(object)
}

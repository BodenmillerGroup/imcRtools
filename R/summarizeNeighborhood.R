#' @rdname summarizeNeighborhood
#' @title Summarizes cell-cell interactions within grouping levels (e.g. images)
#'
#' @description Function to calculate the average number of neighbors B that
#' a cell of type A has using different approaches. 
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object .
#' @param group_by a single character indicating the \code{colData(object)}
#' entry by which interactions are grouped. This is usually the image ID or
#' patient ID.
#' @param label single character specifying the \code{colData(object)} entry
#' which stores the cell labels. These can be cell-types labels or other m
#' metadata. 
#' @param method which cell-cell interaction counting method to use 
#' (see details)
#' @param patch_size if \code{method = "patch"}, a single numeric specifying
#' the minimum number of neighbors of the same type to be considered a patch
#' (see details)
#' @param colPairName single character indicating the \code{colPair(object)}
#' entry containing cell-cell interactions in form of an edge list.
#'
#' @section Counting and summarizing cell-cell interactions:
#' In principle, the \code{summarizeNeighborhood} function counts the number
#' of edges (interactions) between each set of unique entries in 
#' \code{colData(object)[[label]]}. Simplified, it counts for each cell of
#' type A the number of neighbors of type B.
#' 
#' This count can be averaged within each unique entry 
#' \code{colData(object)[[group_by]]} in three different ways:
#' 
#' 1. \code{method = "classic"}: The count is divided by the total number of 
#' cells of type A. The final count can be interpreted as "How many neighbors 
#' of type B does a cell of type A have on average?" 
#' 
#' 2. \code{method = "classic"}: The count is divided by the number of cells
#' of type A that have at least one neighbor of type B. The final count can be 
#' interpreted as "How many many neighbors of type B has a cell of type A on 
#' average, given it has at least one neighbor of type B?"
#' 
#' 3. \code{method = "patch"}: For each cell, the count is binarized to 0 
#' (less than \code{patch_size} neighbors of type B) or 1 (more or equal to 
#' \code{patch_size} neighbors of type B). The binarized counts are averaged 
#' across all cells of type A. The final count can be interpreted as "What 
#' fraction of cells of type A have at least a given number of neighbors of 
#' type B?"
#' 
#' @return a DataFrame containing one row per \code{group_by} entry and unique
#' label entry combination (\code{from_label}, \code{to_label}). The \code{ct}
#' entry stores the interaction count as described in the details. \code{NA}
#' is returned if a certain label is not present in this grouping level.
#'  
#' @examples 
#' library(cytomapper)
#' data(pancreasSCE)
#'
#' pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
#'                                k = 3)
#'                                
#' # Classic style calculation
#' (out <- summarizeNeighborhood(pancreasSCE, 
#'                                 group_by = "ImageNb",
#'                                 label = "CellType", 
#'                                 method = "classic",
#'                                 colPairName = "knn_interaction_graph"))
#'                                 
#' # Histocat style calculation
#' (out <- summarizeNeighborhood(pancreasSCE, 
#'                                 group_by = "ImageNb",
#'                                 label = "CellType", 
#'                                 method = "histocat",
#'                                 colPairName = "knn_interaction_graph"))
#'                                 
#' # Patch style calculation
#' (out <- summarizeNeighborhood(pancreasSCE, 
#'                                 group_by = "ImageNb",
#'                                 label = "CellType", 
#'                                 method = "pathc",
#'                                 patch_size = 3,
#'                                 colPairName = "knn_interaction_graph"))
#'
#' @author Vito Zanotelli
#' @author Jana Fischer
#' @author adapted by Nils Eling (\email{nils.eling@@uzh.ch})
#' 
#' @references
#' \href{https://www.sciencedirect.com/science/article/pii/S2405471217305434?via%3Dihub}{
#' Schulz, D. et al., Simultaneous Multiplexed Imaging of mRNA and Proteins with 
#' Subcellular Resolution in Breast Cancer Tissue Samples by Mass Cytometry., 
#' Cell Systems 2018 6(1):25-36.e5}
#' 
#' \href{https://www.nature.com/articles/nmeth.4391}{
#' Shapiro, D. et al., histoCAT: analysis of cell phenotypes and interactions in 
#' multiplex image cytometry data, Nature Methods 2017 14, p. 873–876}
#'
#' @importFrom data.table setorder
#'
#' @export
summarizeNeighborhood <- function(object, 
                                 group_by,
                                 label,
                                 colPairName,
                                 method = c("classic", "histocat", "patch"),
                                 patch_size = NULL){
    
    # Input check
    method <- match.arg(method)
    .valid.summarizeNeighborhood.input(object, group_by, label, method,
                                       patch_size, colPairName)
    
    cur_label <- colData(object)[[label]]
    cur_table <- .prepare_table(object, group_by, cur_label, colPairName)
    
    # Count interactions
    if (method == "classic") {
        cur_count <- .aggregate_classic(cur_table, object, group_by, label)
    } else if (method == "histocat") {
        cur_count <- .aggregate_histo(cur_table)
    } else if (method == "patch") {
        cur_count <- .aggregate_classic_patch(cur_table, 
                                              patch_size = patch_size,
                                              object, group_by, label)
    }
    
    setorder(cur_count, "group_by", "from_label", "to_label")
    cur_count <- as(cur_count, "DataFrame")
    
    return(cur_count)

}

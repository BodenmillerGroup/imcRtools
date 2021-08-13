#' @rdname neighborhoodPermTest
#' @title Performs perturbations to test if two cell-types interact more or 
#' less frequently than random.
#'
#' @description Cell-cell interactions are summarized in different ways and
#' the resulting count is compared to a distribution of counts arising 
#' from random permutations.
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
#' @param iter single numeric specifying the number of permutations to perform
#' @param p_threshold single numeric indicating the empirical p-value 
#' threshold at which interactions are considered to be significantly 
#' enriched or depleted per group.
#' @param BBPARAM parameters for parallelized processing. 
#'
#' @section Counting and summarizing cell-cell interactions:
#' 
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
#' of type A that have at least one neigbor of type B. The final count can be 
#' interpreted as "How many many neightbours of type B has a cell of type A on 
#' average, given it has at least one neigthbour of type B?"
#' 
#' 3. \code{method = "patch"}: For each cell, the count is binarized to 0 
#' (less than \code{patch_size} neighbors of type B) or 1 (more or equal to 
#' \code{patch_size} neighbors of type B). The binarized counts are averaged 
#' across all cells of type A. The final count can be interpreted as "What 
#' fraction of cells of type A have at least a given number of neighbors of 
#' type B?"
#' 
#' @section Testing for significance:
#' 
#' 
#'  
#' @examples 
#' # TODO
#' 
#' @seealso 
#' \code{\link{summarizeNeighborhood}} for counting (but not testing) cell-cell
#' interactions per grouping level.
#'   
#' \code{\link[BiocParallel]{bpparam}} for the parallelised backend
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
#' @export
neighborhoodPermTest <- function(object, 
                                 group_by,
                                 label,
                                 colPairName,
                                 method = c("classic", "histocat", "patch"),
                                 patch_size = NULL,
                                 iter = 1000,
                                 p_threshold = 0.01,
                                 BBPARAM = SerialParam()){

    # Input check
    .valid.summarizeNeighborhood.input(object, group_by, label, method,
                                       patch_size, colPairName)
    .valid.neighborhoodPermTest.input(iter, p_threshold)
    
    method <- match.arg(method)

    cur_label <- as.factor(colData(object)[[label]])
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
    
    # Permute the labels
    cur_out <- .permute_labels(object, group_by, label, iter, patch_size,
                               colPairName, method, BBPARAM)
    
    cur_out <- .calc_p_vals(cur_count, cur_out, n_perm = iter, 
                            p_thres = p_threshold)
    
    setorder(cur_out, "group_by", "from_label", "to_label")
    
    cur_out <- as(cur_out, "DataFrame")
    
    return(cur_out)
}

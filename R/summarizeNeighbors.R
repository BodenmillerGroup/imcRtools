#' @title Function to summarize the neighborhood from each cell.
#'
#' @description Function to generates a data frame that either summarized the  count of cells of a specified type
#' in the selected neighborhood of each cell or summarizes the marker expression of all cells in the selected neighborhood of each cell.
#'
#' @param object a \code{SingleCellExperiment} object
#' @param colPairName single character indicating the \code{colPair(object)} entry
#' containing the the neighbor information.
#' @param summarize_by character specifying whether the neighborhood should be summarized by cellular
#' features from the \code{colData(object)} or by the marker expression of the neighboring cells.
#' @param group if selecting "celltypes" in \code{summarize_by} a single character specifying the \code{colData(object)} entry containing the
#' cellular features that should be summarized in the neighborhood.
#' Supported entries are of type character.
#' @param subset_row if selecting "expression" in \code{summarize_by} a character vector specifying the entries from \code{rownames(object)}
#' to use for the summary statistics.
#' @param summaryStats if selecting "expression" in \code{summarize_by} then a single character specifying the summary statistics should be provided.
#' Supported entries are "mean", "median", "sd", "var".
#' @param name single character specifying the name of the data frame to be saved in the \code{colData(object)}.
#'
#' @return returns a \code{SingleCellExperiment}
#' containing the data frame in form of a \code{DataFrame} object in
#' \code{colData(object)}.
#'
#'
#' @author Daniel Schulz (\email{daniel.schulz@@uzh.ch})
#'
#' @importFrom dplyr as_tibble group_by count select summarise across
#' @importFrom tidyr pivot_wider
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment assay
#' @export

summarizeNeighbors <- function(object,
                               colPairName,
                               summarize_by = c("celltypes","expression"),
                               group = NULL,
                               assay_type = NULL,
                               subset_row = NULL,
                               summaryStats = NULL,
                               name = NULL){

  summarize_by = match.arg(summarize_by)

  .valid.summarizeNeighbors.input(object, colPairName, summarize_by, group, assay_type, subset_row, summaryStats)


  cur_dat <- as_tibble(colPair(object,colPairName))

  if (summarize_by == "celltypes") {

    cur_out <- cur_dat %>%
      mutate(celltype = colData(object)[[group]][cur_dat$to]) %>%
      group_by(from,.drop=FALSE) %>%
      count(celltype) %>%
      pivot_wider(names_from = celltype,values_from = n, values_fill = 0) %>%
      ungroup() %>%
      select(unique(colData(object)[[group]]))

    # add the dataframe to the object

    name <- ifelse(is.null(name), "summarizedNeighbors", name)

    colData(object)[[name]] <- DataFrame(cur_out)

    return(object)

  }

  else if (summarize_by == "expression") {

    exp_dat <- cbind(cur_dat,as_tibble(t(assay(object,assay_type)))[cur_dat$to,])

    marker_list <- list(mean = mean, median = median, sd = sd , var = var)



    cur_out <- exp_dat %>%
      select(c(from,to,subset_row)) %>%
      group_by(from,.drop=FALSE) %>%
      summarise(across(.cols = all_of(cur_markers),marker_list[summaryStats])) %>%
      select(! from)

    name <- ifelse(is.null(name), paste0(summaryStats,"_summarizedExpression"), name)

    colData(object)[[name]] <- DataFrame(cur_out)

    return(object)

  }

}

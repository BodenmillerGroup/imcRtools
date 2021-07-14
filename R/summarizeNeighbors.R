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
#' Supported entries are of type character. Defaults to celltypes.
#' @param assay_type if selecting "expression" in \code{summarize_by} name of the assay from the \code{SingleCellExperiment} to use.
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
#' @importFrom data.table as.data.table dcast melt
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment colPair colData
#' @importFrom stats median sd var
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


  if (summarize_by == "celltypes") {

    cur_dat <- as.data.table(colPair(object,colPairName))

    cur_dat$celltype <- colData(object)[[group]][cur_dat$to]

    cur_out <- DataFrame(dcast(cur_dat,formula = "from ~ celltype", fun.aggregate = length)[,-1])

    name <- ifelse(is.null(name), "summarizedNeighbors", name)

    colData(object)[[name]] <- DataFrame(cur_out)

    return(object)

  }

  else if (summarize_by == "expression") {

    cur_dat <- as.data.table(colPair(object,colPairName))

    cur_dat <- cbind(cur_dat,t(assay(object,assay_type))[cur_dat$to,])

    stat_list <- list(median = median, sd = sd , var = var)

    cur_dat <- cur_dat[,c("from","to",..subset_row)]
    cur_dat <- melt(cur_dat,measure.vars = subset_row)
    if (summaryStats == "mean"){
      cur_dat[, x := mean(value),by=c("from","variable")]
    }
    if (summaryStats %in% c("median","sd","var")) {
      cur_dat[, x := stat_list[[summaryStats]](value),by=c("from","variable")]
    }
    cur_dat <- unique(cur_dat[,c("from","variable","x")])
    cur_out <- dcast(cur_dat,formula = "from ~ variable",value.var = "x")

    name <- ifelse(is.null(name), paste0(summaryStats,"_summarizedExpression"), name)

    colData(object)[[name]] <- DataFrame(cur_out[,-1])

    return(object)

  }

}

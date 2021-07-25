#' @title Function to aggregate across all neighbors of each cell.
#'
#' @description Function to summarize categorical or expression values of all
#' neighbors of each cell.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param colPairName single character indicating the \code{colPair(object)}
#' entry containing the neighbor information.
#' @param aggregate_by character specifying whether the neighborhood should be
#' summarized by cellular features stored in \code{colData(object)}
#' (\code{aggregate_by = "metdata"}) or by the marker expression of the
#' neighboring cells (\code{aggregate_by = "expression"}).
#' @param count_by for \code{summarize_by = "metadata"}, a single character
#' specifying the \code{colData(object)} entry containing the cellular
#' metadata that should be summarized across each cell's neighborhood.
#' @param proportions single logical indicating whether aggregated metadata
#' should be returned in form of proportions instead of absolute counts.
#' @param assay_type for \code{summarize_by = "expression"}, single character
#' indicating the assay slot to use.
#' @param subset_row for \code{summarize_by = "expression"}, an integer, logical
#' or character vector specifying the features to use. If NULL, defaults to
#' all features.
#' @param statistic for \code{summarize_by = "expression"}, a single character
#' specifying the statistic to be used for summarizing the expression values
#' across all neighboring cells. Supported entries are "mean", "median", "sd",
#' "var". Defaults to "mean" if not specified.
#' @param name single character specifying the name of the data frame to be
#' saved in the \code{colData(object)}. Defaults to "aggregatedNeighbors" when
#' \code{summarize_by = "metadata"} or "{statistic}_aggregatedExpression" when
#' \code{summarize_by = "expression"}.
#'
#' @return returns an object of \code{class(object)} containing the aggregated
#' values in form of a \code{DataFrame} object in
#' \code{colData(object)[[name]]}.
#'   
#' @author Daniel Schulz (\email{daniel.schulz@@uzh.ch})
#'
#' @importFrom data.table as.data.table dcast melt :=
#' @importFrom S4Vectors DataFrame
#' @importFrom utils globalVariables
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment colPair colData
#' @importFrom stats median sd var
#' @export
aggregateNeighbors <- function(object,
                               colPairName,
                               aggregate_by = c("metadata", "expression"),
                               count_by = NULL,
                               proportions = TRUE,
                               assay_type = NULL,
                               subset_row = NULL,
                               statistic = c("mean", "median", "sd", "var"),
                               name = NULL){

    summarize_by = match.arg(aggregate_by)
  
    summaryStats = match.arg(statistic)

    .valid.summarizeNeighbors.input(object, colPairName, aggregate_by, count_by, 
                                    proportions, assay_type, subset_row,
                                    name)

    if (summarize_by == "metadata") {

        cur_dat <- as.data.table(colPair(object, colPairName))
        
        cur_dat[, "celltype" := colData(object)[[count_by]][cur_dat$to]]

        cur_dat <- dcast(cur_dat, formula = "from ~ celltype", 
                               fun.aggregate = length)[,-1]

        if (proportions) {
            .SD <- NULL
            all_col <- names(cur_dat)
            row_sums <- rowSums(cur_dat)
            cur_dat[, (all_col) := lapply(.SD, function(x){x / row_sums}), 
                    .SDcols = all_col]
        }
    
        name <- ifelse(is.null(name), "aggregatedNeighbors", name)

        colData(object)[[name]] <- DataFrame(cur_dat)

        return(object)

    } else {
      
        if (is.null(subset_row)) {
            subset_row <- rownames(object)
        }

        cur_dat <- as.data.table(colPair(object, colPairName))

        cur_dat <- cbind(cur_dat,t(assay(object, assay_type))[cur_dat$to,
                                                              subset_row])
        
        cur_dat <- melt(cur_dat, id.vars = c("from", "to"))
        
        cur_dat <- cur_dat[,eval(parse(text = paste0(statistic, "(value)"))), 
                           by=c("from","variable")]
        
        cur_out <- dcast(cur_dat, formula = "from ~ variable", value.var = "V1")

        name <- ifelse(is.null(name), 
                       paste0(statistic,"_aggregatedExpression"), name)

        colData(object)[[name]] <- DataFrame(cur_out[,-1])

        return(object)
    }
}

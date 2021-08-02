# Helper functions for the neighbourhood permutation test
#' @importFrom data.table as.data.table
.prepare_table <- function(object, img_id, cur_label, colPairName) {
    cur_tab <- as.data.frame(colPair(object, colPairName))
    cur_tab$image_id <- colData(object)[[img_id]][cur_tab$from]
    cur_tab$from_label <-cur_label[cur_tab$from]
    cur_tab$to_label <- cur_label[cur_tab$to]
    return(as.data.table(cur_tab))
}

# Helper functions for the neighbourhood permutation test
#' @importFrom data.table as.data.table
.prepare_table <- function(object, group_by, cur_label, colPairName) {
    cur_tab <- as.data.frame(colPair(object, colPairName))
    cur_tab$group_by <- colData(object)[[group_by]][cur_tab$from]
    cur_tab$from_label <-cur_label[cur_tab$from]
    cur_tab$to_label <- cur_label[cur_tab$to]
    cur_tab$ct <- 1
    return(as.data.table(cur_tab))
}

.aggregate_histo <- function(dat_table) {
    dat_temp <- dat_table[, .(ct=.N), by = c("group_by", "from_label", "to_label", "from")]
    dat_temp <- dat_temp[, .(ct=mean(ct)), by = c("group_by", "from_label", "to_label")]
    return(dat_temp)
}

.aggregate_classic <- function(dat_table){
    dat_temp <- dcast.data.table(dat_table, "group_by + from_label + from ~ to_label",
                                value.var = "ct", fun.aggregate = sum, fill=0) 
    dat_temp <- melt.data.table(dat_temp, id.vars = c("group_by", "from_label", "from"),
                        variable.name = "to_label",
                        value.name = "ct") 
    dat_temp <- dcast.data.table(dat_temp, "group_by + from_label ~ to_label",
                         value.var = "ct",
                         fun.aggregate = mean , fill = 0) 
    dat_temp <- melt.data.table(dat_temp, id.vars = c("group_by", "from_label"),
                        variable.name = "to_label",
                        value.name = "ct")
    return(dat_temp)
}

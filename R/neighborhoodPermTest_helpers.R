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

.aggregate_classic_patch <- function(dat_table, patch_size){
    dat_temp <- dcast.data.table(dat_table, "group_by + from_label + from ~ to_label",
                                value.var = "ct", fun.aggregate = sum, fill=0) 
    dat_temp <- melt.data.table(dat_temp, id.vars = c("group_by", "from_label", "from"),
                        variable.name = "to_label",
                        value.name = "ct")
    
    dat_temp[, ct := patch_size <= ct ]
    
    dat_temp <- dcast.data.table(dat_temp, "group_by + from_label ~ to_label",
                     value.var = "ct",
                     fun.aggregate = mean, fill = 0) 
    dat_temp <- melt.data.table(dat_temp, id.vars = c("group_by", "from_label"),
                                variable.name = "to_label",
                                value.name = "ct")
    return(dat_temp)
}

.permute_labels <- function(object, group_by, cur_label, iter,
                            colPairName, method, BBPARAM) {
    
    cur_out <- bplapply(seq_len(iter), 
                        function(x){
                            cur_lab_table <- data.table(label = cur_label,
                                                        group_by = colData(object)[[group_by]])
                            
                            label_perm <- cur_lab_table[ , .(label=sample(label)), by=group_by]
                            cur_perm <- .prepare_table(object, group_by, label_perm$label, colPairName)
                            
                            if (method == "classic") {
                                cur_perm <- .aggregate_classic(cur_perm)
                            } else if (method == "histocat") {
                                cur_perm <- .aggregate_histo(cur_perm)
                            } else if (method == "patch") {
                                cur_perm <- .aggregate_classic_patch(cur_perm, 
                                                                     patch_size = patch_size)
                            }
                            cur_perm$iter <- x
                            
                            return(cur_perm)
                            
                        }, BPPARAM = BBPARAM)
    
    cur_out <- do.call("rbind", cur_out)
    
    return(cur_out)
}

.calc_p_vals<- function(dat_baseline, dat_perm, n_perm, p_tresh){
    dat_perm <- merge(dat_perm, 
                      dat_baseline[, c("from_label", "to_label", "group_by", "ct")], 
                      by = c("from_label", "to_label", "group_by"),
                      suffixes = c("_perm", "_obs"), all = TRUE)
    
    dat_perm[, ':='(ct_perm = replace(ct_perm, is.na(ct_perm), 0),
                    ct_obs = replace(ct_obs, is.na(ct_obs), 0))]
    
    dat_stat = dat_perm[ , .(p_gt = ifelse(max(ct_obs) == 0, 1, (sum(ct_perm >= ct_obs) + 1) / (n_perm + 1)),
                             p_lt = (n_perm - sum(ct_perm > ct_obs) + 1) / (n_perm + 1)), 
                         by=c("group_by", "from_label", "to_label")]
    
    dat_stat[, direction := p_gt < p_lt]
    dat_stat[, p := p_gt * direction + p_lt * (direction == FALSE)]
    dat_stat[, sig := p < p_tresh]
    dat_stat[, sigval := as.integer(sig) * sign((direction-0.5))]
    
    return(dat_stat)
}

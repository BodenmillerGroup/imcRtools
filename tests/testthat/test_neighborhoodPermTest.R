# Define the original neighbouRhood functions here

IMGNR = 'ImageNumber'
OBJNR = 'ObjectNumber'
OBJNAME = 'ObjectName'
OBJID = 'ObjectID'
GROUP = 'group'
LABEL = 'label'

FIRSTOBJNAME = 'First Object Name'
FIRSTIMAGENUMBER = 'First Image Number'
FIRSTOBJNUMBER = 'First Object Number'
FIRSTOBJID = 'First Object ID'

FIRSTLABEL = 'FirstLabel'

SECONDOBJNAME = 'Second Object Name'
SECONDIMAGENUMBER = 'Second Image Number'
SECONDOBJNUMBER = 'Second Object Number'
SECONDOBJID = 'Second Object ID'
SECONDLABEL = 'SecondLabel'

RELATIONSHIP = 'Relationship'
DEFAULTRELATIONSHIP = 'Neighbors'
DEFAULTOBJNAME = 'cell'

COUNTVAR = 'ct'

prepare_tables <- function(dat_obj, dat_rel, objname=DEFAULTOBJNAME,
                           relationship = DEFAULTRELATIONSHIP,
                           col_imnr=IMGNR,
                           col_objnr=OBJNR,
                           col_label=LABEL,
                           col_group=NULL,
                           col_relationship=RELATIONSHIP,
                           col_objname=OBJNAME
){
    
    cols = c(col_imnr, col_objnr, col_label)
    
    has_objname = F
    if (col_objname %in% colnames(dat_obj)){
        cols = c(cols, col_objname)
        has_objname = T
    }
    
    has_group = F
    if (!is.null(col_group)){
        cols = c(cols, col_group)
        has_group = T
    }
    
    # subset and copy the relevant part of the data
    dat_obj <- copy(dat_obj[,cols, with=F])
    
    # set the objectname if needed
    if (has_objname == F){
        dat_obj[,(col_objname) := objname]
    } else{
        dat_obj = dat_obj[get(col_objname) == objname]
        
    }
    
    # set the group if needed
    if (has_group == F){
        col_group = GROUP
        dat_obj[,(col_group):=get(col_imnr)]
    }
    
    # rename all the columns
    setnames(dat_obj,
             c(col_imnr, col_objnr, col_label, col_group,
               col_objname),
             c(IMGNR, OBJNR, LABEL, GROUP,
               OBJNAME)
    )
    
    dat_obj[, label:=as.factor(label)]
    
    # same for the relationship table
    
    dat_rel = copy(dat_rel[(get(RELATIONSHIP) == relationship) &
                               (get(FIRSTOBJNAME) == objname)&
                               (get(SECONDOBJNAME) == objname), ][, c(FIRSTOBJNAME, FIRSTIMAGENUMBER, FIRSTOBJNUMBER,
                                                                      SECONDOBJNAME,SECONDIMAGENUMBER,  SECONDOBJNUMBER,
                                                                      RELATIONSHIP), with=F])
    
    
    
    # give new ids
    dat_obj[, (OBJID) := 1:.N]
    
    dat_rel <- merge(dat_rel, dat_obj[, c(IMGNR, OBJNR, OBJID, GROUP), with=F],
                     by.x=c(FIRSTIMAGENUMBER, FIRSTOBJNUMBER),
                     by.y=c(IMGNR, OBJNR)
    )
    setnames(dat_rel, OBJID, FIRSTOBJID)
    dat_rel <- merge(dat_rel, dat_obj[, c(IMGNR, OBJNR, OBJID), with=F],
                     by.x=c(SECONDIMAGENUMBER, SECONDOBJNUMBER),
                     by.y=c(IMGNR, OBJNR)
    )
    
    setnames(dat_rel, OBJID, SECONDOBJID)
    dat_rel[, (COUNTVAR) := 1]
    dat_rel <- dat_rel[, c(GROUP, FIRSTOBJID, SECONDOBJID, COUNTVAR), with=F]
    return(list(dat_obj, dat_rel))
}

shuffle_labels <- function(dat_labels){
    return(dat_labels[ , .(label=sample(label), ObjectID=ObjectID), by=group])
}

apply_labels <- function(dat_labels, dat_rel){
    labels = dat_labels[, get(LABEL)]
    objid = dat_labels[, get(OBJID)]
    labvec = rep(labels[1], max(objid))
    labvec[objid] = labels
    
    dat_rel[, (FIRSTLABEL) := labvec[get(FIRSTOBJID)]]
    dat_rel[, (SECONDLABEL) := labvec[get(SECONDOBJID)]]
    dat_rel
}


aggregate_classic<- function(dat_nb){
    dcast.data.table(dat_nb, paste0(GROUP, '+', FIRSTLABEL, '+`', FIRSTOBJID, '`~', SECONDLABEL),
                     value.var = COUNTVAR,
                     fun.aggregate=sum, fill=0) %>%
        melt.data.table(id.vars=c(GROUP, FIRSTLABEL, FIRSTOBJID),
                        variable.name=SECONDLABEL,
                        value.name=COUNTVAR) %>%
        dcast.data.table(paste0(GROUP, '+', FIRSTLABEL, '~', SECONDLABEL),
                         value.var = COUNTVAR,
                         fun.aggregate=mean, fill=0) %>%
        melt.data.table(id.vars=c(GROUP, FIRSTLABEL),
                        variable.name=SECONDLABEL,
                        value.name=COUNTVAR)
}

aggregate_histo <- function(dat_nb){
    dat_temp = dat_nb[, .(ct=.N), by=.(group, FirstLabel, SecondLabel, `First Object ID`)]
    dat_temp[, .(ct=mean(ct)), by=.(group, FirstLabel, SecondLabel)]
}

aggregate_classic_patch<- function(dat_nb, patch_size){
    dat_temp = dcast.data.table(dat_nb, paste0(GROUP, '+', FIRSTLABEL, '+`', FIRSTOBJID, '`~', SECONDLABEL),
                                value.var = COUNTVAR,
                                fun.aggregate=sum, fill=0) %>%
        melt.data.table(id.vars=c(GROUP, FIRSTLABEL, FIRSTOBJID),
                        variable.name=SECONDLABEL,
                        value.name=COUNTVAR)
    
    # Addition patch detection
    dat_temp[, (COUNTVAR) := patch_size <= ct ]
    
    dcast.data.table(dat_temp,paste0(GROUP, '+', FIRSTLABEL, '~', SECONDLABEL),
                     value.var = COUNTVAR,
                     fun.aggregate=mean, fill=0) %>%
        melt.data.table(id.vars=c(GROUP, FIRSTLABEL),
                        variable.name=SECONDLABEL,
                        value.name=COUNTVAR)
}

calc_p_vals<- function(dat_baseline, dat_perm, n_perm, p_tresh=0.01){
    dat_perm <-
        merge(dat_perm, dat_baseline[, c(FIRSTLABEL, SECONDLABEL, GROUP, COUNTVAR), with=F], by=c(FIRSTLABEL, SECONDLABEL, GROUP),
              suffixes = c("_perm", "_obs"),all=T)
    dat_perm[, ':='(ct_perm=replace(ct_perm, is.na(ct_perm), 0),
                    ct_obs=replace(ct_obs, is.na(ct_obs), 0)
    )]
    
    
    dat_stat = dat_perm[ , .(p_gt=ifelse(max(ct_obs)==0, 1,(sum(ct_perm>=ct_obs)+1)/(n_perm+1)),
                             p_lt=(n_perm-sum(ct_perm>ct_obs)+1)/(n_perm+1)) , by=.(group, FirstLabel, SecondLabel)]
    
    dat_stat[, direction := p_gt < p_lt]
    dat_stat[, p := p_gt * direction + p_lt * (direction == F)]
    dat_stat[, sig := p < p_tresh]
    dat_stat[, sigval := as.integer(sig)*sign((direction-0.5))]
    dat_stat
}


test_that("neighbourhoodPermTest function works", {
    # Test neighbouRhood results
    fn_cells <- system.file("extdata/mockData/cpout", "cell.csv", package = "imcRtools")
    fn_relationship <- system.file("extdata/mockData/cpout", "Object relationships.csv", package = "imcRtools")

    dat_cells <- fread(fn_cells)
    dat_relation <- fread(fn_relationship)
    
    set.seed(123)
    cur_label <- sample.int(10, size=nrow(dat_cells), replace = T)
    dat_cells[, label := cur_label]
    dat_cells[, group := ImageNumber]
    
    d <- prepare_tables(dat_cells, dat_relation)
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    histocat_aggregation <- 
    
    
    # Compare to imcRtools
    cur_path <- system.file("extdata/mockData/cpout", package = "imcRtools")
    cur_spe <- read_cpout(cur_path)
    
    cur_colpair <- as.data.frame(colPair(cur_spe, "neighbourhood"))
    cur_colpair <- cur_colpair[order(paste(cur_colpair[,1], cur_colpair[,2])),]
    cur_d <- as.data.frame(d[[2]][,c("First Object ID", "Second Object ID")])
    cur_d <- cur_d[order(paste(cur_d[,1], cur_d[,2])),]
    
    expect_equal(cur_colpair, cur_d, check.attributes = FALSE)
    
    cur_spe$label <- cur_label
    
    cur_table <- .prepare_table(cur_spe, img_id = "sample_id", 
                                cur_label  = as.factor(colData(cur_spe)[["label"]]), 
                                colPairName = "neighbourhood")
    
    labels_applied <- as.data.frame(labels_applied)
    labels_applied <- labels_applied[order(paste(labels_applied[,"First Object ID"], labels_applied[,"Second Object ID"])),]
    cur_table <- as.data.frame(cur_table)
    cur_table <- cur_table[order(paste(cur_table[,"from"], cur_table[,"to"])),]
    
    expect_equal(labels_applied$group, as.numeric(cur_table$image_id))
    expect_equal(labels_applied$`First Object ID`, cur_table$from)
    expect_equal(labels_applied$`Second Object ID`, cur_table$to)
    expect_equal(labels_applied$FirstLabel, cur_table$from_label)
    expect_equal(labels_applied$SecondLabel, cur_table$to_label),
})

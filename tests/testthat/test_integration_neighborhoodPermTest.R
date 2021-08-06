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


test_that("neighbourhoodPermTest gives same results as neighbouRhood", {
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
    
    # Compare to imcRtools
    cur_path <- system.file("extdata/mockData/cpout", package = "imcRtools")
    cur_spe <- read_cpout(cur_path)
    
    cur_colpair <- as.data.frame(colPair(cur_spe, "neighbourhood"))
    cur_colpair <- cur_colpair[order(paste(cur_colpair[,1], cur_colpair[,2])),]
    cur_d <- as.data.frame(d[[2]][,c("First Object ID", "Second Object ID")])
    cur_d <- cur_d[order(paste(cur_d[,1], cur_d[,2])),]
    
    expect_equal(cur_colpair, cur_d, check.attributes = FALSE)
    
    cur_spe$label <- cur_label
    
    cur_table <- .prepare_table(cur_spe, group_by = "sample_id", 
                                cur_label  = as.factor(colData(cur_spe)[["label"]]), 
                                colPairName = "neighbourhood")
    
    labels_applied <- as.data.frame(labels_applied)
    labels_applied <- labels_applied[order(paste(labels_applied[,"First Object ID"], labels_applied[,"Second Object ID"])),]
    cur_table <- as.data.frame(cur_table)
    cur_table <- cur_table[order(paste(cur_table[,"from"], cur_table[,"to"])),]
    
    expect_equal(labels_applied$group, as.numeric(cur_table$group_by))
    expect_equal(labels_applied$`First Object ID`, cur_table$from)
    expect_equal(labels_applied$`Second Object ID`, cur_table$to)
    expect_equal(labels_applied$FirstLabel, cur_table$from_label)
    expect_equal(labels_applied$SecondLabel, cur_table$to_label)
    
    cur_table <- .prepare_table(cur_spe, group_by = "sample_id", 
                                cur_label  = as.factor(colData(cur_spe)[["label"]]), 
                                colPairName = "neighbourhood")

    cur_histo <- .aggregate_histo(cur_table)
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    histocat_aggregation <- aggregate_histo(labels_applied)
    
    cur_histo <- as.data.frame(cur_histo)
    cur_histo <- cur_histo[order(paste(cur_histo[,1], cur_histo[,2], cur_histo[,3])),]
    histocat_aggregation <- as.data.frame(histocat_aggregation)
    histocat_aggregation <- histocat_aggregation[order(paste(histocat_aggregation[,1], histocat_aggregation[,2], histocat_aggregation[,3])),]
    
    expect_equal(cur_histo$group_by, as.character(histocat_aggregation$group))
    expect_equal(cur_histo$from_label, histocat_aggregation$FirstLabel)
    expect_equal(cur_histo$to_label, histocat_aggregation$SecondLabel)
    expect_equal(cur_histo$ct, histocat_aggregation$ct)
    
    classic_aggregation <- aggregate_classic(labels_applied)
    cur_classic <- .aggregate_classic(cur_table)
    
    cur_classic <- as.data.frame(cur_classic)
    cur_classic <- cur_classic[order(paste(cur_classic[,1], cur_classic[,2], cur_classic[,3])),]
    classic_aggregation <- as.data.frame(classic_aggregation)
    classic_aggregation <- classic_aggregation[order(paste(classic_aggregation[,1], classic_aggregation[,2], classic_aggregation[,3])),]
    
    expect_equal(cur_classic$group_by, as.character(classic_aggregation$group))
    expect_equal(cur_classic$from_label, classic_aggregation$FirstLabel)
    expect_equal(cur_classic$to_label, classic_aggregation$SecondLabel)
    expect_equal(cur_classic$ct, classic_aggregation$ct)
    
    patch_aggregation <- aggregate_classic_patch(labels_applied, patch_size = 3)
    cur_patch <- .aggregate_classic_patch(cur_table, patch_size = 3)
    
    cur_patch <- as.data.frame(cur_patch)
    cur_patch <- cur_patch[order(paste(cur_patch[,1], cur_patch[,2], cur_patch[,3])),]
    patch_aggregation <- as.data.frame(patch_aggregation)
    patch_aggregation <- patch_aggregation[order(paste(patch_aggregation[,1], patch_aggregation[,2], patch_aggregation[,3])),]
    
    expect_equal(cur_patch$group_by, as.character(patch_aggregation$group))
    expect_equal(cur_patch$from_label, patch_aggregation$FirstLabel)
    expect_equal(cur_patch$to_label, patch_aggregation$SecondLabel)
    expect_equal(cur_patch$ct, patch_aggregation$ct)
    
    # With cytomapper data
    library(cytomapper)
    data(pancreasSCE)

    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    dat_cells <- as.data.table(colData(pancreasSCE))
    dat_relation <- as.data.table(colPair(pancreasSCE))
    dat_cells[, label := CellType]
    dat_cells[, group := ImageNb]
    dat_cells[, ImageNumber := ImageNb]
    dat_cells[, ObjectNumber := CellNb]
    dat_relation[, Relationship := "Neighbors"]
    dat_relation[, "First Image Number" := colData(pancreasSCE)[["ImageNb"]][from]]
    dat_relation[, "First Object Number" := colData(pancreasSCE)[["CellNb"]][from]]
    dat_relation[, "Second Image Number" := colData(pancreasSCE)[["ImageNb"]][to]]
    dat_relation[, "Second Object Number" := colData(pancreasSCE)[["CellNb"]][to]]
    dat_relation[, "First Object Name" := "cell"]
    dat_relation[, "Second Object Name" := "cell"]
    
    d <- prepare_tables(dat_cells, dat_relation)
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    cur_classic <- aggregate_classic(labels_applied)
})

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
    
    has_objname = FALSE
    if (col_objname %in% colnames(dat_obj)){
        cols = c(cols, col_objname)
        has_objname = TRUE
    }
    
    has_group = FALSE
    if (!is.null(col_group)){
        cols = c(cols, col_group)
        has_group = TRUE
    }
    
    # subset and copy the relevant part of the data
    dat_obj <- copy(dat_obj[,cols, with=FALSE])
    
    # set the objectname if needed
    if (has_objname == FALSE){
        dat_obj[,(col_objname) := objname]
    } else{
        dat_obj = dat_obj[get(col_objname) == objname]
        
    }
    
    # set the group if needed
    if (has_group == FALSE){
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
                                                                      RELATIONSHIP), with=FALSE])
    
    
    
    # give new ids
    dat_obj[, (OBJID) := 1:.N]
    
    dat_rel <- merge(dat_rel, dat_obj[, c(IMGNR, OBJNR, OBJID, GROUP), with=FALSE],
                     by.x=c(FIRSTIMAGENUMBER, FIRSTOBJNUMBER),
                     by.y=c(IMGNR, OBJNR)
    )
    setnames(dat_rel, OBJID, FIRSTOBJID)
    dat_rel <- merge(dat_rel, dat_obj[, c(IMGNR, OBJNR, OBJID), with=FALSE],
                     by.x=c(SECONDIMAGENUMBER, SECONDOBJNUMBER),
                     by.y=c(IMGNR, OBJNR)
    )
    
    setnames(dat_rel, OBJID, SECONDOBJID)
    dat_rel[, (COUNTVAR) := 1]
    dat_rel <- dat_rel[, c(GROUP, FIRSTOBJID, SECONDOBJID, COUNTVAR), with=FALSE]
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
        merge(dat_perm, dat_baseline[, c(FIRSTLABEL, SECONDLABEL, GROUP, COUNTVAR), with=FALSE], by=c(FIRSTLABEL, SECONDLABEL, GROUP),
              suffixes = c("_perm", "_obs"),all=TRUE)
    dat_perm[, ':='(ct_perm=replace(ct_perm, is.na(ct_perm), 0),
                    ct_obs=replace(ct_obs, is.na(ct_obs), 0)
    )]
    
    
    dat_stat = dat_perm[ , .(p_gt=ifelse(max(ct_obs)==0, 1,(sum(ct_perm>=ct_obs)+1)/(n_perm+1)),
                             p_lt=(n_perm-sum(ct_perm>ct_obs)+1)/(n_perm+1)) , by=.(group, FirstLabel, SecondLabel)]
    
    dat_stat[, direction := p_gt < p_lt]
    dat_stat[, p := p_gt * direction + p_lt * (direction == FALSE)]
    dat_stat[, sig := p < p_tresh]
    dat_stat[, sigval := as.integer(sig)*sign((direction-0.5))]
    dat_stat
}


test_that("testInteractions gives same results as neighbouRhood", {
    library(data.table)
    
    # Test neighbouRhood results
    fn_cells <- system.file("extdata/mockData/cpout", "cell.csv", package = "imcRtools")
    fn_relationship <- system.file("extdata/mockData/cpout", "Object_relationships.csv", package = "imcRtools")

    dat_cells <- fread(fn_cells)
    dat_relation <- fread(fn_relationship)
    
    set.seed(1234)
    cur_label <- sample.int(10, size=nrow(dat_cells), replace = TRUE)
    dat_cells[, label := cur_label]
    dat_cells[, group := ImageNumber]
    
    d <- prepare_tables(dat_cells, dat_relation)
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    # Compare to imcRtools
    cur_path <- system.file("extdata/mockData/cpout", package = "imcRtools")
    cur_spe <- read_cpout(cur_path, graph_file = "Object_relationships.csv")
    
    cur_colpair <- as.data.frame(colPair(cur_spe, "neighborhood"))
    cur_colpair <- cur_colpair[order(paste(cur_colpair[,1], cur_colpair[,2])),]
    cur_d <- as.data.frame(d[[2]][,c("First Object ID", "Second Object ID")])
    cur_d <- cur_d[order(paste(cur_d[,1], cur_d[,2])),]
    
    expect_equal(cur_colpair, cur_d, check.attributes = FALSE)
    
    cur_spe$label <- as.factor(cur_label)
    
    cur_table <- .prepare_table(cur_spe, group_by = "sample_id", 
                                cur_label  = as.factor(colData(cur_spe)[["label"]]), 
                                colPairName = "neighborhood")
    
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
                                colPairName = "neighborhood")

    cur_histo <- .aggregate_histo(cur_table, object = cur_spe, group_by = "sample_id", label = "label")
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    histocat_aggregation <- aggregate_histo(labels_applied)
    
    setorder(histocat_aggregation, "group", "FirstLabel", "SecondLabel")
    cur_histo$group_by <- as.numeric(cur_histo$group_by)
    cur_histo <- cur_histo[histocat_aggregation,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(cur_histo, "group_by", "from_label", "to_label")
    
    expect_equal(cur_histo$group_by, histocat_aggregation$group)
    expect_equal(cur_histo$from_label, histocat_aggregation$FirstLabel)
    expect_equal(cur_histo$to_label, histocat_aggregation$SecondLabel)
    expect_equal(cur_histo$ct, histocat_aggregation$ct)
    
    classic_aggregation <- aggregate_classic(labels_applied)
    cur_classic <- .aggregate_classic(cur_table, object = cur_spe, group_by = "sample_id", label = "label")
    
    setorder(classic_aggregation, "group", "FirstLabel", "SecondLabel")
    cur_classic$group_by <- as.numeric(cur_classic$group_by)
    cur_classic <- cur_classic[classic_aggregation,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(cur_classic, "group_by", "from_label", "to_label")
    
    expect_equal(cur_classic$group_by, classic_aggregation$group)
    expect_equal(cur_classic$from_label, classic_aggregation$FirstLabel)
    expect_equal(cur_classic$to_label, classic_aggregation$SecondLabel)
    expect_equal(cur_classic$ct[!is.na(cur_classic$ct)], classic_aggregation$ct[!is.na(cur_classic$ct)])
    expect_true(all(classic_aggregation$ct[is.na(cur_classic$ct)] == 0))
    
    patch_aggregation <- aggregate_classic_patch(labels_applied, patch_size = 3)
    cur_patch <- .aggregate_classic_patch(cur_table, patch_size = 3, object = cur_spe, group_by = "sample_id", label = "label")
    
    setorder(patch_aggregation, "group", "FirstLabel", "SecondLabel")
    cur_patch$group_by <- as.numeric(cur_patch$group_by)
    cur_patch <- cur_patch[patch_aggregation,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(cur_patch, "group_by", "from_label", "to_label")
    
    expect_equal(cur_patch$group_by, patch_aggregation$group)
    expect_equal(cur_patch$from_label, patch_aggregation$FirstLabel)
    expect_equal(cur_patch$to_label, patch_aggregation$SecondLabel)
    expect_equal(cur_patch$ct[!is.na(cur_patch$ct)], patch_aggregation$ct[!is.na(cur_patch$ct)])
    expect_true(all(patch_aggregation$ct[is.na(cur_patch$ct)] == 0))
    
    # With cytomapper data
    ###################################### classic #############################
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
    setorder(cur_classic, group, FirstLabel, SecondLabel)
    
    imcRtools_classic <- countInteractions(pancreasSCE, 
                                               group_by = "ImageNb", 
                                               label = "CellType",
                                               colPairName = "knn_interaction_graph")
    
    imcRtools_classic <- as.data.table(imcRtools_classic)[cur_classic,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_classic, "group_by", "from_label", "to_label")
    
    expect_equal(cur_classic$group, imcRtools_classic$group_by)
    expect_equal(cur_classic$FirstLabel, imcRtools_classic$from_label)
    expect_equal(cur_classic$SecondLabel, imcRtools_classic$to_label)
    expect_equal(cur_classic$ct[!is.na(imcRtools_classic$ct)], imcRtools_classic$ct[!is.na(imcRtools_classic$ct)])

    # Perturbation
    n_perm <- 100
    
    dat_perm <- bplapply(1:n_perm, function(x){
        dat_labels = shuffle_labels(d[[1]])
        apply_labels(dat_labels, d[[2]]) %>%
            aggregate_classic()
    }, BPPARAM = SerialParam(RNGseed = 123))
    dat_perm <- rbindlist(dat_perm, idcol = 'run')
    
    dat_p <- calc_p_vals(cur_classic, dat_perm, n_perm = n_perm, p_tresh = 0.01) 
    setorder(dat_p, group, FirstLabel, SecondLabel)
    dat_p$group <- as.character(dat_p$group)
    
    imcRtools_classic_perm <- testInteractions(pancreasSCE, 
                                                   group_by = "ImageNb", 
                                                   label = "CellType",
                                                   colPairName = "knn_interaction_graph",
                                                   iter = 100,
                                               BPPARAM = SerialParam(RNGseed = 123))
    
    imcRtools_classic_perm <- as.data.table(imcRtools_classic_perm)[dat_p,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_classic_perm, "group_by", "from_label", "to_label")
    
    expect_equal(dat_p$group, imcRtools_classic_perm$group_by)
    expect_equal(as.character(dat_p$FirstLabel), imcRtools_classic_perm$from_label)
    expect_equal(as.character(dat_p$SecondLabel), imcRtools_classic_perm$to_label)
    expect_equal(dat_p$p_gt[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$p_gt[!is.na(imcRtools_classic$ct)])
    expect_equal(dat_p$p_lt[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$p_lt[!is.na(imcRtools_classic$ct)])
    expect_equal(dat_p$direction[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$interaction[!is.na(imcRtools_classic$ct)])
    expect_equal(dat_p$p[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$p[!is.na(imcRtools_classic$ct)])
    expect_equal(dat_p$sig[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$sig[!is.na(imcRtools_classic$ct)])
    expect_equal(dat_p$sigval[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$sigval[!is.na(imcRtools_classic$ct)])
    
    ###################################### histocat #############################
    d <- prepare_tables(dat_cells, dat_relation)
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    cur_histo <- aggregate_histo(labels_applied)
    setorder(cur_histo, group, FirstLabel, SecondLabel)
    
    imcRtools_histo <- countInteractions(pancreasSCE, 
                                               group_by = "ImageNb", 
                                               label = "CellType",
                                               colPairName = "knn_interaction_graph",
                                               method = "histocat")
    
    imcRtools_histo <- as.data.table(imcRtools_histo)[cur_histo,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_histo, "group_by", "from_label", "to_label")
    
    expect_equal(cur_histo$group, imcRtools_histo$group_by)
    expect_equal(cur_histo$FirstLabel, imcRtools_histo$from_label)
    expect_equal(cur_histo$SecondLabel, imcRtools_histo$to_label)
    expect_equal(cur_histo$ct, imcRtools_histo$ct)
    
    # Perturbation
    n_perm <- 100
    
    dat_perm <- bplapply(1:n_perm, function(x){
        dat_labels = shuffle_labels(d[[1]])
        apply_labels(dat_labels, d[[2]]) %>%
            aggregate_histo()
    }, BPPARAM = SerialParam(RNGseed = 123))
    dat_perm <- rbindlist(dat_perm, idcol = 'run')
    
    dat_p <- calc_p_vals(cur_histo, dat_perm, n_perm = n_perm, p_tresh = 0.01) 
    setorder(dat_p, group, FirstLabel, SecondLabel)
    dat_p$group <- as.character(dat_p$group)
    
    imcRtools_histo_perm <- testInteractions(pancreasSCE, 
                                            group_by = "ImageNb", 
                                            label = "CellType",
                                            colPairName = "knn_interaction_graph",
                                            iter = 100,
                                            method = "histocat", 
                                            BPPARAM = SerialParam(RNGseed = 123))
    
    imcRtools_histo_perm <- as.data.table(imcRtools_histo_perm)[dat_p,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_histo_perm, "group_by", "from_label", "to_label")
    
    expect_equal(dat_p$group, imcRtools_histo_perm$group_by)
    expect_equal(as.character(dat_p$FirstLabel), imcRtools_histo_perm$from_label)
    expect_equal(as.character(dat_p$SecondLabel), imcRtools_histo_perm$to_label)
    expect_equal(dat_p$p_gt, imcRtools_histo_perm$p_gt)
    expect_equal(dat_p$p_lt, imcRtools_histo_perm$p_lt)
    expect_equal(dat_p$direction, imcRtools_histo_perm$interaction)
    expect_equal(dat_p$p, imcRtools_histo_perm$p)
    expect_equal(dat_p$sig, imcRtools_histo_perm$sig)
    expect_equal(dat_p$sigval, imcRtools_histo_perm$sigval)
    
    ###################################### patch #############################
    d <- prepare_tables(dat_cells, dat_relation)
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    cur_patch <- aggregate_classic_patch(labels_applied, patch_size = 3)
    setorder(cur_patch, group, FirstLabel, SecondLabel)
    
    imcRtools_patch <- countInteractions(pancreasSCE, 
                                             group_by = "ImageNb", 
                                             label = "CellType",
                                             colPairName = "knn_interaction_graph",
                                             method = "patch",
                                             patch_size = 3)
    
    imcRtools_patch <- as.data.table(imcRtools_patch)[cur_patch,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_patch, "group_by", "from_label", "to_label")
    
    expect_equal(cur_patch$group, imcRtools_patch$group_by)
    expect_equal(cur_patch$FirstLabel, imcRtools_patch$from_label)
    expect_equal(cur_patch$SecondLabel, imcRtools_patch$to_label)
    expect_equal(cur_patch$ct[!is.na(imcRtools_patch$ct)], imcRtools_patch$ct[!is.na(imcRtools_patch$ct)])
    
    # Perturbation
    n_perm <- 100
    
    dat_perm <- bplapply(1:n_perm, function(x){
        dat_labels = shuffle_labels(d[[1]])
        apply_labels(dat_labels, d[[2]]) %>%
            aggregate_classic_patch(patch_size = 3)
    }, BPPARAM = SerialParam(RNGseed = 123))
    dat_perm <- rbindlist(dat_perm, idcol = 'run')
    
    dat_p <- calc_p_vals(cur_patch, dat_perm, n_perm = n_perm, p_tresh = 0.01) 
    setorder(dat_p, group, FirstLabel, SecondLabel)
    dat_p$group <- as.character(dat_p$group)
    
    imcRtools_patch_perm <- testInteractions(pancreasSCE, 
                                                 group_by = "ImageNb", 
                                                 label = "CellType",
                                                 colPairName = "knn_interaction_graph",
                                                 iter = 100,
                                                 method = "patch",
                                                 patch_size = 3,
                                             BPPARAM = SerialParam(RNGseed = 123))
    
    imcRtools_patch_perm <- as.data.table(imcRtools_patch_perm)[dat_p,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_patch_perm, "group_by", "from_label", "to_label")
    
    expect_equal(dat_p$group, imcRtools_patch_perm$group_by)
    expect_equal(as.character(dat_p$FirstLabel), imcRtools_patch_perm$from_label)
    expect_equal(as.character(dat_p$SecondLabel), imcRtools_patch_perm$to_label)
    expect_equal(dat_p$p_gt[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$p_gt[!is.na(imcRtools_patch_perm$ct)])
    expect_equal(dat_p$p_lt[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$p_lt[!is.na(imcRtools_patch_perm$ct)])
    expect_equal(dat_p$direction[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$interaction[!is.na(imcRtools_patch_perm$ct)])
    expect_equal(dat_p$p[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$p[!is.na(imcRtools_patch_perm$ct)])
    expect_equal(dat_p$sig[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$sig[!is.na(imcRtools_patch_perm$ct)])
    expect_equal(dat_p$sigval[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$sigval[!is.na(imcRtools_patch_perm$ct)])
    
    # Corner case
    cur_sce <- pancreasSCE
    cur_sce$CellType[123] <- "test"
    cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                 type = "expansion", threshold = 7)
    
    plotSpatial(cur_sce, node_color_by = "CellType", img_id = "ImageNb", 
                draw_edges = TRUE, colPairName = "expansion_interaction_graph")
    
    dat_cells <- as.data.table(colData(cur_sce))
    dat_relation <- as.data.table(colPair(cur_sce, "expansion_interaction_graph"))
    dat_cells[, label := CellType]
    dat_cells[, group := ImageNb]
    dat_cells[, ImageNumber := ImageNb]
    dat_cells[, ObjectNumber := CellNb]
    dat_relation[, Relationship := "Neighbors"]
    dat_relation[, "First Image Number" := colData(cur_sce)[["ImageNb"]][from]]
    dat_relation[, "First Object Number" := colData(cur_sce)[["CellNb"]][from]]
    dat_relation[, "Second Image Number" := colData(cur_sce)[["ImageNb"]][to]]
    dat_relation[, "Second Object Number" := colData(cur_sce)[["CellNb"]][to]]
    dat_relation[, "First Object Name" := "cell"]
    dat_relation[, "Second Object Name" := "cell"]
    
    d <- prepare_tables(dat_cells, dat_relation)
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    cur_classic <- aggregate_classic(labels_applied)
    setorder(cur_classic, group, FirstLabel, SecondLabel)
    
    imcRtools_classic <- countInteractions(cur_sce, 
                                           group_by = "ImageNb", 
                                           label = "CellType",
                                           colPairName = "expansion_interaction_graph")
    
    imcRtools_classic <- as.data.table(imcRtools_classic)[cur_classic,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_classic, "group_by", "from_label", "to_label")
    
    expect_equal(cur_classic$group, imcRtools_classic$group_by)
    expect_equal(cur_classic$FirstLabel, imcRtools_classic$from_label)
    expect_equal(cur_classic$SecondLabel, imcRtools_classic$to_label)
    #expect_equal(cur_classic$ct[!is.na(imcRtools_classic$ct)], imcRtools_classic$ct[!is.na(imcRtools_classic$ct)])
    
    # Perturbation
    n_perm <- 100
    
    dat_perm <- bplapply(1:n_perm, function(x){
        dat_labels = shuffle_labels(d[[1]])
        apply_labels(dat_labels, d[[2]]) %>%
            aggregate_classic()
    }, BPPARAM = SerialParam(RNGseed = 123))
    dat_perm <- rbindlist(dat_perm, idcol = 'run')
    
    dat_p <- calc_p_vals(cur_classic, dat_perm, n_perm = n_perm, p_tresh = 0.01) 
    setorder(dat_p, group, FirstLabel, SecondLabel)
    dat_p$group <- as.character(dat_p$group)
    
    imcRtools_classic_perm <- testInteractions(cur_sce, 
                                               group_by = "ImageNb", 
                                               label = "CellType",
                                               colPairName = "expansion_interaction_graph",
                                               iter = 100,
                                               BPPARAM = SerialParam(RNGseed = 123))
    
    imcRtools_classic_perm <- as.data.table(imcRtools_classic_perm)[dat_p,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_classic_perm, "group_by", "from_label", "to_label")
    
    expect_equal(dat_p$group, imcRtools_classic_perm$group_by)
    expect_equal(as.character(dat_p$FirstLabel), imcRtools_classic_perm$from_label)
    expect_equal(as.character(dat_p$SecondLabel), imcRtools_classic_perm$to_label)
    #expect_equal(dat_p$p_gt[!is.na(imcRtools_classic_perm$ct)], imcRtools_classic_perm$p_gt[!is.na(imcRtools_classic_perm$ct)])
    #expect_equal(dat_p$p_lt[!is.na(imcRtools_classic_perm$ct)], imcRtools_classic_perm$p_lt[!is.na(imcRtools_classic_perm$ct)])
    expect_equal(dat_p$direction[!is.na(imcRtools_classic_perm$ct)], imcRtools_classic_perm$interaction[!is.na(imcRtools_classic_perm$ct)])
    #expect_equal(dat_p$p[!is.na(imcRtools_classic_perm$ct)], imcRtools_classic_perm$p[!is.na(imcRtools_classic_perm$ct)])
    #expect_equal(dat_p$sig[!is.na(imcRtools_classic_perm$ct)], imcRtools_classic_perm$sig[!is.na(imcRtools_classic_perm$ct)])
    #expect_equal(dat_p$sigval[!is.na(imcRtools_classic_perm$ct)], imcRtools_classic_perm$sigval[!is.na(imcRtools_classic_perm$ct)])
    
    # histocat
    d <- prepare_tables(dat_cells, dat_relation)
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    cur_histo <- aggregate_histo(labels_applied)
    setorder(cur_histo, group, FirstLabel, SecondLabel)
    
    imcRtools_histo <- countInteractions(cur_sce, 
                                           group_by = "ImageNb", 
                                           label = "CellType",
                                           method = "histocat",
                                           colPairName = "expansion_interaction_graph")
    
    imcRtools_histo <- as.data.table(imcRtools_histo)[cur_histo,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_histo, "group_by", "from_label", "to_label")
    
    expect_equal(cur_histo$group, imcRtools_histo$group_by)
    expect_equal(cur_histo$FirstLabel, imcRtools_histo$from_label)
    expect_equal(cur_histo$SecondLabel, imcRtools_histo$to_label)
    expect_equal(cur_histo$ct[!is.na(imcRtools_histo$ct)], imcRtools_histo$ct[!is.na(imcRtools_histo$ct)])
    
    # Perturbation
    n_perm <- 100
    
    dat_perm <- bplapply(1:n_perm, function(x){
        dat_labels = shuffle_labels(d[[1]])
        apply_labels(dat_labels, d[[2]]) %>%
            aggregate_histo()
    }, BPPARAM = SerialParam(RNGseed = 123))
    dat_perm <- rbindlist(dat_perm, idcol = 'run')
    
    dat_p <- calc_p_vals(cur_histo, dat_perm, n_perm = n_perm, p_tresh = 0.01) 
    setorder(dat_p, group, FirstLabel, SecondLabel)
    dat_p$group <- as.character(dat_p$group)
    
    imcRtools_histo_perm <- testInteractions(cur_sce, 
                                               group_by = "ImageNb", 
                                               label = "CellType",
                                               method = "histocat",
                                               colPairName = "expansion_interaction_graph",
                                               iter = 100,
                                               BPPARAM = SerialParam(RNGseed = 123))
    
    imcRtools_histo_perm <- as.data.table(imcRtools_histo_perm)[dat_p,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_histo_perm, "group_by", "from_label", "to_label")
    
    expect_equal(dat_p$group, imcRtools_histo_perm$group_by)
    expect_equal(as.character(dat_p$FirstLabel), imcRtools_histo_perm$from_label)
    expect_equal(as.character(dat_p$SecondLabel), imcRtools_histo_perm$to_label)
    expect_equal(dat_p$p_gt[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$p_gt[!is.na(imcRtools_histo_perm$ct)])
    expect_equal(dat_p$p_lt[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$p_lt[!is.na(imcRtools_histo_perm$ct)])
    expect_equal(dat_p$direction[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$interaction[!is.na(imcRtools_histo_perm$ct)])
    expect_equal(dat_p$p[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$p[!is.na(imcRtools_histo_perm$ct)])
    expect_equal(dat_p$sig[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$sig[!is.na(imcRtools_histo_perm$ct)])
    expect_equal(dat_p$sigval[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$sigval[!is.na(imcRtools_histo_perm$ct)])
    
    ###################################### patch ##############################
    d <- prepare_tables(dat_cells, dat_relation)
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    cur_classic_patch <- aggregate_classic_patch(labels_applied, patch_size = 1)
    setorder(cur_classic_patch, group, FirstLabel, SecondLabel)
    
    imcRtools_classic_patch <- countInteractions(cur_sce, 
                                           group_by = "ImageNb", 
                                           label = "CellType",
                                           method = "patch",
                                           patch_size = 1,
                                           colPairName = "expansion_interaction_graph")
    
    imcRtools_classic_patch <- as.data.table(imcRtools_classic_patch)[cur_classic_patch,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_classic_patch, "group_by", "from_label", "to_label")
    
    expect_equal(cur_classic_patch$group, imcRtools_classic_patch$group_by)
    expect_equal(cur_classic_patch$FirstLabel, imcRtools_classic_patch$from_label)
    expect_equal(cur_classic_patch$SecondLabel, imcRtools_classic_patch$to_label)
    #expect_equal(cur_classic_patch$ct[!is.na(imcRtools_classic_patch$ct)], imcRtools_classic_patch$ct[!is.na(imcRtools_classic_patch$ct)])
    
    # Perturbation
    n_perm <- 100
    
    dat_perm <- bplapply(1:n_perm, function(x){
        dat_labels = shuffle_labels(d[[1]])
        apply_labels(dat_labels, d[[2]]) %>%
            aggregate_classic_patch(patch_size = 1)
    }, BPPARAM = SerialParam(RNGseed = 123))
    dat_perm <- rbindlist(dat_perm, idcol = 'run')
    
    dat_p <- calc_p_vals(cur_classic_patch, dat_perm, n_perm = n_perm, p_tresh = 0.01) 
    setorder(dat_p, group, FirstLabel, SecondLabel)
    dat_p$group <- as.character(dat_p$group)
    
    imcRtools_classic_patch_perm <- testInteractions(cur_sce, 
                                               group_by = "ImageNb", 
                                               label = "CellType",
                                               colPairName = "expansion_interaction_graph",
                                               iter = 100,
                                               method = "patch",
                                               patch_size = 1,
                                               BPPARAM = SerialParam(RNGseed = 123))
    
    imcRtools_classic_patch_perm <- as.data.table(imcRtools_classic_patch_perm)[dat_p,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_classic_patch_perm, "group_by", "from_label", "to_label")
    
    expect_equal(dat_p$group, imcRtools_classic_patch_perm$group_by)
    expect_equal(as.character(dat_p$FirstLabel), imcRtools_classic_patch_perm$from_label)
    expect_equal(as.character(dat_p$SecondLabel), imcRtools_classic_patch_perm$to_label)
    #expect_equal(dat_p$p_gt[!is.na(imcRtools_classic_patch_perm$ct)], imcRtools_classic_patch_perm$p_gt[!is.na(imcRtools_classic_patch_perm$ct)])
    #expect_equal(dat_p$p_lt[!is.na(imcRtools_classic_patch_perm$ct)], imcRtools_classic_patch_perm$p_lt[!is.na(imcRtools_classic_patch_perm$ct)])
    expect_equal(dat_p$direction[!is.na(imcRtools_classic_patch_perm$ct)], imcRtools_classic_patch_perm$interaction[!is.na(imcRtools_classic_patch_perm$ct)])
    #expect_equal(dat_p$p[!is.na(imcRtools_classic_patch_perm$ct)], imcRtools_classic_patch_perm$p[!is.na(imcRtools_classic_patch_perm$ct)])
    expect_equal(dat_p$sig[!is.na(imcRtools_classic_patch_perm$ct)], imcRtools_classic_patch_perm$sig[!is.na(imcRtools_classic_patch_perm$ct)])
    expect_equal(dat_p$sigval[!is.na(imcRtools_classic_patch_perm$ct)], imcRtools_classic_patch_perm$sigval[!is.na(imcRtools_classic_patch_perm$ct)])
    
})

test_that("testInteractions gives same results as neighbouRhood when using different grouping", {
    library(data.table)
    library(cytomapper)
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    dat_cells <- as.data.table(colData(pancreasSCE))
    dat_relation <- as.data.table(colPair(pancreasSCE))
    dat_cells[, label := CellType]
    dat_cells[, group := Pattern]
    dat_cells[, ImageNumber := ImageNb]
    dat_cells[, ObjectNumber := CellNb]
    dat_relation[, Relationship := "Neighbors"]
    dat_relation[, "First Image Number" := colData(pancreasSCE)[["ImageNb"]][from]]
    dat_relation[, "First Object Number" := colData(pancreasSCE)[["CellNb"]][from]]
    dat_relation[, "Second Image Number" := colData(pancreasSCE)[["ImageNb"]][to]]
    dat_relation[, "Second Object Number" := colData(pancreasSCE)[["CellNb"]][to]]
    dat_relation[, "First Object Name" := "cell"]
    dat_relation[, "Second Object Name" := "cell"]
    
    d <- prepare_tables(dat_cells, dat_relation, col_group = "Pattern")
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    cur_classic <- aggregate_classic(labels_applied)
    setorder(cur_classic, group, FirstLabel, SecondLabel)
    
    imcRtools_classic <- countInteractions(pancreasSCE, 
                                           group_by = "Pattern", 
                                           label = "CellType",
                                           colPairName = "knn_interaction_graph")
    
    imcRtools_classic <- as.data.table(imcRtools_classic)[cur_classic,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_classic, "group_by", "from_label", "to_label")
    
    expect_equal(cur_classic$group, imcRtools_classic$group_by)
    expect_equal(cur_classic$FirstLabel, imcRtools_classic$from_label)
    expect_equal(cur_classic$SecondLabel, imcRtools_classic$to_label)
    expect_equal(cur_classic$ct[!is.na(imcRtools_classic$ct)], imcRtools_classic$ct[!is.na(imcRtools_classic$ct)])
    
    # Comment: if the to cell is located outside of the group it's still being counted
    
    # Perturbation
    n_perm <- 100
    
    dat_perm <- bplapply(1:n_perm, function(x){
        dat_labels = shuffle_labels(d[[1]])
        apply_labels(dat_labels, d[[2]]) %>%
            aggregate_classic()
    }, BPPARAM = SerialParam(RNGseed = 123))
    dat_perm <- rbindlist(dat_perm, idcol = 'run')
    
    dat_p <- calc_p_vals(cur_classic, dat_perm, n_perm = n_perm, p_tresh = 0.01) 
    setorder(dat_p, group, FirstLabel, SecondLabel)
    dat_p$group <- as.character(dat_p$group)
    
    imcRtools_classic_perm <- testInteractions(pancreasSCE, 
                                               group_by = "Pattern", 
                                               label = "CellType",
                                               colPairName = "knn_interaction_graph",
                                               iter = 100,
                                               BPPARAM = SerialParam(RNGseed = 123))
    
    imcRtools_classic_perm <- as.data.table(imcRtools_classic_perm)[dat_p,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_classic_perm, "group_by", "from_label", "to_label")
    
    expect_equal(dat_p$group, imcRtools_classic_perm$group_by)
    expect_equal(as.character(dat_p$FirstLabel), imcRtools_classic_perm$from_label)
    expect_equal(as.character(dat_p$SecondLabel), imcRtools_classic_perm$to_label)
    expect_equal(dat_p$p_gt[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$p_gt[!is.na(imcRtools_classic$ct)])
    expect_equal(dat_p$p_lt[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$p_lt[!is.na(imcRtools_classic$ct)])
    expect_equal(dat_p$direction[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$interaction[!is.na(imcRtools_classic$ct)])
    expect_equal(dat_p$p[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$p[!is.na(imcRtools_classic$ct)])
    expect_equal(dat_p$sig[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$sig[!is.na(imcRtools_classic$ct)])
    expect_equal(dat_p$sigval[!is.na(imcRtools_classic$ct)], imcRtools_classic_perm$sigval[!is.na(imcRtools_classic$ct)])
    
    ###################################### histocat #############################
    d <- prepare_tables(dat_cells, dat_relation, col_group = "Pattern")
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    cur_histo <- aggregate_histo(labels_applied)
    setorder(cur_histo, group, FirstLabel, SecondLabel)
    
    imcRtools_histo <- countInteractions(pancreasSCE, 
                                         group_by = "Pattern", 
                                         label = "CellType",
                                         colPairName = "knn_interaction_graph",
                                         method = "histocat")
    
    imcRtools_histo <- as.data.table(imcRtools_histo)[cur_histo,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_histo, "group_by", "from_label", "to_label")
    
    expect_equal(cur_histo$group, imcRtools_histo$group_by)
    expect_equal(cur_histo$FirstLabel, imcRtools_histo$from_label)
    expect_equal(cur_histo$SecondLabel, imcRtools_histo$to_label)
    expect_equal(cur_histo$ct[!is.na(imcRtools_histo$ct)], imcRtools_histo$ct[!is.na(imcRtools_histo$ct)])
    
    # Perturbation
    n_perm <- 100
    
    dat_perm <- bplapply(1:n_perm, function(x){
        dat_labels = shuffle_labels(d[[1]])
        apply_labels(dat_labels, d[[2]]) %>%
            aggregate_histo()
    }, BPPARAM = SerialParam(RNGseed = 123))
    dat_perm <- rbindlist(dat_perm, idcol = 'run')
    
    dat_p <- calc_p_vals(cur_histo, dat_perm, n_perm = n_perm, p_tresh = 0.01) 
    setorder(dat_p, group, FirstLabel, SecondLabel)
    dat_p$group <- as.character(dat_p$group)
    
    imcRtools_histo_perm <- testInteractions(pancreasSCE, 
                                             group_by = "Pattern", 
                                             label = "CellType",
                                             colPairName = "knn_interaction_graph",
                                             iter = 100,
                                             method = "histocat", 
                                             BPPARAM = SerialParam(RNGseed = 123))
    
    imcRtools_histo_perm <- as.data.table(imcRtools_histo_perm)[dat_p,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_histo_perm, "group_by", "from_label", "to_label")
    
    expect_equal(dat_p$group, imcRtools_histo_perm$group_by)
    expect_equal(as.character(dat_p$FirstLabel), imcRtools_histo_perm$from_label)
    expect_equal(as.character(dat_p$SecondLabel), imcRtools_histo_perm$to_label)
    expect_equal(dat_p$p_gt[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$p_gt[!is.na(imcRtools_histo_perm$ct)])
    expect_equal(dat_p$p_lt[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$p_lt[!is.na(imcRtools_histo_perm$ct)])
    expect_equal(dat_p$direction[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$interaction[!is.na(imcRtools_histo_perm$ct)])
    expect_equal(dat_p$p[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$p[!is.na(imcRtools_histo_perm$ct)])
    expect_equal(dat_p$sig[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$sig[!is.na(imcRtools_histo_perm$ct)])
    expect_equal(dat_p$sigval[!is.na(imcRtools_histo_perm$ct)], imcRtools_histo_perm$sigval[!is.na(imcRtools_histo_perm$ct)])
    
    ###################################### patch #############################
    d <- prepare_tables(dat_cells, dat_relation, col_group = "Pattern")
    
    labels_applied <- apply_labels(d[[1]], d[[2]])
    
    cur_patch <- aggregate_classic_patch(labels_applied, patch_size = 3)
    setorder(cur_patch, group, FirstLabel, SecondLabel)
    
    imcRtools_patch <- countInteractions(pancreasSCE, 
                                         group_by = "Pattern", 
                                         label = "CellType",
                                         colPairName = "knn_interaction_graph",
                                         method = "patch",
                                         patch_size = 3)
    
    imcRtools_patch <- as.data.table(imcRtools_patch)[cur_patch,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_patch, "group_by", "from_label", "to_label")
    
    expect_equal(cur_patch$group, imcRtools_patch$group_by)
    expect_equal(cur_patch$FirstLabel, imcRtools_patch$from_label)
    expect_equal(cur_patch$SecondLabel, imcRtools_patch$to_label)
    expect_equal(cur_patch$ct[!is.na(imcRtools_patch$ct)], imcRtools_patch$ct[!is.na(imcRtools_patch$ct)])
    
    # Perturbation
    n_perm <- 100
    
    dat_perm <- bplapply(1:n_perm, function(x){
        dat_labels = shuffle_labels(d[[1]])
        apply_labels(dat_labels, d[[2]]) %>%
            aggregate_classic_patch(patch_size = 3)
    }, BPPARAM = SerialParam(RNGseed = 123))
    dat_perm <- rbindlist(dat_perm, idcol = 'run')
    
    dat_p <- calc_p_vals(cur_patch, dat_perm, n_perm = n_perm, p_tresh = 0.01) 
    setorder(dat_p, group, FirstLabel, SecondLabel)
    dat_p$group <- as.character(dat_p$group)
    
    imcRtools_patch_perm <- testInteractions(pancreasSCE, 
                                             group_by = "Pattern", 
                                             label = "CellType",
                                             colPairName = "knn_interaction_graph",
                                             iter = 100,
                                             method = "patch",
                                             patch_size = 3,
                                             BPPARAM = SerialParam(RNGseed = 123))
    
    imcRtools_patch_perm <- as.data.table(imcRtools_patch_perm)[dat_p,, on = c("group_by==group", "from_label==FirstLabel", "to_label==SecondLabel")]
    setorder(imcRtools_patch_perm, "group_by", "from_label", "to_label")
    
    expect_equal(dat_p$group, imcRtools_patch_perm$group_by)
    expect_equal(as.character(dat_p$FirstLabel), imcRtools_patch_perm$from_label)
    expect_equal(as.character(dat_p$SecondLabel), imcRtools_patch_perm$to_label)
    expect_equal(dat_p$p_gt[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$p_gt[!is.na(imcRtools_patch_perm$ct)])
    expect_equal(dat_p$p_lt[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$p_lt[!is.na(imcRtools_patch_perm$ct)])
    expect_equal(dat_p$direction[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$interaction[!is.na(imcRtools_patch_perm$ct)])
    expect_equal(dat_p$p[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$p[!is.na(imcRtools_patch_perm$ct)])
    expect_equal(dat_p$sig[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$sig[!is.na(imcRtools_patch_perm$ct)])
    expect_equal(dat_p$sigval[!is.na(imcRtools_patch_perm$ct)], imcRtools_patch_perm$sigval[!is.na(imcRtools_patch_perm$ct)])
    
    
})

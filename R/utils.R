############################# Reader helpers ###################################

#### read_steinbock helpers ####

# Helper functions for reading in steinbock data
#' @importFrom vroom vroom
#' @importFrom BiocParallel bplapply
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom readr cols
.steinbock_read_intensities <- function(x, cell_id, return_as, BPPARAM){

    cur_out <-  bplapply(seq_along(x),
                         function(y){
                             cur_int <- vroom(x[y], progress = FALSE,
                                              col_types = cols())
                             cur_counts <- cur_int %>% select(-all_of(cell_id))

                             cur_name <- sub("\\.[^.]*$", "", basename(x[y]))

                             if (return_as == "spe") {
                                 object <- SpatialExperiment(
                                     assays = list(counts = t(as.matrix(cur_counts))),
                                     sample_id = cur_name)
                                 object$ObjectNumber <- cur_int[[cell_id]]
                             } else {
                                 object <- SingleCellExperiment(
                                     assays = list(counts = t(as.matrix(cur_counts))))
                                 object$sample_id <- cur_name
                                 object$ObjectNumber <- cur_int[[cell_id]]
                             }

                             return(object)

                         }, BPPARAM = BPPARAM)

    return(cur_out)
}

.steinbock_read_regionprops <- function(x, cur_path, cell_id, coords,
                                        return_as, BPPARAM){

    cur_out <-  bplapply(x,
                         function(y){
                             cur_sample <- unique(y$sample_id)

                             cur_file <- list.files(cur_path,
                                                    pattern = paste0("^", cur_sample, ".csv", "$"),
                                                    full.names = TRUE)

                             if (length(cur_file) == 0) {
                                 return(y)
                             }

                             cur_props <- vroom(cur_file,
                                                progress = FALSE,
                                                col_types = cols()) %>%
                                 as.data.frame()
                             rownames(cur_props) <- cur_props[[cell_id]]
                             cur_props <- cur_props[as.character(y$ObjectNumber),]

                             if (return_as == "spe") {
                                 spatialCoords(y) <- matrix(c(cur_props[[coords[1]]],
                                                              cur_props[[coords[2]]]),
                                                            ncol = 2, byrow = FALSE,
                                                            dimnames = list(as.character(y$ObjectNumber),
                                                                            c("Pos_X", "Pos_Y")))
                                 colData(y) <- cbind(colData(y),
                                                     cur_props[,!colnames(cur_props) %in%
                                                                   c(cell_id, coords)])
                             } else {
                                 colData(y)$Pos_X <- cur_props[[coords[1]]]
                                 colData(y)$Pos_Y <- cur_props[[coords[2]]]
                                 colData(y) <- cbind(colData(y),
                                                     cur_props[,!colnames(cur_props) %in%
                                                                   c(cell_id, coords)])
                             }

                             return(y)

                         }, BPPARAM = BPPARAM)
    return(cur_out)
}

#' @importFrom S4Vectors SelfHits
#' @importFrom SpatialExperiment spatialCoords<-
.steinbock_read_graphs <- function(x, cur_path, return_as, BPPARAM){

    cur_out <-  bplapply(x,
                         function(y){
                             cur_sample <- unique(y$sample_id)

                             cur_file <- list.files(cur_path,
                                                    pattern = paste0("^", cur_sample, ".csv", "$"),
                                                    full.names = TRUE)

                             if (length(cur_file) == 0) {
                                 return(y)
                             }

                             cur_graphs <- vroom(cur_file,
                                                 progress = FALSE,
                                                 col_types = cols()) %>%
                                 as.data.frame()

                             cur_hits <- SelfHits(from = match(cur_graphs[,1],
                                                                y$ObjectNumber),
                                                  to = match(cur_graphs[,2],
                                                                y$ObjectNumber),
                                                  nnode = ncol(y))

                             colPair(y, "neighborhood") <- cur_hits

                             return(y)
                         }, BPPARAM = BPPARAM)

    return(cur_out)
}

#' @importFrom methods as
.steinbock_add_image_metadata <- function(object, image_file,
                                          extract_imagemetadata_from) {

    cur_img_meta <- vroom(image_file,
                          col_select = all_of(c("image",
                                                extract_imagemetadata_from)),
                          show_col_types = FALSE)
    cur_img_meta$sample_id <- sub("\\.[^.]*$", "", cur_img_meta$image)

    colData(object) <- as(merge(x = colData(object), y = cur_img_meta[,-1],
                                by = "sample_id",
                                sort = FALSE), "DataFrame")

    return(object)
}

.add_panel <- function(x, path, panel, extract_names_from) {

    if (!is.null(panel)) {
        if (file.exists(file.path(path, panel))) {
            cur_panel <- vroom(file.path(path, panel),
                               progress = FALSE,
                               col_types = cols())
        } else if (file.exists(panel)) {
            cur_panel <- vroom(panel,
                               progress = FALSE,
                               col_types = cols())
        } else {
            warning("'panel_file' does not exist.")
            return(x)
        }

        cur_panel <- as.data.frame(cur_panel)

        cur_ind <- match(rownames(x), cur_panel[,extract_names_from])

        cur_panel <- cur_panel[cur_ind,]

        rowData(x) <- cur_panel
    }

    return(x)

}

#### read_cpout helpers ####

# Helper functions for reading in ImcSegmentationPipeline data
#' @importFrom vroom vroom
#' @importFrom dplyr contains mutate across
.cpout_create_object <- function(path, object_file, image_file,
                                 object_feature_file,
                                 intensities, extract_imgid_from,
                                 extract_cellid_from, extract_coords_from,
                                 extract_cellmetadata_from,
                                 extract_sampleid_from,
                                 scale_intensities, extract_scalingfactor_from,
                                 return_as) {

    cur_counts <- vroom(file.path(path, object_file),
                        col_select = c(contains(intensities),
                                       all_of(c(extract_imgid_from,
                                                extract_cellid_from,
                                                extract_coords_from,
                                                extract_cellmetadata_from))),
                        show_col_types = FALSE)

    scaling_factor <- as.data.frame(vroom(file.path(path, image_file),
                                          col_select = all_of(c(extract_imgid_from,
                                                                extract_scalingfactor_from)),
                                          show_col_types = FALSE))
    rownames(scaling_factor) <-
        as.character(scaling_factor[[extract_imgid_from]])

    # Scale counts
    if (scale_intensities) {
        cur_counts <- split(cur_counts, cur_counts[[extract_imgid_from]])
        scaling_factor <- scaling_factor[names(cur_counts),]

        cur_counts <- mapply(function(cells, s_factor) {
            cells %>% mutate(across(contains(intensities),
                                    function(x){x * s_factor}))
        }, cur_counts, as.list(scaling_factor[[extract_scalingfactor_from]]),
        SIMPLIFY = FALSE)

        cur_counts <- do.call(rbind, cur_counts)
    }

    # Order channels
    cur_channels <- colnames(cur_counts)[grepl(intensities,
                                               colnames(cur_counts))]
    cur_channels_id <- as.numeric(str_extract(cur_channels, "[0-9]{1,3}$"))

    cur_channels <- cur_channels[order(cur_channels_id, decreasing = FALSE)]

    if (return_as == "spe") {
        object <- SpatialExperiment(
            assays = list(counts = t(as.matrix(cur_counts[,cur_channels]))),
            sample_id = as.character(cur_counts[[extract_imgid_from]]))
        object$ObjectNumber <- cur_counts[[extract_cellid_from]]
    } else {
        object <- SingleCellExperiment(
            assays = list(counts = t(as.matrix(cur_counts[,cur_channels]))))
        object$sample_id <- as.character(cur_counts[[extract_imgid_from]])
        object$ObjectNumber <- cur_counts[[extract_cellid_from]]
    }

    # Build colData
    if (return_as == "spe") {
        if (!is.null(extract_coords_from)) {
            spatialCoords(object) <- matrix(c(
                cur_counts[[extract_coords_from[1]]],
                cur_counts[[extract_coords_from[2]]]),
                ncol = 2, byrow = FALSE,
                dimnames = list(as.character(object$ObjectNumber),
                                c("Pos_X", "Pos_Y")))
        }

        colData(object) <- cbind(colData(object),
                                 cur_counts[,extract_cellmetadata_from])
    } else {
        if (!is.null(extract_coords_from)) {
            object$Pos_X <- cur_counts[[extract_coords_from[1]]]
            object$Pos_Y <- cur_counts[[extract_coords_from[2]]]
        }

        colData(object) <- cbind(colData(object),
                                 cur_counts[,extract_cellmetadata_from])
    }

    # Set correct metal names
    cur_channels <- vroom(file.path(path, object_feature_file),
                          col_select = c("channel", "channel_id"),
                          show_col_types = FALSE)
    cur_channels <- unique(cur_channels)
    cur_channels <- cur_channels[grepl("[a-zA-Z]{1,2}[0-9]{2,3}",
                                       cur_channels[["channel_id"]]),]
    cur_channels <- cur_channels[order(cur_channels[["channel"]],
                                       decreasing = FALSE),]

    rownames(object) <- cur_channels[["channel_id"]]

    return(object)
}

#' @importFrom methods as
.cpout_add_image_metadata <- function(object, path, image_file,
                                      extract_imgid_from,
                                      extract_imagemetadata_from) {

    cur_img_meta <- vroom(file.path(path, image_file),
                          col_select = all_of(c(extract_imgid_from,
                                                extract_imagemetadata_from)),
                          show_col_types = FALSE)
    colData(object) <- as(merge(x = colData(object), y = cur_img_meta,
                                by.x = "sample_id", by.y = extract_imgid_from,
                                sort = FALSE), "DataFrame")

    return(object)
}

.cpout_add_graph <- function(object, path, graph_file,
                             extract_graphimageid_from,
                             extract_graphcellids_from) {
    cur_graph <- vroom(file.path(path, graph_file),
                       col_select = all_of(c(extract_graphimageid_from,
                                             extract_graphcellids_from)),
                       show_col_types = FALSE)

    cur_graph$firstid <- paste0(cur_graph[[extract_graphimageid_from]], "_",
                                cur_graph[[extract_graphcellids_from[1]]])
    cur_graph$secondid <- paste0(cur_graph[[extract_graphimageid_from]], "_",
                                 cur_graph[[extract_graphcellids_from[2]]])

    cur_objectids <- paste0(object$sample_id, "_", object$ObjectNumber)

    colPair(object, "neighborhood") <- SelfHits(from = match(cur_graph$firstid,
                                                             cur_objectids),
                                                to = match(cur_graph$secondid,
                                                           cur_objectids),
                                                nnode = length(cur_objectids))

    return(object)
}



############################# Spatial helpers ##################################

#### plotSpatial helpers ####
.makeNodes <- function(object, node_color_by, img_id, node_shape_by,
                       node_size_by, assay_type){
    if (is.null(node_color_by)) {
        nodes <- colData(object)[,c(img_id, node_shape_by,
                                    node_size_by),
                                 drop=FALSE]
    } else if (node_color_by %in% names(colData(object))) {
        nodes <- colData(object)[,c(img_id, node_color_by, node_shape_by,
                                    node_size_by),
                                 drop=FALSE]
    } else {
        nodes <- colData(object)[,c(img_id, node_shape_by,
                                    node_size_by),
                                 drop=FALSE]
        nodes <- cbind(nodes, t(assay(object, assay_type)[node_color_by,,
                                                          drop = FALSE]))
    }

    if (!is.null(node_shape_by)) {
        nodes[,node_shape_by] <- as.character(nodes[,node_shape_by])
    }
    
    nodes <- nodes[,unique(colnames(nodes)), drop = FALSE]

    return(nodes)
}

# Function to generate the tidygraph
#' @importFrom S4Vectors isRedundantHit
.generateGraph <- function(object, nodes, colPairName, draw_edges,
                           edge_color_by, edge_width_by, directed){
    if (draw_edges) {

        if (!directed) {
            cur_SH <- colPair(object, colPairName)
            cur_SH <- cur_SH[!isRedundantHit(cur_SH)]
        } else {
            cur_SH <- colPair(object, colPairName)
        }

        edges <- as.data.frame(as(cur_SH, "DataFrame"))

        if (!is.null(edge_color_by) &&
            edge_color_by %in% colnames(colData(object))) {
            edges[,edge_color_by] <-
                colData(object)[[edge_color_by]][edges$from]
        }

        if (!is.null(edge_width_by) &&
            edge_width_by %in% colnames(colData(object))) {
            edges[,edge_width_by] <-
                colData(object)[[edge_width_by]][edges$from]
        }

        cur_graph <- tbl_graph(nodes = nodes,
                               edges = edges,
                               directed = directed)
    } else {
        cur_graph <- tbl_graph(nodes = nodes,
                               directed = directed)
    }

    return(cur_graph)

}

# Function to generate the base plot
.generatePlot <- function(layout, draw_edges, directed, arrow, end_cap,
                          node_color_by, node_size_by, node_shape_by,
                          node_color_fix, node_size_fix, node_shape_fix,
                          edge_color_by, edge_width_by, edge_color_fix,
                          edge_width_fix, nodes_first){

    node_color_by <- if(is.null(node_color_by)) NULL else as.name(node_color_by)
    node_size_by <- if(is.null(node_size_by)) NULL else as.name(node_size_by)
    node_shape_by <- if(is.null(node_shape_by)) NULL else as.name(node_shape_by)
    edge_color_by <- if(is.null(edge_color_by)) NULL else as.name(edge_color_by)
    edge_width_by <- if(is.null(edge_width_by)) NULL else as.name(edge_width_by)

    if (!is.null(node_color_fix)){ node_color_by <- as.character(node_color_fix)
    } else { node_color_by <- node_color_by }
    if (!is.null(node_size_fix)){ node_size_by <- as.character(node_size_fix)
    } else { node_size_by <- node_size_by }
    if (!is.null(node_shape_fix)){ node_shape_by <- as.character(node_shape_fix)
    } else { node_shape_by <- node_shape_by }
    if (!is.null(edge_color_fix)){ edge_color_by <- as.character(edge_color_fix)
    } else { edge_color_by <-  edge_color_by }
    if (!is.null(edge_width_fix)){edge_width_by <- as.character(edge_width_fix)
    } else { edge_width_by <- edge_width_by}

    if (draw_edges) {
        if (!is.null(arrow)) {

            if (is.null(end_cap)) {
                end_cap <- circle(0.1, 'cm')
            }

            if (directed) {
                cur_geom_edge <- geom_edge_fan(aes_(edge_colour = edge_color_by,
                                                    edge_width = edge_width_by),
                                               end_cap = end_cap,
                                               arrow = arrow)
            } else {
                cur_geom_edge <- geom_edge_link(aes_(
                    edge_colour = edge_color_by,
                    edge_width = edge_width_by),
                    end_cap = end_cap,
                    arrow = arrow)
            }
        } else {
            if (directed) {
                cur_geom_edge <- geom_edge_fan0(aes_(
                    edge_colour = edge_color_by,
                    edge_width = edge_width_by))
            } else {
                cur_geom_edge <- geom_edge_link0(aes_(
                    edge_colour = edge_color_by,
                    edge_width = edge_width_by))
            }
        }

        if (nodes_first) {
            p <- ggraph(layout) +
                geom_node_point(aes_(colour = node_color_by,
                                     size = node_size_by,
                                     shape = node_shape_by)) +
                cur_geom_edge
        } else {
            p <- ggraph(layout) +
                cur_geom_edge +
                geom_node_point(aes_(colour = node_color_by,
                                     size = node_size_by,
                                     shape = node_shape_by))
        }
    } else {
        p <- ggraph(layout) +
            geom_node_point(aes_(colour = node_color_by,
                                 size = node_size_by,
                                 shape = node_shape_by))
    }

    return(p)
}

# Post process the plots
#' @importFrom ggplot2 ggtitle scale_x_reverse scale_y_reverse
#' @importFrom viridis scale_color_viridis
.postProcessPlot <- function(p, object, img_id, nrows, ncols, node_color_by,
                             node_color_fix,
                             node_shape_fix, node_size_fix, edge_color_fix,
                             edge_width_fix, scales, flip_x, flip_y){

    if (!is.null(node_color_fix)) {
        names(node_color_fix) <- as.character(node_color_fix)
        p <- p + scale_color_manual(values = node_color_fix,
                                    guide = "none")
    }
    if (!is.null(node_shape_fix)) {
        names(node_shape_fix) <- as.character(node_shape_fix)
        p <- p + scale_shape_manual(values = node_shape_fix,
                                    guide = "none")
    }
    if (!is.null(node_size_fix)) {
        names(node_size_fix) <- as.character(node_size_fix)
        p <- p + scale_size_manual(values = node_size_fix,
                                   guide = "none")
    }
    if (!is.null(edge_color_fix)) {
        names(edge_color_fix) <- as.character(edge_color_fix)
        p <- p + scale_edge_color_manual(values = edge_color_fix,
                                         guide = "none")
    }
    if (!is.null(edge_width_fix)) {
        names(edge_width_fix) <- as.character(edge_width_fix)
        p <- p + scale_edge_width_manual(values = edge_width_fix,
                                         guide = "none")
    }
    if (!is.null(node_color_by) && node_color_by %in% rownames(object)) {
        p <- p + scale_color_viridis()
    }

    if (length(unique(colData(object)[[img_id]])) > 1) {
        p <- p + facet_nodes(img_id, scales = scales,
                             nrow = nrows, ncol = ncols) +
            theme(axis.text = element_text(),
                  panel.background = element_blank())
    } else {
        p <- p + theme(axis.text = element_text(),
                       panel.background = element_blank()) +
            ggtitle(unique(colData(object))[[img_id]])
    }

    if (flip_x) {
        p <- p + scale_x_reverse()
    }

    if (flip_y) {
        p <- p + scale_y_reverse()
    }

    return(p)
}

#### patchDetection helpers ####
# Helper function for the patch detection method
#' @importFrom sf st_multipoint st_cast st_sfc st_distance
#' @importFrom dplyr as_tibble filter sym nest_by summarize
#' @importFrom S4Vectors metadata
.expand_patch <- function(object,
                          name,
                          expand_by,
                          coords,
                          convex,
                          img_id,
                          BPPARAM){

    cur_meta <- metadata(object)
    metadata(object) <- list()
    
    cur_intmeta <- int_metadata(object)

    cur_out <- bplapply(
        unique(colData(object)[[img_id]]),
        function(x){

            cur_obj <- object[,as.character(colData(object)[[img_id]]) == x]

            if (is(cur_obj, "SpatialExperiment")) {
                cur_coords <- spatialCoords(cur_obj)[,coords]
            } else {
                cur_coords <- colData(cur_obj)[,coords]
            }

            cells <- st_multipoint(as.matrix(cur_coords))
            cells_sfc <- st_cast(st_sfc(cells), "POINT")

            if (sum(!is.na(colData(cur_obj)[[name]])) == 0) {
                return(cur_obj)
            }

            data <- polygon <- NULL

            if (is(cur_obj, "SpatialExperiment")) {
                cur_out <- cbind(colData(cur_obj), cur_coords) %>% as_tibble %>%
                    filter(!is.na(!!sym(name))) %>%
                    nest_by(!!sym(name)) %>%
                    summarize(
                        polygon = list(.polygon_function(x = data,
                                                     coords = coords,
                                                     convex = convex)),
                        cells = list(.milieu_function(x = polygon,
                                                  distance = expand_by,
                                                  cells = cells_sfc)))
            } else {
                cur_out <- colData(cur_obj) %>% as_tibble %>%
                    filter(!is.na(!!sym(name))) %>%
                    nest_by(!!sym(name)) %>%
                    summarize(
                        polygon = list(.polygon_function(x = data,
                                                    coords = coords,
                                                    convex = convex)),
                        cells = list(.milieu_function(x = polygon,
                                                    distance = expand_by,
                                                    cells = cells_sfc)))
            }

            # Find cells that are not unique in extended patches
            cur_cells <- do.call(c, cur_out$cells)
            cur_cells <- cur_cells[duplicated(cur_cells)]
            cur_cells <- cur_cells[!is.na(cur_cells)]

            if (length(cur_cells) > 0) {
                cur_dists <- mapply(function(y, patch_name){
                    if (is.na(y)) {return(NULL)}
                    cur_mat <- st_distance(cells_sfc[cur_cells], y)
                    colnames(cur_mat) <- patch_name
                    return(cur_mat)
                }, cur_out$polygon, cur_out[[name]], SIMPLIFY = FALSE)
                cur_dists <- do.call("cbind", cur_dists)
                cur_patch_id <- apply(cur_dists, 1, function(y){
                    return(colnames(cur_dists)[which.min(y)])
                })
            }

            cur_patch <- colData(cur_obj)[[name]]
            for (i in seq_len(nrow(cur_out))) {
                if (all(!is.na(cur_out$cells[[i]]))) {
                    cur_patch[cur_out$cells[[i]]] <- cur_out[[name]][i]
                }
            }

            if (length(cur_cells) > 0) {
                cur_patch[cur_cells] <- cur_patch_id
            }

            cur_obj[[name]] <- cur_patch

            return(cur_obj)

        }, BPPARAM = BPPARAM)

    cur_out <- do.call("cbind", cur_out)
    metadata(cur_out) <- cur_meta
    int_metadata(cur_out) <- cur_intmeta
    
    return(cur_out)


}
#' @importFrom concaveman concaveman
#' @importFrom grDevices chull
#' @importFrom sf st_polygon
.polygon_function <- function(x, coords, convex){
    if (nrow(x) <= 2) {
        return(NA)
    }

    if (convex) {
        hull <- chull(x = x[[coords[1]]], y = x[[coords[2]]])

        # cells that build the border of a patch
        border_cells = x[hull,]
        coordinates = as.matrix(border_cells[,coords])
        coordinates <- rbind(coordinates, coordinates[1,])

        polygon <- st_polygon(list(coordinates))

        return(polygon)
    } else {
        cur_coords <- as.matrix(cbind(x[[coords[1]]], x[[coords[2]]]))
        hull <- data.frame(concaveman(cur_coords, concavity = 1))

        polygon <- st_polygon(list(as.matrix(hull)))

        return(polygon)
    }
}

#' @importFrom sf st_buffer st_sfc st_intersects
.milieu_function <- function(x, distance, cells){

    if (is.na(x)) {
        return(NA)
    }

    polygon_buff <- st_buffer(x, distance)
    polygon_buff_sfc <- st_sfc(polygon_buff)

    intersect_cells <- st_intersects(polygon_buff_sfc, cells)

    return(intersect_cells[[1]])
}

#### testInteractions and countInteractions helpers ####
# Helper functions for the neighbourhood permutation test
#' @importFrom data.table as.data.table
.prepare_table <- function(object, group_by, cur_label, colPairName) {
    cur_tab <- as.data.table(colPair(object, colPairName))
    cur_tab$group_by <- colData(object)[[group_by]][cur_tab$from]
    cur_tab$from_label <- cur_label[cur_tab$from]
    cur_tab$to_label <- cur_label[cur_tab$to]
    cur_tab$ct <- 1

    . <- .N <- NULL

    cur_tab_2 <- data.table(group_by = colData(object)[[group_by]],
                            from_label = cur_label)
    cur_tab_2 <- cur_tab_2[,.(total = .N), by = c("group_by", "from_label")]
    cur_tab <- merge(cur_tab, cur_tab_2, by = c("group_by", "from_label"),
                     sort = FALSE)

    return(cur_tab)
}

#' @importFrom data.table CJ
.aggregate_histo <- function(dat_table, object, group_by, label,
                             check_missing = TRUE) {
    . <- ct <- .N <- NULL
    dat_temp <- dat_table[, .(ct=.N), by = c("group_by", "from_label",
                                             "to_label", "from")]
    dat_temp <- dat_temp[, .(ct=mean(ct)), by = c("group_by", "from_label",
                                                  "to_label")]

    if (check_missing) {
        dat_temp <- dat_temp[CJ(group_by = unique(dat_temp$group_by),
                                from_label = as.factor(levels(dat_temp$from_label)),
                                to_label = as.factor(levels(dat_temp$to_label))),
                             on = c("group_by", "from_label", "to_label")]
        ct <- from_label <- to_label <- NULL
        dat_temp[is.na(dat_temp$ct), ct := 0]

        # Set all cells that are not contained in specific groups to NA
        cur_dat <- unclass(table(colData(object)[[group_by]],
                                 colData(object)[[label]]))
        cur_ind <- which(cur_dat == 0, arr.ind = TRUE)

        if (nrow(cur_ind) > 0) {
            apply(cur_ind, 1 , function(x){
                dat_temp[group_by == rownames(cur_dat)[x[1]] &
                             (from_label == colnames(cur_dat)[x[2]] |
                                  to_label == colnames(cur_dat)[x[2]]),
                         ct := NA]
            })
        }

    }
    return(dat_temp)
}

#' @importFrom data.table dcast.data.table melt.data.table
.aggregate_classic <- function(dat_table, object, group_by, label,
                               check_missing = TRUE){
    dat_temp <- dcast.data.table(dat_table,
                                 "group_by + from_label + total + from ~ to_label",
                                 value.var = "ct", fun.aggregate = sum,
                                 fill = 0)
    dat_temp <- melt.data.table(dat_temp, id.vars = c("group_by", "from_label",
                                                      "from", "total"),
                                variable.name = "to_label",
                                value.name = "ct")

    total <- NULL

    dat_temp[,ct := ct/total]

    dat_temp <- dcast.data.table(dat_temp, "group_by + from_label ~ to_label",
                                 value.var = "ct",
                                 fun.aggregate = sum,
                                 fill = 0, drop = FALSE)

    dat_temp <- melt.data.table(dat_temp, id.vars = c("group_by", "from_label"),
                                variable.name = "to_label",
                                value.name = "ct")

    if (check_missing) {
        dat_temp <- dat_temp[CJ(group_by = unique(dat_table$group_by),
                                from_label = as.factor(levels(dat_table$from_label)),
                                to_label = as.factor(levels(dat_table$to_label))),
                             on = c("group_by", "from_label", "to_label")]
        ct <- from_label <- to_label <- NULL
        dat_temp[is.na(dat_temp$ct), ct := 0]

        # Set all cells that are not contained in specific groups to NA
        cur_dat <- unclass(table(colData(object)[[group_by]],
                                 colData(object)[[label]]))
        cur_ind <- which(cur_dat == 0, arr.ind = TRUE)

        if (nrow(cur_ind) > 0) {
            apply(cur_ind, 1 , function(x){
                dat_temp[group_by == rownames(cur_dat)[x[1]] &
                             (from_label == colnames(cur_dat)[x[2]] |
                                  to_label == colnames(cur_dat)[x[2]]),
                         ct := NA]
            })
        }
    }

    return(dat_temp)
}

.aggregate_classic_patch <- function(dat_table, patch_size, object, group_by,
                                     label, check_missing = TRUE){
    dat_temp <- dcast.data.table(dat_table,
                                 "group_by + from_label + total + from ~ to_label",
                                 value.var = "ct", fun.aggregate = sum,
                                 fill = 0)
    dat_temp <- melt.data.table(dat_temp, id.vars = c("group_by", "from_label",
                                                      "from", "total"),
                                variable.name = "to_label",
                                value.name = "ct")

    total <- NULL

    dat_temp[, ct := patch_size <= ct ]
    dat_temp[, ct := ct/total]

    dat_temp <- dcast.data.table(dat_temp, "group_by + from_label ~ to_label",
                                 value.var = "ct",
                                 fun.aggregate = sum, fill = 0, drop = FALSE)
    dat_temp <- melt.data.table(dat_temp, id.vars = c("group_by", "from_label"),
                                variable.name = "to_label",
                                value.name = "ct")

    if (check_missing) {
        dat_temp <- dat_temp[CJ(group_by = unique(dat_table$group_by),
                                from_label = as.factor(levels(dat_table$from_label)),
                                to_label = as.factor(levels(dat_table$to_label))),
                             on = c("group_by", "from_label", "to_label")]
        ct <- from_label <- to_label <- NULL
        dat_temp[is.na(dat_temp$ct), ct := 0]

        # Set all cells that are not contained in specific groups to NA
        cur_dat <- unclass(table(colData(object)[[group_by]],
                                 colData(object)[[label]]))
        cur_ind <- which(cur_dat == 0, arr.ind = TRUE)

        if (nrow(cur_ind) > 0) {
            apply(cur_ind, 1 , function(x){
                dat_temp[group_by == rownames(cur_dat)[x[1]] &
                             (from_label == colnames(cur_dat)[x[2]] |
                                  to_label == colnames(cur_dat)[x[2]]),
                         ct := NA]
            })
        }
    }

    return(dat_temp)
}

#' @importFrom data.table data.table
.permute_labels <- function(object, group_by, label, iter, patch_size,
                            colPairName, method, BPPARAM) {

    cur_lab_table <- data.table(label = as.factor(colData(object)[[label]]),
                                group_by = colData(object)[[group_by]])

    . <- label <- NULL

    cur_out <- bplapply(seq_len(iter),
                        function(x){

                            label_perm <- cur_lab_table[ ,
                                                         .(label=sample(label)), by=group_by]
                            cur_perm <- .prepare_table(object, group_by,
                                                       label_perm$label, colPairName)

                            if (method == "classic") {
                                cur_perm <- .aggregate_classic(cur_perm, object,
                                                               group_by, label,
                                                               check_missing = FALSE)
                            } else if (method == "histocat") {
                                cur_perm <- .aggregate_histo(cur_perm, object,
                                                             group_by, label,
                                                             check_missing = FALSE)
                            } else if (method == "patch") {
                                cur_perm <- .aggregate_classic_patch(cur_perm,
                                                                     patch_size = patch_size,
                                                                     object, group_by, label,
                                                                     check_missing = FALSE)
                            }
                            cur_perm$iter <- x

                            return(cur_perm)

                        }, BPPARAM = BPPARAM)

    cur_out <- do.call("rbind", cur_out)

    return(cur_out)
}

.calc_p_vals<- function(dat_baseline, dat_perm, n_perm, p_thres){
    dat_perm <- merge(dat_perm,
                      dat_baseline[, c("from_label", "to_label",
                                       "group_by", "ct")],
                      by = c("from_label", "to_label", "group_by"),
                      suffixes = c("_perm", "_obs"), all = TRUE)

    . <- ct_perm <- ct_obs <- p_gt <- p_lt <- NULL
    direction <- sig <- sigval <- p <-  NULL

    dat_perm[, ':='(ct_perm = replace(ct_perm, is.na(ct_perm), 0),
                    ct_obs = replace(ct_obs, is.na(ct_obs), 0))]

    dat_stat <- dat_perm[ , .(ct = mean(ct_obs),
                              p_gt = ifelse(max(ct_obs) == 0, 1,
                                            (sum(ct_perm >= ct_obs) + 1) / (n_perm + 1)),
                              p_lt = (n_perm - sum(ct_perm > ct_obs) + 1) / (n_perm + 1)),
                          by=c("group_by", "from_label", "to_label")]

    dat_stat[, interaction := p_gt < p_lt]
    dat_stat[, p := p_gt * interaction + p_lt * (!interaction)]
    dat_stat[, sig := p < p_thres]
    dat_stat[, sigval := as.integer(sig) * sign(interaction - 0.5)]

    setorder(dat_stat, "group_by", "from_label", "to_label")
    setorder(dat_baseline, "group_by", "from_label", "to_label")

    dat_stat[is.na(dat_baseline$ct),
             c("p_gt", "p_lt", "ct", "interaction",
               "p", "sig", "sigval") := NA]

    return(dat_stat)
}



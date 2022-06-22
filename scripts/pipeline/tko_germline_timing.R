
# tko germline timing


tko_germline_timing <- function(mat_nm) {
    cgraph_id <- mat_nm
    mat <- scdb_mat(mat_nm)
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    # exclude ExE endoderm and Visceral endoderm
    excluded_colors <- c("#F6BFCB", "#7F6874")

    load(sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))

    f <- mat@cell_metadata[colnames(mat@mat), "cell_type"] %in% c("TKO")
    tetra_cls <- colnames(mat@mat)[f]
    tetra_cls <- intersect(tetra_cls, names(cmp_annot$query_cls_col))

    excluded_cls <- tetra_cls[cmp_annot$query_cls_col[tetra_cls] %in% excluded_colors]
    tetra_cls <- setdiff(tetra_cls, excluded_cls)

    tetra_embryos <- unique(mat@cell_metadata[tetra_cls, "embryo"])

    wt10_cls <- intersect(names(mc_wt@mc)[!(mc_wt@colors[mc_wt@mc] %in% excluded_colors)], colnames(mat@mat))
    atlas_time <- mat@cell_metadata[wt10_cls, "transcriptional_rank"]
    names(atlas_time) <- wt10_cls

    data_dir <- sprintf("data/%s", mat_nm)
    if (!dir.exists(data_dir)) {
        dir.create(data_dir)
    }

    data_dir <- sprintf("data/%s/time_match", mat_nm)
    if (!dir.exists(data_dir)) {
        dir.create(data_dir)
    }



    # first get atlas time distribution

    atlas_time_dist <- get_atlas_time_dist(atlas_time = atlas_time, graph_id = cgraph_id)

    # first timing using control cells only
    query_cls <- tetra_cls[(mat@cell_metadata[tetra_cls, "cell_type"] %in% c("TKO")) & (mat@cell_metadata[tetra_cls, "embryo"] %in% tetra_embryos)]
    query_cls_md <- mat@cell_metadata[query_cls, "embryo"]
    names(query_cls_md) <- query_cls

    query_time_dist <- get_query_time_dist(query_cls_md = query_cls_md, atlas_time = atlas_time, graph_id = cgraph_id)

    time_dist_list <- list(
        atlas_time_dist = atlas_time_dist$atlas_time_dist,
        query = query_time_dist$query_time_dist
    )

    save(time_dist_list, file = sprintf("%s/time_dist_control_ko.Rda", data_dir))

    query_emb_summary <- as.data.frame.matrix(table(mat@cell_metadata[tetra_cls, "embryo"], mat@cell_metadata[tetra_cls, "cell_type"]))
    query_emb_summary$embryo <- rownames(query_emb_summary)
    query_emb_summary$best_query <- NA
    query_emb_summary[rownames(time_dist_list$query), "best_query"] <- time_dist_best_match(
        atlas_time_dist = time_dist_list$atlas_time_dist,
        query_time_dist = time_dist_list$query
    )

    write.table(query_emb_summary, file = sprintf("%s/time_match_summary.txt", data_dir), sep = "\t", row.names = F)
}

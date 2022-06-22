# tko tetraploid timing


tetra_timing <- function(mat_nm) {
    cgraph_id <- mat_nm
    mat <- scdb_mat(mat_nm)
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    excluded_colors <- c("#F6BFCB", "#7F6874")

    load(sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))

    f <- mat@cell_metadata[colnames(mat@mat), "cell_type"] %in% c("KO", "control")
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
    query_cls <- tetra_cls[(mat@cell_metadata[tetra_cls, "cell_type"] %in% c("control", "KO")) & (mat@cell_metadata[tetra_cls, "embryo"] %in% tetra_embryos)]
    query_cls_md <- mat@cell_metadata[query_cls, "embryo"]
    names(query_cls_md) <- query_cls

    query_time_dist <- get_query_time_dist(query_cls_md = query_cls_md, atlas_time = atlas_time, graph_id = cgraph_id)

    time_dist_list <- list(
        atlas_time_dist = atlas_time_dist$atlas_time_dist,
        query = query_time_dist$query_time_dist
    )

    save(time_dist_list, file = sprintf("%s/time_dist_control_ko.Rda", data_dir))

    chim_emb_summary <- as.data.frame.matrix(table(mat@cell_metadata[tetra_cls, "embryo"], mat@cell_metadata[tetra_cls, "cell_type"]))
    chim_emb_summary$embryo <- rownames(chim_emb_summary)
    chim_emb_summary$best_query <- NA
    chim_emb_summary[rownames(time_dist_list$query), "best_query"] <- time_dist_best_match(
        atlas_time_dist = time_dist_list$atlas_time_dist,
        query_time_dist = time_dist_list$query
    )

    write.table(chim_emb_summary, file = sprintf("%s/time_match_summary.txt", data_dir), sep = "\t", row.names = F)
}


tetraploid_plot_cor_time_dist_per_emb <- function(mat_nm) {
    rank_to_time <- read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt", stringsAsFactors = F, h = T, sep = "\t")
    dev_time <- rank_to_time$developmental_time
    data_dir <- sprintf("data/%s/time_match", mat_nm)
    load(sprintf("%s/time_dist_control_ko.Rda", data_dir))

    tetra_time_summary <- read.table(file = sprintf("%s/time_match_summary.txt", data_dir), sep = "\t", h = T, stringsAsFactors = F)

    ctrl_embryos <- tetra_time_summary$embryo[is.na(tetra_time_summary$best_rank_ko)]
    ko_embryos <- tetra_time_summary$embryo[!is.na(tetra_time_summary$best_rank_ko)]

    query_dist_all <- time_dist_list$query
    atlas_time_dist <- time_dist_list$atlas_time_dist

    query_ref_cor <- tgs_cor(
        t(as.matrix(as.data.frame.matrix(query_dist_all))),
        t(as.matrix(as.data.frame.matrix(atlas_time_dist)))
    )

    tetra_embryos <- tetra_time_summary$embryo

    fig_dir <- sprintf("%s/%s/time_match", .scfigs_base, mat_nm)
    if (!dir.exists(sprintf("%s/%s", .scfigs_base, mat_nm))) {
        dir.create(sprintf("%s/%s", .scfigs_base, mat_nm))
    }
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    for (tetraploid in tetra_embryos) {
        x <- dev_time
        # y = rollmean(query_ref_cor[tetraploid,],k = 17,fill = c(query_ref_cor[tetraploid,1],0,query_ref_cor[tetraploid,ncol(query_ref_cor)]))

        smooth_spline <- smooth.spline(x = dev_time, y = query_ref_cor[tetraploid, ], spar = 0.8)
        y <- smooth_spline$y
        png(sprintf("%s/time_dist_cor_%s.png", fig_dir, tetraploid))
        plot(
            x = dev_time, y = query_ref_cor[tetraploid, ], pch = 19, main = tetraploid, ylim = c(min(query_ref_cor), 1), xlim = c(min(dev_time), max(dev_time)),
            xlab = "developmental time", cex = 0.7, ylab = "Correlation time distributions"
        )
        lines(x = x, y = y)
        abline(v = x[which.max(y)], lty = "dashed")
        dev.off()
    }
}

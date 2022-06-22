
dko_chim_timing <- function(mat_nm) {
    n_cls_min <- 20

    cgraph_id <- mat_nm
    mat <- scdb_mat(mat_nm)
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    # excluding Erythroid 1,2 and visceral and ExE endoderm
    excluded_colors <- mc_wt@color_key$color[c(22, 23, 28, 29)]

    load(sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))

    f <- mat@cell_metadata[colnames(mat@mat), "cell_type"] %in% c("DKO12", "DKO13", "DKO23", "host")
    chim_cls <- colnames(mat@mat)[f]
    chim_cls <- intersect(chim_cls, names(cmp_annot$query_cls_col))

    excluded_cls <- chim_cls[cmp_annot$query_cls_col[chim_cls] %in% excluded_colors]
    chim_cls <- setdiff(chim_cls, excluded_cls)

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


    #
    cell_type_list <- list(ko = c("DKO12", "DKO13", "DKO23"), host = c("host"))

    time_dist_host_ko <- list(atlas_time_dist = atlas_time_dist$atlas_time_dist)

    for (tag in names(cell_type_list)) {
        cell_type <- cell_type_list[[tag]]
        query_cls <- chim_cls[mat@cell_metadata[chim_cls, "cell_type"] %in% cell_type]
        n_cls_per_chimera <- table(mat@cell_metadata[query_cls, "embryo"])
        chimeras_f <- names(n_cls_per_chimera)[n_cls_per_chimera >= n_cls_min]
        query_cls <- query_cls[mat@cell_metadata[query_cls, "embryo"] %in% chimeras_f]
        query_cls_md <- mat@cell_metadata[query_cls, "embryo"]

        names(query_cls_md) <- query_cls

        query_time_dist <- get_query_time_dist(query_cls_md = query_cls_md, atlas_time = atlas_time, graph_id = cgraph_id)

        time_dist_host_ko[[tag]] <- query_time_dist$query_time_dist
    }

    save(time_dist_host_ko, file = sprintf("%s/time_dist_host_ko.Rda", data_dir))

    chim_emb_summary <- as.data.frame.matrix(table(mat@cell_metadata[chim_cls, "embryo"], mat@cell_metadata[chim_cls, "cell_type"]))
    chim_emb_summary$embryo <- rownames(chim_emb_summary)
    chim_emb_summary$best_rank_host <- NA
    chim_emb_summary[rownames(time_dist_host_ko$host), "best_rank_host"] <- time_dist_best_match(
        atlas_time_dist = time_dist_host_ko$atlas_time_dist,
        query_time_dist = time_dist_host_ko$host
    )
    chim_emb_summary$best_rank_ko <- NA
    chim_emb_summary[rownames(time_dist_host_ko$ko), "best_rank_ko"] <- time_dist_best_match(
        atlas_time_dist = time_dist_host_ko$atlas_time_dist,
        query_time_dist = time_dist_host_ko$ko
    )

    write.table(chim_emb_summary, file = sprintf("%s/time_match_summary.txt", data_dir), sep = "\t", row.names = F)
}


dko_chimera_plot_cor_time_dist_per_emb <- function(mat_nm) {
    n_cls_min <- 19
    rank_to_time <- read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt", stringsAsFactors = F, h = T, sep = "\t")
    dev_time <- rank_to_time$developmental_time
    data_dir <- sprintf("data/%s/time_match", mat_nm)
    load(sprintf("%s/time_dist_host_ko.Rda", data_dir))

    host_dist_all <- time_dist_host_ko$host
    ko_dist_all <- time_dist_host_ko$ko
    atlas_time_dist <- time_dist_host_ko$atlas_time_dist



    host_ref_cor <- tgs_cor(
        t(as.matrix(as.data.frame.matrix(host_dist_all))),
        t(as.matrix(as.data.frame.matrix(atlas_time_dist)))
    )
    ko_ref_cor <- tgs_cor(
        t(as.matrix(as.data.frame.matrix(ko_dist_all))),
        t(as.matrix(as.data.frame.matrix(atlas_time_dist)))
    )

    chim_time_summary <- read.table(file = sprintf("%s/time_match_summary.txt", data_dir), sep = "\t", h = T, stringsAsFactors = F)

    f <- (chim_time_summary$host > n_cls_min)

    chim_embryos <- chim_time_summary$embryo[f]

    if (!dir.exists(sprintf("figs/%s", mat_nm))) {
        dir.create(sprintf("figs/%s", mat_nm))
    }

    fig_dir <- sprintf("figs/%s/time_match", mat_nm)
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    for (chimera in chim_embryos) {
        smooth_spline_host <- smooth.spline(x = dev_time, y = host_ref_cor[chimera, ], spar = 0.8)
        smooth_spline_ko <- smooth.spline(x = dev_time, y = ko_ref_cor[chimera, ], spar = 0.8)

        png(sprintf("%s/time_dist_cor_%s.png", fig_dir, chimera))
        plot(
            x = dev_time, y = host_ref_cor[chimera, ], pch = 19, main = chimera, ylim = c(min(host_ref_cor, ko_ref_cor), 1),
            xlab = "developmental time", cex = 0.7, ylab = "Correlation time distributions"
        )
        points(x = dev_time, y = ko_ref_cor[chimera, ], pch = 19, col = "coral2", cex = 0.7)
        lines(x = dev_time, y = smooth_spline_host$y)
        abline(v = dev_time[which.max(smooth_spline_host$y)], lty = "dashed")
        lines(x = dev_time, y = smooth_spline_ko$y, col = "coral2")
        abline(v = dev_time[which.max(smooth_spline_ko$y)], lty = "dashed", col = "coral2")
        legend(x = "topleft", legend = c("KO", "host"), pch = 19, col = c("coral2", "black"))
        dev.off()
    }
}

plot_host_vs_ko_rank <- function(mat_nm) {
    n_cls_min <- 19
    rank_to_time <- read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt", stringsAsFactors = F, h = T, sep = "\t")
    dev_time <- rank_to_time$developmental_time
    data_dir <- sprintf("data/%s/time_match", mat_nm)

    chim_emb_summary <- read.table(file = sprintf("%s/time_match_summary.txt", data_dir), sep = "\t", h = T, stringsAsFactors = F)

    f <- (chim_emb_summary$host > n_cls_min)

    rank_min <- min(c(chim_emb_summary$best_rank_ko[f], chim_emb_summary$best_rank_host[f])) - 1
    rank_max <- max(c(chim_emb_summary$best_rank_ko[f], chim_emb_summary$best_rank_host[f])) + 1
    time_min <- min(c(dev_time[chim_emb_summary$best_rank_ko[f]], dev_time[chim_emb_summary$best_rank_host[f]])) - 0.1
    time_max <- max(c(dev_time[chim_emb_summary$best_rank_ko[f]], dev_time[chim_emb_summary$best_rank_host[f]])) + 0.1

    col_emb <- rep("black", length(chim_emb_summary$embryo))


    fig_dir <- sprintf("%s/%s/time_match", .scfigs_base, mat_nm)




    png(sprintf("%s/best_rank_ko_vs_host.png", fig_dir))
    plot(chim_emb_summary$best_rank_ko[f], chim_emb_summary$best_rank_host[f],
        pch = 19, xlim = c(rank_min, rank_max), ylim = c(rank_min, rank_max), main = "KO vs host",
        xlab = "Rank KO cells", ylab = "Rank host cells", col = col_emb, cex = 2, cex.lab = 2
    )
    abline(a = 0, b = 1, lty = "dashed")
    dev.off()

    png(sprintf("%s/best_time_ko_vs_host.png", fig_dir))
    plot(rank_to_time$developmental_time[chim_emb_summary$best_rank_ko[f]],
        rank_to_time$developmental_time[chim_emb_summary$best_rank_host[f]],
        pch = 19,
        xlim = c(time_min, time_max), ylim = c(time_min, time_max), main = "KO vs host",
        xlab = "Time KO cells", ylab = "Time host cells", col = col_emb, cex = 2, cex.lab = 2
    )
    abline(a = 0, b = 1, lty = "dashed")
    dev.off()
}





tko_chim_timing <- function() {
    n_cls_min <- 20
    mat_nm <- "tko_chim_wt10"
    cgraph_id <- mat_nm
    mat <- scdb_mat(mat_nm)
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    # excluding Erythroid 1,2 and visceral and ExE endoderm
    excluded_colors <- mc_wt@color_key$color[c(22, 23, 28, 29)]

    load(sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))

    f <- mat@cell_metadata[colnames(mat@mat), "cell_type"] %in% c("KO", "host", "control")
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
    cell_type_list <- list(ko = c("KO"), host = c("host"), host_control = c("host", "control"), control = c("control"))

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
    chim_emb_summary$best_rank_host_control <- NA
    chim_emb_summary[rownames(time_dist_host_ko$host_control), "best_rank_host_control"] <- time_dist_best_match(
        atlas_time_dist = time_dist_host_ko$atlas_time_dist,
        query_time_dist = time_dist_host_ko$host_control
    )

    chim_emb_summary$best_rank_control <- NA
    chim_emb_summary[rownames(time_dist_host_ko$control), "best_rank_control"] <- time_dist_best_match(
        atlas_time_dist = time_dist_host_ko$atlas_time_dist,
        query_time_dist = time_dist_host_ko$control
    )

    write.table(chim_emb_summary, file = sprintf("%s/time_match_summary.txt", data_dir), sep = "\t", row.names = F)
}

plot_cumulative_distribution_ko_host_per_chimera <- function(mat_nm) {
    rank_to_time <- read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt", stringsAsFactors = F, h = T, sep = "\t")
    dev_time <- rank_to_time$developmental_time

    fig_dir <- sprintf("%s/%s/time_match", .scfigs_base, mat_nm)
    fig_scale <- 1.4
    lwd <- 15
    x_pos_text <- 6.7
    xlim_min <- min(dev_time)
    xlim_max <- max(dev_time)

    data_dir <- sprintf("data/%s/time_match", mat_nm)
    load(sprintf("%s/time_dist_host_ko.Rda", data_dir))

    host_dist_all <- time_dist_host_ko$host_control
    ko_dist_all <- time_dist_host_ko$ko

    chim_time_summary <- read.table(file = sprintf("%s/time_match_summary.txt", data_dir), sep = "\t", h = T, stringsAsFactors = F)
    f <- (chim_time_summary$host + chim_time_summary$control > 5) & (chim_time_summary$KO > 5)
    chim_embryos <- chim_time_summary$embryo[f]



    for (chimera in chim_embryos) {
        host_dens <- cumsum(host_dist_all[chimera, ]) / sum(host_dist_all[chimera, ])
        ko_dens <- cumsum(ko_dist_all[chimera, ]) / sum(ko_dist_all[chimera, ])

        # ks_ko_host_chim = ks.test(x = rep(c(1:153),host_dist_all[chimera,]),y = rep(c(1:153),ko_dist_all[chimera,]))
        ks_ko_host_chim <- ks.test(x = rep(dev_time, host_dist_all[chimera, ]), y = rep(dev_time, ko_dist_all[chimera, ]))

        N_host <- sum(host_dist_all[chimera, ])
        N_ko <- sum(ko_dist_all[chimera, ])

        png(sprintf("%s/cumulative_time_dist_%s.png", fig_dir, chimera), w = 500 * fig_scale, h = 400 * fig_scale)
        plot(dev_time, host_dens,
            type = "l", main = chimera, lwd = lwd,
            xlab = "Developmental time", ylab = "", ylim = c(0, 1), col = "gray30", cex.main = 3, cex.lab = 2, cex.axis = 2
        )

        lines(dev_time, ko_dens, type = "l", lwd = lwd, col = "#83c26d")

        text(
            x = rep(x_pos_text, 4), y = c(0.9, 0.75, 0.6, 0.45),
            labels = c(
                sprintf("D = %1.1e", ks_ko_host_chim$statistic), sprintf("p < %1.1e", ks_ko_host_chim$p.value),
                sprintf("N host = %d", N_host), sprintf("N KO = %d", N_ko)
            ), cex = 1.2
        )
        dev.off()
    }
}

plot_ko_vs_host_best_ranks <- function(mat_nm) {
    n_cls_min <- 19
    rank_to_time <- read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt", stringsAsFactors = F, h = T, sep = "\t")
    dev_time <- rank_to_time$developmental_time
    data_dir <- sprintf("data/%s/time_match", mat_nm)

    chim_emb_summary <- read.table(file = sprintf("%s/time_match_summary.txt", data_dir), sep = "\t", h = T, stringsAsFactors = F)

    f <- (chim_emb_summary$host + chim_emb_summary$control > n_cls_min) & (chim_emb_summary$KO > n_cls_min)

    rank_min <- min(c(chim_emb_summary$best_rank_ko[f], chim_emb_summary$best_rank_host_control[f])) - 1
    rank_max <- max(c(chim_emb_summary$best_rank_ko[f], chim_emb_summary$best_rank_host_control[f])) + 1
    time_min <- min(c(dev_time[chim_emb_summary$best_rank_ko[f]], dev_time[chim_emb_summary$best_rank_host_control[f]])) - 0.1
    time_max <- max(c(dev_time[chim_emb_summary$best_rank_ko[f]], dev_time[chim_emb_summary$best_rank_host_control[f]])) + 0.1

    col_emb <- rep("black", length(chim_emb_summary$embryo))


    fig_dir <- sprintf("%s/%s/time_match", .scfigs_base, mat_nm)




    png(sprintf("%s/best_rank_ko_vs_host_control.png", fig_dir))
    plot(chim_emb_summary$best_rank_ko[f], chim_emb_summary$best_rank_host_control[f],
        pch = 19, xlim = c(rank_min, rank_max), ylim = c(rank_min, rank_max), main = "KO vs host/control",
        xlab = "Rank KO cells", ylab = "Rank host/control cells", col = col_emb, cex = 2, cex.lab = 2
    )
    abline(a = 0, b = 1, lty = "dashed")
    dev.off()

    png(sprintf("%s/best_time_ko_vs_host_control.png", fig_dir))
    plot(rank_to_time$developmental_time[chim_emb_summary$best_rank_ko[f]],
        rank_to_time$developmental_time[chim_emb_summary$best_rank_host_control[f]],
        pch = 19,
        xlim = c(time_min, time_max), ylim = c(time_min, time_max), main = "KO vs host/control",
        xlab = "Time KO cells", ylab = "Time host/control cells", col = col_emb, cex = 2, cex.lab = 2
    )
    abline(a = 0, b = 1, lty = "dashed")
    dev.off()
}

tko_chimera_plot_cor_time_dist_per_emb <- function(mat_nm) {
    n_cls_min <- 19
    rank_to_time <- read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt", stringsAsFactors = F, h = T, sep = "\t")
    dev_time <- rank_to_time$developmental_time
    data_dir <- sprintf("data/%s/time_match", mat_nm)
    load(sprintf("%s/time_dist_host_ko.Rda", data_dir))

    host_dist_all <- time_dist_host_ko$host_control
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

    f <- (chim_time_summary$host + chim_time_summary$control > n_cls_min) & (chim_time_summary$KO > n_cls_min)
    chim_embryos <- chim_time_summary$embryo[f]

    fig_dir <- sprintf("%s/%s/time_match", .scfigs_base, mat_nm)

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

source("scripts/paper_figures/fig4.R")
source("scripts/paper_figures/fig2.R")

# Figure s5

generate_figure_s5_plots <- function() {
    if (!dir.exists("figs/paper_figs/fig_s5")) {
        dir.create("figs/paper_figs/fig_s5")
    }

    # fig s5c
    host_vs_dko_age("DKO12", "dko12_chim_wt10", T)
    host_vs_dko_age("DKO13", "dko13_chim_wt10", T)
    host_vs_dko_age("DKO23", "dko23_chim_wt10", T)

    # fig s5d
    plot_dko_chimera_dotplots(dko_genotype = "dko12")
    plot_dko_chimera_dotplots(dko_genotype = "dko13")
    plot_dko_chimera_dotplots(dko_genotype = "dko23")

    dko_dotplots_new_lineages(gt = "dko12", plot_pdf = T)
    dko_dotplots_new_lineages(gt = "dko13", plot_pdf = T)
    dko_dotplots_new_lineages(gt = "dko23", plot_pdf = T)
}


host_vs_dko_age <- function(genotype, mat_nm, plot_pdf = F) {
    n_cls_min <- 20

    rank_to_time <- read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt", stringsAsFactors = F, h = T, sep = "\t")
    dev_time <- rank_to_time$developmental_time

    fig_dir <- "figs/paper_figs/fig_s5"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    chimera_age <- read.table(sprintf("data/%s/time_match/time_match_summary.txt", mat_nm), sep = "\t", stringsAsFactors = F, h = T)

    f <- (chimera_age$host >= n_cls_min) & (chimera_age[, 1] >= n_cls_min)

    time_min <- 7.5
    time_max <- 8.2

    if (plot_pdf) {
        pdf(sprintf("%s/%s_best_time_ko_vs_host.pdf", fig_dir, genotype), useDingbats = F)
    } else {
        png(sprintf("%s/%s_best_time_ko_vs_host.png", fig_dir, genotype))
    }
    plot(rank_to_time$developmental_time[chimera_age$best_rank_ko[f]],
        rank_to_time$developmental_time[chimera_age$best_rank_host[f]],
        pch = 19,
        xlim = c(time_min, time_max), ylim = c(time_min, time_max), main = paste0(genotype, " vs host"),
        xlab = "Time DKO cells", ylab = "Time host cells", cex = 4, cex.lab = 1
    )
    abline(a = 0, b = 1, lty = "dashed")
    dev.off()
}

plot_dko_chimera_dotplots <- function(dko_genotype, plot_pdf = T) {
    ko_color <- "indianred3"
    host_color <- "gray30"
    mat_nm <- paste0(dko_genotype, "_chim_wt10")
    chim_freq <- chimera_dotplot_frequencies(mat_nm = mat_nm, minimal_number_of_cells = 20)

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    col_to_ct <- mc_wt@color_key$group
    col_to_ct[22] <- "Blood"
    names(col_to_ct) <- mc_wt@color_key$color


    ko_color <- "indianred3"
    host_color <- "gray30"

    highlighted_colors <- mc_wt@color_key$color[c(2, 3, 5, 6, 8, 12, 13, 14, 15, 17, 18, 19, 20, 22, 24, 27)]

    fig_dir <- "figs/paper_figs/fig_s5"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    fig_dir <- "figs/paper_figs/fig_s5/cell_type_dot_plots"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    fig_dir <- paste(fig_dir, dko_genotype, sep = "/")
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    tko_freq_n <- chim_freq$tko
    host_freq_n <- chim_freq$host
    wt_freq_n <- chim_freq$wt

    chim_freq_for_p_value_calculation <- chimera_dotplot_frequencies(mat_nm = mat_nm, minimal_number_of_cells = 90, downsample_number_of_cells = 90)
    tko_freq_ds <- chim_freq_for_p_value_calculation$tko
    host_freq_ds <- chim_freq_for_p_value_calculation$host
    wt_freq_ds <- chim_freq_for_p_value_calculation$wt

    tko_to_wt_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = tko_freq_ds[, ct_col], y = wt_freq_ds[, ct_col])
        return(p_val$p.value)
    })

    tko_to_host_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = tko_freq_ds[, ct_col], y = host_freq_ds[, ct_col])
        return(p_val$p.value)
    })

    host_to_wt_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = host_freq_ds[, ct_col], y = wt_freq_ds[, ct_col])
        return(p_val$p.value)
    })

    q_val_tko <- qvalue(p = tko_to_wt_p_values, pi0 = 1)
    q_val_host <- qvalue(p = host_to_wt_p_values, pi0 = 1)
    q_val_tko_to_host <- qvalue(p = tko_to_host_p_values, pi0 = 1)

    q_to_signif <- function(v) {
        v_signif <- sapply(v, function(x) {
            if (x >= 0.05) {
                a <- "ns"
            } else {
                a <- "*"
            }
            return(a)
        })
        return(v_signif)
    }


    genotype_color <- c("DKO" = ko_color, "Host" = host_color, "WT" = "gray70")
    names(genotype_color)[1] <- toupper(dko_genotype)



    stat_comparison <- data.frame(
        group1 = c(rep(toupper(dko_genotype), length(highlighted_colors)), rep("Host", length(highlighted_colors)), rep(toupper(dko_genotype), length(highlighted_colors))),
        group2 = c(rep("Host", length(highlighted_colors)), rep("WT", length(highlighted_colors)), rep("WT", length(highlighted_colors))),
        cell_type = c(col_to_ct[names(q_val_tko_to_host$qvalues)], col_to_ct[names(q_val_host$qvalues)], col_to_ct[names(q_val_tko$qvalues)]),
        cell_type_color = c(names(q_val_tko_to_host$qvalues), names(q_val_host$qvalues), names(q_val_tko$qvalues)),
        q.val = c(q_val_tko_to_host$qvalues, q_val_host$qvalues, q_val_tko$qvalues),
        q.signif = c(q_to_signif(q_val_tko_to_host$qvalues), q_to_signif(q_val_host$qvalues), q_to_signif(q_val_tko$qvalues)), stringsAsFactors = F
    )

    plot_list <- list()

    col_to_ct[20] <- "Haematoendothelial prog."
    col_to_ct[14] <- "Later. & interm. mesoderm"

    for (ct_col in highlighted_colors) {
        main_tag <- gsub("/", "_", col_to_ct[ct_col])

        df_plot_points <- data.frame(
            genotype = factor(x = c(rep(toupper(dko_genotype), nrow(tko_freq_n)), rep("Host", nrow(host_freq_n)), rep("WT", nrow(wt_freq_n))), levels = c(toupper(dko_genotype), "Host", "WT")),
            freq = c(tko_freq_n[, ct_col], host_freq_n[, ct_col], wt_freq_n[, ct_col])
        )

        my_comparisons <- list(c(toupper(dko_genotype), "Host"), c(toupper(dko_genotype), "WT"), c("Host", "WT"))

        stat.test <- compare_means(data = df_plot_points, formula = freq ~ genotype)

        stat_f <- stat_comparison[stat_comparison$cell_type_color == ct_col, ]


        p <- ggplot(data = df_plot_points, aes(x = genotype, y = freq)) +
            geom_dotplot(aes(fill = genotype), dotsize = 1.3, binaxis = "y", stackdir = "center", show.legend = F) +
            stat_pvalue_manual(stat_f, y.position = max(df_plot_points$freq) * 1.1, step.increase = 0.1, label = "q.signif") +
            scale_fill_manual(values = genotype_color) +
            ggtitle(label = main_tag) +
            theme(plot.title = element_text(hjust = 0.5, size = 10)) +
            ylab("") +
            ylim(0, max(df_plot_points$freq) * 1.4) +
            xlab("")
        # theme(axis.text.x = element_text(size=14))
        # stat_compare_means(label = "p.signif",comparisons = my_comparisons) +

        plot_list[[ct_col]] <- p
        if (plot_pdf) {
            ggsave(filename = sprintf("%s/2N_%s.pdf", fig_dir, main_tag), width = 3, height = 2.3, plot = p)
        } else {
            ggsave(filename = sprintf("%s/2N_%s.png", fig_dir, main_tag), width = 3, height = 2.3, plot = p)
        }
    }

    p_all <- grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)

    if (plot_pdf) {
        ggsave(filename = sprintf("%s/2N_all_cell_types.pdf", fig_dir), width = 8.5, height = 6.5, plot = p_all)
    } else {
        ggsave(filename = sprintf("%s/2N_all_cell_types.png", fig_dir), width = 8.5, height = 6.5, plot = p_all)
    }
}



dko_dotplots_new <- function(gt = "dko12", plot_pdf = T) {
    ko_color <- "indianred3"
    host_color <- "gray30"
    mat_nm <- paste0(gt, "_chim_wt10")
    dko_freq <- dko_dotplot_frequencies(mat_nm = mat_nm)
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    col_to_ct <- mc_wt@color_key$group
    col_to_ct[22] <- "Blood"
    names(col_to_ct) <- mc_wt@color_key$color


    ko_color <- "indianred3"
    host_color <- "gray30"

    highlighted_colors <- mc_wt@color_key$color[c(2, 3, 5, 6, 8, 12, 13, 14, 15, 17, 18, 19, 20, 22, 24, 27)]

    if (!dir.exists("figs/paper_figs")) {
        dir.create("figs/paper_figs")
    }

    fig_dir <- "figs/paper_figs/fig_s5"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    fig_dir <- "figs/paper_figs/fig_s5/cell_type_dot_plots"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    fig_dir <- sprintf("figs/paper_figs/fig_s5/cell_type_dot_plots/%s", gt)
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }


    ko_freq <- dko_freq$KO
    ko_freq_n <- ko_freq / rowSums(ko_freq)
    host_freq <- dko_freq$host
    host_freq_n <- host_freq / rowSums(host_freq)
    wt_freq <- dko_freq$wt
    wt_freq_n <- wt_freq / rowSums(wt_freq)

    ko_to_wt_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = ko_freq_n[, ct_col], y = wt_freq_n[, ct_col])
        return(p_val$p.value)
    })

    ko_to_host_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = ko_freq_n[, ct_col], y = host_freq_n[, ct_col])
        return(p_val$p.value)
    })

    host_to_wt_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = host_freq_n[, ct_col], y = wt_freq_n[, ct_col])
        return(p_val$p.value)
    })

    q_val_ko <- qvalue(p = ko_to_wt_p_values, pi0 = 1)
    q_val_host <- qvalue(p = host_to_wt_p_values, pi0 = 1)
    q_val_ko_to_host <- qvalue(p = ko_to_host_p_values, pi0 = 1)

    q_to_signif <- function(v) {
        v_signif <- sapply(v, function(x) {
            if (x >= 0.05) {
                a <- "ns"
            } else {
                a <- "*"
            }
            return(a)
        })
        return(v_signif)
    }

    genotype_color <- c("DKO" = ko_color, "Host" = host_color, "WT" = "gray70")


    stat_comparison <- data.frame(
        group1 = c(rep("DKO", length(highlighted_colors)), rep("Host", length(highlighted_colors)), rep("DKO", length(highlighted_colors))),
        group2 = c(rep("Host", length(highlighted_colors)), rep("WT", length(highlighted_colors)), rep("WT", length(highlighted_colors))),
        cell_type = c(col_to_ct[names(q_val_ko_to_host$qvalues)], col_to_ct[names(q_val_host$qvalues)], col_to_ct[names(q_val_ko$qvalues)]),
        cell_type_color = c(names(q_val_ko_to_host$qvalues), names(q_val_host$qvalues), names(q_val_ko$qvalues)),
        q.val = c(q_val_ko_to_host$qvalues, q_val_host$qvalues, q_val_ko$qvalues),
        q.signif = c(q_to_signif(q_val_ko_to_host$qvalues), q_to_signif(q_val_host$qvalues), q_to_signif(q_val_ko$qvalues)), stringsAsFactors = F
    )

    for (ct_col in highlighted_colors) {
        main_tag <- gsub("/", "_", col_to_ct[ct_col])

        df_plot_points <- data.frame(
            genotype = factor(x = c(rep("DKO", nrow(ko_freq_n)), rep("Host", nrow(host_freq_n)), rep("WT", nrow(wt_freq_n))), levels = c("DKO", "Host", "WT")),
            freq = c(ko_freq_n[, ct_col], host_freq_n[, ct_col], wt_freq_n[, ct_col])
        )

        my_comparisons <- list(c("DKO", "Host"), c("DKO", "WT"), c("Host", "WT"))

        stat.test <- compare_means(data = df_plot_points, formula = freq ~ genotype)

        stat_f <- stat_comparison[stat_comparison$cell_type_color == ct_col, ]


        p <- ggplot(data = df_plot_points, aes(x = genotype, y = freq)) +
            geom_dotplot(aes(fill = genotype), dotsize = 1.3, binaxis = "y", stackdir = "center", show.legend = F) +
            stat_pvalue_manual(stat_f, y.position = max(df_plot_points$freq) * 1.1, step.increase = 0.1, label = "q.signif") +
            scale_fill_manual(values = genotype_color) +
            ggtitle(label = main_tag) +
            theme(plot.title = element_text(hjust = 0.5, size = 10)) +
            ylab("") +
            ylim(0, max(df_plot_points$freq) * 1.4) +
            xlab("")
        # theme(axis.text.x = element_text(size=14))
        # stat_compare_means(label = "p.signif",comparisons = my_comparisons) +


        if (plot_pdf) {
            ggsave(filename = sprintf("%s/2N_%s.pdf", fig_dir, main_tag), width = 3, height = 2.3, plot = p)
        } else {
            ggsave(filename = sprintf("%s/2N_%s.png", fig_dir, main_tag), width = 3, height = 2.3, plot = p)
        }
    }
}

dko_dotplot_frequencies_lineages <- function(mat_nm) {


    # time window Et7.5 - Et8.1
    t_f <- c(125:153)

    mat_chim <- scdb_mat(mat_nm)

    load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))
    ko_query_cls_col <- cmp_annot$query_cls_col[mat_chim@cell_metadata[names(cmp_annot$query_cls_col), "cell_type"] %in% c("DKO12", "DKO13", "DKO23")]
    ko_query_cls <- names(ko_query_cls_col)
    host_query_cls_col <- cmp_annot$query_cls_col[mat_chim@cell_metadata[names(cmp_annot$query_cls_col), "cell_type"] %in% c("control", "host")]
    host_query_cls <- names(host_query_cls_col)

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")
    col_to_ct <- mc_wt@color_key$group
    names(col_to_ct) <- mc_wt@color_key$color
    # visceral and extraembryonic endoderm are excluded
    excluded_colors <- c("#F6BFCB", "#7F6874")
    included_colors <- setdiff(mc_wt@color_key$color, excluded_colors)

    ct_to_lineage <- rep("other", 27)
    ct_to_lineage[3:6] <- "Ectoderm"
    ct_to_lineage[25:27] <- "Endoderm"
    ct_to_lineage[c(8, 12, 13, 14, 15, 16)] <- "Embryonic mesoderm"
    ct_to_lineage[17:19] <- "Extraembryonic mesoderm"
    ct_to_lineage[21:23] <- "Blood"
    names(ct_to_lineage) <- mc_wt@color_key$color[c(1:27)]


    # Blood progenitors, Erythroid 1, Erythroid 2 are aggregated into one cell type blood
    old_to_new_col <- mc_wt@color_key$color
    old_to_new_col[c(21, 22, 23)] <- "#C72228"
    names(old_to_new_col) <- mc_wt@color_key$color

    df_chim <- read.table(sprintf("data/%s/time_match/time_match_summary.txt", mat_nm), sep = "\t", h = T, stringsAsFactors = F)
    rownames(df_chim) <- df_chim$embryo

    f <- (df_chim$host > 20) & (df_chim[, 1] > 20)
    df_chim <- df_chim[f, ]

    embryos <- df_chim$embryo
    embryos <- embryos[order(df_chim[embryos, "best_rank_host"])]

    f_t_q <- df_chim[embryos, "best_rank_host"] %in% t_f

    ko_query_cls_f <- ko_query_cls[(mat_chim@cell_metadata[ko_query_cls, "cell_type"] %in% c("DKO12", "DKO13", "DKO23")) & (ko_query_cls_col[ko_query_cls] %in% included_colors)]
    host_query_cls_f <- host_query_cls[(mat_chim@cell_metadata[host_query_cls, "cell_type"] %in% c("host")) & (host_query_cls_col[host_query_cls] %in% included_colors)]
    wt10_cls <- names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] %in% included_colors]

    emb_wt_age <- unique(mat_chim@cell_metadata[wt10_cls, c("transcriptional_rank", "age_group")])
    emb_wt_age <- emb_wt_age[order(emb_wt_age$transcriptional_rank), ]

    wt10_emb_vs_ct <- table(
        mat_chim@cell_metadata[wt10_cls, "transcriptional_rank"],
        ct_to_lineage[mc_wt@colors[mc_wt@mc[wt10_cls]]]
    )
    wt10_emb_vs_ct <- wt10_emb_vs_ct[t_f, ]
    ko_emb_vs_ct <- table(
        mat_chim@cell_metadata[ko_query_cls_f, "embryo"],
        ct_to_lineage[ko_query_cls_col[ko_query_cls_f]]
    )
    ko_emb_vs_ct <- ko_emb_vs_ct[embryos[f_t_q], ]
    host_emb_vs_ct <- table(
        mat_chim@cell_metadata[host_query_cls_f, "embryo"],
        ct_to_lineage[host_query_cls_col[host_query_cls_f]]
    )
    host_emb_vs_ct <- host_emb_vs_ct[embryos[f_t_q], ]


    return(list(wt = wt10_emb_vs_ct, KO = ko_emb_vs_ct, host = host_emb_vs_ct))
}



dko_dotplot_frequencies <- function(mat_nm) {


    # time window Et7.75 - Et8.1
    t_f <- c(125:153)
    mat_chim <- scdb_mat(mat_nm)

    load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))
    ko_query_cls_col <- cmp_annot$query_cls_col[mat_chim@cell_metadata[names(cmp_annot$query_cls_col), "cell_type"] %in% c("DKO12", "DKO13", "DKO23")]
    ko_query_cls <- names(ko_query_cls_col)
    host_query_cls_col <- cmp_annot$query_cls_col[mat_chim@cell_metadata[names(cmp_annot$query_cls_col), "cell_type"] %in% c("host")]
    host_query_cls <- names(host_query_cls_col)

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")
    col_to_ct <- mc_wt@color_key$group
    names(col_to_ct) <- mc_wt@color_key$color
    # visceral and extraembryonic endoderm are excluded
    excluded_colors <- c("#F6BFCB", "#7F6874")
    included_colors <- setdiff(mc_wt@color_key$color, excluded_colors)

    # Blood progenitors, Erythroid 1, Erythroid 2 are aggregated into one cell type blood
    old_to_new_col <- mc_wt@color_key$color
    old_to_new_col[c(21, 22, 23)] <- "#C72228"
    names(old_to_new_col) <- mc_wt@color_key$color

    df_tetra <- read.table(sprintf("data/%s/time_match/time_match_summary.txt", mat_nm), sep = "\t", h = T, stringsAsFactors = F)
    rownames(df_tetra) <- df_tetra$embryo

    f <- (df_tetra$host > 20) & (df_tetra[, 1] > 20)
    df_tetra <- df_tetra[f, ]

    embryos <- df_tetra$embryo
    embryos <- embryos[order(df_tetra[embryos, "best_rank_host"])]

    f_t_q <- df_tetra[embryos, "best_rank_host"] %in% t_f

    ko_query_cls_f <- ko_query_cls[(mat_chim@cell_metadata[ko_query_cls, "cell_type"] %in% c("DKO12", "DKO13", "DKO23")) & (ko_query_cls_col[ko_query_cls] %in% included_colors)]
    host_query_cls_f <- host_query_cls[(mat_chim@cell_metadata[host_query_cls, "cell_type"] %in% c("host")) & (host_query_cls_col[host_query_cls] %in% included_colors)]
    wt10_cls <- names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] %in% included_colors]

    emb_wt_age <- unique(mat_chim@cell_metadata[wt10_cls, c("transcriptional_rank", "age_group")])
    emb_wt_age <- emb_wt_age[order(emb_wt_age$transcriptional_rank), ]

    wt10_emb_vs_ct <- table(
        mat_chim@cell_metadata[wt10_cls, "transcriptional_rank"],
        factor(x = old_to_new_col[mc_wt@colors[mc_wt@mc[wt10_cls]]], levels = mc_wt@color_key$color[c(1:27)[-c(21, 23)]])
    )
    wt10_emb_vs_ct <- wt10_emb_vs_ct[t_f, ]
    ko_emb_vs_ct <- table(
        mat_chim@cell_metadata[ko_query_cls_f, "embryo"],
        factor(x = old_to_new_col[ko_query_cls_col[ko_query_cls_f]], levels = mc_wt@color_key$color[c(1:27)[-c(21, 23)]])
    )
    ko_emb_vs_ct <- ko_emb_vs_ct[embryos[f_t_q], ]
    host_emb_vs_ct <- table(
        mat_chim@cell_metadata[host_query_cls_f, "embryo"],
        factor(x = old_to_new_col[host_query_cls_col[host_query_cls_f]], levels = mc_wt@color_key$color[c(1:27)[-c(21, 23)]])
    )
    host_emb_vs_ct <- host_emb_vs_ct[embryos[f_t_q], ]


    return(list(wt = wt10_emb_vs_ct, KO = ko_emb_vs_ct, host = host_emb_vs_ct))
}




dko_dotplots_new_lineages <- function(gt = "dko12", plot_pdf = T) {
    ko_color <- "indianred3"
    host_color <- "gray30"
    mat_nm <- paste0(gt, "_chim_wt10")
    dko_freq <- dko_dotplot_frequencies_lineages(mat_nm)

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    highlighted_colors <- c("Endoderm", "Ectoderm", "Embryonic mesoderm", "Extraembryonic mesoderm", "Blood")


    if (!dir.exists("figs/paper_figs")) {
        dir.create("figs/paper_figs")
    }

    fig_dir <- "figs/paper_figs/fig_s5"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    fig_dir <- "figs/paper_figs/fig_s5/cell_type_dot_plots"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    fig_dir <- sprintf("figs/paper_figs/fig_s5/cell_type_dot_plots/%s", gt)
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }


    ko_freq <- dko_freq$KO
    ko_freq_n <- ko_freq / rowSums(ko_freq)
    host_freq <- dko_freq$host
    host_freq_n <- host_freq / rowSums(host_freq)
    wt_freq <- dko_freq$wt
    wt_freq_n <- wt_freq / rowSums(wt_freq)

    ko_to_wt_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = ko_freq_n[, ct_col], y = wt_freq_n[, ct_col])
        return(p_val$p.value)
    })

    ko_to_host_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = ko_freq_n[, ct_col], y = host_freq_n[, ct_col])
        return(p_val$p.value)
    })

    host_to_wt_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = host_freq_n[, ct_col], y = wt_freq_n[, ct_col])
        return(p_val$p.value)
    })

    q_val_ko <- qvalue(p = ko_to_wt_p_values, pi0 = 1)
    q_val_host <- qvalue(p = host_to_wt_p_values, pi0 = 1)
    q_val_ko_to_host <- qvalue(p = ko_to_host_p_values, pi0 = 1)

    q_to_signif <- function(v) {
        v_signif <- sapply(v, function(x) {
            if (x >= 0.05) {
                a <- "ns"
            } else {
                a <- "*"
            }
            return(a)
        })
        return(v_signif)
    }

    genotype_color <- c("DKO" = ko_color, "Host" = host_color, "WT" = "gray70")


    stat_comparison <- data.frame(
        group1 = c(rep("DKO", length(highlighted_colors)), rep("Host", length(highlighted_colors)), rep("DKO", length(highlighted_colors))),
        group2 = c(rep("Host", length(highlighted_colors)), rep("WT", length(highlighted_colors)), rep("WT", length(highlighted_colors))),
        cell_type = c(names(q_val_ko_to_host$qvalues), names(q_val_host$qvalues), names(q_val_ko$qvalues)),
        q.val = c(q_val_ko_to_host$qvalues, q_val_host$qvalues, q_val_ko$qvalues),
        q.signif = c(q_to_signif(q_val_ko_to_host$qvalues), q_to_signif(q_val_host$qvalues), q_to_signif(q_val_ko$qvalues)), stringsAsFactors = F
    )

    for (ct_col in highlighted_colors) {
        main_tag <- ct_col

        df_plot_points <- data.frame(
            genotype = factor(x = c(rep("DKO", nrow(ko_freq_n)), rep("Host", nrow(host_freq_n)), rep("WT", nrow(wt_freq_n))), levels = c("DKO", "Host", "WT")),
            freq = c(ko_freq_n[, ct_col], host_freq_n[, ct_col], wt_freq_n[, ct_col])
        )

        my_comparisons <- list(c("DKO", "Host"), c("DKO", "WT"), c("Host", "WT"))

        stat.test <- compare_means(data = df_plot_points, formula = freq ~ genotype)

        stat_f <- stat_comparison[stat_comparison$cell_type == ct_col, ]


        p <- ggplot(data = df_plot_points, aes(x = genotype, y = freq)) +
            geom_dotplot(aes(fill = genotype), dotsize = 1.3, binaxis = "y", stackdir = "center", show.legend = F) +
            stat_pvalue_manual(stat_f, y.position = max(df_plot_points$freq) * 1.1, step.increase = 0.1, label = "q.signif") +
            scale_fill_manual(values = genotype_color) +
            ggtitle(label = main_tag) +
            theme(plot.title = element_text(hjust = 0.5, size = 10)) +
            ylab("") +
            ylim(0, max(df_plot_points$freq) * 1.4) +
            xlab("")
        # theme(axis.text.x = element_text(size=14))
        # stat_compare_means(label = "p.signif",comparisons = my_comparisons) +


        if (plot_pdf) {
            ggsave(filename = sprintf("%s/2N_%s.pdf", fig_dir, main_tag), width = 3, height = 2.3, plot = p)
        } else {
            ggsave(filename = sprintf("%s/2N_%s.png", fig_dir, main_tag), width = 3, height = 2.3, plot = p)
        }
    }
}

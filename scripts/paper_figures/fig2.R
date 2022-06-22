

generate_figure2_plots <- function() {
    fig_dir <- "figs/paper_figs/fig2"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }


    # fig 2c
    atlas_projection_tko_chimera()

    # fig 2d
    host_vs_ko_age(T)

    # fig 2e
    tko_barplot_ct_frequency(mat_nm = "tko_chim_wt10", ko_type = "KO", plot_pdf = T)
    tko_barplot_ct_frequency(mat_nm = "tko_chim_wt10", ko_type = c("control", "host"), plot_pdf = T, tag = "control_host")

    # fig 2f
    plot_chimera_dotplots()
    plot_tetraploid_dotplots()
}

atlas_projection_tko_chimera <- function(plot_pdf = F) {
    mat_chim <- scdb_mat("tko_chim_wt10")

    gset <- scdb_gset("tko_chim_wt10")
    feat_genes <- names(gset@gene_set)



    ko_cls <- colnames(mat_chim@mat)[mat_chim@cell_metadata[colnames(mat_chim@mat), "cell_type"] == "KO"]
    host_cls <- colnames(mat_chim@mat)[mat_chim@cell_metadata[colnames(mat_chim@mat), "cell_type"] == "host"]

    fig_dir <- "figs/paper_figs/fig2"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    mat_query <- mat_chim@mat[, ko_cls]

    fn <- "figs/paper_figs/fig2/atlas_projection_tko_chim_ko_cls_new_12.png"
    w <- 1000
    h <- 1000
    if (plot_pdf) {
        fn <- gsub(pattern = ".png", replacement = ".pdf", x = fn)
        w <- 1000 / 72
        h <- 1000 / 72
    }
    atlas_proj_on_wt10(mat_query = mat_query, feat_genes = feat_genes, fn = fn, cex_points = 1.2, w = w, h = h, plot_pdf = plot_pdf, plot_gray_background = F)

    mat_query <- mat_chim@mat[, host_cls]

    fn <- "figs/paper_figs/fig2/atlas_projection_tko_chim_host_control_cls_new_12.png"

    atlas_proj_on_wt10(mat_query = mat_query, feat_genes = feat_genes, fn = fn, cex_points = 1.2, w = w, h = h, plot_pdf = plot_pdf, plot_gray_background = F)
}



host_vs_ko_age <- function(plot_pdf = F) {
    n_cls_min <- 20

    rank_to_time <- read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt", stringsAsFactors = F, h = T, sep = "\t")
    dev_time <- rank_to_time$developmental_time

    fig_dir <- "figs/paper_figs/fig2"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    chimera_age <- read.table("data/tko_chim_wt10/time_match/time_match_summary.txt", sep = "\t", stringsAsFactors = F, h = T)

    f <- (chimera_age$control + chimera_age$host >= n_cls_min) & (chimera_age$KO >= n_cls_min)

    time_min <- 6.9
    time_max <- 8.2

    if (plot_pdf) {
        pdf(sprintf("%s/best_time_ko_vs_host_control.pdf", fig_dir), useDingbats = F)
    } else {
        png(sprintf("%s/best_time_ko_vs_host_control.png", fig_dir))
    }
    plot(rank_to_time$developmental_time[chimera_age$best_rank_ko[f]],
        rank_to_time$developmental_time[chimera_age$best_rank_host_control[f]],
        pch = 19,
        xlim = c(time_min, time_max), ylim = c(time_min, time_max), main = "KO vs host/control",
        xlab = "Time KO cells", ylab = "Time host/control cells", cex = 4, cex.lab = 1
    )
    abline(a = 0, b = 1, lty = "dashed")
    dev.off()
}

tko_barplot_ct_frequency <- function(mat_nm, ko_type = "KO", plot_pdf = F, tag = "KO") {
    n_cls_min <- 19
    fig_dir <- "figs/paper_figs/fig2"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }


    df_chim <- read.table(sprintf("data/%s/time_match/time_match_summary.txt", mat_nm), sep = "\t", stringsAsFactors = F, h = T)
    rownames(df_chim) <- df_chim$embryo
    mat <- scdb_mat(mat_nm)
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")
    col_to_rank <- c(1:nrow(mc_wt@color_key))
    names(col_to_rank) <- mc_wt@color_key$color
    col_to_ct <- mc_wt@color_key$group
    names(col_to_ct) <- mc_wt@color_key$color

    excluded_colors <- c("#F6BFCB", "#7F6874")
    included_colors <- setdiff(unique(mc_wt@color_key$color), excluded_colors)



    chim_embryos <- df_chim$embryo[(df_chim$KO > n_cls_min) & (df_chim$control + df_chim$host > n_cls_min)]
    chim_embryos <- chim_embryos[order(df_chim[chim_embryos, "best_rank_host_control"])]

    tmp <- matrix(0, nrow = length(chim_embryos), ncol = length(included_colors))
    rownames(tmp) <- chim_embryos
    colnames(tmp) <- included_colors

    load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))

    query_cls_col <- cmp_annot$query_cls_col
    query_cls <- names(query_cls_col)[!(query_cls_col %in% excluded_colors)]
    query_cls <- query_cls[mat@cell_metadata[query_cls, "embryo"] %in% chim_embryos]

    modified_cell_type_levels <- c(
        mc_wt@color_key$group[1:7], c("space1"),
        mc_wt@color_key$group[c(9, 10, 11)], c("space2"),
        mc_wt@color_key$group[c(8, 12, 13, 14, 15, 16)], c("space3"),
        mc_wt@color_key$group[c(17, 18, 19)], c("space4"),
        mc_wt@color_key$group[c(20:23)], c("space5"),
        mc_wt@color_key$group[c(24:27)]
    )

    modified_colors <- c(
        mc_wt@color_key$color[1:7], c("white"),
        mc_wt@color_key$color[c(9, 10, 11)], c("white"),
        mc_wt@color_key$color[c(8, 12, 13, 14, 15, 16)], c("white"),
        mc_wt@color_key$color[c(17, 18, 19)], c("white"),
        mc_wt@color_key$color[c(20:23)], c("white"),
        mc_wt@color_key$color[c(24:27)]
    )

    filtered_cls <- query_cls[mat@cell_metadata[query_cls, "cell_type"] %in% ko_type]
    filtered_vs_ct <- table(factor(x = mat@cell_metadata[filtered_cls, "embryo"], levels = chim_embryos), factor(x = col_to_ct[query_cls_col[filtered_cls]], levels = modified_cell_type_levels))
    filtered_vs_ct_n <- filtered_vs_ct / rowSums(filtered_vs_ct)
    filtered_vs_ct_n[is.na(filtered_vs_ct_n)] <- 0

    filtered_vs_ct_n[1:2, c("space1", "space5")] <- 0.04
    filtered_vs_ct_n[3:nrow(filtered_vs_ct_n), c("space1", "space3", "space4", "space5")] <- 0.02

    filtered_vs_ct_n <- filtered_vs_ct_n / rowSums(filtered_vs_ct_n)

    if (plot_pdf) {
        pdf(sprintf("%s/barplot_ct_freq_%s.pdf", fig_dir, tag), w = 12, h = 7.5, useDingbats = F)
        barplot(t(filtered_vs_ct_n), col = modified_colors, las = 2, axes = F, axisnames = F, border = NA)
        dev.off()
    } else {
        png(sprintf("%s/barplot_ct_freq_%s.png", fig_dir, tag), w = 1250, h = 750)
        barplot(t(filtered_vs_ct_n), col = modified_colors, las = 2, axes = F, axisnames = F, border = NA)
        dev.off()
    }
}

plot_chimera_dotplots <- function(plot_pdf = T, included_transcriptional_ranks = NULL, highlighted_colors = NULL, minimal_number_of_cells_for_p_value = 100) {
    ko_color <- "indianred3"
    host_color <- "gray30"
    mat_nm <- "tko_chim_wt10"

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    if (is.null(included_transcriptional_ranks)) {
        included_transcriptional_ranks <- c(125:153)
    }
    if (is.null(highlighted_colors)) {
        highlighted_colors <- mc_wt@color_key$color[c(2, 3, 5, 6, 8, 12, 13, 14, 15, 17, 18, 19, 20, 22, 24, 27)]
    }

    chim_freq <- chimera_dotplot_frequencies(mat_nm = mat_nm, minimal_number_of_cells = 20, included_transcriptional_ranks = included_transcriptional_ranks)

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    col_to_ct <- mc_wt@color_key$group
    col_to_ct[22] <- "Blood"
    names(col_to_ct) <- mc_wt@color_key$color


    ko_color <- "indianred3"
    host_color <- "gray30"

    mat_nm <- "tko_chim_wt10"



    fig_dir <- "figs/paper_figs/fig2"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    fig_dir <- "figs/paper_figs/fig2/cell_type_dot_plots"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    tko_freq_n <- chim_freq$tko
    host_freq_n <- chim_freq$host
    wt_freq_n <- chim_freq$wt

    chim_freq_for_p_value_calculation <- chimera_dotplot_frequencies(mat_nm = mat_nm, minimal_number_of_cells = minimal_number_of_cells_for_p_value, downsample_number_of_cells = minimal_number_of_cells_for_p_value)
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

    if (length(highlighted_colors) > 1) {
        q_val_tko <- qvalue(p = tko_to_wt_p_values, pi0 = 1)
        q_val_host <- qvalue(p = host_to_wt_p_values, pi0 = 1)
        q_val_tko_to_host <- qvalue(p = tko_to_host_p_values, pi0 = 1)
    } else {
        q_val_tko <- list(qvalues = tko_to_wt_p_values)
        q_val_host <- list(qvalues = host_to_wt_p_values)
        q_val_tko_to_host <- list(qvalues = tko_to_host_p_values)
    }


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

    genotype_color <- c("TKO" = ko_color, "Host/Control" = host_color, "WT" = "gray70")



    stat_comparison <- data.frame(
        group1 = c(rep("TKO", length(highlighted_colors)), rep("Host/Control", length(highlighted_colors)), rep("TKO", length(highlighted_colors))),
        group2 = c(rep("Host/Control", length(highlighted_colors)), rep("WT", length(highlighted_colors)), rep("WT", length(highlighted_colors))),
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
            genotype = factor(x = c(rep("TKO", nrow(tko_freq_n)), rep("Host/Control", nrow(host_freq_n)), rep("WT", nrow(wt_freq_n))), levels = c("TKO", "Host/Control", "WT")),
            freq = c(tko_freq_n[, ct_col], host_freq_n[, ct_col], wt_freq_n[, ct_col])
        )

        my_comparisons <- list(c("TKO", "Host/Control"), c("TKO", "WT"), c("Host/Control", "WT"))

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


aggregate_blood_subtypes_into_one_type <- function(color_vector) {

    # I replace Blood progenitors color and Erythroid 2 color by Erythroid 1 color #C72228
    color_vector_names <- names(color_vector)
    color_vector <- as.character(color_vector)
    color_vector[color_vector %in% c("#c9a997", "#EF4E22")] <- "#C72228"
    names(color_vector) <- color_vector_names

    return(color_vector)
}

downsample_cells_indexed_by_metadata <- function(cells, cells_metadata, n_downsample, seed = NULL) {
    n_cells_per_metadata <- table(cells_metadata)
    included_metadata_levels <- names(n_cells_per_metadata)[n_cells_per_metadata >= n_downsample]

    f <- cells_metadata %in% included_metadata_levels
    cells <- cells[f]
    cells_metadata <- cells_metadata[f]

    cells_ds <- tapply(cells, cells_metadata, function(v) {
        if (!is.null(seed)) {
            set.seed(seed)
        }
        v_ds <- sample(v, size = n_downsample)
        return(v_ds)
    })
    cells_ds <- unlist(cells_ds)
    return(cells_ds)
}

chimera_dotplot_frequencies <- function(mat_nm, minimal_number_of_cells = 20, downsample_number_of_cells = NULL, included_transcriptional_ranks = c(125:153)) {
    mat_chim <- scdb_mat(mat_nm)
    # time window Et7.75 - Et8.1

    df_chim <- read.table(sprintf("data/%s/time_match/time_match_summary.txt", mat_nm), sep = "\t", h = T, stringsAsFactors = F)
    rownames(df_chim) <- df_chim$embryo
    if (mat_nm == "tko_chim_wt10") {
        f <- (df_chim$control + df_chim$host >= minimal_number_of_cells) & (df_chim$KO >= minimal_number_of_cells)
    } else {
        f <- (df_chim$host >= minimal_number_of_cells) & (df_chim[, 1] >= minimal_number_of_cells)
    }
    df_chim <- df_chim[f, ]
    if (mat_nm == "tko_chim_wt10") {
        f <- df_chim[, "best_rank_host_control"] %in% included_transcriptional_ranks
    } else {
        f <- df_chim[, "best_rank_host"] %in% included_transcriptional_ranks
    }
    df_chim <- df_chim[f, ]
    included_chimeras <- df_chim$embryo

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    # visceral and extraembryonic endoderm are excluded
    included_colors <- mc_wt@color_key$color[1:27]


    load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))

    # Modify blood subtypes color to one color
    query_cells_color <- aggregate_blood_subtypes_into_one_type(cmp_annot$query_cls_col)
    wt_cells_color <- mc_wt@colors[mc_wt@mc]
    names(wt_cells_color) <- names(mc_wt@mc)
    wt_cells_color <- aggregate_blood_subtypes_into_one_type(wt_cells_color)
    # remove Blood progenitors and Erythroid 2 from color levels
    included_colors <- included_colors[-c(21, 23)]

    query_cells_color <- query_cells_color[query_cells_color %in% included_colors]
    wt_cells_color <- wt_cells_color[wt_cells_color %in% included_colors]

    ko_cells <- names(query_cells_color)[mat_chim@cell_metadata[names(query_cells_color), "cell_type"] %in% c("KO", "DKO12", "DKO13", "DKO23")]
    host_cells <- names(query_cells_color)[mat_chim@cell_metadata[names(query_cells_color), "cell_type"] %in% c("control", "host")]
    wt_cells <- names(wt_cells_color)

    # downsample cells to common number per embryo
    if (!is.null(downsample_number_of_cells)) {
        if (minimal_number_of_cells < downsample_number_of_cells) {
            stop("minimal_number_of_cells smaller than downsample_number_of_cells")
        }
        ko_cells <- downsample_cells_indexed_by_metadata(cells = ko_cells, cells_metadata = mat_chim@cell_metadata[ko_cells, "embryo"], n_downsample = downsample_number_of_cells, seed = 123)
        host_cells <- downsample_cells_indexed_by_metadata(cells = host_cells, cells_metadata = mat_chim@cell_metadata[host_cells, "embryo"], n_downsample = downsample_number_of_cells, seed = 123)
        wt_cells <- downsample_cells_indexed_by_metadata(cells = wt_cells, cells_metadata = mat_chim@cell_metadata[wt_cells, "transcriptional_rank"], n_downsample = downsample_number_of_cells, seed = 123)
    }


    # compute two way tables
    ko_emb_vs_ct <- compute_two_way_table(
        values_row = mat_chim@cell_metadata[ko_cells, "embryo"],
        values_col = query_cells_color[ko_cells],
        included_levels_row = included_chimeras,
        included_levels_col = included_colors, normalize_rows = T
    )

    host_emb_vs_ct <- compute_two_way_table(
        values_row = mat_chim@cell_metadata[host_cells, "embryo"],
        values_col = query_cells_color[host_cells],
        included_levels_row = included_chimeras,
        included_levels_col = included_colors, normalize_rows = T
    )

    wt10_emb_vs_ct <- compute_two_way_table(
        values_row = mat_chim@cell_metadata[wt_cells, "transcriptional_rank"],
        values_col = wt_cells_color[wt_cells],
        included_levels_row = included_transcriptional_ranks,
        included_levels_col = included_colors, normalize_rows = F
    )

    f_wt <- rowSums(wt10_emb_vs_ct) > 0
    wt10_emb_vs_ct <- wt10_emb_vs_ct[f_wt, ]
    wt10_emb_vs_ct <- wt10_emb_vs_ct / rowSums(wt10_emb_vs_ct)

    return(list(wt = wt10_emb_vs_ct, tko = ko_emb_vs_ct, host = host_emb_vs_ct))
}

tetraploid_dotplot_frequencies <- function(mat_nm, minimal_number_of_cells = 20, downsample_number_of_cells = NULL, included_transcriptional_ranks = c(125:153)) {
    mat <- scdb_mat(mat_nm)
    # time window Et7.75 - Et8.1

    df_tetra <- read.table(sprintf("data/%s/time_match/time_match_summary.txt", mat_nm), sep = "\t", h = T, stringsAsFactors = F)
    rownames(df_tetra) <- df_tetra$embryo
    f <- (df_tetra[, 1] >= minimal_number_of_cells)
    df_tetra <- df_tetra[f, ]
    f <- df_tetra[, "best_query"] %in% included_transcriptional_ranks
    df_tetra <- df_tetra[f, ]
    included_chimeras <- df_tetra$embryo[]

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    # visceral and extraembryonic endoderm are excluded
    included_colors <- mc_wt@color_key$color[1:27]


    load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))

    # Modify blood subtypes color to one color
    query_cells_color <- aggregate_blood_subtypes_into_one_type(cmp_annot$query_cls_col)
    wt_cells_color <- mc_wt@colors[mc_wt@mc]
    names(wt_cells_color) <- names(mc_wt@mc)
    wt_cells_color <- aggregate_blood_subtypes_into_one_type(wt_cells_color)
    # remove Blood progenitors and Erythroid 2 from color levels
    included_colors <- included_colors[-c(21, 23)]

    query_cells_color <- query_cells_color[query_cells_color %in% included_colors]
    wt_cells_color <- wt_cells_color[wt_cells_color %in% included_colors]

    query_cells <- names(query_cells_color)[mat@cell_metadata[names(query_cells_color), "cell_type"] %in% c("KO", "control")]
    wt_cells <- names(wt_cells_color)

    # downsample cells to common number per embryo
    if (!is.null(downsample_number_of_cells)) {
        if (minimal_number_of_cells < downsample_number_of_cells) {
            stop("minimal_number_of_cells smaller than downsample_number_of_cells")
        }
        query_cells <- downsample_cells_indexed_by_metadata(cells = query_cells, cells_metadata = mat@cell_metadata[query_cells, "embryo"], n_downsample = minimal_number_of_cells, seed = 123)
        wt_cells <- downsample_cells_indexed_by_metadata(cells = wt_cells, cells_metadata = mat@cell_metadata[wt_cells, "transcriptional_rank"], n_downsample = minimal_number_of_cells, seed = 123)
    }


    # compute two way tables
    query_emb_vs_ct <- compute_two_way_table(
        values_row = mat@cell_metadata[query_cells, "embryo"],
        values_col = query_cells_color[query_cells],
        included_levels_row = included_chimeras,
        included_levels_col = included_colors, normalize_rows = T
    )

    wt10_emb_vs_ct <- compute_two_way_table(
        values_row = mat@cell_metadata[wt_cells, "transcriptional_rank"],
        values_col = wt_cells_color[wt_cells],
        included_levels_row = included_transcriptional_ranks,
        included_levels_col = included_colors, normalize_rows = F
    )

    f_wt <- rowSums(wt10_emb_vs_ct) > 0
    wt10_emb_vs_ct <- wt10_emb_vs_ct[f_wt, ]
    wt10_emb_vs_ct <- wt10_emb_vs_ct / rowSums(wt10_emb_vs_ct)


    return(list(wt = wt10_emb_vs_ct, query = query_emb_vs_ct))
}




compute_two_way_table <- function(values_row, values_col, included_levels_row = NULL, included_levels_col = NULL, normalize_rows = F) {
    if (length(values_row) != length(values_col)) {
        stop("values_row and values_col don't have the same length")
    }

    if (!is.null(included_levels_row)) {
        f <- values_row %in% included_levels_row
        values_row <- values_row[f]
        values_row <- factor(x = values_row, levels = included_levels_row)
        values_col <- values_col[f]
    }
    if (!is.null(included_levels_col)) {
        f <- values_col %in% included_levels_col
        values_row <- values_row[f]
        values_col <- values_col[f]
        values_col <- factor(x = values_col, levels = included_levels_col)
    }
    row_vs_col_freq <- table(values_row, values_col)

    if (normalize_rows) {
        row_vs_col_freq <- row_vs_col_freq / rowSums(row_vs_col_freq)
    }
    return(row_vs_col_freq)
}

plot_tetraploid_dotplots <- function(plot_pdf = T, included_transcriptional_ranks = NULL, highlighted_colors = NULL) {
    minimal_number_of_cells <- 250
    ko_color <- "indianred3"
    host_color <- "gray30"

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    if (is.null(included_transcriptional_ranks)) {
        included_transcriptional_ranks <- c(125:153)
    }
    if (is.null(highlighted_colors)) {
        highlighted_colors <- mc_wt@color_key$color[c(2, 3, 5, 6, 8, 12, 13, 14, 15, 17, 18, 19, 20, 22, 24, 27)]
    }

    tko_tetra_freq <- tetraploid_dotplot_frequencies(mat_nm = "tko_tetra_wt10", minimal_number_of_cells = 20, included_transcriptional_ranks = included_transcriptional_ranks)
    control_tetra_freq <- tetraploid_dotplot_frequencies(mat_nm = "control_tetra_all_wt10", minimal_number_of_cells = 20)

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    col_to_ct <- mc_wt@color_key$group
    col_to_ct[22] <- "Blood"
    names(col_to_ct) <- mc_wt@color_key$color


    ko_color <- "indianred3"
    host_color <- "gray30"

    fig_dir <- "figs/paper_figs/fig2"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    fig_dir <- "figs/paper_figs/fig2/cell_type_dot_plots_4N"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    tko_freq_n <- tko_tetra_freq$query
    ctrl_freq_n <- control_tetra_freq$query
    wt_freq_n <- tko_tetra_freq$wt

    tko_tetra_freq_ds <- tetraploid_dotplot_frequencies(mat_nm = "tko_tetra_wt10", minimal_number_of_cells = minimal_number_of_cells, downsample_number_of_cells = minimal_number_of_cells)
    control_tetra_freq_ds <- tetraploid_dotplot_frequencies(mat_nm = "control_tetra_all_wt10", minimal_number_of_cells = minimal_number_of_cells, downsample_number_of_cells = minimal_number_of_cells)

    tko_freq_ds <- tko_tetra_freq_ds$query
    ctrl_freq_ds <- control_tetra_freq_ds$query
    wt_freq_ds <- tko_tetra_freq_ds$wt


    tko_to_wt_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = tko_freq_ds[, ct_col], y = wt_freq_ds[, ct_col])
        return(p_val$p.value)
    })

    tko_to_ctrl_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = tko_freq_ds[, ct_col], y = ctrl_freq_ds[, ct_col])
        return(p_val$p.value)
    })

    ctrl_to_wt_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = ctrl_freq_ds[, ct_col], y = wt_freq_ds[, ct_col])
        return(p_val$p.value)
    })



    if (length(highlighted_colors) > 1) {
        q_val_tko <- qvalue(p = tko_to_wt_p_values, pi0 = 1)
        q_val_ctrl <- qvalue(p = ctrl_to_wt_p_values, pi0 = 1)
        q_val_tko_to_ctrl <- qvalue(p = tko_to_ctrl_p_values, pi0 = 1)
    } else {
        q_val_tko <- list(qvalues = tko_to_wt_p_values)
        q_val_ctrl <- list(qvalues = ctrl_to_wt_p_values)
        q_val_tko_to_ctrl <- list(qvalues = tko_to_ctrl_p_values)
    }



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

    genotype_color <- c("TKO" = ko_color, "Ctrl" = host_color, "WT" = "gray70")


    stat_comparison <- data.frame(
        group1 = c(rep("TKO", length(highlighted_colors)), rep("Ctrl", length(highlighted_colors)), rep("TKO", length(highlighted_colors))),
        group2 = c(rep("Ctrl", length(highlighted_colors)), rep("WT", length(highlighted_colors)), rep("WT", length(highlighted_colors))),
        cell_type = c(col_to_ct[names(q_val_tko_to_ctrl$qvalues)], col_to_ct[names(q_val_ctrl$qvalues)], col_to_ct[names(q_val_tko$qvalues)]),
        cell_type_color = c(names(q_val_tko_to_ctrl$qvalues), names(q_val_ctrl$qvalues), names(q_val_tko$qvalues)),
        q.val = c(q_val_tko_to_ctrl$qvalues, q_val_ctrl$qvalues, q_val_tko$qvalues),
        q.signif = c(q_to_signif(q_val_tko_to_ctrl$qvalues), q_to_signif(q_val_ctrl$qvalues), q_to_signif(q_val_tko$qvalues)), stringsAsFactors = F
    )

    plot_list <- list()

    col_to_ct[20] <- "Haematoendothelial prog."
    col_to_ct[14] <- "Later. & interm. mesoderm"

    for (ct_col in highlighted_colors) {
        main_tag <- gsub("/", "_", col_to_ct[ct_col])

        df_plot_points <- data.frame(
            genotype = factor(x = c(rep("TKO", nrow(tko_freq_n)), rep("Ctrl", nrow(ctrl_freq_n)), rep("WT", nrow(wt_freq_n))), levels = c("TKO", "Ctrl", "WT")),
            freq = c(tko_freq_n[, ct_col], ctrl_freq_n[, ct_col], wt_freq_n[, ct_col])
        )

        my_comparisons <- list(c("TKO", "Ctrl"), c("TKO", "WT"), c("Ctrl", "WT"))

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
            ggsave(filename = sprintf("%s/4N_%s.pdf", fig_dir, main_tag), width = 3, height = 2.3, plot = p)
        } else {
            ggsave(filename = sprintf("%s/4N_%s.png", fig_dir, main_tag), width = 3, height = 2.3, plot = p)
        }
    }


    p_all <- grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)

    if (plot_pdf) {
        ggsave(filename = sprintf("%s/4N_all_cell_types.pdf", fig_dir), width = 8.5, height = 6.5, plot = p_all)
    } else {
        ggsave(filename = sprintf("%s/4N_all_cell_types.png", fig_dir), width = 8.5, height = 6.5, plot = p_all)
    }
}

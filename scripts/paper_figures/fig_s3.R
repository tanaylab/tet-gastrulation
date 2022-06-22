source("scripts/paper_figures/fig2.R")
# fig s3

generate_figure_s3_plots <- function() {
    fig_dir <- "figs/paper_figs/fig_s3"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    # fig s3b
    ctrl_host_barplot_ct_frequency(ko_type = "host", plot_pdf = T, tag = "host")
    ctrl_host_barplot_ct_frequency(ko_type = "control", plot_pdf = T, tag = "control")

    # fig s3c
    fig_s3_ctrl_host_frequency_plots()

    # fig s3d
    fig_s3_ctrl_vs_host_timing()

    # fig s3e
    cell_type_comparisons_s3()

    # fig s3g
    dosage_effect_tko_in_chimeras_dotplots()

    # fig s3H
    fig_s3_chimera_dotplots_new_lineages(T)

    # fig s3i
    # see plots fig 2f

    # fig s3g
    # old dotplot fig2
}

ctrl_host_barplot_ct_frequency <- function(ko_type = "host", plot_pdf = F, tag = "host") {
    mat_nm <- "tko_chim_wt10"
    n_cls_min <- 19
    fig_dir <- "figs/paper_figs/fig_s3"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }


    df_chim <- read.table(sprintf("data/%s/time_match/time_match_summary.txt", mat_nm), sep = "\t", stringsAsFactors = F, h = T)
    rownames(df_chim) <- df_chim$embryo
    mat <- scdb_mat(mat_nm)
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")
    col_to_rank <- c(1:nrow(mc_wt@color_key))
    names(col_to_rank) <- mc_wt@color_key$color

    excluded_colors <- c("#F6BFCB", "#7F6874")
    included_colors <- setdiff(unique(mc_wt@color_key$color), excluded_colors)



    chim_embryos <- df_chim$embryo[(df_chim$host > n_cls_min) & (df_chim$control > n_cls_min)]
    chim_embryos <- chim_embryos[order(df_chim[chim_embryos, "best_rank_host_control"])]

    tmp <- matrix(0, nrow = length(chim_embryos), ncol = length(included_colors))
    rownames(tmp) <- chim_embryos
    colnames(tmp) <- included_colors

    load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))

    query_cls_col <- cmp_annot$query_cls_col
    query_cls <- names(query_cls_col)[!(query_cls_col %in% excluded_colors)]
    query_cls <- query_cls[mat@cell_metadata[query_cls, "embryo"] %in% chim_embryos]

    filtered_cls <- query_cls[mat@cell_metadata[query_cls, "cell_type"] %in% ko_type]
    filtered_vs_ct <- table(factor(x = mat@cell_metadata[filtered_cls, "embryo"], levels = chim_embryos), factor(x = query_cls_col[filtered_cls], levels = mc_wt@color_key$color[1:27]))
    filtered_vs_ct_n <- filtered_vs_ct / rowSums(filtered_vs_ct)
    filtered_vs_ct_n[is.na(filtered_vs_ct_n)] <- 0

    if (plot_pdf) {
        pdf(sprintf("%s/barplot_ct_freq_%s.pdf", fig_dir, tag), w = 12, h = 7.5, useDingbats = F)
        barplot(t(filtered_vs_ct_n), col = colnames(filtered_vs_ct_n), las = 2, axes = F, axisnames = F)
        dev.off()
    } else {
        png(sprintf("%s/barplot_ct_freq_%s.png", fig_dir, tag), w = 1250, h = 750)
        barplot(t(filtered_vs_ct_n), col = colnames(filtered_vs_ct_n), las = 2, axes = F, axisnames = F)
        dev.off()
    }
}


cell_type_comparisons_s3 <- function() {
    tko_chim <- tko_bulk_comparison_s3(
        mat_id = "tko_chim_wt10",
        included_cell_types = c(
            "Epiblast", "Early nascent mesoderm", "Amnion/Chorion", "ExE mesoderm",
            "Primitive streak"
        ),
        genotype = "KO",
        xlab_plot = "TKO",
        tag_fn = "tko_2N"
    )

    ctrl_chim <- tko_bulk_comparison_s3(
        mat_id = "tko_chim_wt10",
        included_cell_types = c(
            "Epiblast", "Early nascent mesoderm", "Amnion/Chorion", "ExE mesoderm",
            "Primitive streak"
        ),
        genotype = c("control", "host"),
        xlab_plot = "Ctrl/Host",
        tag_fn = "ctrl_host_2N"
    )


    a <- arrangeGrob(tko_chim[[1]], ctrl_chim[[1]],
        tko_chim[[5]], ctrl_chim[[5]],
        tko_chim[[2]], ctrl_chim[[2]],
        tko_chim[[3]], ctrl_chim[[3]],
        tko_chim[[4]], ctrl_chim[[4]],
        ncol = 2, nrow = 5
    )


    a <- arrangeGrob(tko_chim[[1]], tko_chim[[5]], tko_chim[[2]], tko_chim[[3]], tko_chim[[4]],
        ctrl_chim[[1]], ctrl_chim[[5]], ctrl_chim[[2]], ctrl_chim[[3]], ctrl_chim[[4]],
        ncol = 5, nrow = 2
    )

    ggsave(filename = "figs/paper_figs/fig_s3/tko_chim_vs_wt_per_ct.png", plot = a, w = 24, h = 10)
    ggsave(filename = "figs/paper_figs/fig_s3/tko_chim_vs_wt_per_ct.pdf", plot = a, w = 24, h = 10)
}


tko_bulk_comparison_s3 <- function(mat_id, included_cell_types, genotype, plot_pdf = F, xlab_plot = "", tag_fn = "tko_4N") {
    reg <- 5e-5
    mat <- scdb_mat(mat_id)

    f_tko <- mat@cell_metadata[colnames(mat@mat), "cell_type"] %in% genotype

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")
    mat_wt <- scdb_mat("sing_emb_wt10")

    egc <- t(tgs_matrix_tapply(x = mat_wt@mat[, names(mc_wt@mc)], mc_wt@mc, sum))
    egc <- t(t(egc) / colSums(egc))

    ct_to_col <- mc_wt@color_key$color
    names(ct_to_col) <- mc_wt@color_key$group
    col_to_ct <- mc_wt@color_key$group
    names(col_to_ct) <- mc_wt@color_key$color

    legc <- log2(egc + 1e-5)

    gset <- scdb_gset("sing_emb_wt10")
    bad_genes <- read.table("data/tet_tko.bad_genes.txt", sep = "\t", stringsAsFactors = F)$x
    bad_genes <- c("Igf2", "AK145379;H19", "Tet1", "Tet2", "Tet3")
    feat_genes <- setdiff(names(gset@gene_set), bad_genes)

    sc_mc_cor <- tgs_cor(as.matrix(mat@mat[feat_genes, f_tko]), legc[feat_genes, ])

    best_ref <- apply(sc_mc_cor, 1, which.max)

    umi_query <- t(tgs_matrix_tapply(mat@mat[, f_tko], col_to_ct[mc_wt@colors[best_ref]], sum))

    tot_umis_query <- colSums(umi_query)
    egc_query <- t(t(umi_query) / colSums(umi_query))
    legc_query <- log2(egc_query + reg)

    egc_ref <- t(tgs_matrix_tapply(egc[, best_ref], col_to_ct[mc_wt@colors[best_ref]], sum))
    egc_ref <- t(t(egc_ref) / colSums(egc_ref))
    legc_ref <- log2(egc_ref + reg)

    bad_genes <- read.table("data/tet_tko.bad_genes.txt", sep = "\t", stringsAsFactors = F)$x
    bad_genes <- c(bad_genes, c("Tet1", "Tet2", "Tet3"))
    genes_f <- setdiff(rownames(mat@mat), bad_genes)

    filtered_genes <- rownames(egc_query)[pmax(apply(egc_query, 1, max), apply(egc_ref, 1, max)) > 1e-5]

    lfc_query_ref <- legc_query[filtered_genes, ] - legc_ref[filtered_genes, ]


    qvalue_list <- list()

    for (ct in included_cell_types) {
        p_val_vector <- sapply(filtered_genes, function(gene) {
            test_result <- chisq.test(x = c(umi_query[gene, ct], tot_umis_query[ct] - umi_query[gene, ct]), p = c(egc_ref[gene, ct], 1 - egc_ref[gene, ct]))
            return(test_result$p.value)
        })

        p_val_vector[is.na(p_val_vector)] <- 1

        q_value_vector <- qvalue(p_val_vector, pi0 = 1)$qvalues

        qvalue_list[[ct]] <- q_value_vector
    }

    plot_ls <- list()
    for (ct in included_cell_types) {
        df_plot <- data.frame(
            gene = filtered_genes,
            expression_tko = legc_query[filtered_genes, ct],
            expression_wt = legc_ref[filtered_genes, ct]
        )

        qvalues_plot <- qvalue_list[[ct]]

        gene_color <- ifelse((qvalues_plot >= 1e-3) | (abs(lfc_query_ref[, ct]) < log2(1.5)), ct_to_col[ct], "gray30")

        f <- rank(-abs(legc_query[filtered_genes, ct] - legc_ref[filtered_genes, ct])) < 4
        p <- ggplot(df_plot, aes(y = expression_wt, x = expression_tko)) +
            geom_abline(slope = 1, intercept = 0, color = "gray") +
            geom_abline(slope = 1, intercept = -1, color = "gray", linetype = "dashed") +
            geom_abline(slope = 1, intercept = 1, color = "gray", linetype = "dashed") +
            geom_point(color = gene_color, size = 1) +
            ylab("WT") +
            xlab(xlab_plot) +
            theme(axis.title = element_text(size = 20)) #+
        # ggtitle(label = ct) + theme(plot.title = element_text(hjust = 0.5))
        plot_ls[[ct]] <- p
    }





    return(plot_ls)
}





fig_s3_ctrl_host_frequency_plots <- function() {
    fn_prefix <- ""
    plot_wt <- F
    plot_ctrl <- T
    plot_host <- T
    plot_pdf <- T

    ctrl_color <- "cornflowerblue"
    host_color <- "gray30"

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    s3_chimera_plot_frequency_group_of_cell_types(
        highlighted_colors = mc_wt@color_key$color[2:6], main_tag = paste0(fn_prefix, "Ectoderm"), ctrl_color = ctrl_color, host_color = host_color,
        plot_pdf = plot_pdf, plot_wt = plot_wt, plot_host = plot_host, plot_ctrl = plot_ctrl
    )
    s3_chimera_plot_frequency_group_of_cell_types(
        highlighted_colors = mc_wt@color_key$color[26:27], main_tag = paste0(fn_prefix, "Endoderm"), ctrl_color = ctrl_color, host_color = host_color,
        plot_pdf = plot_pdf, plot_wt = plot_wt, plot_host = plot_host, plot_ctrl = plot_ctrl
    )
    s3_chimera_plot_frequency_group_of_cell_types(
        highlighted_colors = mc_wt@color_key$color[c(8, 12, 13, 14, 15, 16)],
        main_tag = paste0(fn_prefix, "Embryonic mesoderm"), ctrl_color = ctrl_color, host_color = host_color,
        plot_pdf = plot_pdf, plot_wt = plot_wt, plot_host = plot_host, plot_ctrl = plot_ctrl
    )
    s3_chimera_plot_frequency_group_of_cell_types(
        highlighted_colors = mc_wt@color_key$color[c(17, 18, 19)],
        main_tag = paste0(fn_prefix, "Extraembryonic mesoderm"), ctrl_color = ctrl_color, host_color = host_color,
        plot_pdf = plot_pdf, plot_wt = plot_wt, plot_host = plot_host, plot_ctrl = plot_ctrl
    )
    s3_chimera_plot_frequency_group_of_cell_types(
        highlighted_colors = mc_wt@color_key$color[c(21:23)],
        main_tag = paste0(fn_prefix, "Blood"), ctrl_color = ctrl_color, host_color = host_color,
        plot_pdf = plot_pdf, plot_wt = plot_wt, plot_host = plot_host, plot_ctrl = plot_ctrl
    )
}

s3_chimera_plot_frequency_group_of_cell_types <- function(highlighted_colors, main_tag = "", plot_pdf = F, ctrl_color = "cornflowerblue", host_color = "gray30", cex_points = 3, plot_wt = F, plot_ctrl = T, plot_host = F) {
    mat_nm <- "tko_chim_wt10"

    rank_to_time <- read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt", stringsAsFactors = F, h = T, sep = "\t")
    dev_time <- rank_to_time$developmental_time
    age_field_ctrl <- "best_rank_host_control"
    age_field_host <- "best_rank_host_control"

    if (!dir.exists("figs/paper_figs")) {
        dir.create("figs/paper_figs")
    }

    fig_dir <- "figs/paper_figs/fig_s3"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    fig_dir <- "figs/paper_figs/fig_s3/lineage_frequency_control_host"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    roll_width <- 5

    mat_chim <- scdb_mat(mat_nm)

    load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))
    ctrl_query_cls_col <- cmp_annot$query_cls_col[mat_chim@cell_metadata[names(cmp_annot$query_cls_col), "cell_type"] == "control"]
    ctrl_query_cls <- names(ctrl_query_cls_col)
    host_query_cls_col <- cmp_annot$query_cls_col[mat_chim@cell_metadata[names(cmp_annot$query_cls_col), "cell_type"] %in% c("host")]
    host_query_cls <- names(host_query_cls_col)

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")
    col_to_ct <- mc_wt@color_key$group
    names(col_to_ct) <- mc_wt@color_key$color
    excluded_colors <- c("#F6BFCB", "#7F6874")
    included_colors <- setdiff(mc_wt@color_key$color, excluded_colors)

    df_chim <- read.table(sprintf("data/%s/time_match/time_match_summary.txt", mat_nm), sep = "\t", h = T, stringsAsFactors = F)
    rownames(df_chim) <- df_chim$embryo

    f <- (df_chim$host > 19) & (df_chim$control > 19)
    df_chim <- df_chim[f, ]

    embryos <- df_chim$embryo
    embryos <- embryos[order(df_chim[embryos, "best_rank_host_control"])]

    ctrl_query_cls_f <- ctrl_query_cls[(mat_chim@cell_metadata[ctrl_query_cls, "cell_type"] %in% c("control")) & (ctrl_query_cls_col[ctrl_query_cls] %in% included_colors)]
    host_query_cls_f <- host_query_cls[(mat_chim@cell_metadata[host_query_cls, "cell_type"] %in% c("host")) & (host_query_cls_col[host_query_cls] %in% included_colors)]
    wt10_cls <- names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] %in% included_colors]

    emb_wt_age <- unique(mat_chim@cell_metadata[wt10_cls, c("transcriptional_rank", "age_group")])
    emb_wt_age <- emb_wt_age[order(emb_wt_age$transcriptional_rank), ]

    wt10_emb_vs_ct <- table(mat_chim@cell_metadata[wt10_cls, "transcriptional_rank"], mc_wt@colors[mc_wt@mc[wt10_cls]])
    ctrl_emb_vs_ct <- table(mat_chim@cell_metadata[ctrl_query_cls_f, "embryo"], ctrl_query_cls_col[ctrl_query_cls_f])
    ctrl_emb_vs_ct <- ctrl_emb_vs_ct[embryos, ]
    host_emb_vs_ct <- table(mat_chim@cell_metadata[host_query_cls_f, "embryo"], host_query_cls_col[host_query_cls_f])
    host_emb_vs_ct <- host_emb_vs_ct[embryos, ]

    if (length(highlighted_colors) > 1) {
        wt_fr <- rowSums(wt10_emb_vs_ct[, intersect(highlighted_colors, colnames(wt10_emb_vs_ct))]) / rowSums(wt10_emb_vs_ct)
        ctrl_fr <- rowSums(ctrl_emb_vs_ct[, intersect(highlighted_colors, colnames(ctrl_emb_vs_ct))]) / rowSums(ctrl_emb_vs_ct)
        host_fr <- rowSums(host_emb_vs_ct[, intersect(highlighted_colors, colnames(host_emb_vs_ct))]) / rowSums(host_emb_vs_ct)
    } else {
        wt_fr <- wt10_emb_vs_ct[, intersect(highlighted_colors, colnames(wt10_emb_vs_ct))] / rowSums(wt10_emb_vs_ct)
        ctrl_fr <- ctrl_emb_vs_ct[, intersect(highlighted_colors, colnames(ctrl_emb_vs_ct))] / rowSums(ctrl_emb_vs_ct)
        host_fr <- host_emb_vs_ct[, intersect(highlighted_colors, colnames(host_emb_vs_ct))] / rowSums(host_emb_vs_ct)
    }

    # next calculate moving  90% moving average window for every cell type
    n <- length(wt_fr)
    freq_ct <- c(wt_fr, wt_fr[(n - roll_width + 1):n])
    mov_mean <- rollmean(x = freq_ct, k = 2 * roll_width + 1)

    mov_sd <- rollapply(data = freq_ct, width = 2 * roll_width + 1, sd)

    upper_sd <- mov_mean + 2 * mov_sd
    lower_sd <- pmax(mov_mean - 2 * mov_sd, 0)

    names(upper_sd) <- c((roll_width + 1):153)
    names(lower_sd) <- c((roll_width + 1):153)
    names(mov_mean) <- c((roll_width + 1):153)

    # min_rank_plot = min(df_chim[tetra_embryos,c("best_query")]) - 5
    min_rank_plot <- 10
    xlim_min <- dev_time[min_rank_plot]
    xlim_max <- dev_time[153] + 0.05

    # x_ranks = c(min_rank_plot:(153-roll_width))
    x_ranks <- c(min_rank_plot:(153))
    upper_freq <- upper_sd[as.character(c(min_rank_plot:(153)))]
    lower_freq <- lower_sd[as.character(c(min_rank_plot:(153)))]


    ylim_max <- max(ctrl_fr, host_fr, wt_fr)

    # next calculate p value between query and reference

    min_age_rank_test <- 109

    f_host <- df_chim[names(host_fr), age_field_host] >= min_age_rank_test
    f_ctrl <- df_chim[names(ctrl_fr), age_field_ctrl] >= min_age_rank_test
    w_test <- wilcox.test(x = ctrl_fr[f_ctrl], y = host_fr[f_host])


    if (plot_pdf) {
        pdf(sprintf("%s/%s.pdf", fig_dir, main_tag), w = 7, h = 5.6, useDingbats = F)
        par(mar = c(6, 6, 6, 4))

        plot(NULL,
            ylim = c(-0.05 * ylim_max, 1.2 * ylim_max),
            main = main_tag, ylab = "Fraction", xlab = "Developmental Time", cex.main = 3,
            cex.lab = 2, cex.axis = 2, xaxs = "i", yaxs = "i", xlim = c(xlim_min, xlim_max)
        )
        # plot(x = dev_time[df_chim[names(ko_fr),age_field_ko]],y = ko_fr,ylim = c(-0.05*ylim_max,1.2*ylim_max),
        #     pch = 19,cex = 0.3,main = main_tag,ylab = "Fraction",xlab = "Developmental Time",cex.main = 3,
        #     cex.lab = 2,cex.axis = 2,xaxs = 'i',yaxs= 'i',xlim = c(xlim_min,xlim_max))




        polygon(x = c(dev_time[x_ranks], dev_time[rev(x_ranks)]), y = c(upper_freq, rev(lower_freq)), col = "gray80", border = NA)
        lines(x = dev_time[c((roll_width + 1):153)], y = mov_mean, lwd = 3)
        if (plot_wt) {
            points(x = dev_time[as.numeric(names(wt_fr))], y = wt_fr, pch = 19, cex = 0.3)
        }
        if (plot_host) {
            points(x = dev_time[df_chim[names(host_fr), age_field_host]], y = host_fr, pch = 19, cex = cex_points, col = host_color)
        }
        if (plot_ctrl) {
            points(x = dev_time[df_chim[names(ctrl_fr), age_field_ctrl]], y = ctrl_fr, pch = 19, cex = cex_points, col = ctrl_color)
        }

        text(x = 7, y = 0.9 * max(host_fr, ctrl_fr, wt_fr), labels = sprintf("p = %.1e", w_test$p.value), cex = 2)
        # legend(x = "topleft",legend = c("wt","ko","control),pch = 19,col = c("black","#83c26d","cornflowerblue"))

        dev.off()
    } else {
        png(sprintf("%s/%s.png", fig_dir, main_tag), w = 600, h = 450)
        par(mar = c(6, 6, 6, 4))
        plot(NULL,
            ylim = c(-0.05 * ylim_max, 1.2 * ylim_max),
            main = main_tag, ylab = "Fraction", xlab = "Developmental Time", cex.main = 3,
            cex.lab = 2, cex.axis = 2, xaxs = "i", yaxs = "i", xlim = c(xlim_min, xlim_max)
        )
        # plot(x = dev_time[df_chim[names(ko_fr),age_field_ko]],y = ko_fr,ylim = c(-0.05*ylim_max,1.2*ylim_max),
        #     pch = 19,cex = 0.3,main = main_tag,ylab = "Fraction",xlab = "Developmental Time",cex.main = 2,
        #     cex.lab = 2,xaxs = 'i',yaxs= 'i',xlim = c(xlim_min,xlim_max))


        polygon(x = c(dev_time[x_ranks], dev_time[rev(x_ranks)]), y = c(upper_freq, rev(lower_freq)), col = "gray80", border = NA)
        lines(x = dev_time[c((roll_width + 1):153)], y = mov_mean, lwd = 3)
        if (plot_wt) {
            points(x = dev_time[as.numeric(names(wt_fr))], y = wt_fr, pch = 19, cex = 0.3)
        }
        if (plot_ctrl) {
            points(x = dev_time[df_chim[names(ctrl_fr), age_field_ctrl]], y = ctrl_fr, pch = 19, cex = cex_points, col = ctrl_color)
        }
        if (plot_host) {
            points(x = dev_time[df_chim[names(host_fr), age_field_host]], y = host_fr, pch = 19, cex = cex_points, col = host_color)
        }
        text(x = 7, y = 0.9 * max(host_fr, ctrl_fr, wt_fr), labels = sprintf("p = %.1e", w_test$p.value), cex = 2)
        # legend(x = "topleft",legend = c("wt","ko","control"),pch = 19,col = c("black","#83c26d","cornflowerblue"))

        dev.off()
    }
}



fig_s3_ctrl_vs_host_timing <- function() {
    rank_to_time <- read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt", stringsAsFactors = F, h = T, sep = "\t")
    dev_time <- rank_to_time$developmental_time

    df_chim <- read.table(file = "data/tko_chim_wt10/time_match/time_match_summary.txt", sep = "\t", stringsAsFactors = F, h = T)

    f <- df_chim$host > 19 & df_chim$control > 19

    df_plot <- df_chim[f, c("best_rank_host", "best_rank_control")]
    df_plot$best_rank_host <- dev_time[df_plot$best_rank_host]
    df_plot$best_rank_control <- dev_time[df_plot$best_rank_control]

    p <- ggplot(data = df_plot, aes(x = best_rank_host, y = best_rank_control)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_point(size = 4) +
        xlab(expression(Host ~ E[t])) +
        ylab(expression(Ctrl ~ E[t])) +
        theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20)) +
        xlim(7, 8.1) +
        ylim(7, 8.1)

    ggsave(filename = "figs/paper_figs/fig_s3/timing_control_vs_host.pdf", plot = p)
}


dosage_effect_tko_in_chimeras_dotplots <- function(plot_pdf = F) {
    mat <- scdb_mat("tko_chim_wt10")
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    load("data/tko_chim_wt10/color_annotation/cmp_annot.Rda")
    wt10_age <- read.table("data/wt10_transcriptional_rank_developmental_time.txt", sep = "\t", h = T)

    tko_developmental_time <- read.table("data/tko_chim_wt10/time_match/time_match_summary.txt", sep = "\t", h = T)
    chimera_transcriptional_rank <- ifelse(is.na(tko_developmental_time$best_rank_host_control), tko_developmental_time$best_rank_ko, tko_developmental_time$best_rank_host_control)
    names(chimera_transcriptional_rank) <- tko_developmental_time$embryo

    embryo_df <- read.csv("data/KO fraction per chimera for assessment of dosage effects.csv", stringsAsFactors = F)
    f <- embryo_df$embryo.ID %in% tko_developmental_time$embryo
    embryo_df <- embryo_df[f, ]
    embryos_f <- embryo_df$embryo.ID[chimera_transcriptional_rank[embryo_df$embryo.ID] > 124]

    cells_filtered <- names(cmp_annot$query_cls_col)
    cells_filtered <- cells_filtered[mat@cell_metadata[cells_filtered, "embryo"] %in% embryos_f]
    cells_filtered <- cells_filtered[cmp_annot$query_cls_col[cells_filtered] %in% mc_wt@color_key$color[1:27]]

    tko_cells <- cells_filtered[mat@cell_metadata[cells_filtered, "cell_type"] == "KO"]
    host_cells <- cells_filtered[mat@cell_metadata[cells_filtered, "cell_type"] == "host"]
    control_cells <- cells_filtered[mat@cell_metadata[cells_filtered, "cell_type"] == "control"]

    embryos_f <- intersect(
        unique(mat@cell_metadata[c(host_cells, control_cells), "embryo"]),
        unique(mat@cell_metadata[c(tko_cells), "embryo"])
    )

    embryo_genotype_fraction_df <- as.matrix(embryo_df[, c("TKO", "Control", "Host")])
    rownames(embryo_genotype_fraction_df) <- embryo_df$embryo.ID
    embryo_genotype_fraction_df <- embryo_genotype_fraction_df[embryos_f, ]

    ko_cell_type_fraction <- table(
        factor(x = mat@cell_metadata[tko_cells, "embryo"], levels = embryos_f),
        factor(cmp_annot$query_cls_col[tko_cells], levels = mc_wt@color_key$color[1:27])
    )
    host_cell_type_fraction <- table(
        factor(x = mat@cell_metadata[host_cells, "embryo"], levels = embryos_f),
        factor(cmp_annot$query_cls_col[host_cells], levels = mc_wt@color_key$color[1:27])
    )
    control_cell_type_fraction <- table(
        factor(x = mat@cell_metadata[control_cells, "embryo"], levels = embryos_f),
        factor(cmp_annot$query_cls_col[control_cells], levels = mc_wt@color_key$color[1:27])
    )

    host_control_cell_type_fraction <- host_cell_type_fraction + control_cell_type_fraction
    host_control_cell_type_fraction <- host_control_cell_type_fraction / rowSums(host_control_cell_type_fraction)
    host_cell_type_fraction <- host_cell_type_fraction / rowSums(host_cell_type_fraction)
    ko_cell_type_fraction <- ko_cell_type_fraction / rowSums(ko_cell_type_fraction)
    control_cell_type_fraction <- control_cell_type_fraction / rowSums(control_cell_type_fraction)

    germ_layer_groups <- read.csv("data/germ_layer_groups.csv", stringsAsFactors = F)

    ko_group_fraction <- t(tgs_matrix_tapply(x = as.matrix(as.data.frame.matrix(ko_cell_type_fraction)), germ_layer_groups$group[1:27], sum))
    host_group_fraction <- t(tgs_matrix_tapply(x = as.matrix(as.data.frame.matrix(host_cell_type_fraction)), germ_layer_groups$group[1:27], sum))
    control_group_fraction <- t(tgs_matrix_tapply(x = as.matrix(as.data.frame.matrix(control_cell_type_fraction)), germ_layer_groups$group[1:27], sum))
    host_control_group_fraction <- t(tgs_matrix_tapply(x = as.matrix(as.data.frame.matrix(host_control_cell_type_fraction)), germ_layer_groups$group[1:27], sum))

    chimeras_with_high_tko_content <- rownames(embryo_genotype_fraction_df)[embryo_genotype_fraction_df[, "TKO"] > 0.1]
    chimeras_with_low_tko_content <- rownames(embryo_genotype_fraction_df)[embryo_genotype_fraction_df[, "TKO"] <= 0.1]


    all_groups <- c("Embryonic ectoderm", "Embryonic mesoderm", "Non-embryonic mesoderm", "Blood", "Embryonic endoderm")


    tko_low_to_high_p_values <- sapply(all_groups, function(ct_group) {
        p_val <- wilcox.test(
            x = ko_group_fraction[chimeras_with_high_tko_content, ct_group],
            y = ko_group_fraction[chimeras_with_low_tko_content, ct_group]
        )
        return(p_val$p.value)
    })

    host_control_low_to_high_p_values <- sapply(all_groups, function(ct_group) {
        p_val <- wilcox.test(
            x = host_control_group_fraction[chimeras_with_high_tko_content, ct_group],
            y = host_control_group_fraction[chimeras_with_low_tko_content, ct_group]
        )
        return(p_val$p.value)
    })


    q_val_tko <- qvalue(p = tko_low_to_high_p_values, pi0 = 1)
    q_val_host_control <- qvalue(p = host_control_low_to_high_p_values, pi0 = 1)

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




    ko_color <- "indianred3"
    host_color <- "gray30"

    fig_dir <- "figs/paper_figs/fig_s3/dosage_effect"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    genotype_color <- c("TKO" = ko_color, "Host/Control" = host_color, "WT" = "gray70")



    stat_comparison <- data.frame(
        group1 = c(rep("TKO Low", length(all_groups)), rep("Host/Control Low", length(all_groups))),
        group2 = c(rep("TKO High", length(all_groups)), rep("Host/Control High", length(all_groups))),
        cell_type_group = rep(all_groups, 2),
        q.val = c(q_val_tko$qvalues, q_val_host_control$qvalues),
        q.signif = c(q_to_signif(q_val_tko$qvalues), q_to_signif(q_val_host_control$qvalues)), stringsAsFactors = F
    )

    plot_list <- list()


    for (ct_group in all_groups) {
        main_tag <- ct_group


        df_plot_points <- data.frame(
            genotype = factor(x = c(
                rep("TKO", length(embryos_f)),
                rep("Host/Control", length(embryos_f))
            ), levels = c("TKO", "Host/Control")),
            genotype_dosage = factor(
                x = c(
                    rep("TKO Low", length(chimeras_with_low_tko_content)), rep("TKO High", length(chimeras_with_high_tko_content)),
                    rep("Host/Control Low", length(chimeras_with_low_tko_content)), rep("Host/Control High", length(chimeras_with_high_tko_content))
                ),
                levels = c("TKO Low", "TKO High", "Host/Control Low", "Host/Control High")
            ),
            freq = c(
                ko_group_fraction[chimeras_with_low_tko_content, ct_group],
                ko_group_fraction[chimeras_with_high_tko_content, ct_group],
                host_control_group_fraction[chimeras_with_low_tko_content, ct_group],
                host_control_group_fraction[chimeras_with_high_tko_content, ct_group]
            )
        )


        my_comparisons <- list(c("TKO Low", "TKO High"), c("Host/Control Low", "Host/Control High"))

        stat.test <- compare_means(data = df_plot_points, formula = freq ~ genotype_dosage)

        stat_f <- stat_comparison[stat_comparison$cell_type_group == ct_group, ]


        p <- ggplot(data = df_plot_points, aes(x = genotype_dosage, y = freq)) +
            geom_dotplot(aes(fill = genotype), dotsize = 2.3, binaxis = "y", stackdir = "center", show.legend = F) +
            stat_pvalue_manual(stat_f, y.position = max(df_plot_points$freq) * 1.1, step.increase = 0.1, label = "q.signif") +
            scale_fill_manual(values = genotype_color) +
            scale_x_discrete(labels = c("TKO \n Low", "TKO \n High", "TKO \n Low", "TKO \n High")) +
            ggtitle(label = main_tag) +
            theme(plot.title = element_text(hjust = 0.5, size = 10)) +
            ylab("") +
            ylim(0, max(df_plot_points$freq) * 1.4) +
            xlab("")
        # theme(axis.text.x = element_text(size=14))
        # stat_compare_means(label = "p.signif",comparisons = my_comparisons) +

        plot_list[[ct_group]] <- p
    }





    p <- ggplot(data = df_plot_points, aes(x = genotype_dosage, y = freq)) +
        geom_dotplot(aes(fill = genotype), dotsize = 2.3, binaxis = "y", stackdir = "center", show.legend = T) +
        stat_pvalue_manual(stat_f, y.position = max(df_plot_points$freq) * 1.1, step.increase = 0.1, label = "q.signif") +
        scale_fill_manual(name = "Genotype", values = genotype_color) +
        scale_x_discrete(labels = c("TKO \n Low", "TKO \n High", "TKO \n Low", "TKO \n High")) +
        ggtitle(label = main_tag) +
        theme(plot.title = element_text(hjust = 0.5, size = 10)) +
        ylab("") +
        ylim(0, max(df_plot_points$freq) * 1.4) +
        xlab("")

    p_legend <- cowplot::get_legend(p)

    plot_list[["legend"]] <- p_legend

    p_all <- grid.arrange(grobs = plot_list, ncol = 3, nrow = 2)

    if (plot_pdf) {
        ggsave(filename = sprintf("%s/all_cell_type_groups.pdf", fig_dir), width = 6.5, height = 4, plot = p_all)
    } else {
        ggsave(filename = sprintf("%s/all_cell_type_group.png", fig_dir), width = 6.5, height = 4, plot = p_all)
    }
}






fig_s3_chimera_dotplots_new_lineages <- function(plot_pdf = T) {
    ko_color <- "indianred3"
    host_color <- "gray30"
    chim_freq <- chimera_dotplot_frequencies_lineages()

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")


    ko_color <- "indianred3"
    host_color <- "gray30"

    mat_nm <- "tko_chim_wt10"

    highlighted_colors <- c("Endoderm", "Ectoderm", "Embryonic mesoderm", "Extraembryonic mesoderm", "Blood")

    if (!dir.exists("figs/paper_figs")) {
        dir.create("figs/paper_figs")
    }

    fig_dir <- "figs/paper_figs/fig_s3"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    fig_dir <- "figs/paper_figs/fig_s3/chimera_dot_plots_per_lineage"
    if (!dir.exists(fig_dir)) {
        dir.create(fig_dir)
    }

    tko_freq <- chim_freq$tko
    tko_freq_n <- tko_freq / rowSums(tko_freq)
    host_freq <- chim_freq$host
    host_freq_n <- host_freq / rowSums(host_freq)
    wt_freq <- chim_freq$wt
    wt_freq_n <- wt_freq / rowSums(wt_freq)



    tko_to_wt_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = tko_freq_n[, ct_col], y = wt_freq_n[, ct_col])
        return(p_val$p.value)
    })

    tko_to_host_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = tko_freq_n[, ct_col], y = host_freq_n[, ct_col])
        return(p_val$p.value)
    })

    host_to_wt_p_values <- sapply(highlighted_colors, function(ct_col) {
        p_val <- wilcox.test(x = host_freq_n[, ct_col], y = wt_freq_n[, ct_col])
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

    genotype_color <- c("TKO" = ko_color, "Host/Control" = host_color, "WT" = "gray70")


    stat_comparison <- data.frame(
        group1 = c(rep("TKO", length(highlighted_colors)), rep("Host/Control", length(highlighted_colors)), rep("TKO", length(highlighted_colors))),
        group2 = c(rep("Host/Control", length(highlighted_colors)), rep("WT", length(highlighted_colors)), rep("WT", length(highlighted_colors))),
        cell_type = c(names(q_val_tko_to_host$qvalues), names(q_val_host$qvalues), names(q_val_tko$qvalues)),
        q.val = c(q_val_tko_to_host$qvalues, q_val_host$qvalues, q_val_tko$qvalues),
        q.signif = c(q_to_signif(q_val_tko_to_host$qvalues), q_to_signif(q_val_host$qvalues), q_to_signif(q_val_tko$qvalues)), stringsAsFactors = F
    )

    for (ct_col in highlighted_colors) {
        main_tag <- ct_col

        df_plot_points <- data.frame(
            genotype = factor(x = c(rep("TKO", nrow(tko_freq_n)), rep("Host/Control", nrow(host_freq_n)), rep("WT", nrow(wt_freq_n))), levels = c("TKO", "Host/Control", "WT")),
            freq = c(tko_freq_n[, ct_col], host_freq_n[, ct_col], wt_freq_n[, ct_col])
        )

        my_comparisons <- list(c("TKO", "Host/Control"), c("TKO", "WT"), c("Host/Control", "WT"))

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


chimera_dotplot_frequencies_lineages <- function() {


    # time window Et7.75 - Et8.1
    t_f <- c(125:153)
    # time window Et7.5 - Et8.1
    # t_f = c(109:153)


    mat_nm <- "tko_chim_wt10"

    age_field_ko <- "best_rank_host_control"
    age_field_host <- "best_rank_host_control"

    mat_chim <- scdb_mat(mat_nm)

    load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda", mat_nm))
    ko_query_cls_col <- cmp_annot$query_cls_col[mat_chim@cell_metadata[names(cmp_annot$query_cls_col), "cell_type"] == "KO"]
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
    ct_to_lineage[2:6] <- "Ectoderm"
    ct_to_lineage[26:27] <- "Endoderm"
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

    f <- (df_chim$control + df_chim$host > 20) & (df_chim$KO > 20)
    df_chim <- df_chim[f, ]

    embryos <- df_chim$embryo
    embryos <- embryos[order(df_chim[embryos, "best_rank_host_control"])]

    f_t_q <- df_chim[embryos, "best_rank_host_control"] %in% t_f

    ko_query_cls_f <- ko_query_cls[(mat_chim@cell_metadata[ko_query_cls, "cell_type"] %in% c("KO")) & (ko_query_cls_col[ko_query_cls] %in% included_colors)]
    host_query_cls_f <- host_query_cls[(mat_chim@cell_metadata[host_query_cls, "cell_type"] %in% c("control", "host")) & (host_query_cls_col[host_query_cls] %in% included_colors)]
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


    return(list(wt = wt10_emb_vs_ct, tko = ko_emb_vs_ct, host = host_emb_vs_ct))
}

fig_s3_dotplots_definitive_endoderm <- function() {
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    # Time range E7.25 to E7.75
    plot_chimera_dotplots(plot_pdf = T, included_transcriptional_ranks = c(93:124), highlighted_colors = mc_wt@color_key$color[26], minimal_number_of_cells_for_p_value = 90)
    plot_tetraploid_dotplots(plot_pdf = T, included_transcriptional_ranks = c(93:124), highlighted_colors = mc_wt@color_key$color[26])
}

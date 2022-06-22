source("scripts/paper_figures/fig3.R")

# Fig s4

generate_figure_s4_plots <- function() {
    if (!dir.exists("figs/paper_figs/fig_s4")) {
        dir.create("figs/paper_figs/fig_s4")
    }

    # fig s4a
    epiblast_score_genes_correlation_heatmap()

    # fig s4b,c
    plot_epiblast_gene_expression_vs_time(genes_f = c("Pou3f1", "Sox11", "Tdgf1", "Nodal"), plot_pdf = T, fig_dir = "figs/paper_figs/fig_s4/epiblast_expression_vs_time")

    # fig s4d
    nascent_mesoderm_score_correlation_heatmap()


    # fig s4e,f
    plot_nascent_mesoderm_gene_expression_vs_time(genes_f = c("Fgf3", "Jag1", "Hand2", "Pcdh8", "Pcdh19", "Pitx2", "Cfc1"), plot_pdf = T, fig_dir = "figs/paper_figs/fig_s4/nascent_mesoderm_expression_vs_time")

    plot_heatmaps_epiblast_and_nascent_mesoderm_score_genes()
}


epiblast_score_genes_correlation_heatmap <- function() {
    mc <- scdb_mc("sing_emb_wt10_recolored")

    legc <- log2(mc@e_gc + 1e-5)

    gg_cor <- tgs_cor(t(legc[epiblast_genes(), ]))

    breaks <- seq(0.5, 1, length.out = 101)
    # shades = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(100))
    shades <- colorRampPalette(c("white", "red", "yellow"))(101)
    pheatmap::pheatmap(gg_cor, treeheight_row = 0, treeheight_col = 0, breaks = breaks, col = shades, filename = "figs/paper_figs/fig_s4/epiblast_genes_correlation_heatmap.pdf", w = 4, h = 4)
}






nascent_mesoderm_score_correlation_heatmap <- function() {
    mc <- scdb_mc("sing_emb_wt10_recolored")

    legc <- log2(mc@e_gc + 1e-5)

    gg_cor <- tgs_cor(t(legc[early_nascent_mesoderm_genes(), ]))

    breaks <- seq(0, 1, length.out = 101)
    # shades = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(100))
    shades <- colorRampPalette(c("white", "red", "yellow"))(101)
    pheatmap::pheatmap(gg_cor, treeheight_row = 0, treeheight_col = 0, breaks = breaks, col = shades, filename = "figs/paper_figs/fig_s4/early_nm_genes_correlation_heatmap.pdf", w = 4, h = 4)
}



plot_heatmaps_epiblast_and_nascent_mesoderm_score_genes <- function() {
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    cell_type_color <- mc_wt@color_key$color[1]

    genes_plot <- epiblast_genes()
    plot_heatmap_projection_on_atlas(cell_type_color = cell_type_color, filename = "figs/paper_figs/fig_s4/heatmap_epiblast_score_genes.pdf", genes_plot = genes_plot)

    cell_type_color <- mc_wt@color_key$color[10]

    genes_plot <- early_nascent_mesoderm_genes()
    plot_heatmap_projection_on_atlas(cell_type_color = cell_type_color, filename = "figs/paper_figs/fig_s4/heatmap_early_nascent_mesoderm_score_genes.pdf", genes_plot = genes_plot)
}


plot_heatmap_projection_on_atlas <- function(cell_type_color, filename, genes_plot) {
    mat_wt <- scdb_mat("sing_emb_wt10")
    mc_wt <- scdb_mc("sing_emb_wt10_recolored")

    mc_ds <- t(tgs_matrix_tapply(mat_wt@mat[, names(mc_wt@mc)], mc_wt@mc, sum))
    f <- colSums(mc_ds) >= 200000
    mc_ds_f <- scm_downsamp(umis = mc_ds[, f], n = 200000)
    mc_ds[, colnames(mc_ds_f)] <- as.matrix(mc_ds_f)

    egc <- t(t(mc_ds) / colSums(mc_ds))

    gset <- scdb_gset("sing_emb_wt10")
    bad_genes <- c("Tet1", "Tet2", "Tet3")


    feat_genes <- setdiff(names(gset@gene_set), bad_genes)

    tko_chim <- scdb_mat("tko_chim_wt10")
    tko_tetra <- scdb_mat("tko_tetra_wt10")
    control_tetra_all <- scdb_mat("control_tetra_all_wt10")

    load("data/tko_chim_wt10/tko_chim_md.Rda")
    load("data/tko_tetra_wt10/tko_tetra_md.Rda")
    load("data/control_tetra_all_wt10/control_tetra_all_md.Rda")

    f_chim <- (tko_chim_md$cell_type == "KO") & (tko_chim_md$ct_color == cell_type_color)
    f_tetra <- (tko_tetra_md$cell_type == "KO") & (tko_tetra_md$ct_color == cell_type_color)
    f_host <- (tko_chim_md$cell_type %in% c("host")) & (tko_chim_md$ct_color == cell_type_color)
    f_ctrl_chim <- (tko_chim_md$cell_type %in% c("control")) & (tko_chim_md$ct_color == cell_type_color)
    f_ctrl_tetra <- (control_tetra_all_md$cell_type == "control") & (control_tetra_all_md$ct_color == cell_type_color)

    q_mat <- cbind(
        as.matrix(tko_chim@mat[, tko_chim_md$cell[f_chim]]),
        as.matrix(tko_tetra@mat[, tko_tetra_md$cell[f_tetra]]),
        as.matrix(tko_chim@mat[, tko_chim_md$cell[f_host]]),
        as.matrix(tko_chim@mat[, tko_chim_md$cell[f_ctrl_chim]]),
        as.matrix(control_tetra_all@mat[, control_tetra_all_md$cell[f_ctrl_tetra]])
    )

    ref_mat <- cmp_proj(query_mat = q_mat, egc_ref = egc, feat_genes = feat_genes)


    md_cls <- rbind(
        tko_chim_md[f_chim, ],
        tko_tetra_md[f_tetra, ],
        tko_chim_md[f_host, ],
        tko_chim_md[f_ctrl_chim, ],
        control_tetra_all_md[f_ctrl_tetra, ]
    )

    md_cls$assay <- c(
        rep("2N", sum(f_chim)),
        rep("4N", sum(f_tetra)),
        rep("2N", sum(f_host)),
        rep("2N", sum(f_ctrl_chim)),
        rep("4N", sum(f_ctrl_tetra))
    )

    md_cls$assay_clone <- paste(md_cls$assay, md_cls$clone, sep = "_")


    egc_q <- t(tgs_matrix_tapply(q_mat, md_cls$assay_clone, sum))
    egc_q <- t(t(egc_q) / colSums(egc_q))

    egc_ref <- t(tgs_matrix_tapply(ref_mat, md_cls$assay_clone, sum))
    egc_ref <- t(t(egc_ref) / colSums(egc_ref))


    legc_clone <- log2(egc_q + 1e-5)
    legc_ref <- log2(egc_ref + 1e-5)
    lfp <- legc_clone - legc_ref

    included_clones <- c("2N_host", "2N_Ctrl1", "4N_Ctrl1", "2N_TKO23", "2N_TKO29", "4N_TKO26", "4N_TKO29")

    lfp_plot <- lfp[genes_plot, included_clones]

    shades <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(100))
    breaks <- seq(-1.5, 1.5, length.out = 101)

    pheatmap::pheatmap(pmin(pmax(lfp_plot, -1.5), 1.5),
        col = shades, breaks = breaks,
        cluster_cols = F, cluster_rows = F,
        filename = filename, w = 3, h = 4
    )
}


cmp_proj <- function(query_mat, egc_ref, feat_genes) {
    legc_ref <- log2(egc_ref[feat_genes, ] + 1e-5)
    q_r_cor <- tgs_cor(query_mat[feat_genes, ], legc_ref)

    best_ref <- apply(q_r_cor, 1, which.max)

    ref_mat <- egc_ref[, best_ref]
    return(ref_mat)
}

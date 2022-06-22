
# Figure S6


generate_figure_s6_plots <- function() {
    if (!dir.exists("figs/paper_figs/fig_s6")) {
        dir.create("figs/paper_figs/fig_s6")
    }

    # Fig S6I
    plot_heatmap_marker_genes_node_and_foregut_lineage()
}


plot_heatmap_marker_genes_node_and_foregut_lineage <- function() {
    if (!dir.exists("figs/paper_figs/fig_s6")) {
        dir.create("figs/paper_figs/fig_s6")
    }

    mc_node <- c(367, 108, 366)
    mc_gut <- c(377, 382, 385)

    node_genes <- c("Foxj1;Rnf157", "Noto", "T")
    gut_genes <- c("Sox17", "Spink3")
    pit_genes <- c("Tppp3", "Fam183b")
    notochord_genes <- c("Uhrf1", "Cer1", "Id3")
    crown_genes <- c("Shh", "Tbx6", "Cdx1")

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")
    col_to_ct <- mc_wt@color_key$group
    names(col_to_ct) <- mc_wt@color_key$color
    ct_to_col <- mc_wt@color_key$color
    names(ct_to_col) <- mc_wt@color_key$group

    legc <- log2(mc_wt@e_gc + 1e-5)

    legc_gut <- rowMeans(legc[, mc_gut])
    lfc_node_gut <- rowMeans(legc[, mc_node]) - rowMeans(legc[, mc_gut])

    genes_all <- c(node_genes, gut_genes, notochord_genes, pit_genes, crown_genes)
    genes_all <- c(
        "T", "Aldh1a2", "Cdx1", "Cer1", "Dand5;SP1", "Dnaja4", "Dnajb13", "Dynlrb2", "Fam183b", "Foxa2", "Foxj1;Rnf157", "Fst", "Gal",
        "Hoxb1", "Id3", "Igfbp5", "Krt18", "Krt8", "Mixl1", "Nog", "Noto", "Rspo3", "Shh", "Slc16a1", "Slc2a3", "Sox17", "Sp5",
        "Spink3", "Tbx6", "Tppp3", "Wnt3a"
    )

    annotation_col <- data.frame("Cell type" = col_to_ct[mc_wt@colors[c(mc_node, mc_gut)]], stringsAsFactors = F)
    colnames(annotation_col) <- "Cell type"
    rownames(annotation_col) <- c(mc_node, mc_gut)
    annotation_colors <- list("Cell type" = ct_to_col[c(24, 27)])
    shades <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, "BuPu"))(101)
    pheatmap::pheatmap(legc[genes_all, c(mc_node, mc_gut)],
        color = shades, cluster_cols = F, treeheight_row = 0,
        annotation_col = annotation_col, annotation_colors = annotation_colors,
        filename = "figs/paper_figs/fig_s6/hm_marker_genes.pdf"
    )
}

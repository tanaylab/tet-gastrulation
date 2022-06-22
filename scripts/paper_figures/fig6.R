

generate_figure6_plots <- function() {
    if (!dir.exists("figs/paper_figs/fig6")) {
        dir.create("figs/paper_figs/fig6")
    }

    # fig 6d
    trajectories_epigenetic_factors()

    # fig 6e
    cell_cycle_analysis()
}


trajectories_epigenetic_factors <- function() {
    if (!dir.exists("figs/paper_figs/fig6/time_trajectories")) {
        dir.create("figs/paper_figs/fig6/time_trajectories")
    }

    age_groups_time <- read.table("data/wt10_age_groups_developmental_time.txt", sep = "\t", stringsAsFactors = F)

    reg <- 1e-5
    mc_node <- c(108, 366, 367)
    mc_gut <- c(377, 382, 385)
    mct <- scdb_mctnetwork("sing_emb_wt10")
    mat <- scdb_mat("sing_emb_wt10")
    mc <- scdb_mc("sing_emb_wt10_recolored")

    egc <- t(tgs_matrix_tapply(mat@mat[, names(mc@mc)], mc@mc, sum))
    egc <- t(t(egc) / colSums(egc))

    genes_f <- c("Tet1", "Tet2", "Tet3", "Uhrf1", "Dnmt1", "Dnmt3a", "Dnmt3b")

    df_egc_all <- data.frame()

    mc_type <- c(rep(mc@color_key$group[24], length(mc_node)), rep(mc@color_key$group[27], length(mc_gut)))
    names(mc_type) <- c(mc_node, mc_gut)

    for (mc_f in c(mc_node)) {
        mc_p <- rep(0, ncol(mc@e_gc))
        mc_p[mc_f] <- 1
        prop <- mctnetwork_propogate_from_t(mct = mct, t = 13, mc_p = mc_p)

        egc_t <- egc %*% prop$probs
        egc_t <- log2(egc_t + reg)

        df_egc <- as.data.frame(t(egc_t[genes_f, ]))
        df_egc$time <- age_groups_time$developmental_time
        df_egc$cell_type <- NA
        df_egc$cell_type <- rep(mc_type[as.character(mc_f)], 13)
        df_egc$mc <- rep(mc_f, 13)

        df_egc_all <- rbind(df_egc_all, df_egc)
    }

    mc_gut_bulk <- which(mc@colors == mc@color_key$color[27])
    mc_p <- rep(0, ncol(mc@e_gc))
    mc_p[mc_gut_bulk] <- mct@mc_t[mc_gut_bulk, 13]
    mc_p <- mc_p / sum(mc_p)
    prop <- mctnetwork_propogate_from_t(mct = mct, t = 13, mc_p = mc_p)

    egc_t <- egc %*% prop$probs
    egc_t <- log2(egc_t + reg)

    df_egc <- as.data.frame(t(egc_t[genes_f, ]))
    df_egc$time <- age_groups_time$developmental_time
    df_egc$cell_type <- NA
    df_egc$cell_type <- rep(mc@color_key$group[27], 13)
    df_egc$mc <- rep("Foregut Bulk", 13)

    df_egc_all <- rbind(df_egc_all, df_egc)

    mc_color <- c(rep(mc@color_key$color[24], length(mc_node)), c(mc@color_key$color[27]))
    names(mc_color) <- c(mc_node, "Foregut Bulk")

    df_plot <- pivot_longer(data = df_egc_all, cols = all_of(genes_f), names_to = "gene", values_to = "expression")
    df_plot$mc <- as.character(df_plot$mc)

    linetype_values <- c("solid", "dashed", "dotted", "solid")
    names(linetype_values) <- c(mc_node, "Foregut Bulk")

    plot_list <- list()

    for (gene in genes_f) {
        df_f <- df_plot[df_plot$gene == gene, ]
        p <- ggplot(data = df_f, aes(x = time, y = expression, color = mc, linetype = mc)) +
            geom_line(size = 1.5) +
            scale_color_manual(values = mc_color, labels = c("Crown cell", "Notochord", "Ciliary node", "Foregut Bulk")) +
            scale_linetype_manual(values = linetype_values, labels = c("Crown cell", "Notochord", "Ciliary node", "Foregut Bulk")) +
            ggtitle(label = gene) +
            theme(legend.position = "none", plot.title = element_text(size = 20, hjust = 0.5)) +
            xlab("Time") +
            ylab("Expression")


        ggsave(filename = sprintf("figs/paper_figs/fig6/time_trajectories/%s.png", gene), height = 4, width = 5)

        plot_list[[gene]] <- p
    }

    gene <- "Dnmt1"

    df_f <- df_plot[df_plot$gene == gene, ]
    p <- ggplot(data = df_f, aes(x = time, y = expression, color = mc, linetype = mc)) +
        geom_line(size = 1.5) +
        scale_color_manual(values = mc_color, labels = c("Crown cell", "Notochord", "Ciliary node", "Foregut Bulk")) +
        scale_linetype_manual(values = linetype_values, labels = c("Crown cell", "Notochord", "Ciliary node", "Foregut Bulk")) +
        ggtitle(label = gene) +
        theme(legend.key.size = unit(1.5, "cm")) +
        labs(color = "Cell type", linetype = "Cell type") +
        theme(legend.background = element_rect(color = "black"))

    p_leg <- as_ggplot(get_legend(p))
    ggsave(plot = p_leg, filename = "figs/paper_figs/fig6/time_trajectories/legend_trajectories.pdf")



    p_all <- arrangeGrob(plot_list[["Dnmt1"]], plot_list[["Uhrf1"]],
        plot_list[["Tet1"]], plot_list[["Tet2"]],
        plot_list[["Tet3"]], p_leg,
        ncol = 2, nrow = 3,
        layout_matrix = matrix(c(1, 2, 3, 4, 5, 6), ncol = 2, nrow = 3, byrow = T)
    )

    ggsave(plot = p_all, filename = "figs/paper_figs/fig6/time_trajectories/genes_combined.pdf", width = 7, height = 12)
}


cell_cycle_analysis <- function() {
    if (!dir.exists("figs/paper_figs/fig6/cell_cycle_analysis")) {
        dir.create("figs/paper_figs/fig6/cell_cycle_analysis")
    }

    m_genes <- c("Mki67", "Cenpf", "Top2a", "Smc4;SMC4", "Ube2c", "Ccnb1", "Cdk1", "Arl6ip1", "Ankrd11", "Hmmr;IHABP", "Cenpa;Cenp-a", "Tpx2", "Aurka", "Kif4", "Kif2c", "Bub1b", "Ccna2", "Kif23", "Kif20a", "Sgol2", "Smc2", "Kif11", "Cdca2", "Incenp", "Cenpe")
    s_genes <- c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2", "Uhrf1", "Orc6", "Tipin") # Npm1

    mat <- scdb_mat("sing_emb_wt10")
    mc <- scdb_mc("sing_emb_wt10_recolored")

    s_genes <- intersect(rownames(mc@mc_fp), s_genes)
    m_genes <- intersect(rownames(mc@mc_fp), m_genes)

    tot <- colSums(mat@mat[, names(mc@mc)])
    s_tot <- colSums(mat@mat[s_genes, names(mc@mc)])
    m_tot <- colSums(mat@mat[m_genes, names(mc@mc)])

    s_score <- s_tot / tot
    m_score <- m_tot / tot

    m_s_score <- 0.5 * m_score / mean(m_score) + 0.5 * s_score / mean(s_score)

    mc_node <- which(mc@colors == mc@color_key$color[24])
    mc_gut <- c(377, 382, 385)

    sc_mc <- mc@mc
    sc_mc[!mc@mc %in% c(mc_node, mc_gut)] <- "bulk"


    f <- !mc@mc %in% c(mc_node, mc_gut)


    m_s_score <- 0.5 * m_score / mean(m_score) + 0.5 * s_score / mean(s_score)


    for (mc_f in mc_node) {
        score_dens <- density(m_s_score[names(mc@mc)[f]], from = 0)
        score_dens_mc_f <- density(m_s_score[names(mc@mc)[mc@mc == mc_f]], from = 0, to = 2)

        df_plot <- data.frame(score = c(score_dens$x, score_dens_mc_f$x), y = c(score_dens$y, score_dens_mc_f$y), type = c(rep("Embryo Bulk", length(score_dens$x)), rep("Notochord", length(score_dens_mc_f$x))))

        fill_colors <- c(mc@color_key$color[24], "gray50")
        names(fill_colors) <- c("Notochord", "Embryo Bulk")

        p <- ggplot(df_plot) +
            geom_area(aes(x = score, y = y, fill = type), alpha = 0.4) +
            scale_fill_manual(values = fill_colors) +
            geom_line(aes(x = score, y = y, color = type), size = 1.5) +
            scale_color_manual(values = fill_colors) +
            ylab("") +
            xlab("M/S Phase Score") +
            ggtitle(sprintf("MC #%d", mc_f)) +
            theme(plot.title = element_text(size = 20, hjust = 0.5)) +
            ylim(0, 3.2) +
            xlim(0, 2) +
            theme(legend.position = "none")

        ggsave(plot = p, sprintf("figs/paper_figs/fig6/cell_cycle_analysis/distribution_m_s_score_notochord_mc%d.pdf", mc_f), width = 10, height = 7)
    }
}



generate_ko_mat_chimera <- function() {
    load("data/tko_chim_wt10/color_annotation/cmp_annot.Rda")

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")
    mat <- scdb_mat("tko_chim")
    mat_f <- scdb_mat("tko_chim_wt10")

    cls_f <- names(cmp_annot$query_cls_col)[!(cmp_annot$query_cls_col %in% mc_wt@color_key$color[c(28, 29)])]
    cls_f2 <- colnames(mat@mat)[mat@cell_metadata[colnames(mat@mat), "cell_type"] == "KO"]

    cls_f <- intersect(cls_f, cls_f2)

    mat_new <- scm_ignore_cells(scmat = mat, ig_cells = cls_f, reverse = T)
    scdb_add_mat(id = "tko_chim_ko", mat = mat_new)
}

generate_host_mat <- function() {
    mat_all <- scdb_mat("tko_chimera")
    mat <- scdb_mat("tko_chim_wt10")
    load("data/tko_chim_wt10/color_annotation/cmp_annot.Rda")

    cls_host <- colnames(mat@mat)[mat@cell_metadata[colnames(mat@mat), "cell_type"] == "host"]
    cls_host <- intersect(colnames(mat_all@mat), cls_host)

    mat_host1 <- scm_ignore_cells(scmat = mat_all, ig_cells = cls_host, reverse = T)

    mat2 <- scdb_mat("control_chim")

    cls_host2 <- colnames(mat2@mat)[mat2@cell_metadata[colnames(mat2@mat), "cell_type"] == "host"]
    mat_host2 <- scm_ignore_cells(scmat = mat2, ig_cells = cls_host2, reverse = T)

    source("scripts/pipeline/merge_umi_mat_with_wt10_umi_mat.r")
    mat_new <- merge_two_scmat(mat1 = mat_host1, mat2 = mat_host2)
    scdb_add_mat("tko_control_host", mat = mat_new)
}



generate_ko_mat_tetraploid <- function() {
    load("data/tko_tetra_wt10/color_annotation/cmp_annot.Rda")

    mc_wt <- scdb_mc("sing_emb_wt10_recolored")
    mat <- scdb_mat("tko_tetra")
    mat_f <- scdb_mat("tko_tetra_wt10")

    cls_f <- names(cmp_annot$query_cls_col)[!(cmp_annot$query_cls_col %in% mc_wt@color_key$color[c(28, 29)])]
    cls_f2 <- colnames(mat@mat)[mat@cell_metadata[colnames(mat@mat), "cell_type"] == "KO"]

    cls_f <- intersect(cls_f, cls_f2)

    mat_new <- scm_ignore_cells(scmat = mat, ig_cells = cls_f, reverse = T)
    scdb_add_mat(id = "tko_tetra_ko", mat = mat_new)
}

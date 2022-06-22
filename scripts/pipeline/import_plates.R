



import_plates <- function(mat_nm, metadata, base_dir = "data/umi.tables/") {
    metadata$Amp.Batch.ID <- metadata$plate
    metadata$Seq.Batch.ID <- metadata$plate
    metadata$Batch.Set.ID <- metadata$plate

    write.table(x = metadata, file = paste("config/key_", mat_nm, ".txt", sep = ""), quote = F, sep = "\t", row.names = F)

    mcell_import_multi_mars(
        mat_nm = mat_nm,
        dataset_table_fn = paste("config/key_", mat_nm, ".txt", sep = ""),
        base_dir = base_dir,
        patch_cell_name = F,
        force = TRUE
    )

    mat <- scdb_mat(mat_nm)
    nms <- c(rownames(mat@mat), rownames(mat@ignore_gmat))
    bad_genes <- c(grep("^mt\\-", nms, v = T), "Neat1", grep("ERCC", nms, v = T), "Atpase6", "Xist", "Malat1", "Cytb", "AK018753", "AK140265", "AK163440", "DQ539915")

    mat <- scm_ignore_genes(scmat = mat, ig_genes = bad_genes)

    # ignore small and large cells

    small_cells <- colnames(mat@mat)[colSums(mat@mat) < 2000]
    large_cells <- colnames(mat@mat)[colSums(mat@mat) > 12000]

    mat <- scm_ignore_cells(scmat = mat, ig_cells = c(small_cells, large_cells))

    # read metadata
    mat_md <- mat@cell_metadata
    mat_md$cell <- rownames(mat@cell_metadata)

    f <- colnames(mat_md) %in% c("molecule", "spike_count")
    mat_md <- mat_md[, !f]

    metadata_cells <- read.table(file = "data/metadata_cells_scrna_sequencing.txt", h = T, sep = "\t", stringsAsFactors = F)

    metadata_cells <- metadata_cells[colnames(metadata_cells) != "plate"]

    mat_md <- left_join(mat_md, metadata_cells, by = "cell")

    rownames(mat_md) <- mat_md$cell

    mat@cell_metadata <- mat_md

    # ignore empty wells and wells with duplicate cells



    new_ignore_cells <- unique(c(mat_md$cell[mat_md$embryo %in% c("empty", "duplicate")], mat@ignore_cells))

    mat <- scm_ignore_cells(scmat = mat, ig_cells = new_ignore_cells)

    scdb_add_mat(id = mat_nm, mat = mat)

    file.remove(paste("config/key_", mat_nm, ".txt", sep = ""))
}

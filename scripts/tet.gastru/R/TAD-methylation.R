#' Get 'background' methylation in TADs
#'
#' @param tad_intervs intervals of the TADs
#' @param meth_dips,k4me3_peaks,k27me3_peaks,atac_peaks intervals of functional annotations (see notebook)
#'
#' @export
get_cpg_meth_tads <- function(tad_intervs, meth_dips, k4me3_peaks, k27me3_peaks, atac_peaks) {
    gvtrack.create("tor", "Encode.esd3.replichip.rep2", "avg")
    cpg_meth_tads <- gextract_meth(c(comb_tracks, "e8_5", "Zhang_Nature_Genetics_2017.Ect_mCG", "Zhang_Nature_Genetics_2017.End_mCG", "Zhang_Nature_Genetics_2017.Mes_mCG", "Zhang_Nature_Genetics_2017.E65Epi_mCG"), names = c(comb_names, "e8.5", "ecto", "endo", "meso", "epi"), intervals = tad_intervs, join_intervals = TRUE, extract_meth_calls = TRUE)

    # Set of intervals that is clean of dips,k4me3 and k27me3
    bg_intervs1 <- gintervals.diff(gintervals.all(), meth_dips) %>%
        gintervals.diff(k4me3_peaks) %>%
        gintervals.diff(k27me3_peaks) %>%
        gintervals.diff(atac_peaks)

    cpg_meth_tads <- cpg_meth_tads %>%
        gintervals.neighbors1(bg_intervs1) %>%
        filter(dist == 0)
    cpg_meth_tads <- cpg_meth_tads %>%
        mutate(e7.5.meth = ecto.meth + meso.meth + endo.meth, e7.5.cov = ecto.cov + meso.cov + endo.cov) %>%
        group_by(chrom1, start1, end1, id) %>%
        summarise_at(vars(ends_with("cov"), ends_with("meth")), sum, na.rm = TRUE) %>%
        ungroup() %>%
        mutate(
            ctrl = ctrl.meth / ctrl.cov,
            tet_tko = tet_tko.meth / tet_tko.cov,
            e8.5 = e8.5.meth / e8.5.cov,
            e7.5 = e7.5.meth / e7.5.cov,
            epi = epi.meth / epi.cov
        ) %>%
        rename(chrom = chrom1, start = start1, end = end1)

    tads_meth <- gextract.left_join("tor", intervals = cpg_meth_tads, iterator = cpg_meth_tads)
    return(tads_meth)
}

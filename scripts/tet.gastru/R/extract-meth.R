
#' Extract methylation at an intervals set
#'
#' @param intervals intervals to extract methylation at
#' @param iterator iterator of gextract_meth call
#' @param min_cov minimal coverage
#'
#' @export
extract_meth_at_context <- function(intervals, iterator = intervals, min_cov = NULL) {
    gextract_meth(c(comb_tracks, "e8_5", "Zhang_Nature_Genetics_2017.Ect_mCG", "Zhang_Nature_Genetics_2017.End_mCG", "Zhang_Nature_Genetics_2017.Mes_mCG", "Zhang_Nature_Genetics_2017.E65Epi_mCG"), names = c(comb_names, "e8.5", "ecto", "endo", "meso", "epi"), extract_meth_calls = TRUE, intervals = intervals, iterator = iterator, min_cov = min_cov) %>%
        mutate(e7.5.meth = ecto.meth + endo.meth + meso.meth, e7.5.cov = ecto.cov + endo.cov + meso.cov, e7.5 = e7.5.meth / e7.5.cov)
}

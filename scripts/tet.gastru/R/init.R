#' Initalize analysis
#'
#' @param use_cache Use cache if available
#'
#' @export
init <- function(use_cache = TRUE) {
    purrr::walk(c("dplyr", "purrr", "tidyr", "ggplot2", "tgstat", "misha", "misha.ext", "here", "patchwork", "tgutil", "glue"), library, character.only = TRUE, warn.conflicts = FALSE)
    db_dir <- here::here("db")
    if (!dir.exists(db_dir)) {
        cli_abort("db directory does not exist. Download it using {.code download_misha_db()}")
    }
    gsetroot(db_dir)

    options(tgutil.verbose = FALSE)
    options(gmax.data.size = 1e9)
    options(gmultitasking = FALSE)

    if (dir.exists(here("output"))) {
        dir.create(here("output"), showWarnings = FALSE)
    }

    tracks <<- c(
        "TET.control_s1",
        "TET.control_s2",
        "TET.tet_tko_s1",
        "TET.tet_tko_s2"
    )

    md <<- tibble(track = tracks, name = c("ctrl1", "ctrl2", "tet_tko1", "tet_tko2"))

    comb_tracks <<- c("TET.control", "TET.tet_tko")
    comb_names <<- c("ctrl", "tet_tko")
    comb_md <<- tibble(track = comb_tracks, name = comb_names)

    if (use_cache) {
        options(tgutil.cache = TRUE)
    } else {
        options(tgutil.cache = FALSE)
    }
}

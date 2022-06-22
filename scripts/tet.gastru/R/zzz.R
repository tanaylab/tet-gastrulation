.onLoad <- function(libname, pkgname) {
    set.seed(17)

    ggplot2::theme_set(tgppt::theme_arial(8))

    options(gmax.data.size = 1e9, gmultitasking = FALSE)
}

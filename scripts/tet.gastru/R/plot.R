#' Plot a smooth scatter of methylation data
#'
#' @param df Dataframe with methylation data
#' @param col1 Column name for the first dataset
#' @param col2 Column name for the second dataset
#' @param min_cov Minmial coverage for a site to be included in the plot
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param colors Vector of colors for the scatter
#'
#' @export
plot_smooth_cpg_scatter <- function(df, col1, col2, min_cov, xlab = col1, ylab = col2, colors = c("white", "lightblue", "blue", "darkblue", "yellow", "gold", "orange", "red", "darkred")) {
    df <- df %>%
        filter(!!sym(glue("{col1}.cov")) >= min_cov, !!sym(glue("{col2}.cov")) >= min_cov) %>%
        as_tibble()
    shades <- colorRampPalette(colors)
    smoothScatter(df[[col1]], df[[col2]], colramp = shades, pch = 19, cex = 0.1, xlab = xlab, ylab = ylab)
    title(glue("{scales::comma(nrow(df))} covered CpGs"))
}

#' Plot a legend for smooth scatter
#'
#' @param colors Vector of colors
#' @param n_colors Number of colors in the legend
#'
#' @examples
#' plot_smooth_legend(c("white", "lightblue", "blue", "darkblue", "yellow", "gold", "orange", "red", "darkred"), 1000)
#'
#' @export
plot_smooth_legend <- function(colors, n_colors = 1000) {
    color <- colorRampPalette(colors)(n_colors)
    plot(NA, xlim = c(0, n_colors), ylim = c(0, n_colors + 1), type = "n", ann = FALSE, axes = FALSE)
    rect(0, 1:n_colors, 50, 2:(n_colors + 1), border = NA, col = color)
}

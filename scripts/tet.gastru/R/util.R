#' Add vertival x axis labels in ggplot
#'
#' @export
vertical_labs <- function() {
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

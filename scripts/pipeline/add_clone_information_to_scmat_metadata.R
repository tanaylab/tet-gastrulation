

add_clone_type_information_to_scmat <- function(mat) {
    md <- mat@cell_metadata

    md$clone_type <- "unclear"
    cls_host <- colnames(mat@mat)[md[colnames(mat@mat), "cell_type"] == "host"]
    md[cls_host, "clone_type"] <- "host"

    for (clone_type in c("Ctrl1", "Ctrl2", "Ctrl3")) {
        cls_f <- colnames(mat@mat)[grep(clone_type, md[colnames(mat@mat), "embryo"])]
        cls_f <- cls_f[md[cls_f, "cell_type"] == "control"]
        md[cls_f, "clone_type"] <- clone_type
    }

    for (clone_type in c("TKO23", "TKO26", "TKO29")) {
        cls_f <- colnames(mat@mat)[grep(clone_type, md[colnames(mat@mat), "embryo"])]
        cls_ko <- cls_f[md[cls_f, "cell_type"] == "KO"]
        md[cls_ko, "clone_type"] <- clone_type

        cls_control <- cls_f[md[cls_f, "cell_type"] == "control"]
        if (clone_type == "TKO29") {
            md[cls_control, "clone_type"] <- "Ctrl1"
        }
        if (clone_type == "TKO26") {
            md[cls_control, "clone_type"] <- "Ctrl2"
        }
    }

    mat@cell_metadata <- md
    return(mat)
}

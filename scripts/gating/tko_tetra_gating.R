
arsinh <- function(x, a = 1, b = 1) {
    return(log(x + sqrt(x^2 / a^2 + b / a)))
}

calc_a_b_from_x_y <- function(x1, y1, x2, y2) {
    b <- (y2 - y1) / (x2 - x1)
    a <- (y1 * x2 - x1 * y2) / (x2 - x1)
    return(c(a, b))
}

experiment <- "Tetraploid complementation assay"

mat_nm <- "tko_tetra"
mat <- scdb_mat(mat_nm)

if (!dir.exists(sprintf("%s/gating", .scfigs_base))) {
    dir.create(sprintf("%s/gating", .scfigs_base))
}

if (!dir.exists(sprintf("%s/gating/%s", .scfigs_base, mat_nm))) {
    dir.create(sprintf("%s/gating/%s", .scfigs_base, mat_nm))
}
fig_dir <- sprintf("%s/gating/%s", .scfigs_base, mat_nm)

cls_type <- rep("unclear", length(colnames(mat@mat)))
names(cls_type) <- colnames(mat@mat)

cls_clone <- rep("unclear", length(colnames(mat@mat)))
names(cls_type) <- colnames(mat@mat)

cls <- colnames(mat@mat)

cls_mcherry <- cls[is.na(mat@cell_metadata[cls, "tdtomato_a"]) & !is.na(mat@cell_metadata[cls, "mcherry_a"])]
cls_tomato <- cls[!is.na(mat@cell_metadata[cls, "tdtomato_a"])]

# color cells by wt9 atlas color

cls_color <- proj_on_wt9_atlas(mat_query = mat)

f <- !(cls_color %in% c("#F6BFCB", "#7F6874", "#1A1A1A", "#989898"))

cls_color[f] <- "gray80"

cls_color[cls_color == "#989898"] <- "indianred3"

# parameters needed for the transformation
a_x <- 1
b_x <- 1000
a_y <- 1
b_y <- 1000




# gfp_tr = arsinh(mat@cell_metadata[cls,"gfp_a"], a= a_x,b = b_x)
tdtomato_tr <- arsinh(mat@cell_metadata[cls, "tdtomato_a"], a = a_y, b = b_y)
mcherry_tr <- arsinh(mat@cell_metadata[cls, "mcherry_a"], a = a_y, b = b_y)
tdtomato_tr[is.na(tdtomato_tr)] <- mcherry_tr[is.na(tdtomato_tr)]
# tdtomato_tr[!is.na(mcherry_tr)] = mcherry_tr[!is.na(mcherry_tr)]
ssc_a <- mat@cell_metadata[cls, "ssc_a"]
fsc_a <- mat@cell_metadata[cls, "fsc_a"]

names(tdtomato_tr) <- cls
names(ssc_a) <- cls
names(fsc_a) <- cls
names(mcherry_tr) <- cls

lim_tdtomato <- c(min(tdtomato_tr), max(tdtomato_tr))
lim_fsc <- c(min(fsc_a), max(fsc_a))
lim_ssc <- c(min(ssc_a), max(ssc_a))

png(sprintf("%s/tdtomato_vs_fsc_all.png", fig_dir))
plot(tdtomato_tr, fsc_a, pch = 19, cex = 0.5, xlab = "tdTomtato", ylab = "FSC-A", col = cls_color, main = paste("All tetraploid TKO embryos", sep = " "), xlim = lim_tdtomato, ylim = lim_fsc)
points(tdtomato_tr[!f], fsc_a[!f], cex = 0.5, pch = 19, col = cls_color[!f])
dev.off()


png(sprintf("%s/tdtomato_vs_ssc_all.png", fig_dir))
plot(tdtomato_tr, ssc_a, pch = 19, cex = 0.5, xlab = "tdTomato", ylab = "SSC-A", col = cls_color, main = paste("All tetraploid TKO embryos", sep = " "), xlim = lim_tdtomato, ylim = lim_ssc)
points(tdtomato_tr[!f], ssc_a[!f], cex = 0.5, pch = 19, col = cls_color[!f])
dev.off()

png(sprintf("%s/fsc_vs_ssc_all.png", fig_dir))
plot(fsc_a, ssc_a, pch = 19, cex = 0.5, xlab = "FSC-A", ylab = "SSC-A", col = cls_color, main = paste("All tetraploid TKO embryos", sep = " "), xlim = lim_fsc, ylim = lim_ssc)
points(fsc_a[!f], ssc_a[!f], cex = 0.5, pch = 19, col = cls_color[!f])
dev.off()


# plot by clone
for (clone in c("TKO26", "TKO29")) {
    cls_f <- cls[grep(clone, mat@cell_metadata[cls, "embryo"])]

    ff <- !f & (cls %in% cls_f)

    png(sprintf("%s/tdtomato_vs_fsc_%s.png", fig_dir, clone))
    plot(tdtomato_tr[cls_f], fsc_a[cls_f], pch = 19, cex = 0.5, xlab = "tdTomtato", ylab = "FSC-A", col = cls_color[cls_f], main = paste(clone, " tetraploid embryos", sep = " "), xlim = lim_tdtomato, ylim = lim_fsc)
    points(tdtomato_tr[ff], fsc_a[ff], cex = 0.5, pch = 19, col = cls_color[ff])
    dev.off()

    png(sprintf("%s/tdtomato_vs_ssc_%s.png", fig_dir, clone))
    plot(tdtomato_tr[cls_f], ssc_a[cls_f], pch = 19, cex = 0.5, xlab = "tdTomato", ylab = "SSC-A", col = cls_color[cls_f], main = paste(clone, " tetraploid embryos", sep = " "), xlim = lim_tdtomato, ylim = lim_ssc)
    points(tdtomato_tr[ff], ssc_a[ff], cex = 0.5, pch = 19, col = cls_color[ff])
    dev.off()

    png(sprintf("%s/fsc_vs_ssc_all_%s.png", fig_dir, clone))
    plot(fsc_a[cls_f], ssc_a[cls_f], pch = 19, cex = 0.5, xlab = "FSC-A", ylab = "SSC-A", col = cls_color[cls_f], main = paste(clone, " tetraploid embryos", sep = " "), xlim = lim_fsc, ylim = lim_ssc)
    points(fsc_a[ff], ssc_a[ff], cex = 0.5, pch = 19, col = cls_color[ff])
    dev.off()
}


# plot by mcherry vs tdtomato

cls_f <- cls_tomato

ff <- !f & (cls %in% cls_f)
a <- -750000
b <- 120000
clone <- "tdtomato"
main_tag <- "Tetraploid TKO cells without sorting date 19/10/20"

png(sprintf("%s/tdtomato_vs_fsc_%s.png", fig_dir, clone))
plot(tdtomato_tr[cls_f], fsc_a[cls_f], pch = 19, cex = 0.5, xlab = "tdTomtato", ylab = "FSC-A", col = cls_color[cls_f], main = main_tag, xlim = lim_tdtomato, ylim = lim_fsc)
abline(a = a, b = b)
points(tdtomato_tr[ff], fsc_a[ff], cex = 0.5, pch = 19, col = cls_color[ff])
dev.off()

png(sprintf("%s/tdtomato_vs_ssc_%s.png", fig_dir, clone))
plot(tdtomato_tr[cls_f], ssc_a[cls_f], pch = 19, cex = 0.5, xlab = "tdTomato", ylab = "SSC-A", col = cls_color[cls_f], main = main_tag, xlim = lim_tdtomato, ylim = lim_ssc)
points(tdtomato_tr[ff], ssc_a[ff], cex = 0.5, pch = 19, col = cls_color[ff])
dev.off()

f_ko <- (fsc_a[cls_f] < a + b * tdtomato_tr[cls_f]) & (cls_color[cls_f] == "gray80")
ko_cls <- cls_f[f_ko]
cls_type[ko_cls] <- "KO"

clone <- "mcherry"

cls_f <- cls_mcherry

ff <- !f & (cls %in% cls_f)

main_tag <- "Tetraploid TKO cells from sorting date 19/10/20"

png(sprintf("%s/mcherry_vs_fsc_%s.png", fig_dir, clone))
plot(mcherry_tr[cls_f], fsc_a[cls_f], pch = 19, cex = 0.5, xlab = "mcherry", ylab = "FSC-A", col = cls_color[cls_f], main = main_tag, xlim = lim_tdtomato, ylim = lim_fsc)
points(mcherry_tr[ff], fsc_a[ff], cex = 0.5, pch = 19, col = cls_color[ff])
abline(v = 5.5)
dev.off()

png(sprintf("%s/mcherry_vs_ssc_%s.png", fig_dir, clone))
plot(mcherry_tr[cls_f], ssc_a[cls_f], pch = 19, cex = 0.5, xlab = "mcherry", ylab = "SSC-A", col = cls_color[cls_f], main = main_tag, xlim = lim_tdtomato, ylim = lim_ssc)
points(mcherry_tr[ff], ssc_a[ff], cex = 0.5, pch = 19, col = cls_color[ff])
dev.off()

f_ko <- (tdtomato_tr[cls_f] > 5.5) & (cls_color[cls_f] == "gray80")
ko_cls <- cls_f[f_ko]
cls_type[ko_cls] <- "KO"


if (1) {
    df_gating <- data.frame(cell = colnames(mat@mat), cell_type = cls_type, stringsAsFactors = F)
    mat <- scdb_mat(mat_nm)
    md <- mat@cell_metadata

    f <- colnames(mat@cell_metadata) %in% c("cell_type", "clone_type", "cell_genotype")
    md <- md[, !f]


    md <- left_join(md, df_gating, by = "cell")
    rownames(md) <- md$cell
    mat@cell_metadata <- md

    mat <- add_clone_type_information_to_scmat(mat = mat)

    scdb_add_mat(id = mat_nm, mat = mat)
}

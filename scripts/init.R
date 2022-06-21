
packages = c(
  "Matrix",
  "metacell",
  "tgstat",
  "tgutil",
  "dplyr",
  "tidyr",
  "slam",
  "sparsesvd",
  "qlcMatrix",
  "zoo",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "qvalue",
  "gridExtra"
)

for (pkg in packages) {
  library(pkg,character.only = T)
}

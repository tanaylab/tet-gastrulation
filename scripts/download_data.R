# script for downloading and initializing the single-cell RNA-seq database for the metacell package
options(timeout = 1e9)

download_raw_umi_tables <- function() {
    download.file("https://tet-gastrulation.s3.eu-west-1.amazonaws.com/tet_umi_tables.tar.gz", "tet_umi_tables.tar.gz")

    system("tar -xzvf tet_umi_tables.tar.gz")

    file.remove("tet_umi_tables.tar.gz")
}

download_embflow_scrna_db <- function() {
    download.file("https://tet-gastrulation.s3.eu-west-1.amazonaws.com/scrna_db_embflow.tar.gz", "scrna_db_embflow.tar.gz")

    system("tar -xzvf scrna_db_embflow.tar.gz")

    file.remove("scrna_db_embflow.tar.gz")
}

download_tet_data <- function() {
    download.file("https://tet-gastrulation.s3.eu-west-1.amazonaws.com/tet_data.tar.gz", "tet_data.tar.gz")

    system("tar -xzvf tet_data.tar.gz")

    file.remove("tet_data.tar.gz")
}

download_tet_scrna_db <- function() {
    download.file("https://tet-gastrulation.s3.eu-west-1.amazonaws.com/scrna_db_tet.tar.gz", "scrna_db_tet.tar.gz")

    if (!dir.exists("scrna_db")) {
        dir.create("scrna_db")
    }

    system("tar -xzvf scrna_db_tet.tar.gz")

    file.remove("scrna_db_tet.tar.gz")
}

download_misha_db <- function() {
    download.file("https://tet-gastrulation.s3.eu-west-1.amazonaws.com/misha_db.tar.gz", "misha_db.tar.gz")

    system("tar -xzvf misha_db.tar.gz")

    file.remove("misha_db.tar.gz")
}

download_methylation_files <- function() {
    download.file("https://tet-gastrulation.s3.eu-west-1.amazonaws.com/methylation_data.tar.gz", "methylation_data.tar.gz")

    system("tar -xzvf methylation_data.tar.gz")

    file.remove("methylation_data.tar.gz")
}

download_minimal_scrna_data <- function() {
    download_raw_umi_tables()
    download_embflow_scrna_db()
}

download_full_data <- function() {
    download_tet_data()
    download_raw_umi_tables()
    download_embflow_scrna_db()
    download_tet_scrna_db()
    download_misha_db()
    download_methylation_files()
}


if (!dir.exists("figs")) {
    dir.create("figs")
}

if (!dir.exists("figs/paper_figs")) {
    dir.create("figs/paper_figs")
}

if (!dir.exists("output")) {
    dir.create("output")
}

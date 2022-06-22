# script for downloading and initializing the single-cell RNA-seq database for the metacell package

download_raw_umi_tables = function() {

  download.file("http://abc","tet_umi_tables.tar.gz")

  if (!dir.exists("data")) {
    dir.create("data")
  }
  system("tar -xzvf tet_umi_tables.tar.gz --directory data/")

  file.remove("tet_umi_tables.tar.gz")



  file.copy("umi_tables")

}

download_embflow_scrna_db = function() {

  download.file("http://abc","scrna_db_embflow.tar.gz")

  if (!dir.exists("scrna_db")) {
    dir.create("scrna_db")
  }

  system("tar -xzvf scrna_db_embflow.tar.gz --directory scrna_db/")

  file.remove("scrna_db_embflow.tar.gz")

}

download_tet_data = function() {

  download.file("http://abc","tet_data.tar.gz")

  if (!dir.exists("data")) {
    dir.create("data")
  }

  system("tar -xzvf tet_data.tar.gz --directory data/")

  file.remove("tet_data.tar.gz")


}

download_tet_scrna_db = function() {

  download.file("http://abc","scrna_db_tet.tar.gz")

  if (!dir.exists("scrna_db")) {
    dir.create("scrna_db")
  }

  system("tar -xzvf scrna_db_tet.tar.gz --directory scrna_db/")

  file.remove("scrna_db_tet.tar.gz")

}

download_minimal_data = function() {

  download_raw_umi_tables()
  download_embflow_scrna_db()

}

download_full_data = function() {

  download_tet_data()
  download_embflow_scrna_db()
  download_tet_scrna_db()

}


if(!dir.exists("figs")) {
  dir.create("figs")
} 
if(!dir.exists("figs/paper_figs")) {
  dir.create("figs/paper_figs")
}


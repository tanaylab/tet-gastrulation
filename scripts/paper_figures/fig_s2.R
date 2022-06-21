

generate_figure_s2_plots = function() {
  
  fig_dir = "figs/paper_figs/fig_s2"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  # fig S2a
  number_of_umis(T)
  
  # fig S2b
  number_of_cells_per_embryo(T)
  
  # fig S2c
  cell_type_comparisons_s2()
  
  # fig S2d
  # see figure 2f
  
  # fig S2f
  tko_tetraploid_vs_wt_time_distribution(T)

  # fig S2I
  plot_tko_germline_barplot_cell_type_composition()
  
  # fig S2J
  plot_tko_germline_nascent_mesosderm_vs_exe_meso_fraction()
}


number_of_umis = function(plot_pdf = F) {
  
  mat1 = scdb_mat("tko_tetra")
  mat2 = scdb_mat("control_tetra_all")
  
  f1 = mat1@cell_metadata[colnames(mat1@ignore_cmat),"embryo"] != "empty"
  f2 = mat2@cell_metadata[colnames(mat2@ignore_cmat),"embryo"] != "empty"
  
  n_umis = c(colSums(mat1@mat),colSums(mat1@ignore_cmat[,f1]))
  
  df = data.frame(umis = n_umis)
  
  p = ggplot(data = df,aes(x = umis)) +
    geom_histogram(binwidth = 500,fill = "gray",color = "black") + xlim(0,20000) + xlab("Number of UMIs / Cell") +
    ggtitle(label = "Tet TKO cells") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  if(plot_pdf) {
    fn = "figs/paper_figs/fig_s2/tetra_tko_number_of_umis_per_cell.pdf"
  } else {
    fn = "figs/paper_figs/fig_s2/tetra_tko_number_of_umis_per_cell.png"
  }
  ggsave(filename = fn,plot = p)
  
  
  
  
  n_umis = c(colSums(mat2@mat),colSums(mat2@ignore_cmat[,f2]))
  
  df = data.frame(umis = n_umis)
  
  p = ggplot(data = df,aes(x = umis)) +
    geom_histogram(binwidth = 500,fill = "gray",color = "black") + xlim(0,20000) + xlab("Number of UMIs / Cell") +
    ggtitle(label = "Ctrl cells") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  if(plot_pdf) {
    fn = "figs/paper_figs/fig_s2/tetra_ctrl_number_of_umis_per_cell.pdf"
  } else {
    fn = "figs/paper_figs/fig_s2/tetra_ctrl_number_of_umis_per_cell.png"
  }
  ggsave(filename = fn,plot = p)
  
  
}

number_of_cells_per_embryo = function(plot_pdf = F) {
  
  mat1 = scdb_mat("tko_tetra")
  mat2 = scdb_mat("control_tetra_all")
  
  n1 = table(mat1@cell_metadata[colnames(mat1@mat),"embryo"],mat1@cell_metadata[colnames(mat1@mat),"cell_type"])
  n2 = table(mat2@cell_metadata[colnames(mat2@mat),"embryo"],mat2@cell_metadata[colnames(mat2@mat),"cell_type"])
  
  all_embryos = bind_rows(as.data.frame(n1,stringsAsFactors = F),
                          as.data.frame(n2,stringsAsFactors = F))
  
  colnames(all_embryos) = c("embryo","Genotype","n_cells")
  
  f = all_embryos$Genotype != "unclear"
  all_embryos = all_embryos[f,]
  f_tko = all_embryos$Genotype == "KO"
  all_embryos[f_tko,"Genotype"] = "TKO"
  f = all_embryos$n_cells > 0
  all_embryos = all_embryos[f,]
  
  mat1_time = read.table("data/tko_tetra_wt10/time_match/time_match_summary.txt",sep = "\t",stringsAsFactors = F,h = T)
  
  mat2_time = read.table("data/control_tetra_all_wt10/time_match/time_match_summary.txt",sep = "\t",stringsAsFactors = F,h = T)
  
  transcriptional_time_all= bind_rows(rename(mat1_time[,c("embryo","best_query")],"transcriptional_rank" = "best_query"),
                                      rename(mat2_time[,c("embryo","best_query")],"transcriptional_rank" = "best_query"))
  
  rank_dev_times = read.table('data/wt10_transcriptional_rank_developmental_time.txt',sep = '\t',stringsAsFactors = F,h = T)
  
  transcriptional_time_all = left_join(transcriptional_time_all,rank_dev_times,by = "transcriptional_rank")
  
  all_embryos = left_join(all_embryos,transcriptional_time_all,by = "embryo")
  
  type_col = c("TKO" = "indianred3","control" = "gray30","host" = "black","DKO12" = "khaki","DKO13" = "lightgreen","DKO23" = "lightpink")
  
  
  f_tko = all_embryos$Genotype == "TKO"
  
  p = ggplot(data = all_embryos[f_tko,],aes(x = developmental_time,y = n_cells,color = Genotype)) +
    geom_point(size = 4) + scale_color_manual(values = type_col) + xlab("Developmental time") +
    ylab("Number of Cells") +
    theme(legend.position = "none") +
    ggtitle(label = "TKO cells per embryo") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  if(plot_pdf) {
    fn = "figs/paper_figs/fig_s2/tko_tetra_number_of_cells_per_embryo.pdf"
  } else {
    fn = "figs/paper_figs/fig_s2/tko_tetra_number_of_cells_per_embryo.png"
  }
  
  ggsave(filename = fn,plot = p)
  
  
  
  f_ctrl = all_embryos$Genotype == "control"
  
  p = ggplot(data = all_embryos[f_ctrl,],aes(x = developmental_time,y = n_cells,color = Genotype)) +
    geom_point(size = 4) + scale_color_manual(values = type_col) + xlab("Developmental time") +
    ylab("Number of Cells") +
    theme(legend.position = "none") +
    ggtitle(label = "Ctrl cells per embryo") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  if(plot_pdf) {
    fn = "figs/paper_figs/fig_s2/ctrl_tetra_number_of_cells_per_embryo.pdf"
  } else {
    fn = "figs/paper_figs/fig_s2/ctrl_tetra_number_of_cells_per_embryo.png"
  }
  
  ggsave(filename = fn,plot = p)
  
  
}


cell_type_comparisons_s2 = function() {
  
  if(!dir.exists("figs/paper_figs/fig_s2/cell_type_annotation")) {
    dir.create("figs/paper_figs/fig_s2/cell_type_annotation")
  }
  
  tko_tetra = tko_bulk_comparison(mat_id = "tko_tetra_wt10",
                                  included_cell_types = c("Epiblast","Early nascent mesoderm","Amnion/Chorion","ExE mesoderm",
                                                          "Primitive streak"),
                                  genotype = "KO",
                                  xlab_plot = "TKO",
                                  tag_fn = "tko_4N")
  
  ctrl_tetra = tko_bulk_comparison(mat_id = "control_tetra_all_wt10",
                                   included_cell_types = c("Epiblast","Early nascent mesoderm","Amnion/Chorion","ExE mesoderm",
                                                           "Primitive streak"),
                                   genotype = "control",
                                   xlab_plot = "Ctrl",
                                   tag_fn = "ctrl_4N")
  
  
  a = arrangeGrob(tko_tetra[[1]],tko_tetra[[5]],tko_tetra[[2]],tko_tetra[[3]],tko_tetra[[4]],
                  ctrl_tetra[[1]],ctrl_tetra[[5]],ctrl_tetra[[2]],ctrl_tetra[[3]],ctrl_tetra[[4]],
                  ncol = 5,nrow = 2)
  
  ggsave(filename = "figs/paper_figs/fig_s2/tko_tetra_vs_wt_per_ct.png",plot = a,w = 24,h = 10)
  ggsave(filename = "figs/paper_figs/fig_s2/tko_tetra_vs_wt_per_ct.pdf",plot = a,w = 12,h = 5) 
  
  
}

tko_bulk_comparison = function(mat_id,included_cell_types,genotype,plot_pdf = F,xlab_plot = "",tag_fn = "tko_4N") {
  
  reg = 5e-5
  reg_chi_square = 5e-5
  mat = scdb_mat(mat_id)
  
  f_tko = mat@cell_metadata[colnames(mat@mat),"cell_type"] %in% genotype
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  mat_wt = scdb_mat("sing_emb_wt10")
  
  egc = t(tgs_matrix_tapply(x = mat_wt@mat[,names(mc_wt@mc)],mc_wt@mc,sum))
  egc = t(t(egc)/colSums(egc))
  
  ct_to_col = mc_wt@color_key$color
  names(ct_to_col) = mc_wt@color_key$group
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color
  
  legc = log2(egc + 1e-5)
  
  gset = scdb_gset("sing_emb_wt10")
  bad_genes = read.table("data/tet_tko.bad_genes.txt",sep = "\t",stringsAsFactors = F)$x
  bad_genes = c("Igf2","AK145379;H19","Tet1","Tet2","Tet3")
  feat_genes = setdiff(names(gset@gene_set),bad_genes)
  
  sc_mc_cor = tgs_cor(as.matrix(mat@mat[feat_genes,f_tko]),legc[feat_genes,])
  
  best_ref = apply(sc_mc_cor,1,which.max)
  
  umi_query = t(tgs_matrix_tapply(mat@mat[,f_tko],col_to_ct[mc_wt@colors[best_ref]],sum))
  
  
  tot_umis_query = colSums(umi_query)
  egc_query = t(t(umi_query)/colSums(umi_query))
  legc_query= log2(egc_query + reg)
  
  egc_ref = t(tgs_matrix_tapply(egc[,best_ref],col_to_ct[mc_wt@colors[best_ref]],sum))
  egc_ref = t(t(egc_ref)/colSums(egc_ref))
  legc_ref = log2(egc_ref + reg)
  
  bad_genes = read.table("data/tet_tko.bad_genes.txt",sep = "\t",stringsAsFactors = F)$x
  bad_genes = c(bad_genes,c("Tet1","Tet2","Tet3"))
  genes_f = setdiff(rownames(mat@mat),bad_genes)
  

  filtered_genes = rownames(egc_query)[pmax(apply(egc_query,1,max),apply(egc_ref,1,max)) > 1e-5]

  #chi_square_stat = (egc_ref[filtered_genes,] - egc_query[filtered_genes,])^2/(egc_ref[filtered_genes,] + reg_chi_square)
  #chi_square_stat = t(t(chi_square_stat)*n_umis_query)
  
  qvalue_list = list()
  
  for (ct in included_cell_types) {
    
    p_val_vector = sapply(filtered_genes,function(gene) {
      
      test_result = chisq.test(x = c(umi_query[gene,ct],tot_umis_query[ct] - umi_query[gene,ct]),p = c(egc_ref[gene,ct],1- egc_ref[gene,ct]))
      return(test_result$p.value)
    })
    
    p_val_vector[is.na(p_val_vector)] = 1
    
    q_value_vector = qvalue(p_val_vector,pi0 = 1)$qvalues
    
    qvalue_list[[ct]] = q_value_vector
  }
  
  lfc_query_ref = legc_query[filtered_genes,] - legc_ref[filtered_genes,]

    plot_ls = list()
  for (ct in included_cell_types) {
    
    df_plot = data.frame(gene = filtered_genes,
                         expression_tko = legc_query[filtered_genes,ct],
                         expression_wt = legc_ref[filtered_genes,ct])
    
    qvalues_plot = qvalue_list[[ct]]
    
    gene_color = ifelse( (qvalues_plot >= 1e-3 )| ( abs(lfc_query_ref[,ct]) < log2(1.5) ),ct_to_col[ct],"gray30")
    
    
    f = rank(-abs(legc_query[filtered_genes,ct] - legc_ref[filtered_genes,ct])) < 4
    p = ggplot(df_plot,aes(x = expression_wt,y = expression_tko)) + 
      geom_abline(slope = 1,intercept = 0,color = "gray") +
      geom_abline(slope = 1,intercept = -1,color = "gray",linetype= "dashed") +
      geom_abline(slope = 1,intercept = 1,color= "gray",linetype = "dashed") +
      geom_point(color = gene_color,size = 1) + xlab("WT") + ylab(xlab_plot) +
      theme(axis.title = element_text(size = 20))
    #+
    # ggtitle(label = ct) + theme(plot.title = element_text(hjust = 0.5))
    plot_ls[[ct]] = p
    
    ct_name = gsub(pattern = "/",replacement = "_",x = ct)
    
    if(plot_pdf) {
      fn = sprintf("figs/paper_figs/fig_s2/cell_type_annotation/%s_%s.pdf",tag_fn,ct_name)
    } else {
      fn = sprintf("figs/paper_figs/fig_s2/cell_type_annotation/%s_%s.png",tag_fn,ct_name)
    }
    
    ggsave(fn,plot = p,w = 3.5,h= 3.5)
  }
  
  
  
  
  
  return(plot_ls)
}


tko_tetraploid_vs_wt_time_distribution = function(plot_pdf = F) {
  
  
  tko_embryos = c("4N_TKO29_12","4N_TKO26_13","4N_TKO26_12")
  
  tetra_timing = read.table('data/tko_tetra_wt10/time_match/time_match_summary.txt',sep = "\t",h = T,stringsAsFactors = F)
  f = tetra_timing$embryo %in% tko_embryos
  
  timing_f = tetra_timing[f,]
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  excluded_colors = c("#F6BFCB","#7F6874")
  mat_wt = scdb_mat("sing_emb_wt10")
  
  wt_age = read.table("data/wt10_transcriptional_rank_developmental_time.txt",sep = "\t",stringsAsFactors = F,h = T)
  
  mat_nm_tko = "tko_tetra_wt10"
  mat_nm_ctrl = "control_tetra_all_wt10"
  
  mat_tko = scdb_mat(mat_nm_tko)
  mat_ctrl = scdb_mat(mat_nm_ctrl)
  
  load(sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm_tko))
  tko_annot = cmp_annot
  load(sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm_ctrl))
  ctrl_annot = cmp_annot
  
  tko_cls = intersect(colnames(mat_tko@mat)[mat_tko@cell_metadata[colnames(mat_tko@mat),"cell_type"] == "KO"],names(tko_annot$query_cls_col))
  tko_cls = tko_cls[! tko_annot$query_cls_col[tko_cls] %in% excluded_colors]
  tko_cls = tko_cls[mat_tko@cell_metadata[tko_cls,"embryo"] %in% tko_embryos]
  
  wt10_cls = intersect(names(mc_wt@mc)[ !(mc_wt@colors[mc_wt@mc] %in% excluded_colors) ],colnames(mat_tko@mat))
  atlas_time = mat_tko@cell_metadata[wt10_cls,"transcriptional_rank"]
  names(atlas_time) = wt10_cls
  
  
  wt_cls_f = wt10_cls[mat_tko@cell_metadata[wt10_cls,"transcriptional_rank"] %in% timing_f$best_query]
  atlas_time_dist = get_atlas_time_dist(atlas_time = atlas_time,graph_id = mat_nm_tko)
  
  atlas_time_f = atlas_time_dist$atlas_time_match[wt_cls_f]
  # first timing using control cells only
  query_cls_md_tko = mat_tko@cell_metadata[tko_cls,"embryo"]
  names(query_cls_md_tko) = tko_cls
  
  query_time_dist_tko = get_query_time_dist(query_cls_md = query_cls_md_tko,atlas_time = atlas_time,graph_id = mat_nm_tko)
  
  
  df_plot = data.frame(age = c(wt_age$developmental_time[query_time_dist_tko$query_time_match],wt_age$developmental_time[atlas_time_f]),
                       genotype = c(rep("TKO",length(tko_cls)),rep("WT",length(wt_cls_f))),
                       embryo = c(query_cls_md_tko[names(query_time_dist_tko$query_time_match)],mat_tko@cell_metadata[wt_cls_f,"transcriptional_rank"]),stringsAsFactors = F)
  
  #col_values = c(rep("indianred3",length(tko_E75_embryos$embryo)),rep("gray30",length(ctrl_E75_embryos$embryo)))
  #names(col_values) = c(tko_E75_embryos$embryo,ctrl_E75_embryos$embryo)
  
  
  embryos_f = c(tko_embryos,timing_f$best_query)
  
  col_values = c(rep("indianred3",3),rep('gray30',3))
  names(col_values) = c(tko_embryos,timing_f$best_query)
  
  f = df_plot$embryo %in% tko_embryos
  
  p = ggplot(data = df_plot[f,],aes(x = age,col = embryo)) +
    geom_density(size = 0.5) + scale_color_manual(values = col_values) +
    theme(legend.position = "none") +
    xlab("Transcriptional Time") +
    ylab("")
  
  if(plot_pdf) {
    fn = "figs/paper_figs/fig_s2/tko_cls_time_distribution.pdf"
  } else {
    fn = "figs/paper_figs/fig_s2/tko_cls_time_distribution.png"
  }
  
  ggsave(filename = fn,plot = p,w = 5,h = 2)
  
  f = df_plot$embryo %in% as.character(timing_f$best_query)
  
  p = ggplot(data = df_plot[f,],aes(x = age,col = embryo)) +
    geom_density(size = 0.5) + scale_color_manual(values = col_values) +
    theme(legend.position = "none") +
    xlab("Transcriptional Time") +
    ylab("")
  
  if(plot_pdf) {
    fn = "figs/paper_figs/fig_s2/wt_cls_time_distribution.pdf"
  } else {
    fn = "figs/paper_figs/fig_s2/wt_cls_time_distribution.png"
  }
  
  ggsave(filename = fn,plot = p,w = 5,h = 2)
}



plot_tko_germline_barplot_cell_type_composition = function() {
  
  mat_nm = "tko_germline_wt10"
  df_embryos = read.table(sprintf("data/%s/time_match/time_match_summary.txt",mat_nm),sep = "\t",stringsAsFactors = F,h= T)
  rownames(df_embryos) = df_embryos$embryo
  mat = scdb_mat(mat_nm)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  col_to_rank = c(1:nrow(mc_wt@color_key))
  names(col_to_rank) = mc_wt@color_key$color
  
  excluded_colors = c("#F6BFCB","#7F6874") 
  included_colors = setdiff(unique(mc_wt@color_key$color),excluded_colors)
  
  tko_germline_embryos = df_embryos$embryo
  
  tko_germline_embryos = tko_germline_embryos[order(df_embryos[tko_germline_embryos,"best_query"])]
  
  ko_type = "TKO"
  tmp = matrix(0,nrow = length(tko_germline_embryos),ncol = length(included_colors))
  rownames(tmp) = tko_germline_embryos
  colnames(tmp) = included_colors
  
  load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm))
  
  query_cls_col = cmp_annot$query_cls_col
  query_cls = names(query_cls_col)[!(query_cls_col %in% excluded_colors)]
  query_cls = query_cls[mat@cell_metadata[query_cls,"embryo"] %in% tko_germline_embryos]
  
  filtered_cls = query_cls[mat@cell_metadata[query_cls,"cell_type"] %in%  ko_type]
  filtered_vs_ct = table(factor(x = mat@cell_metadata[filtered_cls,"embryo"],
                                levels = tko_germline_embryos),
                         factor(x = query_cls_col[filtered_cls],
                                levels = mc_wt@color_key$color[1:27]))
  
  filtered_vs_ct_n = filtered_vs_ct/rowSums(filtered_vs_ct)
  filtered_vs_ct_n[is.na(filtered_vs_ct_n)] = 0
  
  pdf("figs/paper_figs/fig_s2/tko_germline_barplot.pdf",useDingbats = F)
  barplot(t(filtered_vs_ct_n),col = colnames(filtered_vs_ct_n),las = 2,axes = F,axisnames = T)  
  dev.off()
  
}





plot_tko_germline_nascent_mesosderm_vs_exe_meso_fraction = function() {
  
  mat_nm_tko_germline = "tko_germline_wt10"
  df_germline_embryos = read.table(sprintf("data/%s/time_match/time_match_summary.txt",mat_nm_tko_germline),sep = "\t",stringsAsFactors = F,h= T)
  rownames(df_germline_embryos) = df_germline_embryos$embryo
  mat_tko_germline = scdb_mat(mat_nm_tko_germline)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  mat_wt = scdb_mat("sing_emb_wt10")
  included_colors = mc_wt@color_key$color
  
  #calculate embryo vs cell type frequency table for TKO germline embryos
  load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm_tko_germline))
  tko_germline_annotation = cmp_annot
  tko_germline_cells = names(tko_germline_annotation$query_cls_col)[mat_tko_germline@cell_metadata[names(tko_germline_annotation$query_cls_col),"embryo"] %in% df_germline_embryos$embryo]  
  
  tko_germline_embryo_vs_cell_type = compute_two_way_table(values_row = mat_tko_germline@cell_metadata[tko_germline_cells,"embryo"],
                                                           values_col = tko_germline_annotation$query_cls_col[tko_germline_cells],
                                                           included_levels_col = mc_wt@color_key$color[1:27],normalize_rows = T)
  
  
  #calculate embryo vs cell type frequency table for TKO tetraploid complementation assay embryos
  mat_nm_tko_tetra = "tko_tetra_wt10"
  mat_tko_tetra = scdb_mat(mat_nm_tko_tetra)
  df_tetra = read.table(sprintf("data/%s/time_match/time_match_summary.txt",mat_nm_tko_tetra),sep = "\t",stringsAsFactors = F,h= T)
  tetraploid_embryos = df_tetra$embryo[!is.na(df_tetra$best_query)]
  load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm_tko_tetra))
  tko_tetra_annotation = cmp_annot
  
  tko_tetra_cells = names(tko_tetra_annotation$query_cls_col)[mat_tko_tetra@cell_metadata[names(tko_tetra_annotation$query_cls_col),"cell_type"] == "KO"]
  tko_tetra_cells = tko_tetra_cells[mat_tko_tetra@cell_metadata[tko_tetra_cells,"embryo"] %in% tetraploid_embryos]
  
  tko_tetra_embryo_vs_cell_type = compute_two_way_table(values_row = mat_tko_tetra@cell_metadata[tko_tetra_cells,"embryo"],
                                                        values_col = tko_tetra_annotation$query_cls_col[tko_tetra_cells],
                                                        included_levels_col = mc_wt@color_key$color[1:27],normalize_rows = T)
  
  
  #calculate embryo vs cell type frequency table for WT embryos
  wt_embryo_vs_cell_type = compute_two_way_table(values_row = mat_wt@cell_metadata[names(mc_wt@mc),"transcriptional_rank"],
                                                 values_col = mc_wt@colors[mc_wt@mc],
                                                 included_levels_col = mc_wt@color_key$color[1:27],normalize_rows = T)
  
  
  
  exe_meso_cell_types = c(17,18,19,21,22,23)
  emb_meso_cell_types = c(8,12,13,14,15,16)
  tko_germline_emb_meso = rowSums(tko_germline_embryo_vs_cell_type[,emb_meso_cell_types])
  tko_germline_exe_meso = rowSums(tko_germline_embryo_vs_cell_type[,exe_meso_cell_types ])
  tko_tetra_emb_meso = rowSums(tko_tetra_embryo_vs_cell_type[,emb_meso_cell_types])
  tko_tetra_exe_meso = rowSums(tko_tetra_embryo_vs_cell_type[,exe_meso_cell_types ])
  wt_emb_meso = rowSums(wt_embryo_vs_cell_type[,emb_meso_cell_types])
  wt_exe_meso = rowSums(wt_embryo_vs_cell_type[,exe_meso_cell_types ])
  
  
  tko_germline_ecto = rowSums(tko_germline_embryo_vs_cell_type[,c(2,3,4,5,6)])
  tko_germline_endo = rowSums(tko_germline_embryo_vs_cell_type[,26:27])
  tko_tetra_ecto = rowSums(tko_tetra_embryo_vs_cell_type[,c(2,3,4,5,6)])
  tko_tetra_endo = rowSums(tko_tetra_embryo_vs_cell_type[,26:27])
  wt_ecto = rowSums(wt_embryo_vs_cell_type[,c(2,3,4,5,6)])
  wt_endo = rowSums(wt_embryo_vs_cell_type[,26:27])
  
  
  pdf("figs/paper_figs/fig_s2/tko_germline_embryonic_vs_nonembryonic_mesoderm.pdf",useDingbats = F)
  plot(x = wt_exe_meso,y = wt_emb_meso,pch = 19,cex = 4,col = 'gray30',xlim = c(0,0.85),ylim = c(0,0.45),
       xlab = "Nonembryonic mesoderm",ylab = "Embryonic mesoderm")
  
  points(x = tko_tetra_exe_meso ,y = tko_tetra_emb_meso,pch = 19,col = "indianred3",cex = 4)
  points(x = tko_germline_exe_meso,y = tko_germline_emb_meso,pch = 19,col = "cornflowerblue",cex = 4)
  dev.off()
  
  pdf("figs/paper_figs/fig_s2/tko_germline_ectoderm_vs_endoderm.pdf",useDingbats = F)
  plot(x = wt_endo,y = wt_ecto,pch = 19,cex = 4,col = 'gray30',xlim = c(0,0.125),ylim = c(0,0.6),
       xlab = "Embryonic endoderm",ylab = "Embryonic ectoderm")
  
  points(x = tko_tetra_endo ,y = tko_tetra_ecto,pch = 19,col = "indianred3",cex = 4)
  points(x = tko_germline_endo,y = tko_germline_ecto,pch = 19,col = "cornflowerblue",cex = 4)
  dev.off()
  
}



library("metacell")
library(dplyr)
library("Matrix")
library(ggplot2)
#scdb_init("scrna_db",force_reinit = T)
scfigs_init("figs")


generate_figure4_plots = function() {
  
  if(!dir.exists("figs/paper_figs/fig4")) {
    dir.create("figs/paper_figs/fig4")
  }
  
  # fig 4b
  plot_pdf = T
  dko_barplot_ct_frequency(mat_nm = "dko12_chim_wt10",ko_type = "DKO12",tag = "dko12_ko",show_axisnames = T,fig_width = 6,plot_pdf = plot_pdf)
  dko_barplot_ct_frequency(mat_nm = "dko12_chim_wt10",ko_type = "host",tag = "dko12_host",show_axisnames = T,fig_width = 6,plot_pdf = plot_pdf)
  dko_barplot_ct_frequency(mat_nm = "dko13_chim_wt10",ko_type = "DKO13",tag = "dko13_ko",show_axisnames = T,fig_width = 4,plot_pdf = plot_pdf)
  dko_barplot_ct_frequency(mat_nm = "dko13_chim_wt10",ko_type = "host",tag = "dko13_host",show_axisnames = T,fig_width = 4,plot_pdf = plot_pdf)
  dko_barplot_ct_frequency(mat_nm = "dko23_chim_wt10",ko_type = "DKO23",tag = "dko23_ko",show_axisnames = T,fig_width = 3,plot_pdf = plot_pdf)
  dko_barplot_ct_frequency(mat_nm = "dko23_chim_wt10",ko_type = "host",tag = "dko23_host",show_axisnames = T,fig_width = 3,plot_pdf = plot_pdf)
  
  # fig 4c
  heatmap_cell_types_query_vs_wt("TKO",T)
  heatmap_cell_types_query_vs_wt("host",T)
  heatmap_cell_types_query_vs_wt("DKO12",T)
  heatmap_cell_types_query_vs_wt("DKO13",T)
  heatmap_cell_types_query_vs_wt("DKO23",T)
  
  # fig 4d, 4e
  scatter_plots_dko()
  
  
}


dko_barplot_ct_frequency = function(mat_nm,ko_type,tag,plot_pdf = F,show_axisnames = F,fig_width = 12,fig_height= 7.5) {
  
  n_cls_min = 19
  fig_dir = "figs/paper_figs/fig4"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  fig_dir = "figs/paper_figs/fig4/dko_ct_frequencies"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  df_chim = read.table(sprintf("data/%s/time_match/time_match_summary.txt",mat_nm),sep = "\t",stringsAsFactors = F,h= T)
  rownames(df_chim) = df_chim$embryo
  mat = scdb_mat(mat_nm)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  col_to_rank = c(1:nrow(mc_wt@color_key))
  names(col_to_rank) = mc_wt@color_key$color
  
  excluded_colors = c("#F6BFCB","#7F6874") 
  included_colors = setdiff(unique(mc_wt@color_key$color),excluded_colors)
  
  
  
  chim_embryos = df_chim$embryo[(df_chim$host > n_cls_min)]
  chim_embryos = chim_embryos[order(df_chim[chim_embryos,"best_rank_host"])]
  
  tmp = matrix(0,nrow = length(chim_embryos),ncol = length(included_colors))
  rownames(tmp) = chim_embryos
  colnames(tmp) = included_colors
  
  load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm))
  
  query_cls_col = cmp_annot$query_cls_col
  query_cls = names(query_cls_col)[!(query_cls_col %in% excluded_colors)]
  query_cls = query_cls[mat@cell_metadata[query_cls,"embryo"] %in% chim_embryos]
  
  filtered_cls = query_cls[mat@cell_metadata[query_cls,"cell_type"] %in%  ko_type]
  filtered_vs_ct = table(factor(x = mat@cell_metadata[filtered_cls,"embryo"],levels = chim_embryos),factor(x = query_cls_col[filtered_cls],levels = mc_wt@color_key$color[1:27]))
  filtered_vs_ct_n = filtered_vs_ct/rowSums(filtered_vs_ct)
  filtered_vs_ct_n[is.na(filtered_vs_ct_n)] = 0
  
  if(plot_pdf) {
    pdf(sprintf("%s/barplot_ct_freq_%s.pdf",fig_dir,tag),w = fig_width,h= fig_height,useDingbats = F)
    barplot(t(filtered_vs_ct_n),col = colnames(filtered_vs_ct_n),las = 2,axes = F,axisnames = F)
    dev.off()
  } else {
    png(sprintf("%s/barplot_ct_freq_%s.png",fig_dir,tag),w = fig_width*100,h= fig_height*100)
    barplot(t(filtered_vs_ct_n),col = colnames(filtered_vs_ct_n),las = 2,axes = F,axisnames = show_axisnames)
    dev.off()
  }
  
  
  
  
}

heatmap_cell_types_query_vs_wt = function(genotype,plot_pdf = F) {
  
  gset_wt = scdb_gset("sing_emb_wt10")
  feat_genes = names(gset_wt@gene_set)
  mat_wt = scdb_mat("sing_emb_wt10")
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color
  egc_wt = t(tgs_matrix_tapply(mat_wt@mat[,names(mc_wt@mc)],col_to_ct[mc_wt@colors[mc_wt@mc]],sum))
  egc_wt = t(t(egc_wt)/colSums(egc_wt))
  
  
  
  ct_list = c(2,3,4,6,8,13,14,15,17,18,27)
  cell_types_plot = mc_wt@color_key$group[ct_list]
  
  cmp = cmp_egc_per_type_without_tetraploid(included_cts = cell_types_plot,n_downsample = 50)
  
  f_query = cmp$ko_type[colnames(cmp$egc_type)] == genotype
  
  egc_query = t(tgs_matrix_tapply(cmp$egc_type[,f_query],cmp$type_ct[f_query],mean))
  
  query_wt_cor = tgs_cor(egc_query[feat_genes,cell_types_plot],egc_wt[feat_genes,cell_types_plot])
  
  ct_annotation = data.frame(cell_type = cell_types_plot)
  rownames(ct_annotation) = cell_types_plot
  ct_to_col = mc_wt@color_key$color[ct_list]
  names(ct_to_col) = mc_wt@color_key$group[ct_list]
  annotation_colors = list(cell_type = ct_to_col)
  
  if(plot_pdf) {
    fn = sprintf("figs/paper_figs/fig4/%s_ct_similarity.pdf",genotype)
  } else {
    fn = sprintf("figs/paper_figs/fig4/%s_ct_similarity.png",genotype)
  }
  
  breaks = seq(0.55,1,length.out = 100)
  
  shades = colorRampPalette(RColorBrewer::brewer.pal(9,"PuBu"))(101)
  pheatmap::pheatmap(query_wt_cor,color = shades,cluster_rows = F,cluster_cols = F,annotation_row = ct_annotation,annotation_col = ct_annotation,annotation_colors = annotation_colors,
                     width = 3,height = 3,filename = fn,annotation_legend = F,show_rownames = F,show_colnames = F,
                     annotation_names_row = F,annotation_names_col = F,breaks = breaks)
  
  
}






scatter_plots_dko = function(plot_pdf = F) {
  
  mat_wt = scdb_mat("sing_emb_wt10")
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color 
  ct_to_col = mc_wt@color_key$color
  names(ct_to_col) = mc_wt@color_key$group
  
  egc_wt = t(tgs_matrix_tapply(mat_wt@mat[,names(mc_wt@mc)],col_to_ct[mc_wt@colors[mc_wt@mc]],sum))
  egc_wt = t(t(egc_wt)/colSums(egc_wt))
  
  #cell_types_plot = mc_wt@color_key$group[c(2,3,4,6,8,13,14,15,17,18,27)]
  #cell_types_plot = mc_wt@color_key$group[c(2,3,4,6,8,17,18)]
  # only those cell types that have more than two TKO clones
  cell_types_plot = mc_wt@color_key$group[c(2,3,4,6,14,17,18)]
  
  cmp = cmp_egc_per_type_without_tetraploid(included_cts = cell_types_plot,n_downsample = 50,downsample_umis = 190000)
  #cmp = cmp_egc_per_type_without_tetraploid_with_foxc(included_cts = cell_types_plot,n_downsample = 50,downsample_umis = 190000)
  
  reg = 6e-5
  lfp_thr = 0.8
  min_egc = 1e-4
  fig_dir = "figs/paper_figs/fig4/mean_differential_expression"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  # select differentially expressed genes
  lfp_type = log2(cmp$egc_type + reg)
  
  legc_min = apply(lfp_type,1,min)
  legc_max = apply(lfp_type,1,max)
  
  for (i in 1:ncol(cmp$egc_type)) {
    
    ct = cmp$type_ct[i]
    
    lfp_type[,i] = log2(cmp$egc_type[,i] + reg) - log2(egc_wt[,ct] + reg)
    
  }
  
  
  
  f_min = apply(cbind(cmp$egc_type,egc_wt[,cell_types_plot]),1,max) > 1e-4
  
  f_diff = legc_max - legc_min > 1.6
  
  f_diff_high = apply(lfp_type,1,max) - apply(lfp_type,1,median) > lfp_thr
  f_diff_low = apply(lfp_type,1,min) - apply(lfp_type,1,median) <  -lfp_thr
  
  genes_f = rownames(cmp$egc_type)[f_diff & f_min]
  #genes_f = rownames(cmp$egc_type)[(f_diff_high & f_min) | (f_diff_low & f_min)]
  #print(length(genes_f))
  
  genes_cts = lapply(cell_types_plot,function(ct) {
    
    return(cmp_ct_score_genes(ct_col = ct_to_col[ct],threshold = 1))
  })
  genes_all = union(genes_f,unlist(genes_cts))
  bad_genes = read.table("data/tet_tko.bad_genes.txt",sep = "\t",h = T,stringsAsFactors = F)$x
  bad_genes = c(bad_genes,c("Igf2","AK145379;H19","Polg","Slc25a4","Peg10","Igf2as","AK086477;Sh3glb1","Grb10","Nnat;Peg5"))
  genes_all = setdiff(genes_all,bad_genes)
  
  f_column = cmp$type_ct %in% cell_types_plot 
  
  df_query = as.data.frame.table(as.matrix(log2(cmp$egc_type[genes_all,f_column] + reg)),stringsAsFactors = F)
  df_wt = as.data.frame.table(as.matrix(log2(egc_wt[genes_all,cell_types_plot] + reg)),stringsAsFactors = F)
  colnames(df_query) = c("gene","type","legc_query")
  colnames(df_wt) = c("gene","cell_type","legc_wt")
  
  df_query$cell_type = cmp$type_ct[df_query$type]
  df_query$ko_type = cmp$ko_type[df_query$type]
  df_query$clone_assay = cmp$clone_assay[df_query$type]
  
  clone_assay_labels = c("Ctrl1 Chimera" = "Control 1",
                         "Ctrl2 Chimera" = "Control 2",
                         "DKO12_10 Chimera" = "DKO1/2#10",
                         "DKO12_25 Chimera" = "DKO1/2#25",
                         "DKO13_3 Chimera" = "DKO1/3",
                         "DKO23_5 Chimera" = "DKO2/3",
                         "host Chimera" = "Host",
                         "TKO23 Chimera" = "TKO3",
                         "TKO26 Chimera" = "TKO1",
                         "TKO29 Chimera" = "TKO2")
  
  #clone_assay_labels_mod = rep("",length(clone_assay_labels))
  #names(clone_assay_labels_mod) = names(clone_assay_labels)
  #clone_assay_labels = clone_assay_labels_mod
  
  df_all = left_join(x = df_query,y = df_wt,by = c("gene" = "gene","cell_type" = "cell_type"))
  
  type_col = c("TKO" = "indianred3","Ctrl" = "deepskyblue4","host" = "black","DKO12" = "khaki","DKO13" = "lightgreen","DKO23" = "lightpink","Feature Gene" = "gray50","Foxc DKO" = "cyan")
  type_col = c("TKO" = "gray50","Ctrl" = "gray50","host" = "gray50","DKO12" = "gray50","DKO13" = "gray50","DKO23" = "gray50","Feature Gene" = "gray50","Foxc DKO" = "gray50")
  
  df_text = data.frame(mean_lfp = paste0("av. lfc = ",round(colMeans(abs(lfp_type[genes_all,])),2)),
                       clone_assay = cmp$clone_assay,
                       x = rep(-13,ncol(lfp_type)),
                       y = rep(-8,ncol(lfp_type)),
                       cell_type = cmp$type_ct[colnames(lfp_type)],stringsAsFactors = F)
  
  for (ct in cell_types_plot) {
    print(ct)
    ct_genes = cmp_ct_score_genes(ct_col = ct_to_col[ct],threshold = 1)
    
    df_f = filter(df_all,df_all$cell_type == ct)
    
    df_f$ct_gene = df_f$ko_type
    df_f$ct_gene[df_f$gene %in% ct_genes] = "Feature Gene"
    
    type_levels = unique(as.character(df_f$type))
    
    
    nrow = ceiling(length(type_levels)/4)
    
    #df_f$ko_type = factor(x = df_f$ko_type,levels = c("TKO","DKO12","Ctrl","DKO13","DKO23","host"))
    
    
    df_text_f = df_text[df_text$cell_type == ct,]
    
    p = ggplot(data = df_f,aes(x = legc_wt,y = legc_query)) + 
      geom_abline(intercept = 0,slope = 1,linetype = "solid",color = "gray") +
      geom_abline(intercept = -1,slope = 1,linetype = "dashed",color = "gray") +
      geom_abline(intercept = 1,slope = 1,linetype = "dashed",color = "gray") +
      geom_point(aes(color = ct_gene),show.legend = F) +
      xlab("") +
      ylab(paste0("")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      facet_wrap(clone_assay ~ .,
                 nrow = ceiling(length(type_levels)/4),
                 labeller = labeller(clone_assay = clone_assay_labels)) +
      scale_color_manual(values = type_col) +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
    #geom_text(data = df_text_f,mapping = aes(x = x,y = y,label = mean_lfp),color = "black",hjust =0)
    #geom_point(data = filter(df_f,df_f$ct_gene == "Feature Gene"),aes(color = ct_gene),show.legend = F) +
    
    ct_name = gsub("/", " ",ct)
    if(plot_pdf) {
      fn = sprintf("%s/%s_combined.pdf",fig_dir,ct_name)
    } else {
      fn = sprintf("%s/%s_combined.png",fig_dir,ct_name)
    }
    ggsave(filename = fn,plot = p,height = ceiling(length(type_levels)/4)*3.8,width = 12)
    
  }
  
  lfp_type = lfp_type[,f_column]
  
  type_mean_lfp = colMeans(abs(lfp_type[genes_all,]))
  
  ct_to_rank = c(1:nrow(mc_wt@color_key))
  names(ct_to_rank) = mc_wt@color_key$group
  
  df_hm = data.frame(mean_dist = type_mean_lfp,cell_type = cmp$type_ct[colnames(lfp_type)],clone_assay = cmp$clone_assay[colnames(lfp_type)],ko_type = cmp$ko_type[colnames(lfp_type)])
  
  #df_hm$cell_type[df_hm$cell_type == "Lateral & intermediate mesoderm"] = "Lat. & interm. mesoderm"
  type_col = c("TKO" = "indianred3","Ctrl" = "deepskyblue4","host" = "black","DKO12" = "khaki","DKO13" = "lightgreen","DKO23" = "lightpink","Feature Gene" = "gray50","Foxc DKO" = "cyan")
  
  df_hm = rename(df_hm,Genotype = ko_type)
  
  # manually adjusting a few values for the graphical layout
  #df_hm$mean_dist[df_hm$cell_type == "Surface ectoderm" & df_hm$clone_assay == "DKO12_25 Chimera"] = df_hm$mean_dist[df_hm$cell_type == "Surface ectoderm" & df_hm$clone_assay == "DKO12_25 Chimera"] + 0.002
  #df_hm$mean_dist[df_hm$cell_type == "Surface ectoderm" & df_hm$clone_assay == "DKO12_10 Chimera"] = df_hm$mean_dist[df_hm$cell_type == "Surface ectoderm" & df_hm$clone_assay == "DKO12_10 Chimera"] - 0.0015
  #df_hm$mean_dist[df_hm$cell_type == "Rostral neural plate" & df_hm$clone_assay == "DKO12_10 Chimera"] = df_hm$mean_dist[df_hm$cell_type == "Rostral neural plate" & df_hm$clone_assay == "DKO12_10 Chimera"] - 0.0015
  #df_hm$mean_dist[df_hm$cell_type == "Rostral neural plate" & df_hm$Genotype == "DKO13"] = df_hm$mean_dist[df_hm$cell_type == "Rostral neural plate" & df_hm$Genotype == "DKO13"] - 0.002
  #df_hm$mean_dist[df_hm$cell_type == "ExE mesoderm" & df_hm$clone_assay == "DKO12_25 Chimera"] = df_hm$mean_dist[df_hm$cell_type == "ExE mesoderm" & df_hm$clone_assay == "DKO12_25 Chimera"] + 0.003
  #df_hm$mean_dist[df_hm$cell_type == "Amnion/Chorion" & df_hm$Genotype == "host"] = df_hm$mean_dist[df_hm$cell_type == "Amnion/Chorion" & df_hm$Genotype == "host"] + 0.002
  
  
  
  p2 = ggplot(data = df_hm,aes(y = cell_type,x = mean_dist,fill = Genotype)) +
    geom_point(shape = 21,size = 2) + scale_fill_manual(values = type_col) +
    xlab("Mean absolute log2 fold change") + ylab("") +
    theme(axis.text.y = element_text(size=13))
  
  if (plot_pdf) {
    fn = sprintf("%s/mean_lfc.pdf",fig_dir)
  } else {
    fn = sprintf("%s/mean_lfc.png",fig_dir)
  }
  
  ggsave(filename = fn,plot = p2,width = 8,height = 4)
  
  
}


cmp_egc_per_type_without_tetraploid = function(included_cts,n_downsample,seed_downsample = 123456,downsample_umis = NULL) {
  
  mat_chim = scdb_mat("tko_chim")
  mat_dko12 = scdb_mat("dko12_chim")
  mat_dko13 = scdb_mat("dko13_chim")
  mat_dko23 = scdb_mat("dko23_chim")
  #chim_cls = intersect(colnames(mat_chim@mat),colnames(mat_chim2@mat))
  #mat_tetra_tko = scdb_mat("tko_tetra")
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  mat_wt = scdb_mat("sing_emb_wt10")
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color
  ct_to_col = mc_wt@color_key$color
  names(ct_to_col) = mc_wt@color_key$group
  
  included_cts = mc_wt@color_key$group[c(1:27)]
  
  mat_all = cbind(mat_chim@mat,mat_dko12@mat,mat_dko13@mat,mat_dko23@mat)
  
  md_all = bind_rows(mat_chim@cell_metadata[colnames(mat_chim@mat),],
                     mat_dko12@cell_metadata[colnames(mat_dko12@mat),],
                     mat_dko13@cell_metadata[colnames(mat_dko13@mat),],
                     mat_dko23@cell_metadata[colnames(mat_dko23@mat),])
  
  rownames(md_all) = md_all$cell
  
  load("data/tko_chim_wt10/color_annotation/cmp_annot.Rda")
  tko_chim_annot = cmp_annot
  load("data/dko12_chim_wt10/color_annotation/cmp_annot.Rda")
  dko12_chim_annot = cmp_annot
  load("data/dko13_chim_wt10/color_annotation/cmp_annot.Rda")
  dko13_chim_annot = cmp_annot
  load("data/dko23_chim_wt10/color_annotation/cmp_annot.Rda")
  dko23_chim_annot = cmp_annot
  
  clone_assay = paste0(md_all$clone_type,c(rep(" Chimera",ncol(mat_chim@mat)),rep(" Chimera",ncol(mat_dko12@mat) + ncol(mat_dko13@mat) +ncol(mat_dko23@mat))))
  names(clone_assay) = c(colnames(mat_chim@mat),colnames(mat_dko12@mat),colnames(mat_dko13@mat),colnames(mat_dko23@mat))
  
  query_cls_col = c(tko_chim_annot$query_cls_col,dko12_chim_annot$query_cls_col,dko13_chim_annot$query_cls_col,dko23_chim_annot$query_cls_col)
  
  cls_a = intersect(names(query_cls_col),colnames(mat_all))
  
  df_cls = data.frame(cell = cls_a,cell_type = col_to_ct[query_cls_col[cls_a]],clone_assay = clone_assay[cls_a],stringsAsFactors = F)
  
  df_cls$clone_cell_type = paste(df_cls$cell_type,df_cls$clone_assay,sep = " ")  
  
  df_cls = left_join(df_cls,rename(mc_wt@color_key[,-3],cell_type = group,cell_type_color = color),by = "cell_type")
  
  
  f_ct = df_cls$cell_type %in% included_cts
  df_cls = df_cls[f_ct,]
  
  n_cls_clone_cell_type = table(df_cls$clone_cell_type)
  
  included_types = names(n_cls_clone_cell_type)[n_cls_clone_cell_type >= n_downsample]
  
  df_cls = df_cls[df_cls$clone_cell_type %in% included_types,]
  
  clone_assay_levels = unique(df_cls$clone_assay) 
  
  # next downsample cells
  set.seed(seed_downsample)
  ind_ds = unlist(tapply(c(1:nrow(df_cls)),df_cls$clone_cell_type,function(v) {
    return(sample(x = v,size = n_downsample,replace = F))
  }))
  df_cls = df_cls[ind_ds,]
  
  egc_type = t(tgs_matrix_tapply(mat_all[,df_cls$cell],df_cls$clone_cell_type,sum))
  n_umis_per_type = colSums(egc_type)
  
  if(!is.null(downsample_umis)) {
    egc_type = scm_downsamp(egc_type,n = downsample_umis)
  }
  egc_type = t(t(egc_type)/colSums(egc_type))
  
  type_ct = rep("unclear",ncol(egc_type))
  for (ct in mc_wt@color_key$group) {
    
    type_ct[grep(ct,colnames(egc_type))] = ct
  }
  names(type_ct) = colnames(egc_type)
  
  clone_assay = rep("unclear",ncol(egc_type))
  for (cl_a in clone_assay_levels) {
    
    clone_assay[grep(cl_a,colnames(egc_type))] = cl_a
  }
  names(clone_assay) = colnames(egc_type)
  
  
  
  ko_type = rep("unclear",ncol(egc_type))
  for (a in c("TKO","Ctrl","host","DKO12","DKO13","DKO23")) {
    ko_type[grep(a,colnames(egc_type))] = a
  }
  names(ko_type) = colnames(egc_type)
  
  experiment_type = rep("unclear",ncol(egc_type))
  for (a in c("Chimera","Tetraploid")) {
    experiment_type[grep(a,colnames(egc_type))] = a
  }
  names(experiment_type) = colnames(egc_type)
  
  return(list(egc_type = egc_type,type_ct = type_ct,ko_type = ko_type,experiment_type = experiment_type,clone_assay = clone_assay))
}


cmp_ct_score_genes = function(ct_col,threshold) {
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  
  lfp_ct = t(tgs_matrix_tapply(log2(mc_wt@mc_fp),mc_wt@colors,mean))
  
  f = lfp_ct[,ct_col] > threshold
  return(rownames(lfp_ct)[f])
}



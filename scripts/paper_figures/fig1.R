

generate_figure1_plots = function() {
  
  fig_dir = "figs/paper_figs/fig1"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  # fig 1c
  atlas_projection_tko_tetraploid()
  atlas_projection_control_tetraploid()
  
  # fig 1d
  barplot_tetraploid_ct_frequency("tko_tetra_wt10",plot_pdf = T)
  barplot_tetraploid_ct_frequency("control_tetra_all_wt10",plot_pdf = T)
  
  # fig 1e
  fig1_frequencies_per_germ_layer(plot_pdf = T)
  
  # fig 1g
  var_vs_mean_tko_ctrl()

}


fig1_frequencies_per_germ_layer = function(plot_pdf = F) {
  
  fn_prefix = ""
  plot_wt = F
  plot_ko = T
  plot_control = T
  
  
  #ko_color = "#a30000"
  ko_color = "indianred3"
  control_color = "gray30"
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  
  plot_frequency_group_of_cell_types(highlighted_colors = mc_wt@color_key$color[2:6],main_tag = paste0(fn_prefix,"Ectoderm"),ko_color = ko_color,control_color = control_color,
                                     plot_pdf = plot_pdf,plot_wt = plot_wt,plot_control = plot_control,plot_ko = plot_ko)
  plot_frequency_group_of_cell_types(highlighted_colors = mc_wt@color_key$color[26:27],main_tag = paste0(fn_prefix,"Endoderm"),ko_color = ko_color,control_color = control_color,
                                     plot_pdf = plot_pdf,plot_wt = plot_wt,plot_control = plot_control,plot_ko = plot_ko)
  plot_frequency_group_of_cell_types(highlighted_colors = mc_wt@color_key$color[c(8,12,13,14,15,16)],
                                     main_tag = paste0(fn_prefix,"Embryonic mesoderm"),ko_color = ko_color,control_color = control_color,
                                     plot_pdf = plot_pdf,plot_wt = plot_wt,plot_control = plot_control,plot_ko = plot_ko)
  plot_frequency_group_of_cell_types(highlighted_colors = mc_wt@color_key$color[c(17,18,19)],
                                     main_tag = paste0(fn_prefix,"Extraembryonic mesoderm"),ko_color = ko_color,control_color = control_color,
                                     plot_pdf = plot_pdf,plot_wt = plot_wt,plot_control = plot_control,plot_ko = plot_ko)

  

}

plot_frequency_group_of_cell_types = function(highlighted_colors,main_tag = "",plot_pdf = F,ko_color = "indianred2",control_color = "gray30",plot_wt = F,plot_control = F,plot_ko = T) {
  
  mat_nm_ko = "tko_tetra_wt10"
  mat_nm_control = "control_tetra_all_wt10"
  
  rank_to_time = read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt",stringsAsFactors = F,h = T,sep = "\t")
  dev_time = rank_to_time$developmental_time
  age_field_ko = "best_query"
  age_field_control = "best_query"
  
  if(!dir.exists("figs/paper_figs")) {
    dir.create("figs/paper_figs")
  }
  if(!dir.exists("figs/paper_figs/fig1")) {
    dir.create("figs/paper_figs/fig1")
  }
  
  fig_dir = "figs/paper_figs/fig1/cell_types"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  roll_width = 4
  
  load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm_ko))
  ko_query_cls_col = cmp_annot$query_cls_col
  ko_query_cls = names(ko_query_cls_col)
  
  load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm_control))
  control_query_cls_col = cmp_annot$query_cls_col
  control_query_cls = names(control_query_cls_col)
  
  mat_ko=  scdb_mat(mat_nm_ko)
  mat_control = scdb_mat(mat_nm_control)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  cgraph_ko = scdb_cgraph(mat_nm_ko)
  cgraph_control = scdb_cgraph(mat_nm_control)
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color
  excluded_colors = c("#F6BFCB","#7F6874")
  included_colors = setdiff(mc_wt@color_key$color,excluded_colors)
  
  df_ko = read.table(sprintf("data/%s/time_match/time_match_summary.txt",mat_nm_ko),sep = "\t",h = T,stringsAsFactors = F)
  rownames(df_ko) = df_ko$embryo
  df_control = read.table(sprintf("data/%s/time_match/time_match_summary.txt",mat_nm_control),sep = "\t",h = T,stringsAsFactors = F)
  rownames(df_control) = df_control$embryo
  
  ko_embryos = df_ko$embryo
  ko_embryos = ko_embryos[order(df_ko[ko_embryos,"best_query"])]
  control_embryos = df_control$embryo
  control_embryos = control_embryos[order(df_control[control_embryos,"best_query"])]
  
  ko_query_cls_f = ko_query_cls[( mat_ko@cell_metadata[ko_query_cls,"cell_type"] %in%  c("KO","control","host") ) & ( ko_query_cls_col[ko_query_cls] %in% included_colors )]
  control_query_cls_f = control_query_cls[( mat_control@cell_metadata[control_query_cls,"cell_type"] %in%  c("KO","control","host") ) & ( control_query_cls_col[control_query_cls] %in% included_colors )]
  wt10_cls = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] %in% included_colors]
  
  emb_wt_age = unique(mat_ko@cell_metadata[wt10_cls,c("transcriptional_rank","age_group")])
  emb_wt_age = emb_wt_age[order(emb_wt_age$transcriptional_rank),]
  
  wt10_emb_vs_ct = table(mat_ko@cell_metadata[wt10_cls,"transcriptional_rank"],mc_wt@colors[mc_wt@mc[wt10_cls]])
  ko_emb_vs_ct = table(mat_ko@cell_metadata[ko_query_cls_f,"embryo"],ko_query_cls_col[ko_query_cls_f])
  ko_emb_vs_ct = ko_emb_vs_ct[ko_embryos,]
  control_emb_vs_ct = table(mat_control@cell_metadata[control_query_cls_f,"embryo"],control_query_cls_col[control_query_cls_f])
  control_emb_vs_ct = control_emb_vs_ct[control_embryos,]
  
  
  wt_fr = rowSums(wt10_emb_vs_ct[,intersect(highlighted_colors,colnames(wt10_emb_vs_ct)),drop = FALSE])/rowSums(wt10_emb_vs_ct)
  ko_fr = rowSums(ko_emb_vs_ct[,intersect(highlighted_colors,colnames(ko_emb_vs_ct)),drop = FALSE])/rowSums(ko_emb_vs_ct)  
  control_fr = rowSums(control_emb_vs_ct[,intersect(highlighted_colors,colnames(control_emb_vs_ct)),drop = FALSE])/rowSums(control_emb_vs_ct)  
 
  # next calculate moving  90% moving average window for every cell type
  n = length(wt_fr)
  freq_ct = c(wt_fr,wt_fr[(n - roll_width + 1):n])
  mov_mean = rollmean(x = freq_ct,k = 2*roll_width+1)
  
  mov_sd =  rollapply(data = freq_ct,width = 2*roll_width+1,sd)
  
  upper_sd = mov_mean + 2*mov_sd
  lower_sd = mov_mean - 2*mov_sd
  
  names(upper_sd) = c((roll_width+1):153)
  names(lower_sd) = c((roll_width+1):153)
  names(mov_mean) = c((roll_width+1):153)
  
  #min_rank_plot = min(df_chim[tetra_embryos,c("best_query")]) - 5
  min_rank_plot = 10
  xlim_min = dev_time[min_rank_plot]
  xlim_max = dev_time[153] + 0.05
  
  #x_ranks = c(min_rank_plot:(153-roll_width))
  x_ranks = c(min_rank_plot:(153))
  upper_freq = upper_sd[as.character(c(min_rank_plot:(153)))]
  lower_freq = lower_sd[as.character(c(min_rank_plot:(153)))]
  
  
  ylim_max = max(ko_fr,control_fr,wt_fr)
  
  # next calculate p value between query and reference
  
  min_age_rank_test = 109
  
  f_control = df_control[names(control_fr),age_field_control] >= min_age_rank_test
  f_ko = df_ko[names(ko_fr),age_field_ko] >= min_age_rank_test
  w_test = wilcox.test(x = ko_fr[f_ko],y = control_fr[f_control])
  
  
  if(plot_pdf) {
    pdf(sprintf("%s/%s.pdf",fig_dir,main_tag),w = 7,h = 5.6,useDingbats = F)
    par(mar = c(6,6,6,4))
    
    plot(NULL,ylim = c(-0.05*ylim_max,1.2*ylim_max),
         main = main_tag,ylab = "Fraction",xlab = "Developmental Time",cex.main = 3,
         cex.lab = 2,cex.axis = 2,xaxs = 'i',yaxs= 'i',xlim = c(xlim_min,xlim_max))
    #plot(x = dev_time[df_ko[names(ko_fr),age_field_ko]],y = ko_fr,ylim = c(-0.05*ylim_max,1.2*ylim_max),
    #     pch = 19,cex = 3,col = control_color,main = main_tag,ylab = "Fraction",xlab = "Developmental Time",cex.main = 3,
    #     cex.lab = 2,cex.axis = 2,xaxs = 'i',yaxs= 'i',xlim = c(xlim_min,xlim_max))
    polygon(x = c(dev_time[x_ranks],dev_time[rev(x_ranks)]),y = c(upper_freq,rev(lower_freq)),col = "gray80",border = NA)
    lines(x = dev_time[c((roll_width + 1):153)],y = mov_mean,lwd = 3)
    

    if (plot_wt) {
      points(x = dev_time[as.numeric(names(wt_fr))],y = wt_fr,pch = 19,cex = 0.3)
    }
    if (plot_control) {
      points(x = dev_time[df_control[names(control_fr),age_field_control]],y = control_fr, pch = 19,cex = 3,col = control_color)
    }
    if (plot_ko) {
      points(x = dev_time[df_ko[names(ko_fr),age_field_ko]],y = ko_fr, pch = 19,cex = 3,col = ko_color)
    }
    text(x = 7,y = 0.9*max(control_fr,ko_fr,wt_fr),labels = sprintf("p = %.1e",w_test$p.value),cex = 2)
    #legend(x = "topleft",legend = c("wt","ko","control),pch = 19,col = c("black","#83c26d","cornflowerblue"))
    
    dev.off()
    
  } else {
    
    png(sprintf("%s/%s.png",fig_dir,main_tag),w = 600,h = 450)
    par(mar = c(6,6,6,4))
    

    plot(NULL,ylim = c(-0.05*ylim_max,1.2*ylim_max),
         main = main_tag,ylab = "Fraction",xlab = "Developmental Time",cex.main = 3,
         cex.lab = 2,cex.axis = 2,xaxs = 'i',yaxs= 'i',xlim = c(xlim_min,xlim_max))
    #plot(x = dev_time[df_ko[names(ko_fr),age_field_ko]],y = ko_fr,ylim = c(-0.05*ylim_max,1.2*ylim_max),
    #     pch = 19,cex = 3,col = control_color,main = main_tag,ylab = "Fraction",xlab = "Developmental Time",cex.main = 2,
    #     cex.lab = 2,xaxs = 'i',yaxs= 'i',xlim = c(xlim_min,xlim_max))
    polygon(x = c(dev_time[x_ranks],dev_time[rev(x_ranks)]),y = c(upper_freq,rev(lower_freq)),col = "gray80",border = NA)
    lines(x = dev_time[c((roll_width + 1):153)],y = mov_mean,lwd = 3)
    
    if (plot_wt) {
      points(x = dev_time[as.numeric(names(wt_fr))],y = wt_fr,pch = 19,cex = 0.3)
    }

    if (plot_control) {
      points(x = dev_time[df_control[names(control_fr),age_field_control]],y = control_fr, pch = 19,cex = 4,col = control_color)
    }
    if (plot_ko) {
      points(x = dev_time[df_ko[names(ko_fr),age_field_ko]],y = ko_fr, pch = 19,cex = 4,col = ko_color)
    }

    
    #legend(x = "topleft",legend = c("wt","ko","control"),pch = 19,col = c("black","#83c26d","cornflowerblue"))
    text(x = 7,y = 0.9*max(control_fr,ko_fr,wt_fr),labels = sprintf("p = %.1e",w_test$p.value),cex = 2)
    
    dev.off()
  }


}


atlas_projection_tko_tetraploid = function(plot_pdf = F) {
  
  mat_tetra = scdb_mat("tko_tetra_ko")
  
  gset = scdb_gset("tko_tetra_ko")
  feat_genes = names(gset@gene_set)
  
  mat_query = mat_tetra@mat
  
  fig_dir = "figs/paper_figs/fig1"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  fn = "figs/paper_figs/fig1/atlas_projection_tko_tetraploid_all_embryos_new_12.png"
  w = 1000
  h = 1000
  if(plot_pdf) {
    fn = gsub(pattern = ".png",replacement = ".pdf",x = fn)
    w = 1000/72
    h = 1000/72
  }
  #before: cex_points = 0.7
  atlas_proj_on_wt10(mat_query = mat_query,feat_genes = feat_genes,fn = fn,cex_points = 1.2,plot_pdf = plot_pdf,w = w,h = h,plot_gray_background = F)

}

atlas_projection_control_tetraploid = function(plot_pdf = F) {
  
  mat_tetra = scdb_mat("control_tetra_all")
  
  gset = scdb_gset("control_tetra_all")
  feat_genes = names(gset@gene_set)
  
  mat_query = mat_tetra@mat
  
  f = mat_tetra@cell_metadata[colnames(mat_tetra@mat),"cell_type"] == "control"
  mat_query = mat_query[,f]

  # filter host and unclear cells from embflow tetraploid control embryos that
  mat_control_embflow = scdb_mat("control_tetra_from_embflow_paper")
  cells_f = intersect(colnames(mat_control_embflow@mat),colnames(mat_query))
  cells_ignore = cells_f[mat_control_embflow@cell_metadata[cells_f,"cell_type"] != "control"]

  f = colnames(mat_query) %in% cells_ignore
  mat_query = mat_query[,!f]

  fig_dir = "figs/paper_figs/fig1"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  fn = "figs/paper_figs/fig1/atlas_projection_control_tetraploid_all_embryos.png"
  w = 1000
  h = 1000
  if(plot_pdf) {
    fn = gsub(pattern = ".png",replacement = ".pdf",x = fn)
    w = 1000/72
    h = 1000/72
  }
  
  atlas_proj_on_wt10(mat_query = mat_query,feat_genes = feat_genes,fn = fn,cex_points = 1.2,plot_pdf = plot_pdf,w = w,h = h,plot_gray_background = F)


}


barplot_tetraploid_ct_frequency = function(mat_nm,plot_pdf = F) {

  fig_dir = "figs/paper_figs/fig1"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  
  df_tetra = read.table(sprintf("data/%s/time_match/time_match_summary.txt",mat_nm),sep = "\t",stringsAsFactors = F,h= T)
  rownames(df_tetra) = df_tetra$embryo
  tetra_embryos = df_tetra$embryo
  mat = scdb_mat(mat_nm)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  col_to_rank = c(1:nrow(mc_wt@color_key))
  names(col_to_rank) = mc_wt@color_key$color
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color
  
  excluded_colors = c("#F6BFCB","#7F6874") 
  included_colors = setdiff(unique(mc_wt@color_key$color),excluded_colors)
  
  tetra_embryos = tetra_embryos[order(df_tetra[tetra_embryos,"best_query"])]
  
  tmp = matrix(0,nrow = length(tetra_embryos),ncol = length(included_colors))
  rownames(tmp) = tetra_embryos
  colnames(tmp) = included_colors
  
  load(file = sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm))
  
  query_cls_col = cmp_annot$query_cls_col
  query_cls = names(query_cls_col)[!(query_cls_col %in% excluded_colors)]
  query_cls = query_cls[mat@cell_metadata[query_cls,"embryo"] %in% tetra_embryos]
  
  modified_cell_type_levels = c(mc_wt@color_key$group[c(1,7,9,10,11,20,24,25)],c("space1"),
                                mc_wt@color_key$group[c(2:6)],c("space2"),
                                mc_wt@color_key$group[c(8,12,13,14,15,16)],c("space3"),
                                mc_wt@color_key$group[c(17,18,19)],c("space4"),
                                mc_wt@color_key$group[c(21:23)],c("space5"),
                                mc_wt@color_key$group[c(26:27)])
  
  modified_colors = c(mc_wt@color_key$color[c(1,7,9,10,11,20,24,25)],c("white"),
                      mc_wt@color_key$color[c(2:6)],c("white"),
                      mc_wt@color_key$color[c(8,12,13,14,15,16)],c("white"),
                      mc_wt@color_key$color[c(17,18,19)],c("white"),
                      mc_wt@color_key$color[c(21:23)],c("white"),
                      mc_wt@color_key$color[c(26:27)])
  
  filtered_cls = query_cls[mat@cell_metadata[query_cls,"cell_type"] %in%  c("KO","control")]
  filtered_vs_ct = table(factor(x = mat@cell_metadata[filtered_cls,"embryo"],levels = tetra_embryos),factor(x = col_to_ct[query_cls_col[filtered_cls]],levels = modified_cell_type_levels))
  filtered_vs_ct_n = filtered_vs_ct/rowSums(filtered_vs_ct)
  filtered_vs_ct_n[is.na(filtered_vs_ct_n)] = 0
  

  filtered_vs_ct_n[,c("space1","space2","space3","space4","space5")] = 0

  filtered_vs_ct_n = filtered_vs_ct_n/rowSums(filtered_vs_ct_n)
  
  if(plot_pdf) {
    pdf(sprintf("%s/barplot_ct_freq_%s.pdf",fig_dir,mat_nm),w = 12,h= 7.5,useDingbats = F)
    barplot(t(filtered_vs_ct_n),col = modified_colors,las = 2,axes = F,axisnames = F,border = NA)
    dev.off()
  } else {
    png(sprintf("%s/barplot_ct_freq_%s.png",fig_dir,mat_nm),w = 1200,h= 750)
    barplot(t(filtered_vs_ct_n),col = modified_colors,las = 2,axes = F,axisnames = F)
    dev.off()
  }

}

var_vs_mean_tko_ctrl = function() {
  
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  excluded_colors = c("#F6BFCB","#7F6874")
  mat_wt = scdb_mat("sing_emb_wt10")
  
  wt_age = read.table("data/wt10_transcriptional_rank_developmental_time.txt",sep = "\t",stringsAsFactors = F,h = T)
  
  tko_4n_time = read.table("data/tko_tetra_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  ctrl_4n_time = read.table("data/control_tetra_all_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  
  tko_embryos = tko_4n_time$embryo
  ctrl_embryos = ctrl_4n_time$embryo
  
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
  
  ctrl_cls = intersect(colnames(mat_ctrl@mat)[mat_ctrl@cell_metadata[colnames(mat_ctrl@mat),"cell_type"] == "control"],names(ctrl_annot$query_cls_col))
  ctrl_cls = ctrl_cls[! ctrl_annot$query_cls_col[ctrl_cls] %in% excluded_colors]
  ctrl_cls = ctrl_cls[mat_ctrl@cell_metadata[ctrl_cls,"embryo"] %in% ctrl_embryos]
  
  wt10_cls = intersect(names(mc_wt@mc)[ !(mc_wt@colors[mc_wt@mc] %in% excluded_colors) ],colnames(mat_tko@mat))
  atlas_time = mat_tko@cell_metadata[wt10_cls,"transcriptional_rank"]
  names(atlas_time) = wt10_cls
  
  
  # first timing using control cells only
  
  atlas_time_dist = get_atlas_time_dist(atlas_time = atlas_time,graph_id = mat_nm_tko)
  
  query_cls_md_tko = mat_tko@cell_metadata[tko_cls,"embryo"]
  names(query_cls_md_tko) = tko_cls
  
  query_time_dist_tko = get_query_time_dist(query_cls_md = query_cls_md_tko,atlas_time = atlas_time,graph_id = mat_nm_tko)
  
  query_cls_md_ctrl = mat_ctrl@cell_metadata[ctrl_cls,"embryo"]
  names(query_cls_md_ctrl) = ctrl_cls
  
  query_time_dist_ctrl = get_query_time_dist(query_cls_md = query_cls_md_ctrl,atlas_time = atlas_time,graph_id = mat_nm_ctrl)
  
  ctrl_dist_n = query_time_dist_ctrl$query_time_dist/rowSums(query_time_dist_ctrl$query_time_dist)
  ko_dist_n = query_time_dist_tko$query_time_dist/rowSums(query_time_dist_tko$query_time_dist)
  wt_dist_n = atlas_time_dist$atlas_time_dist/rowSums(atlas_time_dist$atlas_time_dist)
  
  ctrl_mean = (ctrl_dist_n %*% wt_age$developmental_time)[,1]
  ko_mean = (ko_dist_n %*% wt_age$developmental_time)[,1]
  wt_mean = (wt_dist_n %*% wt_age$developmental_time)[,1]
  
  ctrl_var = apply(query_time_dist_ctrl$query_time_dist,1,function(v) {
    
    time_samples = rep(c(1:153),v)
    
    time_var = var(wt_age$developmental_time[time_samples])
    return(time_var)
  })
  
  ko_var = apply(query_time_dist_tko$query_time_dist,1,function(v) {
    
    time_samples = rep(c(1:153),v)
    
    time_var = var(wt_age$developmental_time[time_samples])
    return(time_var)
  })
  
  wt_var = apply(atlas_time_dist$atlas_time_dist,1,function(v) {
    
    time_samples = rep(c(1:153),v)
    
    time_var = var(wt_age$developmental_time[time_samples])
    return(time_var)
  })

  rollwidth = 8
  wt_var = wt_var[order(wt_mean)]
  wt_mean = wt_mean[order(wt_mean)]

  wt_var_roll_mean = rollmean(wt_var,2*rollwidth + 1)

  wt_var_roll_sd = rollapply(wt_var,2*rollwidth + 1,sd)

  pdf('figs/paper_figs/fig1/4n_embryos_variance_vs_mean_all_4N_embryos.pdf',useDingbats = F)
  plot(x = wt_mean[(rollwidth + 1):(length(wt_mean)-rollwidth)],y = wt_var_roll_mean,type = "l",lwd = 2,
       ylim = c(0,0.16),xlim = c(6.85,7.92),
       xlab = "Mean Et",ylab = "Variance Et")
  polygon(x = c(wt_mean[(rollwidth + 1):(length(wt_mean)-rollwidth)],rev(wt_mean[(rollwidth + 1):(length(wt_mean)-rollwidth)])),
          y = c(wt_var_roll_mean + 2*wt_var_roll_sd,rev(wt_var_roll_mean - 2*wt_var_roll_sd)),col = "gray80",border = NA)
  lines(x = wt_mean[(rollwidth + 1):(length(wt_mean)-rollwidth)],y = wt_var_roll_mean,type = "l",lwd = 2)
  points(x = ctrl_mean,y = ctrl_var,pch = 19,col = "gray30",cex = 3)
  points(x = ko_mean,y = ko_var,pch = 19,col = "indianred3",cex = 3)
  dev.off()



}


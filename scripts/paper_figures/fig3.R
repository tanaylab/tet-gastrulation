# figure 3 plots
library(ggplot2)
library(zoo)
library("metacell")
library("tgstat")
library("Matrix")
#source("scripts/differential_expression/cell_type_scores.R")
#source("scripts/paper_figs/mc_gg_plots.R")

generate_figure3_plots = function() {
  
  if(!dir.exists("figs/paper_figs/fig3")) {
    dir.create("figs/paper_figs/fig3")
  }
  
  # fig 3a
  plot_epiblast_score_distribution (T)
  
  # fig 3b
  # jupyter notebook heatmap_projection_epiblast.ipynb
  
  
  # fig 3c
  plot_epiblast_gene_expression_vs_time(genes_f = c("Dppa4","Gdf3"),T)
  
  # fig 3d
  plot_nascent_mesoderm_score_distribution(T)
  

  # fig 3e
  # jupyter notebook heatmap_projection_nascent_meso.ipynb
  
  # fig 3f
  plot_nascent_mesoderm_gene_expression_vs_time(genes_f = c("Lefty2","Fgf15","Spry4","Dll1"),T)

}

epiblast_genes = function() {
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  legc = log2(mc@e_gc + 1e-5)
  
  gg_cor = tgs_cor(t(legc))
  
  genes_f = colnames(gg_cor)[gg_cor["Utf1",] > 0.87]
  genes_f = setdiff(genes_f,c("Dnmt3b"))
  return(genes_f)
}

early_nascent_mesoderm_genes = function() {
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  egc = t(tgs_matrix_tapply(mat@mat[,names(mc@mc)],mc@mc,sum))
  egc = t(t(egc)/colSums(egc))
  mc_incl = c(449:456)
  
  legc = log2(egc + 5e-5)
  
  lfc = legc - apply(legc,1,median)
  
  lfc_meso = rowMeans(lfc[,mc_incl])
  
  genes_f = names(lfc_meso)[lfc_meso > 1.7]
  genes_f = setdiff(genes_f,c("Dnmt3b","Pim2","Fgf5"))
  return(genes_f)
}


plot_epiblast_gene_expression_vs_time = function(genes_f = c("Dppa4","Gdf3"),plot_pdf = F,fig_dir = "figs/paper_figs/fig3/epiblast_expression_vs_time") {
  
  reg = 1e-5
  n_min_per_embryo = 10
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mat_wt = scdb_mat("sing_emb_wt10")
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  tko_chim = scdb_mat("tko_chim_wt10")
  tko_tetra = scdb_mat("tko_tetra_wt10")
  control_tetra_all = scdb_mat("control_tetra_all_wt10")
  
  load("data/tko_chim_wt10/tko_chim_md.Rda")
  load("data/tko_tetra_wt10/tko_tetra_md.Rda")
  load("data/control_tetra_all_wt10/control_tetra_all_md.Rda")
  time_chim = read.table("data/tko_chim_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  rownames(time_chim) = time_chim$embryo
  time_tetra = read.table("data/tko_tetra_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  rownames(time_tetra) = time_tetra$embryo
  time_control = read.table("data/control_tetra_all_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  rownames(time_control) = time_control$embryo
  
  
  age_group_time = read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt",h = T,sep = "\t",stringsAsFactors = F)
  age_group_time = age_group_time$developmental_time
  
  f_chim = (tko_chim_md$cell_type == "KO") & (tko_chim_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(tko_chim_md$age_group))
  n_cls_tko_chim_per_embryo =  table(tko_chim_md$embryo[f_chim])
  embryos_f = names(n_cls_tko_chim_per_embryo)[n_cls_tko_chim_per_embryo >= n_min_per_embryo]
  f_chim = f_chim & (tko_chim_md$embryo %in% embryos_f)
  e_t_chim = t(tgs_matrix_tapply(tko_chim@mat[,tko_chim_md$cell[f_chim]],tko_chim_md$embryo[f_chim],sum))
  e_t_chim = t(t(e_t_chim)/colSums(e_t_chim))
  
  f_chim_host_control = (tko_chim_md$cell_type %in% c("control","host")) & (tko_chim_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(tko_chim_md$age_group))
  n_cls_control_chim_per_embryo =  table(tko_chim_md$embryo[f_chim_host_control])
  embryos_f = names(n_cls_control_chim_per_embryo)[n_cls_control_chim_per_embryo >= n_min_per_embryo]
  f_chim_host_control = f_chim_host_control & (tko_chim_md$embryo %in% embryos_f)
  e_t_chim_control = t(tgs_matrix_tapply(tko_chim@mat[,tko_chim_md$cell[f_chim_host_control]],tko_chim_md$embryo[f_chim_host_control],sum))
  e_t_chim_control = t(t(e_t_chim_control)/colSums(e_t_chim_control))
  
  f_tetra = (tko_tetra_md$cell_type == "KO") & (tko_tetra_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(tko_tetra_md$age_group))
  n_cls_tko_tetra_per_embryo =  table(tko_tetra_md$embryo[f_tetra])
  embryos_f = names(n_cls_tko_tetra_per_embryo)[n_cls_tko_tetra_per_embryo >= n_min_per_embryo]
  f_tetra = f_tetra & (tko_tetra_md$embryo %in% embryos_f)
  e_t_tetra = t(tgs_matrix_tapply(tko_tetra@mat[,tko_tetra_md$cell[f_tetra]],tko_tetra_md$embryo[f_tetra],sum))
  e_t_tetra = t(t(e_t_tetra)/colSums(e_t_tetra))
  
  f_control = (control_tetra_all_md$cell_type == "control") & (control_tetra_all_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(control_tetra_all_md$age_group))
  n_cls_control_tetra_per_embryo =  table(control_tetra_all_md$embryo[f_control])
  embryos_f = names(n_cls_control_tetra_per_embryo)[n_cls_control_tetra_per_embryo >= n_min_per_embryo]
  f_control = f_control & (control_tetra_all_md$embryo %in% embryos_f)
  e_t_control = t(tgs_matrix_tapply(control_tetra_all@mat[,control_tetra_all_md$cell[f_control]],control_tetra_all_md$embryo[f_control],sum))
  e_t_control = t(t(e_t_control)/colSums(e_t_control))
   
    
  cls_wt = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] == mc_wt@color_key$color[1]]
  n_cls_wt_per_embryo = table(mat_wt@cell_metadata[cls_wt,"transcriptional_rank"])
  ranks_included = as.numeric(names(n_cls_wt_per_embryo)[n_cls_wt_per_embryo > n_min_per_embryo])
  cls_wt = cls_wt[mat_wt@cell_metadata[cls_wt,"transcriptional_rank"] %in% ranks_included]
  e_t_wt = t(tgs_matrix_tapply(mat_wt@mat[,cls_wt],mat_wt@cell_metadata[cls_wt,"transcriptional_rank"],sum))
  e_t_wt = t(t(e_t_wt)/colSums(e_t_wt))
  
  le_t_chim = log2(e_t_chim + reg)
  le_t_chim_control = log2(e_t_chim_control + reg)
  le_t_tetra = log2(e_t_tetra + reg)
  le_t_control = log2(e_t_control + reg)
  le_t_wt = log2(e_t_wt+ reg)
  
  roll_width = 6
  
  for (gene in genes_f) {
    
    freq_wt = le_t_wt[gene,]
    
    
    
    #freq_median = rollmedian(x = freq_wt,k = 2*roll_width+1,fill = c(mean(freq_wt[1:roll_width]),0,mean(freq_wt[(length(freq_wt)-roll_width + 1):length(freq_wt)])))
    freq_median = rollmean(x = freq_wt,k = 2*roll_width+1)
    freq_sd = rollapply(freq_wt,width = 2*roll_width + 1,function(v) {
      v_sd = sd(v)
      return(v_sd)
    })
    
    lower_freq = freq_median - 2*freq_sd
    upper_freq = freq_median + 2*freq_sd
    
    f_chim = time_chim[colnames(le_t_chim),"best_rank_host_control"] %in% as.numeric(names(freq_wt))
    f_tetra = time_tetra[colnames(le_t_tetra),"best_query"] %in% as.numeric(names(freq_wt))
    
    chim_wilcox_test = wilcox.test(x = le_t_chim[gene,f_chim],y = freq_wt[time_chim[colnames(le_t_chim)[f_chim],"best_rank_host_control"]])
    tetra_wilcox_test = wilcox.test(x = le_t_tetra[gene,f_tetra],y = freq_wt[time_tetra[colnames(le_t_tetra)[f_tetra],"best_query"]])
    
   
    if(plot_pdf) {
      pdf(sprintf("%s/epiblast_%s_chimera.pdf",fig_dir,gene),useDingbats = F)
    } else {
      png(sprintf("%s/epiblast_%s_chimera.png",fig_dir,gene))
    }
    #plot(x = age_group_time[as.numeric(colnames(le_t_wt))],y = le_t_wt[gene,],
    #     ylim = c(min(c(le_t_chim[gene,],le_t_tetra[gene,],le_t_wt[gene,])),max(c(le_t_chim[gene,],le_t_tetra[gene,],le_t_wt[gene,]))),
    #     xlim = c(6.4,8.2),
    #     pch = 19,ylab = "Log2 Expression",xlab = "Developmental Time",
    #     main = gene,col = "gray")
    ylim_plot = c(log2(reg),max(upper_freq,le_t_chim[gene,],le_t_tetra[gene,]))
    plot(x = age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],y = freq_median,type = "l",lwd = 2,ylim = ylim_plot,
         ylab = "Log2 Expression",xlab = "Developmental Time",
         main = gene,cex.main = 3)
    polygon(x = c(age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],rev(age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])])),
            y = c(upper_freq,rev(lower_freq)),col = "gray80",border = NA)
    lines(x = age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],
          y = freq_median,type = "l",lwd = 2)
    points(x = age_group_time[as.numeric(colnames(le_t_wt))],y = le_t_wt[gene,],pch = 19,cex = 0.3)
    
    points(x = age_group_time[time_chim[colnames(le_t_chim)[f_chim],"best_rank_host_control"]],y = le_t_chim[gene,f_chim],pch = 19,col = "indianred3",cex = 2)  
    text(x = 6.6,y = -16,labels = sprintf("p = %.1e",chim_wilcox_test$p.value))
    dev.off()
    
    
    
    
    if(plot_pdf) {
      pdf(sprintf("%s/epiblast_%s_tetraploid.pdf",fig_dir,gene),useDingbats = F)
    } else {
      png(sprintf("%s/epiblast_%s_tetraploid.png",fig_dir,gene))
    }
    #plot(x = age_group_time[as.numeric(colnames(le_t_wt))],y = le_t_wt[gene,],
    #     ylim = c(min(c(le_t_chim[gene,],le_t_tetra[gene,],le_t_wt[gene,])),max(c(le_t_chim[gene,],le_t_tetra[gene,],le_t_wt[gene,]))),
    #     xlim = c(6.4,8.2),
    #     pch = 19,ylab = "Log2 Expression",xlab = "Developmental Time",
    #     main = gene,col = "gray")
    ylim_plot = c(log2(reg),max(upper_freq,le_t_chim[gene,],le_t_tetra[gene,]))
    plot(x = age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],y = freq_median,type = "l",lwd = 2,ylim = ylim_plot,
         ylab = "Log2 Expression",xlab = "Developmental Time",
         main = gene,cex.main = 3)
    polygon(x = c(age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],rev(age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])])),
            y = c(upper_freq,rev(lower_freq)),col = "gray80",border = NA)
    lines(x = age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],
          y = freq_median,type = "l",lwd = 2)
    points(x = age_group_time[as.numeric(colnames(le_t_wt))],y = le_t_wt[gene,],pch = 19,cex = 0.3)

    points(x = age_group_time[time_tetra[colnames(le_t_tetra)[f_tetra],"best_query"]],y = le_t_tetra[gene,f_tetra],pch = 17,col = "indianred3",cex = 2)
    abline(v = 7.3,lty = "dashed")
    text(x = 6.6,y = -16,labels = sprintf("p = %.1e",tetra_wilcox_test$p.value))
    dev.off()
  }

}


plot_nascent_mesoderm_gene_expression_vs_time = function(genes_f = c("Dppa4","Gdf3"),plot_pdf = F,fig_dir = "figs/paper_figs/fig3/nascent_mesoderm_expression_vs_time") {
  
  reg = 1e-5
  n_min_per_embryo = 10
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mat_wt = scdb_mat("sing_emb_wt10")
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  tko_chim = scdb_mat("tko_chim_wt10")
  tko_tetra = scdb_mat("tko_tetra_wt10")
  control_tetra_all = scdb_mat("control_tetra_all_wt10")
  
  load("data/tko_chim_wt10/tko_chim_md.Rda")
  load("data/tko_tetra_wt10/tko_tetra_md.Rda")
  load("data/control_tetra_all_wt10/control_tetra_all_md.Rda")
  time_chim = read.table("data/tko_chim_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  rownames(time_chim) = time_chim$embryo
  time_tetra = read.table("data/tko_tetra_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  rownames(time_tetra) = time_tetra$embryo
  time_control = read.table("data/control_tetra_all_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  rownames(time_control) = time_control$embryo
  
  
  age_group_time = read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt",h = T,sep = "\t",stringsAsFactors = F)
  age_group_time = age_group_time$developmental_time
  
  f_chim = (tko_chim_md$cell_type == "KO") & (tko_chim_md$ct_color == mc_wt@color_key$color[10]) & (!is.na(tko_chim_md$age_group))
  n_cls_tko_chim_per_embryo =  table(tko_chim_md$embryo[f_chim])
  embryos_f = names(n_cls_tko_chim_per_embryo)[n_cls_tko_chim_per_embryo >= n_min_per_embryo]
  f_chim = f_chim & (tko_chim_md$embryo %in% embryos_f)
  e_t_chim = t(tgs_matrix_tapply(tko_chim@mat[,tko_chim_md$cell[f_chim]],tko_chim_md$embryo[f_chim],sum))
  e_t_chim = t(t(e_t_chim)/colSums(e_t_chim))
  
  f_tetra = (tko_tetra_md$cell_type == "KO") & (tko_tetra_md$ct_color == mc_wt@color_key$color[10]) & (!is.na(tko_tetra_md$age_group))
  n_cls_tko_tetra_per_embryo =  table(tko_tetra_md$embryo[f_tetra])
  embryos_f = names(n_cls_tko_tetra_per_embryo)[n_cls_tko_tetra_per_embryo >= n_min_per_embryo]
  f_tetra = f_tetra & (tko_tetra_md$embryo %in% embryos_f)
  e_t_tetra = t(tgs_matrix_tapply(tko_tetra@mat[,tko_tetra_md$cell[f_tetra]],tko_tetra_md$embryo[f_tetra],sum))
  e_t_tetra = t(t(e_t_tetra)/colSums(e_t_tetra))
  
  f_control = (control_tetra_all_md$cell_type == "control") & (control_tetra_all_md$ct_color == mc_wt@color_key$color[10]) & (!is.na(control_tetra_all_md$age_group))
  n_cls_control_tetra_per_embryo =  table(control_tetra_all_md$embryo[f_control])
  embryos_f = names(n_cls_control_tetra_per_embryo)[n_cls_control_tetra_per_embryo >= n_min_per_embryo]
  f_control = f_control & (control_tetra_all_md$embryo %in% embryos_f)
  e_t_control = t(tgs_matrix_tapply(control_tetra_all@mat[,control_tetra_all_md$cell[f_control]],control_tetra_all_md$embryo[f_control],sum))
  e_t_control = t(t(e_t_control)/colSums(e_t_control))
  
  
  cls_wt = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] == mc_wt@color_key$color[10]]
  n_cls_wt_per_embryo = table(mat_wt@cell_metadata[cls_wt,"transcriptional_rank"])
  ranks_included = as.numeric(names(n_cls_wt_per_embryo)[n_cls_wt_per_embryo > n_min_per_embryo])
  cls_wt = cls_wt[mat_wt@cell_metadata[cls_wt,"transcriptional_rank"] %in% ranks_included]
  e_t_wt = t(tgs_matrix_tapply(mat_wt@mat[,cls_wt],mat_wt@cell_metadata[cls_wt,"transcriptional_rank"],sum))
  e_t_wt = t(t(e_t_wt)/colSums(e_t_wt))
  
  le_t_chim = log2(e_t_chim + reg)
  le_t_tetra = log2(e_t_tetra + reg)
  le_t_control = log2(e_t_control + reg)
  le_t_wt = log2(e_t_wt+ reg)
  
  roll_width = 6
  
  for (gene in genes_f) {
    
    freq_wt = le_t_wt[gene,]
    
    #freq_median = rollmedian(x = freq_wt,k = 2*roll_width+1,fill = c(mean(freq_wt[1:roll_width]),0,mean(freq_wt[(length(freq_wt)-roll_width + 1):length(freq_wt)])))
    freq_median = rollmean(x = freq_wt,k = 2*roll_width+1)
    freq_sd = rollapply(freq_wt,width = 2*roll_width + 1,function(v) {
      v_sd = sd(v)
      return(v_sd)
    })
    
    lower_freq = freq_median - 2*freq_sd
    upper_freq = freq_median + 2*freq_sd
    
    f_chim = as.character(time_chim[colnames(le_t_chim),"best_rank_host_control"]) %in% names(freq_wt)
    f_tetra = as.character(time_tetra[colnames(le_t_tetra),"best_query"]) %in% names(freq_wt)
    
    chim_wilcox_test = wilcox.test(x = le_t_chim[gene,f_chim],y = freq_wt[as.character(time_chim[colnames(le_t_chim),"best_rank_host_control"])[f_chim]])
    tetra_wilcox_test = wilcox.test(x = le_t_tetra[gene,f_tetra],y = freq_wt[as.character(time_tetra[colnames(le_t_tetra),"best_query"])[f_tetra]])
    
    
    if(plot_pdf) {
      pdf(sprintf("%s/nascent_mesoderm_%s_chimera.pdf",fig_dir,gene),useDingbats = F)
    } else {
      png(sprintf("%s/nascent_mesoderm_%s_chimera.png",fig_dir,gene))
    }
    #plot(x = age_group_time[as.numeric(colnames(le_t_wt))],y = le_t_wt[gene,],
    #     ylim = c(min(c(le_t_chim[gene,],le_t_tetra[gene,],le_t_wt[gene,])),max(c(le_t_chim[gene,],le_t_tetra[gene,],le_t_wt[gene,]))),
    #     xlim = c(6.4,8.2),
    #     pch = 19,ylab = "Log2 Expression",xlab = "Developmental Time",
    #     main = gene,col = "gray")
    ylim_plot = c(log2(reg),max(upper_freq,le_t_chim[gene,],le_t_tetra[gene,]))
    xlim_plot = c(6.9,7.5)
    plot(x = age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],y = freq_median,type = "l",lwd = 2,ylim = ylim_plot,
         ylab = "Log2 Expression",xlab = "Developmental Time",
         main = gene,cex.main = 3,xlim = xlim_plot)
    polygon(x = c(age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],rev(age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])])),
            y = c(upper_freq,rev(lower_freq)),col = "gray80",border = NA)
    lines(x = age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],
          y = freq_median,type = "l",lwd = 2)
    points(x = age_group_time[as.numeric(colnames(le_t_wt))],y = le_t_wt[gene,],pch = 19,cex = 0.3)
    
    points(x = age_group_time[time_chim[colnames(le_t_chim),"best_rank_host_control"]],y = le_t_chim[gene,],pch = 19,col = "indianred3",cex = 2)
    text(x = 7.0,y = -16,labels = sprintf("p = %.1e",chim_wilcox_test$p.value))
    dev.off()
    
    
    if(plot_pdf) {
      pdf(sprintf("%s/nascent_mesoderm_%s_tetraploid.pdf",fig_dir,gene),useDingbats = F)
    } else {
      png(sprintf("%s/nascent_mesoderm_%s_tetraploid.png",fig_dir,gene))
    }
    #plot(x = age_group_time[as.numeric(colnames(le_t_wt))],y = le_t_wt[gene,],
    #     ylim = c(min(c(le_t_chim[gene,],le_t_tetra[gene,],le_t_wt[gene,])),max(c(le_t_chim[gene,],le_t_tetra[gene,],le_t_wt[gene,]))),
    #     xlim = c(6.4,8.2),
    #     pch = 19,ylab = "Log2 Expression",xlab = "Developmental Time",
    #     main = gene,col = "gray")
    ylim_plot = c(log2(reg),max(upper_freq,le_t_chim[gene,],le_t_tetra[gene,]))
    xlim_plot = c(6.9,7.5)
    plot(x = age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],y = freq_median,type = "l",lwd = 2,ylim = ylim_plot,
         ylab = "Log2 Expression",xlab = "Developmental Time",
         main = gene,cex.main = 3,xlim = xlim_plot)
    polygon(x = c(age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],rev(age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])])),
            y = c(upper_freq,rev(lower_freq)),col = "gray80",border = NA)
    lines(x = age_group_time[as.numeric(colnames(le_t_wt)[(roll_width + 1):(ncol(le_t_wt)-roll_width)])],
          y = freq_median,type = "l",lwd = 2)
    points(x = age_group_time[as.numeric(colnames(le_t_wt))],y = le_t_wt[gene,],pch = 19,cex = 0.3)
    points(x = age_group_time[time_tetra[colnames(le_t_tetra),"best_query"]],y = le_t_tetra[gene,],pch = 17,col = "indianred3",cex = 2)
    abline(v = 7.3,lty = "dashed")
    text(x = 7.0,y = -16,labels = sprintf("p = %.1e",tetra_wilcox_test$p.value))
    dev.off()
    
  }

}


cmp_single_cell_score = function(feat_genes,ct_col) {
  
  included_clones = c("TKO23","TKO26","TKO29","Ctrl1","Ctrl2","host")
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  wt_cls = names(mc@mc)[mc@colors[mc@mc] == ct_col]
  wt_cls_score = log2(colSums(mat@mat[feat_genes,wt_cls])/colSums(mat@mat[,wt_cls]) + 1e-3)
  
  tko_chim = scdb_mat("tko_chim_wt10")
  tko_tetra = scdb_mat("tko_tetra_wt10")
  ctrl_tetra = scdb_mat("control_tetra_all_wt10")
  
  load("data/tko_chim_wt10/color_annotation/cmp_annot.Rda")
  tko_chim_annot = cmp_annot
  load("data/tko_tetra_wt10/color_annotation/cmp_annot.Rda")
  tko_tetra_annot = cmp_annot
  load("data/control_tetra_all_wt10/color_annotation/cmp_annot.Rda")
  ctrl_tetra_annot = cmp_annot
  
  f = !is.na(tko_chim@cell_metadata[colnames(tko_chim@mat),"clone_type"])
  query_chim_cls = colnames(tko_chim@mat)[f][tko_chim@cell_metadata[colnames(tko_chim@mat)[f],"clone_type"] %in% included_clones]
  query_chim_cls = intersect(query_chim_cls,names(tko_chim_annot$query_cls_col))
  query_chim_cls = query_chim_cls[tko_chim_annot$query_cls_col[query_chim_cls] == ct_col]
  query_chim_score = log2(colSums(tko_chim@mat[feat_genes,query_chim_cls])/colSums(tko_chim@mat[,query_chim_cls]) + 1e-3)
  
  f = !is.na(ctrl_tetra@cell_metadata[colnames(ctrl_tetra@mat),"clone_type"])
  ctrl_tetra_cls = colnames(ctrl_tetra@mat)[f][ctrl_tetra@cell_metadata[colnames(ctrl_tetra@mat)[f],"clone_type"] %in% c("Ctrl1","Ctrl2")]
  ctrl_tetra_cls = intersect(ctrl_tetra_cls,names(ctrl_tetra_annot$query_cls_col))
  ctrl_tetra_cls = ctrl_tetra_cls[ctrl_tetra_annot$query_cls_col[ctrl_tetra_cls] == ct_col]
  ctrl_tetra_score = log2(colSums(ctrl_tetra@mat[feat_genes,ctrl_tetra_cls])/colSums(ctrl_tetra@mat[,ctrl_tetra_cls]) + 1e-3)
  
  f = !is.na(tko_tetra@cell_metadata[colnames(tko_tetra@mat),"clone_type"])
  tko_tetra_cls = colnames(tko_tetra@mat)[f][tko_tetra@cell_metadata[colnames(tko_tetra@mat)[f],"clone_type"] %in% c("TKO26","TKO29")]
  tko_tetra_cls = intersect(tko_tetra_cls,names(tko_tetra_annot$query_cls_col))
  tko_tetra_cls = tko_tetra_cls[tko_tetra_annot$query_cls_col[tko_tetra_cls] == ct_col]
  tko_tetra_score = log2(colSums(tko_tetra@mat[feat_genes,tko_tetra_cls])/colSums(tko_tetra@mat[,tko_tetra_cls]) + 1e-3)
  
  genotype_nm = c("KO" = "TKO","control" = "Ctrl","host" = "Host")
  
  df_score = data.frame(clone = c(tko_chim@cell_metadata[query_chim_cls,"clone_type"],tko_tetra@cell_metadata[tko_tetra_cls,"clone_type"],ctrl_tetra@cell_metadata[ctrl_tetra_cls,"clone_type"]),
                        experimental_assay = c(rep("2N",length(query_chim_cls)),rep("4N",length(c(tko_tetra_cls,ctrl_tetra_cls)))),
                        sc_names = c(query_chim_cls,tko_tetra_cls,ctrl_tetra_cls),
                        Genotype = genotype_nm[c(tko_chim@cell_metadata[query_chim_cls,"cell_type"],tko_tetra@cell_metadata[tko_tetra_cls,"cell_type"],ctrl_tetra@cell_metadata[ctrl_tetra_cls,"cell_type"])],
                        sc_score = c(query_chim_score,tko_tetra_score,ctrl_tetra_score),stringsAsFactors = F)
  
  df_score$multi_panel = paste(df_score$experimental_assay,df_score$Genotype,sep = " ")
  
  multi_panel_wt = rep(c("2N Ctrl","4N Ctrl","2N Host","2N TKO","4N TKO"),each = length(wt_cls_score))
  
  df_score_wt = data.frame(clone = rep("WT",length(multi_panel_wt)),
                           experimental_assay = rep("WT",length(multi_panel_wt)),
                           sc_names = rep(wt_cls,5),
                           Genotype = rep("WT",length(multi_panel_wt)),
                           sc_score = rep(wt_cls_score,5),
                           multi_panel = multi_panel_wt,stringsAsFactors = F)
  
  
  df_score_all = bind_rows(df_score,df_score_wt)
  df_score_all$multi_panel = factor(x = df_score_all$multi_panel,levels = c("2N Ctrl","4N Ctrl","2N TKO","4N TKO","2N Host"))
  
  return(df_score_all)
}

plot_epiblast_score_distribution = function(plot_pdf = F) {
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  ct_col = mc@color_key$color[1]
  feat_genes = epiblast_genes()
  
  df_score_all = cmp_single_cell_score(feat_genes,ct_col)
  
  cell_type_cols = c("WT" = "gray60","Ctrl" = "blue","Host" = "black","TKO" = "indianred3")
  
  f = df_score_all$Genotype != "WT"
  n_cls = table(df_score_all$multi_panel[f])
  
  df_mean_and_sd = compute_mean_and_sd_of_expression(df_score_all = df_score_all)

  
  write.csv(df_mean_and_sd,file = "figs/paper_figs/fig3/fig3a_mean_and_sd.csv")
  
  df_text = data.frame(multi_panel = names(n_cls), n_cls = paste0("N = ",as.vector(n_cls)), x = rep(-9,length(n_cls)),y = rep(1,length(n_cls)))
  

  p <- ggplot(df_score_all, aes(x = sc_score,color = Genotype)) +
    geom_density(size = 1) +
    facet_wrap(multi_panel ~ .,nrow = 5) +
    scale_color_manual(values = cell_type_cols) +
    geom_text(data = df_text,mapping = aes(x = x,y =y,label = n_cls),col = "black",hjust = 0) + 
    xlab("Epiblast score")
  
  #print(p)
  if (plot_pdf) {
    fn = "figs/paper_figs/fig3/epiblast_score_per_experiment_genotype.pdf"
  } else {
    fn = "figs/paper_figs/fig3/epiblast_score_per_experiment_genotype.png"
  }
  #ggsave(filename = fn,plot = p,width = 4,h = 4)
  ggsave(filename = fn,plot = p,width = 3,h = 10)

}

heatmap_epiblast_differential_expression = function() {
  a = 1
}

heatmap_early_nascent_mesoderm_differential_expression = function() {
  a = 2
}

plot_nascent_mesoderm_score_distribution = function(plot_pdf = F) {
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  ct_col = mc@color_key$color[10]
  feat_genes = early_nascent_mesoderm_genes()
  
  df_score_all = cmp_single_cell_score(feat_genes,ct_col)
  
  cell_type_cols = c("WT" = "gray60","Ctrl" = "blue","Host" = "black","TKO" = "indianred3")
  
  f = df_score_all$Genotype != "WT"
  n_cls = table(df_score_all$multi_panel[f])
  
  df_mean_and_sd = compute_mean_and_sd_of_expression(df_score_all = df_score_all)
  
  
  write.csv(df_mean_and_sd,file = "figs/paper_figs/fig3/fig3d_mean_and_sd.csv")
  
  df_text = data.frame(multi_panel = names(n_cls), n_cls = paste0("N = ",as.vector(n_cls)), x = rep(-9,length(n_cls)),y = rep(1,length(n_cls)))
  
  
  p <- ggplot(df_score_all, aes(x = sc_score,color = Genotype)) +
    geom_density(size = 1) +
    facet_wrap(multi_panel ~ .,nrow = 5) +
    scale_color_manual(values = cell_type_cols) +
    geom_text(data = df_text,mapping = aes(x = x,y =y,label = n_cls),col = "black",hjust = 0) + 
    xlab("Epiblast score")
  
  #print(p)
  if (plot_pdf) {
    fn = "figs/paper_figs/fig3/nascent_mesoderm_score_per_experiment_genotype.pdf"
  } else {
    fn = "figs/paper_figs/fig3/nascent_mesoderm_score_per_experiment_genotype.png"
  }
  #ggsave(filename = fn,plot = p,width = 4,h = 4)
  ggsave(filename = fn,plot = p,width = 3,h = 10)


}

compute_mean_and_sd_of_expression = function(df_score_all) {
  
  
  # Add KS statistics to the plot
  df_ks = df_score_all[df_score_all$Genotype != "WT",]
  
  mean_expression = tapply(df_ks$sc_score,df_ks$multi_panel,mean)
  sd_expression = tapply(df_ks$sc_score,df_ks$multi_panel,sd)
  
  wt_score= df_score_all$sc_score[df_score_all$multi_panel == "2N Host" & df_score_all$Genotype == "WT"]
  
  wt_mean = mean(wt_score)
  wt_sd = sd(wt_score)
  
  sd_expression = c(sd_expression,c("WT" = wt_sd))
  mean_expression = c(mean_expression,c("WT" = wt_mean))
  
  df_out = data.frame(multi_panel = names(mean_expression),
                      mean = mean_expression,
                      sd = sd_expression[names(mean_expression)],stringsAsFactors = F)
  
  return(df_out)
}


if(0) {
  
  plot_epiblast_score_distribution_per_genotype_and_experiment = function(plot_pdf = F) {
    
    
    mat = scdb_mat("sing_emb_wt10")
    mc = scdb_mc("sing_emb_wt10_recolored")
    
    wt_cls = names(mc@mc)[mc@colors[mc@mc] == mc@color_key$color[1]]
    wt_cls_score = log2(colSums(mat@mat[epiblast_genes(),wt_cls])/colSums(mat@mat[,wt_cls]) + 1e-3)
    
    n_top = 500
    top_value = sort(wt_cls_score,decreasing = T)[n_top]
    
    cmp = cmp_epiblast_score()
    df_score = cmp$df_score
    df_clone = cmp$df_clone
    
    tmp = c("control Chimera" = "2N Ctrl","control Tetraploid" = "4N Ctrl","host Chimera" = "2N Host","KO Chimera" = "2N TKO","KO Tetraploid" = "4N TKO")
    experiment_genotype = tmp[paste(df_score$cell_type,df_score$experimental_assay,sep = " ")]
    df_score$multi_panel = factor(x = experiment_genotype,levels = c("2N Ctrl","2N TKO","2N Host","4N Ctrl","4N TKO"))
    
    # add wildtype cells to the score
    experiment_genotype_wt = rep(tmp,each = length(wt_cls_score))
    wt_score = rep(wt_cls_score,length(df_clone$clone_assay))
    
    df_wt_score = data.frame(clone_assay = "WT",experimental_assay = "WT",cell_type = "WT",clone_type = "WT",epiblast_score = wt_score,multi_panel = experiment_genotype_wt,
                             sc_name = names(wt_score),stringsAsFactors = F)
    
    df_score_all = bind_rows(df_score,df_wt_score)
    
    cell_type_cols = c("WT" = "gray60","control" = "blue","host" = "black","KO" = "indianred3")
    
    n_cls = table(df_score$multi_panel)
    n_cls = n_cls[names(n_cls) != "WT"]
    
    df_text = data.frame(multi_panel = names(n_cls), n_cls = paste0("N = ",as.vector(n_cls)), x = rep(-9,length(n_cls)),y = rep(1,length(n_cls)))
    
    df_score_all$multi_panel = factor(x = df_score_all$multi_panel,levels = c("2N Ctrl","4N Ctrl","2N TKO","4N TKO","2N Host"))
    
    df_score_all = rename(df_score_all,Genotype = cell_type)
    
    p <- ggplot(df_score_all, aes(x = epiblast_score,color = Genotype)) +
      geom_density(size = 1) +
      facet_wrap(multi_panel ~ .,nrow = 5) +
      geom_vline(xintercept = top_value,linetype = "dashed",col = "gray60") +
      scale_color_manual(values = cell_type_cols) +
      geom_text(data = df_text,mapping = aes(x = x,y =y,label = n_cls),col = "black",hjust = 0) + 
      xlab("Epiblast score")
    
    #print(p)
    if (plot_pdf) {
      fn = "figs/paper_figs/fig3/epiblast_score_per_experiment_genotype.pdf"
    } else {
      fn = "figs/paper_figs/fig3/epiblast_score_per_experiment_genotype.png"
    }
    #ggsave(filename = fn,plot = p,width = 4,h = 4)
    ggsave(filename = fn,plot = p,width = 3,h = 10)
    
    
  }
  
  
  
  plot_nascent_mesoderm_score_distribution_per_genotype_and_experiment = function(plot_pdf = F) {
    
    mat = scdb_mat("sing_emb_wt10")
    mc = scdb_mc("sing_emb_wt10_recolored")
    
    wt_cls = names(mc@mc)[mc@colors[mc@mc] == mc@color_key$color[10]]
    wt_cls_score = log2(colSums(mat@mat[early_nascent_mesoderm_genes(),wt_cls])/colSums(mat@mat[,wt_cls]) + 1e-3)
    
    n_top = 500
    top_value = sort(wt_cls_score,decreasing = T)[n_top]
    
    cmp = cmp_nascent_mesoderm_score()
    df_score = cmp$df_score
    df_clone = cmp$df_clone
    
    tmp = c("control Chimera" = "2N Ctrl","control Tetraploid" = "4N Ctrl","host Chimera" = "2N Host","KO Chimera" = "2N TKO","KO Tetraploid" = "4N TKO")
    experiment_genotype = tmp[paste(df_score$cell_type,df_score$experimental_assay,sep = " ")]
    df_score$multi_panel = factor(x = experiment_genotype,levels = c("2N Ctrl","2N TKO","2N Host","4N Ctrl","4N TKO"))
    
    # add wildtype cells to the score
    experiment_genotype_wt = rep(tmp,each = length(wt_cls_score))
    wt_score = rep(wt_cls_score,length(df_clone$clone_assay))
    
    df_wt_score = data.frame(clone_assay = "WT",experimental_assay = "WT",cell_type = "WT",clone_type = "WT",nm_score = wt_score,multi_panel = experiment_genotype_wt,
                             sc_name = names(wt_score),stringsAsFactors = F)
    
    df_score_all = bind_rows(df_score,df_wt_score)
    
    cell_type_cols = c("WT" = "gray60","control" = "blue","host" = "black","KO" = "indianred3")
    
    n_cls = table(df_score$multi_panel)
    n_cls = n_cls[names(n_cls) != "WT"]
    
    df_text = data.frame(multi_panel = names(n_cls), n_cls = paste0("N = ",as.vector(n_cls)), x = rep(-9,length(n_cls)),y = rep(0.3,length(n_cls)))
    
    df_score_all$multi_panel = factor(x = df_score_all$multi_panel,levels = c("2N Ctrl","4N Ctrl","2N TKO","4N TKO","2N Host"))
    
    df_score_all = rename(df_score_all,Genotype = cell_type)
    
    p <- ggplot(df_score_all, aes(x = nm_score,color = Genotype)) +
      geom_density(size = 1) +
      facet_wrap(multi_panel ~ .,nrow = 2) +
      geom_vline(xintercept = top_value,linetype = "dashed",col = "gray60") +
      scale_color_manual(values = cell_type_cols) +
      geom_text(data = df_text,mapping = aes(x = x,y =y,label = n_cls),col = "black",hjust = 0) + 
      xlab("Nascent mesoderm score")
    
    #print(p)
    if (plot_pdf) {
      fn = "figs/paper_figs/fig3/nm_score_per_experiment_genotype.pdf"
    } else {
      fn = "figs/paper_figs/fig3/nm_score_per_experiment_genotype.png"
    }
    ggsave(filename = fn,plot = p,width = 7.5,h = 5)
    
    
    
  }
  
  
  
  
  
  epiblast_expression_vs_time_cmp_lfp = function(reg = 1e-5) {
    
    rank_max = 114
    
    n_min_per_embryo = 20
    
    mat_wt = scdb_mat("sing_emb_wt10")
    mc_wt = scdb_mc("sing_emb_wt10_recolored")
    tko_chim = scdb_mat("tko_chim_wt10")
    tko_tetra = scdb_mat("tko_tetra_wt10")
    control_tetra_all = scdb_mat("control_tetra_all_wt10")
    
    load("data/tko_chim_wt10/tko_chim_md.Rda")
    load("data/tko_tetra_wt10/tko_tetra_md.Rda")
    load("data/control_tetra_all_wt10/control_tetra_all_md.Rda")
    time_chim = read.table("data/tko_chim_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
    rownames(time_chim) = time_chim$embryo
    time_tetra = read.table("data/tko_tetra_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
    rownames(time_tetra) = time_tetra$embryo
    time_control = read.table("data/control_tetra_all_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
    rownames(time_control) = time_control$embryo
    
    
    age_group_time = read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt",h = T,sep = "\t",stringsAsFactors = F)
    age_group_time = age_group_time$developmental_time
    
    f_chim = (tko_chim_md$cell_type == "KO") & (tko_chim_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(tko_chim_md$age_group))
    n_cls_tko_chim_per_embryo =  table(tko_chim_md$embryo[f_chim])
    embryos_f = names(n_cls_tko_chim_per_embryo)[n_cls_tko_chim_per_embryo >= n_min_per_embryo]
    f_chim = f_chim & (tko_chim_md$embryo %in% embryos_f)
    e_t_chim = t(tgs_matrix_tapply(tko_chim@mat[,tko_chim_md$cell[f_chim]],tko_chim_md$embryo[f_chim],sum))
    e_t_chim = t(t(e_t_chim)/colSums(e_t_chim))
    
    f_tetra = (tko_tetra_md$cell_type == "KO") & (tko_tetra_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(tko_tetra_md$age_group))
    n_cls_tko_tetra_per_embryo =  table(tko_tetra_md$embryo[f_tetra])
    embryos_f = names(n_cls_tko_tetra_per_embryo)[n_cls_tko_tetra_per_embryo >= n_min_per_embryo]
    f_tetra = f_tetra & (tko_tetra_md$embryo %in% embryos_f)
    e_t_tetra = t(tgs_matrix_tapply(tko_tetra@mat[,tko_tetra_md$cell[f_tetra]],tko_tetra_md$embryo[f_tetra],sum))
    e_t_tetra = t(t(e_t_tetra)/colSums(e_t_tetra))
    
    f_control = (control_tetra_all_md$cell_type == "control") & (control_tetra_all_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(control_tetra_all_md$age_group))
    n_cls_control_tetra_per_embryo =  table(control_tetra_all_md$embryo[f_control])
    embryos_f = names(n_cls_control_tetra_per_embryo)[n_cls_control_tetra_per_embryo >= n_min_per_embryo]
    f_control = f_control & (control_tetra_all_md$embryo %in% embryos_f)
    e_t_control = t(tgs_matrix_tapply(control_tetra_all@mat[,control_tetra_all_md$cell[f_control]],control_tetra_all_md$embryo[f_control],sum))
    e_t_control = t(t(e_t_control)/colSums(e_t_control))
    
    
    cls_wt = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] == mc_wt@color_key$color[1]]
    n_cls_wt_per_embryo = table(mat_wt@cell_metadata[cls_wt,"transcriptional_rank"])
    ranks_included = as.numeric(names(n_cls_wt_per_embryo)[n_cls_wt_per_embryo > n_min_per_embryo])
    cls_wt = cls_wt[mat_wt@cell_metadata[cls_wt,"transcriptional_rank"] %in% ranks_included]
    e_t_wt = t(tgs_matrix_tapply(mat_wt@mat[,cls_wt],mat_wt@cell_metadata[cls_wt,"transcriptional_rank"],sum))
    e_t_wt = t(t(e_t_wt)/colSums(e_t_wt))
    
    le_t_chim = log2(e_t_chim + reg)
    le_t_tetra = log2(e_t_tetra + reg)
    le_t_control = log2(e_t_control + reg)
    
    roll_width = 6
    
    e_t_ref = sapply(X = c(1:rank_max),FUN = function(t) {
      
      f_ind = as.numeric(colnames(e_t_wt)) %in% (t-roll_width):(t+roll_width) 
      
      e_g = rowMeans(e_t_wt[,f_ind])
      return(e_g)
    })
    colnames(e_t_ref) = c(1:rank_max)
    
    le_t_ref = log2(e_t_ref + reg)
    
    f_tetra_t = as.character(time_tetra[colnames(le_t_tetra),"best_query"]) %in% colnames(le_t_ref)
    f_chim_t = as.character(time_chim[colnames(le_t_chim),"best_rank_host_control"]) %in% colnames(le_t_ref)
    f_control_t = as.character(time_control[colnames(le_t_control),"best_query"]) %in% colnames(le_t_ref)
    
    lfc_chim = le_t_chim[,f_chim_t] - le_t_ref[,as.character(time_chim[colnames(le_t_chim),"best_rank_host_control"])[f_chim_t]]
    lfc_tetra = le_t_tetra[,f_tetra_t] - le_t_ref[,as.character(time_tetra[colnames(le_t_tetra),"best_query"])[f_tetra_t]]
    lfc_control = le_t_control[,f_control_t] - le_t_ref[,as.character(time_control[colnames(le_t_control),"best_query"])[f_control_t]]
    
    return(list(lfc_chim = lfc_chim,lfc_tetra = lfc_tetra,lfc_control = lfc_control))
  }
  
  
  old_epiblast_gene_expression_vs_age_group = function(plot_pdf = F) {
    
    fig_dir = "figs/paper_figs/fig3/epiblast_expression_vs_age_group"
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    
    mat_wt = scdb_mat("sing_emb_wt10")
    mc_wt = scdb_mc("sing_emb_wt10_recolored")
    tko_chim = scdb_mat("tko_chim_wt10")
    tko_tetra = scdb_mat("tko_tetra_wt10")
    control_tetra_all = scdb_mat("control_tetra_all_wt10")
    
    load("data/tko_chim_wt10/tko_chim_md.Rda")
    load("data/tko_tetra_wt10/tko_tetra_md.Rda")
    load("data/control_tetra_all_wt10/control_tetra_all_md.Rda")
    
    age_group_time = read.table(file = "data/wt10_age_groups_developmental_time.txt",sep = "\t",stringsAsFactors = F)
    age_group_time = age_group_time$developmental_time
    
    f_chim = (tko_chim_md$cell_type == "KO") & (tko_chim_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(tko_chim_md$age_group))
    e_t_chim = t(tgs_matrix_tapply(tko_chim@mat[,tko_chim_md$cell[f_chim]],tko_chim_md$age_group[f_chim],sum))
    e_t_chim = t(t(e_t_chim)/colSums(e_t_chim))
    
    f_tetra = (tko_tetra_md$cell_type == "KO") & (tko_tetra_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(tko_tetra_md$age_group))
    e_t_tetra = t(tgs_matrix_tapply(tko_tetra@mat[,tko_tetra_md$cell[f_tetra]],tko_tetra_md$age_group[f_tetra],sum))
    e_t_tetra = t(t(e_t_tetra)/colSums(e_t_tetra))
    
    f_control = (control_tetra_all_md$cell_type == "control") & (control_tetra_all_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(control_tetra_all_md$age_group))
    e_t_control = t(tgs_matrix_tapply(control_tetra_all@mat[,control_tetra_all_md$cell[f_control]],control_tetra_all_md$age_group[f_control],sum))
    e_t_control = t(t(e_t_control)/colSums(e_t_control))
    
    cls_wt = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] == mc_wt@color_key$color[1]]
    e_t_wt = t(tgs_matrix_tapply(mat_wt@mat[,cls_wt],mat_wt@cell_metadata[cls_wt,"age_group"],sum))
    e_t_wt = t(t(e_t_wt)/colSums(e_t_wt))
    
    genes_table = read.csv("data/gene_list_fig3.csv")
    genes_f = union(c("Pou3f1","Sox11","Tcf3","Foxp1"),genes_table$Gene.name)
    
    le_t_chim = log2(e_t_chim[,as.character(c(4:5))] + 1e-5)
    le_t_tetra = log2(e_t_tetra[,as.character(c(2:6))] + 1e-5)
    le_t_control = log2(e_t_control[,as.character(c(4:6))] + 1e-5)
    le_t_wt = log2(e_t_wt[,as.character(c(1:6))] + 1e-5)
    
    for (gene in genes_f) {
      
      if(plot_pdf) {
        pdf(sprintf("%s/epiblast_%s.pdf",fig_dir,gene),useDingbats = F)
      } else {
        png(sprintf("%s/epiblast_%s.png",fig_dir,gene))
      }
      plot(x = age_group_time[as.numeric(colnames(le_t_wt))],y = le_t_wt[gene,],type = "l",
           ylim = c(min(c(le_t_chim[gene,],le_t_tetra[gene,],le_t_wt[gene,])),max(c(le_t_chim[gene,],le_t_tetra[gene,],le_t_wt[gene,]))),
           lwd = 3,ylab = "Log2 Expression",xlab = "Developmental Time",
           main = gene,col = "gray")
      points(x = age_group_time[as.numeric(colnames(le_t_chim))],y = le_t_chim[gene,],pch = 19,col = "darkgoldenrod2",cex = 4)  
      points(x = age_group_time[as.numeric(colnames(le_t_tetra))],y = le_t_tetra[gene,],pch = 19,col = "darkgreen",cex = 4)
      dev.off()
    }
    
  }
  

  cmp_epiblast_score = function() {
    
    feat_genes = epiblast_genes()
    
    mc = scdb_mc("sing_emb_wt10_recolored")
    tko_chim = scdb_mat("tko_chim_wt10")
    tko_tetra = scdb_mat("tko_tetra_wt10")
    ctrl_tetra = scdb_mat("control_tetra_all_wt10")
    
    df_clone = data.frame(clone_assay = c("TKO26 Chimera","TKO26 Tetraploid","TKO29 Chimera","TKO29 Tetraploid","TKO23 Chimera","Host Chimera","Ctrl1 Chimera","Ctrl1 Tetraploid","Ctrl2 Chimera","Ctrl2 Tetraploid"),
                          cell_type =  c("KO"   ,"KO"   ,"KO"   ,"KO"   ,"KO"   ,"host","control","control","control","control"),
                          clone_type = c("TKO26","TKO26","TKO29","TKO29","TKO23","host","Ctrl1"  ,"Ctrl1"  ,"Ctrl2"  ,"Ctrl2"),
                          experimental_assay = c("Chimera","Tetraploid","Chimera","Tetraploid","Chimera","Chimera","Chimera","Tetraploid","Chimera","Tetraploid"),
                          mat_nm = c("tko_chim_wt10","tko_tetra_wt10","tko_chim_wt10","tko_tetra_wt10","tko_chim_wt10","tko_chim_wt10","tko_chim_wt10","control_tetra_all_wt10","tko_chim_wt10","control_tetra_all_wt10"),stringsAsFactors = F)
    
    epiblast_score = list()
    
    df_score = data.frame(clone_assay = c(),experimental_assay = c(),cell_type = c(),clone_type = c(),epiblast_score = c(),stringsAsFactors = F)
    
    for (i in 1:nrow(df_clone)) {
      print(i)
      mat_nm = df_clone$mat_nm[i]
      scmat = scdb_mat(mat_nm)
      cell_type = df_clone$cell_type[i]
      clone_type = df_clone$clone_type[i]
      tag = df_clone$clone_assay[i]
      experimental_assay = df_clone$experimental_assay[i]
      
      query_cls = colnames(scmat@mat)[scmat@cell_metadata[colnames(scmat@mat),"cell_type"] == cell_type]
      query_cls = query_cls[scmat@cell_metadata[query_cls,"clone_type"] == clone_type]
      query_score = query_cmp_sc_score_per_clone(mat_nm = mat_nm,cls_f = query_cls,ct_col = mc@color_key$color[1],genes = feat_genes)
      
      epiblast_score[[tag]] = query_score
      
      df_tmp = data.frame(clone_assay = rep(tag,length(query_score)),
                          experimental_assay = rep(experimental_assay,length(query_score)),
                          cell_type = rep(cell_type,length(query_score)),
                          clone_type = rep(clone_type,length(query_score)),
                          epiblast_score = query_score,
                          sc_names = names(query_score),stringsAsFactors = F)
      
      df_score = rbind(df_score,df_tmp)
    }
    
    df_score$multi_panel = df_score$clone_assay
    
    
    return(list(df_score = df_score,df_clone = df_clone))
  }
  
  cmp_nascent_mesoderm_score = function() {
    
    feat_genes = early_nascent_mesoderm_genes()
    
    
    mc = scdb_mc("sing_emb_wt10_recolored")
    tko_chim = scdb_mat("tko_chim_wt10")
    tko_tetra = scdb_mat("tko_tetra_wt10")
    ctrl_tetra = scdb_mat("control_tetra_all_wt10")
    
    df_clone = data.frame(clone_assay = c("TKO26 Chimera","TKO26 Tetraploid","TKO29 Chimera","TKO29 Tetraploid","TKO23 Chimera","Host Chimera","Ctrl1 Chimera","Ctrl1 Tetraploid","Ctrl2 Chimera","Ctrl2 Tetraploid"),
                          cell_type =  c("KO"   ,"KO"   ,"KO"   ,"KO"   ,"KO"   ,"host","control","control","control","control"),
                          clone_type = c("TKO26","TKO26","TKO29","TKO29","TKO23","host","Ctrl1"  ,"Ctrl1"  ,"Ctrl2"  ,"Ctrl2"),
                          experimental_assay = c("Chimera","Tetraploid","Chimera","Tetraploid","Chimera","Chimera","Chimera","Tetraploid","Chimera","Tetraploid"),
                          mat_nm = c("tko_chim_wt10","tko_tetra_wt10","tko_chim_wt10","tko_tetra_wt10","tko_chim_wt10","tko_chim_wt10","tko_chim_wt10","control_tetra_all_wt10","tko_chim_wt10","control_tetra_all_wt10"),stringsAsFactors = F)
    
    nm_score = list()
    
    df_score = data.frame(clone_assay = c(),experimental_assay = c(),cell_type = c(),clone_type = c(),nm_score = c(),stringsAsFactors = F)
    
    for (i in 1:nrow(df_clone)) {
      print(i)
      mat_nm = df_clone$mat_nm[i]
      scmat = scdb_mat(mat_nm)
      cell_type = df_clone$cell_type[i]
      clone_type = df_clone$clone_type[i]
      tag = df_clone$clone_assay[i]
      experimental_assay = df_clone$experimental_assay[i]
      
      query_cls = colnames(scmat@mat)[scmat@cell_metadata[colnames(scmat@mat),"cell_type"] == cell_type]
      query_cls = query_cls[scmat@cell_metadata[query_cls,"clone_type"] == clone_type]
      query_score = query_cmp_sc_score_per_clone(mat_nm = mat_nm,cls_f = query_cls,ct_col = mc@color_key$color[10],genes = feat_genes)
      
      nm_score[[tag]] = query_score
      
      df_tmp = data.frame(clone_assay = rep(tag,length(query_score)),
                          experimental_assay = rep(experimental_assay,length(query_score)),
                          cell_type = rep(cell_type,length(query_score)),
                          clone_type = rep(clone_type,length(query_score)),
                          nm_score = query_score,
                          sc_names = names(query_score),stringsAsFactors = F)
      
      df_score = rbind(df_score,df_tmp)
    }
    
    df_score$multi_panel = df_score$clone_assay
    
    return(list(df_score = df_score,df_clone = df_clone))
  }
  
  
  plot_epiblast_score_distribution_per_clone_and_experiment = function(plot_pdf = F) {
    
    feat_genes = epiblast_genes()
    mat = scdb_mat("sing_emb_wt10")
    mc = scdb_mc("sing_emb_wt10_recolored")
    
    wt_cls = names(mc@mc)[mc@colors[mc@mc] == mc@color_key$color[1]]
    wt_cls_score = log2(colSums(mat@mat[feat_genes,wt_cls])/colSums(mat@mat[,wt_cls]) + 1e-3)
    
    n_top = 500
    top_value = sort(wt_cls_score,decreasing = T)[n_top]
    
    cmp = cmp_epiblast_score()
    df_score = cmp$df_score
    df_clone = cmp$df_clone
    
    # add wildtype cells to the score
    clone_assay_wt = rep(df_clone$clone_assay,each = length(wt_cls_score))
    wt_score = rep(wt_cls_score,length(df_clone$clone_assay))
    
    df_wt_score = data.frame(clone_assay = "WT",experimental_assay = "WT",cell_type = "WT",clone_type = "WT",epiblast_score = wt_score,multi_panel = clone_assay_wt,
                             sc_name = names(wt_score),stringsAsFactors = F)
    
    df_score_all = bind_rows(df_score,df_wt_score)
    
    cell_type_cols = c("WT" = "gray60","control" = "blue","host" = "black","KO" = "indianred3")
    
    n_cls = table(df_score_all$clone_assay)
    n_cls = n_cls[names(n_cls) != "WT"]
    
    df_text = data.frame(multi_panel = names(n_cls), n_cls = paste0("N = ",as.vector(n_cls)), x = rep(-9,length(n_cls)),y = rep(1,length(n_cls)))
    
    p <- ggplot(df_score_all, aes(x = epiblast_score,color = cell_type)) +
      geom_density(size = 1) +
      facet_wrap(multi_panel ~ .,nrow = 2) +
      geom_vline(xintercept = top_value,linetype = "dashed",col = "gray60") +
      ggtitle("Epiblast score distribution among epiblast single cells") + scale_color_manual() +
      scale_color_manual(values = cell_type_cols) +
      geom_text(data = df_text,mapping = aes(x = x,y =y,label = n_cls),col = "black",hjust = 0)
    
    #print(p)
    if (plot_pdf) {
      fn = "figs/paper_figs/fig3/epiblast_score_per_clone_assay.pdf"
    } else {
      fn = "figs/paper_figs/fig3/epiblast_score_per_clone_assay.png"
    }
    ggsave(filename = fn,plot = p,width = 9,h = 4)
    
    
    
  }
  
  
  
  
  plot_nascent_mesoderm_score_distribution_per_clone_and_experiment = function(plot_pdf = F) {
    feat_genes = early_nascent_mesoderm_genes()
    mat = scdb_mat("sing_emb_wt10")
    mc = scdb_mc("sing_emb_wt10_recolored")
    
    wt_cls = names(mc@mc)[mc@colors[mc@mc] == mc@color_key$color[10]]
    wt_cls_score = log2(colSums(mat@mat[feat_genes,wt_cls])/colSums(mat@mat[,wt_cls]) + 1e-3)
    
    n_top = 500
    top_value = sort(wt_cls_score,decreasing = T)[n_top]
    
    cmp = cmp_nascent_mesoderm_score()
    df_score = cmp$df_score
    df_clone = cmp$df_clone
    # add wildtype cells to the score
    clone_assay_wt = rep(df_clone$clone_assay,each = length(wt_cls_score))
    wt_score = rep(wt_cls_score,length(df_clone$clone_assay))
    
    df_wt_score = data.frame(clone_assay = "WT",experimental_assay = "WT",cell_type = "WT",clone_type = "WT",nm_score = wt_score,multi_panel = clone_assay_wt,
                             sc_name = names(wt_score),stringsAsFactors = F)
    
    df_score_all = bind_rows(df_score,df_wt_score)
    
    cell_type_cols = c("WT" = "gray60","control" = "blue","host" = "black","KO" = "indianred3")
    
    n_cls = table(df_score_all$clone_assay)
    n_cls = n_cls[names(n_cls) != "WT"]
    
    df_text = data.frame(multi_panel = names(n_cls), n_cls = paste0("N = ",as.vector(n_cls)), x = rep(-9,length(n_cls)),y = rep(1,length(n_cls)))
    
    p <- ggplot(df_score_all, aes(x = nm_score,color = cell_type)) +
      geom_density(size = 1) +
      facet_wrap(multi_panel ~ .,nrow = 2) +
      geom_vline(xintercept = top_value,linetype = "dashed",col = "gray60") +
      ggtitle("Nascent mesoderm score distribution among nascent mesoderm single cells") + scale_color_manual() +
      scale_color_manual(values = cell_type_cols) +
      geom_text(data = df_text,mapping = aes(x = x,y =y,label = n_cls),col = "black",hjust = 0)
    
    print(p)
    
    if (plot_pdf) {
      fn = "figs/paper_figs/fig3/nm_score_per_clone_assay.pdf"
    } else {
      fn = "figs/paper_figs/fig3/nm_score_per_clone_assay.png"
    }
    ggsave(filename = fn,plot = p,width = 9,h = 4)
    
  }
  
  plot_heatmap_diff_genes_epiblast_score = function(plot_pdf = F) {
    
    mat = scdb_mat("sing_emb_wt10")
    mc = scdb_mc("sing_emb_wt10_recolored")
    
    wt_cls = names(mc@mc)[mc@colors[mc@mc] == mc@color_key$color[1]]
    wt_cls_score = log2(colSums(mat@mat[epiblast_genes(),wt_cls])/colSums(mat@mat[,wt_cls]) + 1e-3)
    
    n_top = 500
    top_value = sort(wt_cls_score,decreasing = T)[n_top]
    
    cmp = cmp_epiblast_score()
    df_score = cmp$df_score
    df_clone = cmp$df_clone
    
    
    df_wt_score_2 = data.frame(clone_assay = "WT",experimental_assay = "WT",cell_type = "WT",clone_type = "WT",epiblast_score = wt_cls_score,sc_names = names(wt_cls_score),multi_panel = "WT",stringsAsFactors = F)
    
    df_score_all_2 = rbind(df_score,df_wt_score_2)
    
    # exclude three conditions
    excluded_clone_assay = c("Ctrl2 Chimera","Ctrl2 Tetraploid","TKO26 Chimera")
    
    df_score_all_2 = filter(df_score_all_2,!(df_score_all_2$clone_assay %in% excluded_clone_assay))
    
    
    top_cells = filter(df_score_all_2,df_score_all_2$epiblast_score >= top_value)
    
    #top_cells$clone_assay[top_cells$clone_assay %in% c("Ctrl1 Chimera","Ctrl1 Tetraploid")] = "Ctrl1 Chimera/Tetraploid"
    
    
    mat_chim = scdb_mat("tko_chimera")
    mat_tetra_tko = scdb_mat("tko_tetra")
    mat_tetra_control = scdb_mat("control_tetra_all")
    
    mat_all = cbind(mat@mat,mat_chim@mat,mat_tetra_tko@mat,mat_tetra_control@mat)
    
    egc_clone = t(tgs_matrix_tapply(mat_all[,top_cells$sc_names],top_cells$clone_assay,sum))
    n_umi_clone = colSums(egc_clone)
    egc_clone = t(t(egc_clone)/colSums(egc_clone))
    
    legc_clone = log2(egc_clone + 5e-5)
    
    sd_egc = sqrt(egc_clone[,"WT"]*1e5)/1e5
    
    lfp_top = legc_clone - legc_clone[,"WT"]
    
    # consistency among TKO clones
    lfp_tetra = rowMeans(lfp_top[,c("TKO26 Tetraploid","TKO29 Tetraploid")])
    lfp_chim = rowMeans(lfp_top[,c("TKO23 Chimera","TKO29 Chimera")])
    lfp_ctrl = rowMeans(lfp_top[,c("Ctrl1 Chimera","Ctrl1 Tetraploid")])
    
    thr_chim = 0.6
    thr_tetra = 0.6
    f_tetra_high = lfp_tetra > thr_tetra
    f_tetra_low = lfp_tetra < -thr_tetra
    f_chim_high = lfp_chim > thr_chim
    f_chim_low = lfp_chim < -thr_chim
    
    lfp_tetra_diff1 = lfp_tetra - lfp_ctrl
    lfp_tetra_diff2 = lfp_tetra - lfp_top[,"Host Chimera"]
    lfp_chim_diff1 = lfp_chim - lfp_ctrl
    lfp_chim_diff2 = lfp_chim - lfp_top[,"Host Chimera"]
    
    
    thr_diff = 0.3
    
    f_tetra_diff1 = abs(lfp_tetra_diff1) > thr_diff
    f_tetra_diff2 = abs(lfp_tetra_diff2) > thr_diff
    
    f_chim_diff1 = abs(lfp_chim_diff1) > thr_diff
    f_chim_diff2 = abs(lfp_chim_diff2) > thr_diff
    
    genes_f_chim = rownames(legc_clone)[(f_chim_high | f_chim_low) & f_chim_diff1 & f_chim_diff2]
    genes_f_tetra = rownames(legc_clone)[(f_tetra_high | f_tetra_low) & f_tetra_diff1 & f_tetra_diff2]
    
    genes_f = unique(c(genes_f_chim,genes_f_tetra))
    
    shades = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(100))
    breaks = seq(-1.5,1.5,length.out = 101)
    
    bad_genes = read.table("data/tet_tko.bad_genes.txt",sep = "\t",h = T,stringsAsFactors = F)$x
    bad_genes = c(bad_genes,c("Igf2","AK145379;H19","Polg","Slc25a4","Peg10","Igf2as","AK086477;Sh3glb1","Grb10","Nnat;Peg5"))
    genes_f = setdiff(genes_f,bad_genes)
    
    km_cl = tglkmeans::TGL_kmeans(df = as.data.frame(lfp_top[genes_f,]),k = 10,id_column = F)
    names(km_cl$cluster) = genes_f 
    
    genes_f = genes_f[order(km_cl$cluster)]
    
    if(plot_pdf) {
      fn = "figs/paper_figs/fig3/epiblast_heatmap.pdf"
    } else {
      fn = "figs/paper_figs/fig3/epiblast_heatmap.png"
    }
    
    lfp_top = lfp_top[,c("Host Chimera","Ctrl1 Chimera","Ctrl1 Tetraploid","TKO23 Chimera","TKO29 Chimera","TKO29 Tetraploid","TKO26 Tetraploid")]
    
    gaps_row = which(diff(km_cl$cluster[genes_f]) > 0)
    pheatmap::pheatmap(pmin(pmax(lfp_top[genes_f,],-1.5),1.5),
                       col = shades,breaks = breaks,
                       cluster_cols = F,cluster_rows = F,
                       filename = fn,w = 7,h = 20,gaps_row = gaps_row)
    
    write.table(x = genes_f,file = 'data/fig3/epiblast_diff_genes.txt',sep = '\t')
    
  } 
  
  
  plot_heatmap_diff_genes_epiblast_time_matched = function() {
    
    reg = 1e-5
    n_min_per_embryo = 10
    fig_dir = "figs/paper_figs/fig3/epiblast_expression_vs_time"
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    
    mat_wt = scdb_mat("sing_emb_wt10")
    mc_wt = scdb_mc("sing_emb_wt10_recolored")
    tko_chim = scdb_mat("tko_chim_wt10")
    tko_tetra = scdb_mat("tko_tetra_wt10")
    control_tetra_all = scdb_mat("control_tetra_all_wt10")
    
    load("data/tko_chim_wt10/tko_chim_md.Rda")
    load("data/tko_tetra_wt10/tko_tetra_md.Rda")
    load("data/control_tetra_all_wt10/control_tetra_all_md.Rda")
    time_chim = read.table("data/tko_chim_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
    rownames(time_chim) = time_chim$embryo
    time_tetra = read.table("data/tko_tetra_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
    rownames(time_tetra) = time_tetra$embryo
    time_control = read.table("data/control_tetra_all_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
    rownames(time_control) = time_control$embryo
    
    
    age_group_time = read.table(file = "data/wt10_transcriptional_rank_developmental_time.txt",h = T,sep = "\t",stringsAsFactors = F)
    age_group_time = age_group_time$developmental_time
    
    f_chim = (tko_chim_md$cell_type == "KO") & (tko_chim_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(tko_chim_md$age_group))
    n_cls_tko_chim_per_embryo =  table(tko_chim_md$embryo[f_chim])
    embryos_f = names(n_cls_tko_chim_per_embryo)[n_cls_tko_chim_per_embryo >= n_min_per_embryo]
    f_chim = f_chim & (tko_chim_md$embryo %in% embryos_f)
    e_t_chim = t(tgs_matrix_tapply(tko_chim@mat[,tko_chim_md$cell[f_chim]],tko_chim_md$embryo[f_chim],sum))
    e_t_chim = t(t(e_t_chim)/colSums(e_t_chim))
    
    f_chim_host_control = (tko_chim_md$cell_type %in% c("control","host")) & (tko_chim_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(tko_chim_md$age_group))
    n_cls_control_chim_per_embryo =  table(tko_chim_md$embryo[f_chim_host_control])
    embryos_f = names(n_cls_control_chim_per_embryo)[n_cls_control_chim_per_embryo >= n_min_per_embryo]
    f_chim_host_control = f_chim_host_control & (tko_chim_md$embryo %in% embryos_f)
    e_t_chim_control = t(tgs_matrix_tapply(tko_chim@mat[,tko_chim_md$cell[f_chim_host_control]],tko_chim_md$embryo[f_chim_host_control],sum))
    e_t_chim_control = t(t(e_t_chim_control)/colSums(e_t_chim_control))
    
    f_tetra = (tko_tetra_md$cell_type == "KO") & (tko_tetra_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(tko_tetra_md$age_group))
    n_cls_tko_tetra_per_embryo =  table(tko_tetra_md$embryo[f_tetra])
    embryos_f = names(n_cls_tko_tetra_per_embryo)[n_cls_tko_tetra_per_embryo >= n_min_per_embryo]
    f_tetra = f_tetra & (tko_tetra_md$embryo %in% embryos_f)
    e_t_tetra = t(tgs_matrix_tapply(tko_tetra@mat[,tko_tetra_md$cell[f_tetra]],tko_tetra_md$embryo[f_tetra],sum))
    e_t_tetra = t(t(e_t_tetra)/colSums(e_t_tetra))
    
    f_control = (control_tetra_all_md$cell_type == "control") & (control_tetra_all_md$ct_color == mc_wt@color_key$color[1]) & (!is.na(control_tetra_all_md$age_group))
    n_cls_control_tetra_per_embryo =  table(control_tetra_all_md$embryo[f_control])
    embryos_f = names(n_cls_control_tetra_per_embryo)[n_cls_control_tetra_per_embryo >= n_min_per_embryo]
    f_control = f_control & (control_tetra_all_md$embryo %in% embryos_f)
    e_t_control = t(tgs_matrix_tapply(control_tetra_all@mat[,control_tetra_all_md$cell[f_control]],control_tetra_all_md$embryo[f_control],sum))
    e_t_control = t(t(e_t_control)/colSums(e_t_control))
    
    
    cls_wt = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] == mc_wt@color_key$color[1]]
    n_cls_wt_per_embryo = table(mat_wt@cell_metadata[cls_wt,"transcriptional_rank"])
    ranks_included = as.numeric(names(n_cls_wt_per_embryo)[n_cls_wt_per_embryo > n_min_per_embryo])
    cls_wt = cls_wt[mat_wt@cell_metadata[cls_wt,"transcriptional_rank"] %in% ranks_included]
    e_t_wt = t(tgs_matrix_tapply(mat_wt@mat[,cls_wt],mat_wt@cell_metadata[cls_wt,"transcriptional_rank"],sum))
    e_t_wt = t(t(e_t_wt)/colSums(e_t_wt))
    
    genes_table = read.csv("data/gene_list_fig3.csv",stringsAsFactors = F)
    diff_genes = read.table(file = 'data/fig3/epiblast_diff_genes.txt',sep = '\t',stringsAsFactors = F)$x
    genes_f = unique(c(c("Pou3f1","Sox11","Tcf3","Foxp1","Tdgf1","Nodal"),genes_table$Gene.name,diff_genes))
    
    le_t_chim = log2(e_t_chim + reg)
    le_t_chim_control = log2(e_t_chim_control + reg)
    le_t_tetra = log2(e_t_tetra + reg)
    le_t_control = log2(e_t_control + reg)
    le_t_wt = log2(e_t_wt+ reg)
    
    roll_width = 6
    
    
    
  }
  
  
  plot_epiblast_score_single_genes = function(plot_pdf = F) {
    
    
    all_genes = c("Pou3f1","Sox11","Foxp1","Tcf3","Nodal")
    fig_dir = "figs/paper_figs/fig3/genes_epiblast_score"
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    for (g in all_genes) {
      print(g)
      mc_gg_chim_tetra_without_control(genes1 = epiblast_genes(),genes2 = g,xlab = "Epiblast score",ylab = g,fig_dir = fig_dir,main_tag = g,plot_pdf = plot_pdf,fig_width = 6,fig_height = 3.5)
    }
    
    
  }
  
  
  plot_heatmap_diff_genes_nascent_mesoderm_score = function(plot_pdf = F) {
    
    
    mat = scdb_mat("sing_emb_wt10")
    mc = scdb_mc("sing_emb_wt10_recolored")
    
    wt_cls = names(mc@mc)[mc@colors[mc@mc] == mc@color_key$color[10]]
    wt_cls_score = log2(colSums(mat@mat[early_nascent_mesoderm_genes(),wt_cls])/colSums(mat@mat[,wt_cls]) + 1e-3)
    
    n_top = 500
    top_value = sort(wt_cls_score,decreasing = T)[n_top]
    
    cmp = cmp_nascent_mesoderm_score()
    df_score = cmp$df_score
    df_clone = cmp$df_clone
    
    
    
    df_wt_score_2 = data.frame(clone_assay = "WT",experimental_assay = "WT",cell_type = "WT",clone_type = "WT",nm_score = wt_cls_score,sc_names = names(wt_cls_score),multi_panel = "WT",stringsAsFactors = F)
    
    df_score_all_2 = rbind(df_score,df_wt_score_2)
    
    # exclude three conditions
    excluded_clone_assay = c("Ctrl1 Chimera","Ctrl2 Chimera","Ctrl2 Tetraploid")
    
    df_score_all_2 = filter(df_score_all_2,!(df_score_all_2$clone_assay %in% excluded_clone_assay))
    
    
    top_cells = filter(df_score_all_2,df_score_all_2$nm_score >= top_value)
    
    
    mat_chim = scdb_mat("tko_chimera")
    mat_tetra_tko = scdb_mat("tko_tetra")
    mat_tetra_control = scdb_mat("control_tetra_all")
    
    mat_all = cbind(mat@mat,mat_chim@mat,mat_tetra_tko@mat,mat_tetra_control@mat)
    
    egc_clone = t(tgs_matrix_tapply(mat_all[,top_cells$sc_names],top_cells$clone_assay,sum))
    n_umi_clone = colSums(egc_clone)
    egc_clone = t(t(egc_clone)/colSums(egc_clone))
    
    legc_clone = log2(egc_clone + 5e-5)
    
    sd_egc = sqrt(egc_clone[,"WT"]*1e5)/1e5
    
    lfp_top = legc_clone - legc_clone[,"WT"]
    
    lfp_tetra = rowMeans(lfp_top[,c("TKO26 Tetraploid","TKO29 Tetraploid")])
    lfp_chim = rowMeans(lfp_top[,c("TKO23 Chimera","TKO26 Chimera","TKO29 Chimera")])
    
    thr_chim = 0.6
    thr_tetra = 0.6
    
    f_tetra_high = lfp_tetra > thr_tetra
    f_tetra_low = lfp_tetra < -thr_tetra
    f_chim_high = lfp_chim > thr_chim
    f_chim_low = lfp_chim < -thr_chim
    
    lfp_tetra_diff1 = lfp_tetra - lfp_top[,"Ctrl1 Tetraploid"]
    lfp_tetra_diff2 = lfp_tetra - lfp_top[,"Host Chimera"]
    lfp_chim_diff1 = lfp_chim - lfp_top[,"Ctrl1 Tetraploid"]
    lfp_chim_diff2 = lfp_chim - lfp_top[,"Host Chimera"]
    
    thr_diff = 0.3
    
    f_tetra_diff1 = abs(lfp_tetra_diff1) > thr_diff
    f_tetra_diff2 = abs(lfp_tetra_diff2) > thr_diff
    
    f_chim_diff1 = abs(lfp_chim_diff1) > thr_diff
    f_chim_diff2 = abs(lfp_chim_diff2) > thr_diff
    
    genes_f_chim = rownames(legc_clone)[(f_chim_high | f_chim_low) & f_chim_diff1 & f_chim_diff2]
    genes_f_tetra = rownames(legc_clone)[(f_tetra_high | f_tetra_low) & f_tetra_diff1 & f_tetra_diff2]
    
    genes_f = unique(c(genes_f_chim,genes_f_tetra))
    
    shades = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(100))
    breaks = seq(-1.5,1.5,length.out = 101)
    
    bad_genes = read.table("data/tet_tko.bad_genes.txt",sep = "\t",h = T,stringsAsFactors = F)$x
    bad_genes = c(bad_genes,c("Igf2","AK145379;H19","Polg","Slc25a4","Peg10","Igf2as","AK086477;Sh3glb1","Grb10","Nnat;Peg5"))
    genes_f = setdiff(genes_f,bad_genes)
    
    km_cl = tglkmeans::TGL_kmeans(df = as.data.frame(lfp_top[genes_f,]),k = 10,id_column = F)
    names(km_cl$cluster) = genes_f 
    
    genes_f = genes_f[order(km_cl$cluster)]
    
    if(plot_pdf) {
      fn = "figs/paper_figs/fig3/nm_heatmap.pdf"
    } else {
      fn = "figs/paper_figs/fig3/nm_heatmap.png"
    }
    
    gaps_row = which(diff(km_cl$cluster[genes_f]) > 0)
    pheatmap::pheatmap(pmin(pmax(lfp_top[genes_f,c("WT","Ctrl1 Tetraploid","Host Chimera","TKO23 Chimera","TKO26 Chimera","TKO29 Chimera","TKO26 Tetraploid","TKO29 Tetraploid")],-1.5),1.5),
                       col = shades,breaks = breaks,
                       cluster_cols = F,cluster_rows = F,
                       filename = fn,w = 5,h = 15,gaps_row = gaps_row)
    
    write.table(x = genes_f,file = 'data/fig3/nascent_mesoderm_diff_genes.txt',sep = '\t')
    
  }
  
  plot_nascent_mesoderm_score_single_genes = function(plot_pdf = T) {
    
    mc_wt = scdb_mc("sing_emb_wt10_recolored")
    
    
    all_genes = c("Lefty2","Fgf15","Fgf4","Nodal","Cyp26a1","Aldh1a2","Hand1","Fzd3")
    fig_dir = "figs/paper_figs/fig3/genes_nm_score"
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    for (g in all_genes) {
      print(g)
      mc_gg_chim_tetra_without_control(genes1 = early_nascent_mesoderm_genes(),
                                       genes2 = g,
                                       xlab = "Nascent mesoderm score",
                                       ylab = g,
                                       fig_dir = fig_dir,
                                       main_tag = g,
                                       plot_pdf = plot_pdf,
                                       included_cts = mc_wt@color_key$color[c(8,c(10:23))],fig_width = 6,fig_height = 3.5)
    }
    
    
    
  }
  
  epiblast_score_distributions  = function() {
    
    mat = scdb_mat("sing_emb_wt10")
    mc = scdb_mc("sing_emb_wt10_recolored")
    
    wt_cls = names(mc@mc)[mc@colors[mc@mc] == mc@color_key$color[1]]
    wt_cls_score = log2(colSums(mat@mat[epiblast_genes(),wt_cls])/colSums(mat@mat[,wt_cls]) + 1e-3)
    
    n_top = 500
    top_value = sort(wt_cls_score,decreasing = T)[n_top]
    
    tko_chim = scdb_mat("tko_chim_wt10")
    tko_tetra = scdb_mat("tko_tetra_wt10")
    ctrl_tetra = scdb_mat("control_tetra_all_wt10")
    
    df_clone = data.frame(clone_assay = c("TKO26 Chimera","TKO26 Tetraploid","TKO29 Chimera","TKO29 Tetraploid","TKO23 Chimera","Host Chimera","Ctrl1 Chimera","Ctrl1 Tetraploid","Ctrl2 Chimera","Ctrl2 Tetraploid"),
                          cell_type =  c("KO"   ,"KO"   ,"KO"   ,"KO"   ,"KO"   ,"host","control","control","control","control"),
                          clone_type = c("TKO26","TKO26","TKO29","TKO29","TKO23","host","Ctrl1"  ,"Ctrl1"  ,"Ctrl2"  ,"Ctrl2"),
                          experimental_assay = c("Chimera","Tetraploid","Chimera","Tetraploid","Chimera","Chimera","Chimera","Tetraploid","Chimera","Tetraploid"),
                          mat_nm = c("tko_chim_wt10","tko_tetra_wt10","tko_chim_wt10","tko_tetra_wt10","tko_chim_wt10","tko_chim_wt10","tko_chim_wt10","control_tetra_all_wt10","tko_chim_wt10","control_tetra_all_wt10"),stringsAsFactors = F)
    
    epiblast_score = list()
    
    df_score = data.frame(clone_assay = c(),experimental_assay = c(),cell_type = c(),clone_type = c(),epiblast_score = c(),stringsAsFactors = F)
    
    for (i in 1:nrow(df_clone)) {
      print(i)
      mat_nm = df_clone$mat_nm[i]
      scmat = scdb_mat(mat_nm)
      cell_type = df_clone$cell_type[i]
      clone_type = df_clone$clone_type[i]
      tag = df_clone$clone_assay[i]
      experimental_assay = df_clone$experimental_assay[i]
      
      query_cls = colnames(scmat@mat)[scmat@cell_metadata[colnames(scmat@mat),"cell_type"] == cell_type]
      query_cls = query_cls[scmat@cell_metadata[query_cls,"clone_type"] == clone_type]
      query_score = query_cmp_sc_score_per_clone(mat_nm = mat_nm,cls_f = query_cls,ct_col = mc@color_key$color[1],genes = epiblast_genes())
      
      epiblast_score[[tag]] = query_score
      
      df_tmp = data.frame(clone_assay = rep(tag,length(query_score)),
                          experimental_assay = rep(experimental_assay,length(query_score)),
                          cell_type = rep(cell_type,length(query_score)),
                          clone_type = rep(clone_type,length(query_score)),
                          epiblast_score = query_score,
                          sc_names = names(query_score),stringsAsFactors = F)
      
      df_score = rbind(df_score,df_tmp)
    }
    
    df_score$multi_panel = df_score$clone_assay
    
    
    # add wildtype cells to the score
    clone_assay_wt = rep(df_clone$clone_assay,each = length(wt_cls_score))
    wt_score = rep(wt_cls_score,length(df_clone$clone_assay))
    
    df_wt_score = data.frame(clone_assay = "WT",experimental_assay = "WT",cell_type = "WT",clone_type = "WT",epiblast_score = wt_score,multi_panel = clone_assay_wt,
                             sc_name = names(wt_score),stringsAsFactors = F)
    
    df_score_all = bind_rows(df_score,df_wt_score)
    
    cell_type_cols = c("WT" = "gray60","control" = "blue","host" = "black","KO" = "indianred3")
    
    n_cls = table(df_score_all$clone_assay)
    n_cls = n_cls[names(n_cls) != "WT"]
    
    df_text = data.frame(multi_panel = names(n_cls), n_cls = paste0("N = ",as.vector(n_cls)), x = rep(-9,length(n_cls)),y = rep(1,length(n_cls)))
    
    p <- ggplot(df_score_all, aes(x = epiblast_score,color = cell_type)) +
      geom_density(size = 1) +
      facet_wrap(multi_panel ~ .,nrow = 2) +
      geom_vline(xintercept = top_value,linetype = "dashed",col = "gray60") +
      ggtitle("Epiblast score distribution among epiblast single cells") + scale_color_manual() +
      scale_color_manual(values = cell_type_cols) +
      geom_text(data = df_text,mapping = aes(x = x,y =y,label = n_cls),col = "black",hjust = 0)
    
    print(p)
    
    ggsave(filename = "figs/paper_figs/fig3/epiblast_score_per_clone_assay.png",plot = p,width = 9,h = 4)
    
    
    df_wt_score_2 = data.frame(clone_assay = "WT",experimental_assay = "WT",cell_type = "WT",clone_type = "WT",epiblast_score = wt_cls_score,sc_names = names(wt_cls_score),multi_panel = "WT",stringsAsFactors = F)
    
    df_score_all_2 = rbind(df_score,df_wt_score_2)
    
    # exclude three conditions
    excluded_clone_assay = c("Ctrl2 Chimera","Ctrl2 Tetraploid","TKO26 Chimera")
    
    df_score_all_2 = filter(df_score_all_2,!(df_score_all_2$clone_assay %in% excluded_clone_assay))
    
    
    top_cells = filter(df_score_all_2,df_score_all_2$epiblast_score >= top_value)
    
    top_cells$clone_assay[top_cells$clone_assay %in% c("Ctrl1 Chimera","Ctrl1 Tetraploid")] = "Ctrl1 Chimera/Tetraploid"
    
    
    mat_chim = scdb_mat("tko_chimera")
    mat_tetra_tko = scdb_mat("tko_tetra")
    mat_tetra_control = scdb_mat("control_tetra_all")
    
    mat_all = cbind(mat@mat,mat_chim@mat,mat_tetra_tko@mat,mat_tetra_control@mat)
    
    egc_clone = t(tgs_matrix_tapply(mat_all[,top_cells$sc_names],top_cells$clone_assay,sum))
    n_umi_clone = colSums(egc_clone)
    egc_clone = t(t(egc_clone)/colSums(egc_clone))
    
    legc_clone = log2(egc_clone + 5e-5)
    
    sd_egc = sqrt(egc_clone[,"WT"]*1e5)/1e5
    
    lfp_top = legc_clone - legc_clone[,"WT"]
    
    lfp_min = apply(lfp_top[,c("TKO23 Chimera","TKO26 Tetraploid","TKO29 Chimera","TKO29 Tetraploid")],1,min) 
    lfp_max = apply(lfp_top[,c("TKO23 Chimera","TKO26 Tetraploid","TKO29 Chimera","TKO29 Tetraploid")],1,max) 
    
    # consistency among TKO clones
    f_high = lfp_min > 0.35
    f_low = lfp_max < -0.35
    
    lfp_tko_mean = rowMeans(lfp_top[,c("TKO23 Chimera","TKO26 Tetraploid","TKO29 Chimera","TKO29 Tetraploid")])
    
    f_mean = abs(lfp_tko_mean) > 0.4
    
    lfp_diff = lfp_tko_mean - lfp_top[,"Ctrl1 Chimera/Tetraploid"]
    f_diff = abs(lfp_diff) > 0.2
    lfp_diff2 = lfp_tko_mean - lfp_top[,"Host Chimera"]
    f_diff2 = abs(lfp_diff2) > 0.2
    genes_high = names(lfp_min)[f_high & f_mean & f_diff & f_diff2]
    genes_low = names(lfp_max)[f_low & f_mean & f_diff & f_diff2]
    
    
    
    genes_f = c(genes_high,genes_low)
    
    shades = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(100))
    breaks = seq(-1.5,1.5,length.out = 101)
    
    bad_genes = read.table("data/tet_tko.bad_genes.txt",sep = "\t",h = T,stringsAsFactors = F)$x
    genes_f = setdiff(genes_f,bad_genes)
    
    pheatmap::pheatmap(pmin(pmax(lfp_top[genes_f,],-1.5),1.5),
                       col = shades,breaks = breaks,
                       cluster_cols = F,treeheight_row = 0,
                       filename = "figs/paper_figs/fig3/epiblast_heatmap.png",w = 5,h = 15)
    
    
    gg_cor = tgs_cor(t(log2(mc@e_gc + 1e-5)))
    
    additional_genes = c("Foxd3","Dnmt3b","Nodal","Lefty2","Dkk1")
    
    all_genes = c(genes_f,additional_genes)
    fig_dir = "figs/paper_figs/fig3/genes"
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    for (g in all_genes) {
      print(g)
      mc_gg_chim_tetra_without_control(genes1 = epiblast_genes(),genes2 = g,xlab = "Epiblast score",ylab = g,fig_dir = fig_dir,main_tag = g)
    }
    
    fig_dir = "figs/paper_figs/fig3/genes2"
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    for (g in all_genes) {
      print(g)
      mc_gg_chim_tetra_without_control_2(genes1 = epiblast_genes(),genes2 = g,xlab = "Epiblast score",ylab = g,fig_dir = fig_dir,main_tag = g)
    }
  }
  
  
  
  query_cmp_sc_score_per_clone = function(mat_nm,cls_f,ct_col,genes) {
    
    mat = scdb_mat(mat_nm)
    
    load(sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm))
    
    cls_annot = cmp_annot$query_cls_col
    
    cls_f = cls_f[cls_f %in% names(cls_annot)]
    cls_f = cls_f[cls_annot[cls_f] == ct_col]
    
    # epiblast
    score_dist = log2(colSums(mat@mat[genes,cls_f])/colSums(mat@mat[,cls_f]) + 1e-3)
    
    return(score_dist)  
  }
  
  
  nascent_mesoderm_differential_expression = function() {
    
    mat = scdb_mat("sing_emb_wt10")
    mc = scdb_mc("sing_emb_wt10_recolored")
    
    wt_cls = names(mc@mc)[mc@colors[mc@mc] == mc@color_key$color[10]]
    wt_cls_score = log2(colSums(mat@mat[early_nascent_mesoderm_genes(),wt_cls])/colSums(mat@mat[,wt_cls]) + 1e-3)
    
    n_top = 500
    top_value = sort(wt_cls_score,decreasing = T)[n_top]
    
    tko_chim = scdb_mat("tko_chim_wt10")
    tko_tetra = scdb_mat("tko_tetra_wt10")
    ctrl_tetra = scdb_mat("control_tetra_all_wt10")
    
    df_clone = data.frame(clone_assay = c("TKO26 Chimera","TKO26 Tetraploid","TKO29 Chimera","TKO29 Tetraploid","TKO23 Chimera","Host Chimera","Ctrl1 Chimera","Ctrl1 Tetraploid","Ctrl2 Chimera","Ctrl2 Tetraploid"),
                          cell_type =  c("KO"   ,"KO"   ,"KO"   ,"KO"   ,"KO"   ,"host","control","control","control","control"),
                          clone_type = c("TKO26","TKO26","TKO29","TKO29","TKO23","host","Ctrl1"  ,"Ctrl1"  ,"Ctrl2"  ,"Ctrl2"),
                          experimental_assay = c("Chimera","Tetraploid","Chimera","Tetraploid","Chimera","Chimera","Chimera","Tetraploid","Chimera","Tetraploid"),
                          mat_nm = c("tko_chim_wt10","tko_tetra_wt10","tko_chim_wt10","tko_tetra_wt10","tko_chim_wt10","tko_chim_wt10","tko_chim_wt10","control_tetra_all_wt10","tko_chim_wt10","control_tetra_all_wt10"),stringsAsFactors = F)
    
    nm_score = list()
    
    df_score = data.frame(clone_assay = c(),experimental_assay = c(),cell_type = c(),clone_type = c(),nm_score = c(),stringsAsFactors = F)
    
    for (i in 1:nrow(df_clone)) {
      print(i)
      mat_nm = df_clone$mat_nm[i]
      scmat = scdb_mat(mat_nm)
      cell_type = df_clone$cell_type[i]
      clone_type = df_clone$clone_type[i]
      tag = df_clone$clone_assay[i]
      experimental_assay = df_clone$experimental_assay[i]
      
      query_cls = colnames(scmat@mat)[scmat@cell_metadata[colnames(scmat@mat),"cell_type"] == cell_type]
      query_cls = query_cls[scmat@cell_metadata[query_cls,"clone_type"] == clone_type]
      query_score = query_cmp_sc_score_per_clone(mat_nm = mat_nm,cls_f = query_cls,ct_col = mc@color_key$color[10],genes = early_nascent_mesoderm_genes())
      
      nm_score[[tag]] = query_score
      
      df_tmp = data.frame(clone_assay = rep(tag,length(query_score)),
                          experimental_assay = rep(experimental_assay,length(query_score)),
                          cell_type = rep(cell_type,length(query_score)),
                          clone_type = rep(clone_type,length(query_score)),
                          nm_score = query_score,
                          sc_names = names(query_score),stringsAsFactors = F)
      
      df_score = rbind(df_score,df_tmp)
    }
    
    df_score$multi_panel = df_score$clone_assay
    
    
    # add wildtype cells to the score
    clone_assay_wt = rep(df_clone$clone_assay,each = length(wt_cls_score))
    wt_score = rep(wt_cls_score,length(df_clone$clone_assay))
    
    df_wt_score = data.frame(clone_assay = "WT",experimental_assay = "WT",cell_type = "WT",clone_type = "WT",nm_score = wt_score,multi_panel = clone_assay_wt,
                             sc_name = names(wt_score),stringsAsFactors = F)
    
    df_score_all = bind_rows(df_score,df_wt_score)
    
    cell_type_cols = c("WT" = "gray60","control" = "blue","host" = "black","KO" = "indianred3")
    
    n_cls = table(df_score_all$clone_assay)
    n_cls = n_cls[names(n_cls) != "WT"]
    
    df_text = data.frame(multi_panel = names(n_cls), n_cls = paste0("N = ",as.vector(n_cls)), x = rep(-9,length(n_cls)),y = rep(1,length(n_cls)))
    
    p <- ggplot(df_score_all, aes(x = nm_score,color = cell_type)) +
      geom_density(size = 1) +
      facet_wrap(multi_panel ~ .,nrow = 2) +
      geom_vline(xintercept = top_value,linetype = "dashed",col = "gray60") +
      ggtitle("Nascent mesoderm score distribution among epiblast single cells") + scale_color_manual() +
      scale_color_manual(values = cell_type_cols) +
      geom_text(data = df_text,mapping = aes(x = x,y =y,label = n_cls),col = "black",hjust = 0)
    
    print(p)
    
    ggsave(filename = "figs/paper_figs/fig3/nm_score_per_clone_assay.png",plot = p,width = 9,h = 4)
    
    
    df_wt_score_2 = data.frame(clone_assay = "WT",experimental_assay = "WT",cell_type = "WT",clone_type = "WT",nm_score = wt_cls_score,sc_names = names(wt_cls_score),multi_panel = "WT",stringsAsFactors = F)
    
    df_score_all_2 = rbind(df_score,df_wt_score_2)
    
    # exclude three conditions
    excluded_clone_assay = c("Ctrl1 Chimera","Ctrl2 Chimera","Ctrl2 Tetraploid")
    
    df_score_all_2 = filter(df_score_all_2,!(df_score_all_2$clone_assay %in% excluded_clone_assay))
    
    
    top_cells = filter(df_score_all_2,df_score_all_2$nm_score >= top_value)
    
    
    mat_chim = scdb_mat("tko_chimera")
    mat_tetra_tko = scdb_mat("tko_tetra")
    mat_tetra_control = scdb_mat("control_tetra_all")
    
    mat_all = cbind(mat@mat,mat_chim@mat,mat_tetra_tko@mat,mat_tetra_control@mat)
    
    egc_clone = t(tgs_matrix_tapply(mat_all[,top_cells$sc_names],top_cells$clone_assay,sum))
    n_umi_clone = colSums(egc_clone)
    egc_clone = t(t(egc_clone)/colSums(egc_clone))
    
    legc_clone = log2(egc_clone + 5e-5)
    
    sd_egc = sqrt(egc_clone[,"WT"]*1e5)/1e5
    
    lfp_top = legc_clone - legc_clone[,"WT"]
    
    lfp_min_chim = apply(lfp_top[,c("TKO23 Chimera","TKO26 Chimera","TKO29 Chimera")],1,min) 
    lfp_max_chim = apply(lfp_top[,c("TKO23 Chimera","TKO26 Chimera","TKO29 Chimera")],1,max) 
    lfp_min_tetra = apply(lfp_top[,c("TKO26 Tetraploid","TKO29 Tetraploid")],1,min) 
    lfp_max_tetra = apply(lfp_top[,c("TKO26 Tetraploid","TKO29 Tetraploid")],1,max)
    
    # consistency among TKO clones
    f_high_chim = lfp_min_chim > 0.5
    f_low_chim = lfp_max_chim < -0.5
    f_high_tetra = lfp_min_tetra > 0.5
    f_low_tetra = lfp_max_tetra < -0.5
    
    lfp_tko_mean = rowMeans(lfp_top[,c("TKO23 Chimera","TKO26 Tetraploid","TKO29 Chimera","TKO29 Tetraploid")])
    
    f_mean = abs(lfp_tko_mean) > 0.4
    
    lfp_diff = lfp_tko_mean - lfp_top[,"Ctrl1 Tetraploid"]
    f_diff = abs(lfp_diff) > 0.3
    lfp_diff2 = lfp_tko_mean - lfp_top[,"Host Chimera"]
    f_diff2 = abs(lfp_diff2) > 0.3
    genes_high_chim = names(lfp_min_chim)[f_high_chim &  f_diff & f_diff2]
    genes_low_chim = names(lfp_max_chim)[f_low_chim &  f_diff & f_diff2]
    genes_high_tetra = names(lfp_min_tetra)[f_high_tetra &  f_diff & f_diff2]
    genes_low_tetra = names(lfp_max_tetra)[f_low_tetra &  f_diff & f_diff2]
    genes_high = union(genes_high_chim,genes_high_tetra)
    genes_low = union(genes_low_chim,genes_low_tetra)
    
    genes_f = c(genes_high,genes_low)
    
    shades = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu"))(100))
    breaks = seq(-1.5,1.5,length.out = 101)
    
    bad_genes = read.table("data/tet_tko.bad_genes.txt",sep = "\t",h = T,stringsAsFactors = F)$x
    genes_f = setdiff(genes_f,bad_genes)
    
    pheatmap::pheatmap(pmin(pmax(lfp_top[genes_f,c("WT","Ctrl1 Tetraploid","Host Chimera","TKO23 Chimera","TKO26 Chimera","TKO29 Chimera","TKO26 Tetraploid","TKO29 Tetraploid")],-1.5),1.5),
                       col = shades,breaks = breaks,
                       cluster_cols = F,treeheight_row = 0,
                       filename = "figs/paper_figs/fig3/nm_heatmap.png",w = 5,h = 20)
    
    
    all_genes = c("Hand1","Msx1","Msx2","Tbx3","Foxc1","Foxc2","Tbx6","Dll1","Dll3","Twist1","Gata4","Gata6","Lhx1","Snai1","Lefty2","Eomes","Tdgf1",
                  "Aldh1a2")
    fig_dir = "figs/paper_figs/fig3/nm_later_genes"
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    for (g in all_genes) {
      print(g)
      mc_gg_nascent_meso_score(genes2 = g,ylab = g,fig_dir = fig_dir,main_tag = g)
    }
    
    fig_dir = "figs/paper_figs/fig3/nm_score_genes"
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    all_genes = early_nascent_mesoderm_genes()
    for (g in all_genes) {
      print(g)
      mc_gg_nascent_meso_score(genes2 = g,ylab = g,fig_dir = fig_dir,main_tag = g)
    }
    for (g in all_genes) {
      print(g)
      mc_gg_chim_tetra_with_control(genes1 = epiblast_genes(),genes2 = g,xlab = "Epiblast score",
                                    ylab = g,fig_dir = fig_dir,main_tag = g)
    }
    fig_dir = "figs/paper_figs/fig3/nm_genes3"
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    all_genes = genes_f
    for (g in all_genes) {
      print(g)
      mc_gg_nascent_meso_score(genes2 = g,ylab = g,fig_dir = fig_dir,main_tag = g)
    }
    
    mc_gg_chim_tetra_with_control(genes1 = epiblast_genes(),genes2 = early_nascent_mesoderm_genes(),xlab = "Epiblast score",ylab = "Nascent mesoderm score",fig_dir = "figs/paper_figs/fig3",main_tag = "Epiblast score vs nascent mesoderm score")
    
    
  }


}




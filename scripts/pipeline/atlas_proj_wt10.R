

atlas_proj_on_wt10 = function(mat_query,fn,feat_genes, selected_metacells_for_plot = NULL,cex_points = 0.5,w = 1000,h = 1000, plot_pdf = F,main_tag = "",cex.main = 1,plot_gray_background = T) {
  
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  mc2d_wt = scdb_mc2d("sing_emb_wt10_recolored")
  
  legc_wt = log2(mc_wt@e_gc + 1e-5)
  
  query_ref_cor = tgs_cor(as.matrix(log2(mat_query[feat_genes,] + 1)),legc_wt[feat_genes,])
  
  best_ref = apply(query_ref_cor,1,which.max)
  
 

  
  all_metacells = 1:ncol(legc_wt)
  
  if (!is.null(selected_metacells_for_plot)) {
    
    all_metacells = selected_metacells_for_plot
    best_ref = best_ref[best_ref %in% selected_metacells_for_plot]
    
  }
  
  n = length(best_ref)
  
  xrange = 0.02 * (max(mc2d_wt@mc_x[all_metacells]) - min(mc2d_wt@mc_x[all_metacells]))
  yrange = 0.02 * (max(mc2d_wt@mc_y[all_metacells]) - min(mc2d_wt@mc_y[all_metacells]))
  ref_x = mc2d_wt@mc_x[best_ref] + rnorm(n, 0, xrange)
  ref_y = mc2d_wt@mc_y[best_ref] + rnorm(n, 0, yrange)
  xlim = c(min(mc2d_wt@mc_x[all_metacells]), max(mc2d_wt@mc_x[all_metacells]))
  ylim = c(min(mc2d_wt@mc_y[all_metacells]), max(mc2d_wt@mc_y[all_metacells]))
  if(plot_pdf) {
    pdf(file = fn,width = w,height = h,useDingbats = F)
  } else {
    png(filename = fn,width = w,height = h)
  }
  
  
  all_cells = names(mc_wt@mc)[mc_wt@mc %in% all_metacells]
  
  # before: gray95
  if(plot_gray_background) {
    plot(mc2d_wt@sc_x[all_cells],mc2d_wt@sc_y[all_cells], col = "gray90",pch = 19,cex = cex_points,xaxt = 'n',yaxt = 'n',xlab = "",ylab = "",axes = F,main = main_tag,
         cex.main = cex.main,ylim = ylim, xlim = xlim)
    points(ref_x, ref_y, pch = 19, col = mc_wt@colors[best_ref],cex = cex_points)
  } else {
    plot(x = ref_x, y = ref_y, pch = 19, col = mc_wt@colors[best_ref], 
           ylim = ylim, xlim = xlim,cex = cex_points,xaxt = 'n',yaxt = 'n',xlab = "",ylab = "",axes = F,main = main_tag,cex.main = cex.main)
  }
 
  dev.off()
  
  
}

plot_wt_manifold = function() {
  
  w = 1000/72
  h = 1000/72
  cex_points = 1
  fn = "figs/paper_figs/fig1/wt_atlas_2d_map.pdf"
  
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  mc2d_wt = scdb_mc2d("sing_emb_wt10_recolored")    
  

  xrange = 0.02 * (max(mc2d_wt@mc_x) - min(mc2d_wt@mc_x))
  yrange = 0.02 * (max(mc2d_wt@mc_y) - min(mc2d_wt@mc_y))

  xlim = c(min(mc2d_wt@mc_x), max(mc2d_wt@mc_x))
  ylim = c(min(mc2d_wt@mc_y), max(mc2d_wt@mc_y))
  
  pdf(file = fn,width = w,height = h,useDingbats = F)
  plot(mc2d_wt@sc_x[names(mc_wt@mc)],mc2d_wt@sc_y[names(mc_wt@mc)], col = mc_wt@colors[mc_wt@mc],pch = 19,cex = cex_points,xaxt = 'n',yaxt = 'n',xlab = "",ylab = "",axes = F,
       ylim = ylim, xlim = xlim)
  
  dev.off()
  
  fn = "figs/paper_figs/fig1/wt_atlas_2d_map_only_metacells.pdf"
  pdf(file = fn,width = w,height = h,useDingbats = F)
  cex_points = 1
  plot(mc2d_wt@mc_x,mc2d_wt@mc_y, col = mc_wt@colors,pch = 19,cex = cex_points,xaxt = 'n',yaxt = 'n',xlab = "",ylab = "",axes = F,
       ylim = ylim, xlim = xlim)
  
  segments_from_mc = mc2d_wt@graph$mc1
  segments_to_mc = mc2d_wt@graph$mc2
  
  segments(x0 = mc2d_wt@mc_x[segments_from_mc],y0 = mc2d_wt@mc_y[segments_from_mc],
           x1 = mc2d_wt@mc_x[segments_to_mc],y1 = mc2d_wt@mc_y[segments_to_mc] )
  
  points(mc2d_wt@mc_x,mc2d_wt@mc_y, bg = mc_wt@colors,pch = 21,cex = 4)
  dev.off()
  
}




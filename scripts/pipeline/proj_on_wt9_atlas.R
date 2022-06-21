

proj_on_wt9_atlas = function(mat_query) {
  
  mc_wt = scdb_mc("sing_emb_wt9_bs500f")
  gset = scdb_gset("sing_emb_wt9")
  feat_genes = names(gset@gene_set)
  
  egc_type = t(tgs_matrix_tapply(mc_wt@e_gc[feat_genes,],mc_wt@colors,mean))
  egc_type = egc_type[,colnames(egc_type) != "gray"]
  rownames(egc_type) = feat_genes
  
  legc = log2(egc_type + 1e-5)
  
  query_ref_cor = tgs_cor(as.matrix(mat_query@mat[feat_genes,]),egc_type)
  
  best_ref_type = colnames(query_ref_cor)[apply(query_ref_cor,1,which.max)]
  
  names(best_ref_type) = colnames(mat_query@mat)
  
  return(best_ref_type) 
}
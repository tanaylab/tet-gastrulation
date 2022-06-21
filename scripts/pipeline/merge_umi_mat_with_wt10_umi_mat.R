
library("Matrix")

merge_umi_mat_with_wt10 = function(scmat,min_umi = 2000) {
  
  mat1 = scmat
  mat2 = scdb_mat("sing_emb_wt10")
  
  mat1_n_umi = colSums(mat1@mat)
  mat1_cls_f = colnames(mat1@mat)[mat1_n_umi >= min_umi]
  mat2_n_umi = colSums(mat2@mat)
  mat2_cls_f = colnames(mat2@mat)[mat2_n_umi >= min_umi]
  
  mat1 = scm_ignore_cells(scmat = mat1,ig_cells = mat1_cls_f,reverse = T)
  mat2 = scm_ignore_cells(scmat = mat2,ig_cells = mat2_cls_f,reverse = T)
  
  mat1@cell_metadata$cell_type = as.character(mat1@cell_metadata$cell_type)
  
  ignored_cls = union(mat1@ignore_cells,mat2@ignore_cells)
  ignored_genes = union(mat1@ignore_genes,mat2@ignore_genes)
  
  mat_ls = c(mat1,mat2)
  
  mat_all = rbind(cbind(mat1@mat,mat1@ignore_cmat),cbind(mat1@ignore_gmat,mat1@ignore_gcmat))
  md_all = mat1@cell_metadata
  
  
  i = 2
  mat = mat_ls[[i]]
  mat_tmp = rbind(cbind(mat@mat,mat@ignore_cmat),cbind(mat@ignore_gmat,mat@ignore_gcmat))
  mat_tmp = mat_tmp[rownames(mat_all),]
  mat_all = cbind(mat_all,mat_tmp)
  md_all = bind_rows(md_all,mat@cell_metadata)
  
  rownames(md_all) = md_all$cell
  md_all[colnames(mat2@mat),"cell_type"] = "wt10"
  
  
  
  mat_new = scm_new_matrix(mat = mat_all, stat_type =  "umi",cell_metadata = md_all)
  mat_new = scm_ignore_cells(scmat = mat_new,ig_cells = ignored_cls)
  mat_new = scm_ignore_genes(scmat = mat_new,ig_genes = ignored_genes)
  
  
  
  return(mat_new)
}

merge_two_scmat = function(mat1,mat2,min_umi = 2000) {
  
  mat1_n_umi = colSums(mat1@mat)
  mat1_cls_f = colnames(mat1@mat)[mat1_n_umi >= min_umi]
  mat2_n_umi = colSums(mat2@mat)
  mat2_cls_f = colnames(mat2@mat)[mat2_n_umi >= min_umi]
  
  mat1 = scm_ignore_cells(scmat = mat1,ig_cells = mat1_cls_f,reverse = T)
  mat2 = scm_ignore_cells(scmat = mat2,ig_cells = mat2_cls_f,reverse = T)
  
  mat1@cell_metadata$cell_type = as.character(mat1@cell_metadata$cell_type)
  
  ignored_cls = union(mat1@ignore_cells,mat2@ignore_cells)
  ignored_genes = union(mat1@ignore_genes,mat2@ignore_genes)
  
  if((length(mat1@ignore_cells) == 0) & (length(mat1@ignore_genes) == 0)) {
    
    mat_all = mat1@mat
  }
  if((length(mat1@ignore_cells) > 0) & (length(mat1@ignore_genes) == 0)) {
    
    mat_all = cbind(mat1@mat,mat1@ignore_cmat)
  }
  if((length(mat1@ignore_cells) == 0) & (length(mat1@ignore_genes) > 0)) {
   
    mat_all = rbind(mat1@mat,mat1@ignore_gmat)
  }
  if((length(mat1@ignore_cells) > 0) & (length(mat1@ignore_genes) > 0)) {
    
    mat_all = rbind(cbind(mat1@mat,mat1@ignore_cmat),cbind(mat1@ignore_gmat,mat1@ignore_gcmat))
  }
  
  md_all = mat1@cell_metadata
  
  if((length(mat2@ignore_cells) == 0) & (length(mat2@ignore_genes) == 0)) {
    
    mat_tmp = mat2@mat
  }
  if((length(mat2@ignore_cells) > 0) & (length(mat2@ignore_genes) == 0)) {
    
    mat_tmp = cbind(mat2@mat,mat2@ignore_cmat)
  }
  if((length(mat2@ignore_cells) == 0) & (length(mat2@ignore_genes) > 0)) {
    
    mat_tmp = rbind(mat2@mat,mat2@ignore_gmat)
  }
  if((length(mat2@ignore_cells) > 0) & (length(mat2@ignore_genes) > 0)) {
    
    mat_tmp = rbind(cbind(mat2@mat,mat2@ignore_cmat),cbind(mat2@ignore_gmat,mat2@ignore_gcmat))
  }
  
  mat_tmp = mat_tmp[rownames(mat_all),]
  mat_all = cbind(mat_all,mat_tmp)
  md_all = bind_rows(md_all,mat2@cell_metadata)
  
  rownames(md_all) = md_all$cell
  
  mat_new = scm_new_matrix(mat = mat_all, stat_type =  "umi",cell_metadata = md_all)
  mat_new = scm_ignore_cells(scmat = mat_new,ig_cells = ignored_cls)
  mat_new = scm_ignore_genes(scmat = mat_new,ig_genes = ignored_genes)
  
  return(mat_new)
}
library("Matrix")
# removing all the other cells from the matrix object
# adding clone type information


preprocessing_tetraploid_control_mat_from_embflow_paper = function() {
  
  mat = scdb_mat("control_tetra_from_embflow_paper")
  
  cls_f = colnames(mat@mat)[colSums(mat@mat) >= 2000]
  
  mat_new = scm_ignore_cells(scmat = mat,ig_cells = cls_f,reverse = T)
  
  mat_new2 = scm_new_matrix(mat = mat_new@mat,cell_metadata = mat@cell_metadata[cls_f,],stat_type = mat@stat_type)
  
  mat_new2@ignore_genes = mat_new@ignore_genes
  mat_new2@ignore_gmat = mat_new@ignore_gmat
  
  mat_new2@cell_metadata$clone_type = "Ctrl_foxc"
  
  scdb_add_mat(id = "control_tetra_from_embflow_paper_f",mat = mat_new2)
  
}
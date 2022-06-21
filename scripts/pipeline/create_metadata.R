

create_tko_chimera_metadata = function() {
  
  mat = scdb_mat("tko_chim_wt10")
  load("data/tko_chim_wt10/color_annotation/cmp_annot.Rda")
  df_time = read.table("data/tko_chim_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  dev_time = read.table("data/wt10_transcriptional_rank_developmental_time.txt",sep = "\t",h = T,stringsAsFactors = F)
  mat_wt = scdb_mat("sing_emb_wt10")
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  age_group_rank = unique(mat_wt@cell_metadata[names(mc_wt@mc),c("transcriptional_rank","age_group")])
  
  cls = intersect(colnames(mat@mat),names(cmp_annot$query_cls_col))
  
  md = mat@cell_metadata[cls,c("cell","embryo","cell_type","clone_type")]
  
  df_color_ct = data.frame(cell = names(cmp_annot$query_cls_col),ct_color = cmp_annot$query_cls_col,stringsAsFactors = F)
  
  md = left_join(md,df_color_ct,by = "cell")
  
  chimera_time = df_time[,c("embryo","best_rank_host_control")]
  colnames(chimera_time)[2] = "transcriptional_rank"
  chimera_time = left_join(chimera_time,dev_time,by = "transcriptional_rank") 
  chimera_time = chimera_time %>% left_join(age_group_rank,by = "transcriptional_rank")
  
  md = md %>% left_join(chimera_time,by = "embryo")
  
  color_key = mc_wt@color_key[,c(1,2)]
  colnames(color_key) = c("infered_wt_cell_type","ct_color")
  
  md = md %>% left_join(color_key,by = "ct_color")
  
  tko_chim_md = md
  save(tko_chim_md,file = "data/tko_chim_wt10/tko_chim_md.Rda")
}

create_tko_tetraploid_metadata = function() {
  
  mat = scdb_mat("tko_tetra_wt10")
  load("data/tko_tetra_wt10/color_annotation/cmp_annot.Rda")
  df_time = read.table("data/tko_tetra_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  dev_time = read.table("data/wt10_transcriptional_rank_developmental_time.txt",sep = "\t",h = T,stringsAsFactors = F)
  mat_wt = scdb_mat("sing_emb_wt10")
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  age_group_rank = unique(mat_wt@cell_metadata[names(mc_wt@mc),c("transcriptional_rank","age_group")])
  
  cls = intersect(colnames(mat@mat),names(cmp_annot$query_cls_col))
  
  md = mat@cell_metadata[cls,c("cell","embryo","cell_type","clone_type")]
  
  df_color_ct = data.frame(cell = names(cmp_annot$query_cls_col),ct_color = cmp_annot$query_cls_col,stringsAsFactors = F)
  
  md = left_join(md,df_color_ct,by = "cell")
  
  tetraploid_time = df_time[,c("embryo","best_query")]
  colnames(tetraploid_time)[2] = "transcriptional_rank"
  tetraploid_time = left_join(tetraploid_time,dev_time,by = "transcriptional_rank") 
  tetraploid_time = tetraploid_time %>% left_join(age_group_rank,by = "transcriptional_rank")
  
  
  md = md %>% left_join(tetraploid_time,by = "embryo")
  
  color_key = mc_wt@color_key[,c(1,2)]
  colnames(color_key) = c("infered_wt_cell_type","ct_color")
  
  md = md %>% left_join(color_key,by = "ct_color")
  
  tko_tetra_md = md
  save(tko_tetra_md,file = "data/tko_tetra_wt10/tko_tetra_md.Rda")
}

create_control_tetraploid_metadata = function() {
  
  mat = scdb_mat("control_tetra_all_wt10")
  load("data/control_tetra_all_wt10/color_annotation/cmp_annot.Rda")
  df_time = read.table("data/control_tetra_all_wt10/time_match/time_match_summary.txt",sep = "\t",h = T,stringsAsFactors = F)
  dev_time = read.table("data/wt10_transcriptional_rank_developmental_time.txt",sep = "\t",h = T,stringsAsFactors = F)
  mat_wt = scdb_mat("sing_emb_wt10")
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  age_group_rank = unique(mat_wt@cell_metadata[names(mc_wt@mc),c("transcriptional_rank","age_group")])
  
  cls = intersect(colnames(mat@mat),names(cmp_annot$query_cls_col))
  
  md = mat@cell_metadata[cls,c("cell","embryo","cell_type","clone_type")]
  
  df_color_ct = data.frame(cell = names(cmp_annot$query_cls_col),ct_color = cmp_annot$query_cls_col,stringsAsFactors = F)
  
  md = left_join(md,df_color_ct,by = "cell")
  
  tetraploid_time = df_time[,c("embryo","best_query")]
  colnames(tetraploid_time)[2] = "transcriptional_rank"
  tetraploid_time = left_join(tetraploid_time,dev_time,by = "transcriptional_rank") 
  tetraploid_time = tetraploid_time %>% left_join(age_group_rank,by = "transcriptional_rank")
  
  
  md = md %>% left_join(tetraploid_time,by = "embryo")
  
  color_key = mc_wt@color_key[,c(1,2)]
  colnames(color_key) = c("infered_wt_cell_type","ct_color")
  
  md = md %>% left_join(color_key,by = "ct_color")
  
  control_tetra_all_md = md
  save(control_tetra_all_md,file = "data/control_tetra_all_wt10/control_tetra_all_md.Rda")
}
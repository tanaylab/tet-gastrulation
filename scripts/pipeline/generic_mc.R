

generate_mc = function(name,
		recompute=F, Knn = 100, Knn_core=30,
		min_mc_sz=20,
		T_vm = 0.1,
		add_bad_genes = NULL)
{
	mat_nm = name
	mc_nm = name
	mcf_nm = name
	coc_nm = name

	db_fn = .scdb_base
	if(recompute | !file.exists(sprintf("%s/mc.%s.Rda", db_fn, mcf_nm))) {
		mcell_add_gene_stat(mat_nm, mat_nm)

		mcell_gset_filter_varmean(mat_nm, mat_nm, T_vm=T_vm, force_new=T)
		mcell_gset_filter_cov(mat_nm, mat_nm, T_tot=50, T_top3=3)

		gset = scdb_gset(mat_nm)
		nms = names(gset@gene_set)
	  #bad gene that will be removed from list of genes that helps to mark metacell 
		bad_g = c(grep("^Rpl",nms,v=T),grep("^Gm",nms,v=T),grep("Rps",nms,v=T))
		if(!is.null(add_bad_genes)) {
			bad_g = c(bad_g, add_bad_genes)
		}
		gset_f = gset_new_restrict_nms(gset=gset, bad_g, inverse=T, "feat filt")
		scdb_add_gset(mat_nm, gset_f)

		mcell_add_cgraph_from_mat_bknn(mat_id=mat_nm, 
				gset_id = mat_nm, 
				graph_id=mat_nm,
				K=Knn,
				dsamp=T)

		mcell_coclust_from_graph_resamp(coc_id = coc_nm, graph_id = mat_nm,min_mc_size=15, p_resamp=0.75, n_resamp=500)

		mcell_mc_from_coclust_balanced(mc_nm, coc_nm, mat_nm, K=Knn_core, min_mc_size=min_mc_sz, alpha=2)

		mcell_mc_split_filt(new_mc_id=mcf_nm, mc_nm, mat_nm,T_lfc=3, plot_mats=F)
	}

  #list of genes that help to simply color the metacells and to lable them at the bottom of the heatmap 
	mcell_gset_from_mc_markers(gset_id=mcf_nm, mc_id=mcf_nm)

	wt_atlas = mcell_gen_atlas(mat_id = "sing_emb_wt10",
	                             mc_id = "sing_emb_wt10_recolored",
	                             gset_id = "sing_emb_wt10",
	                             mc2d_id = "sing_emb_wt10_recolored")
	
	#cmp = mcell_proj_on_atlas(mat_id = mat_nm,mc_id = mcf_nm,atlas = wt_atlas,fig_cmp_dir = paste0("figs/atlas_projection_",mat_nm),recolor_mc_id = mcf_nm,ten2mars = F)
}


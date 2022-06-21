
arsinh = function(x,a = 1,b = 1) {
  
  return(log(x + sqrt(x^2/a^2 + b/a)))
}

calc_a_b_from_x_y = function(x1,y1,x2,y2) {

  b = (y2 - y1)/(x2 - x1)
  a = (y1*x2 - x1*y2)/(x2 - x1) 
  return(c(a,b))
}

apply_vertical_gate = function(cls_type,v_l,v_h,tag_l,tag_h,x_val,y_val,cls_f,clone,sort_date,cls_color,xlab = "GFP",ylab = "tdTomato",xlim = c(-0.5,12),
                               ylim = c(1,11)) {
  
  f_low = x_val[cls_f] < v_l
  f_high = x_val[cls_f] > v_h
  
  f = !(cls_color %in% c("#F6BFCB","#7F6874","#1A1A1A","#989898"))
  ff = !f & (cls %in% cls_f)
  
  cls_type[cls_f[f_low]] = tag_l
  cls_type[cls_f[f_high]] = tag_h
  
  gate_color = rep("gray70",length(cls_f))
  gate_color[f_low] = "blue"
  gate_color[f_high] = "coral3"
  
  sort_date_nm = gsub("/","",sort_date)
  
  if(0) {
    png(sprintf("%s/%s_vs_%s_%s_%s.png",fig_dir,xlab,ylab,clone,sort_date_nm))
    plot(x_val[cls_f],y_val[cls_f],pch = 19,cex = 0.5,ylab = ylab,xlab = xlab,col = cls_color[cls_f],main = paste(clone, sort_date," chimera embryos",sep = " "),
         xlim =xlim,ylim = ylim)
    abline(v = v_l,lty = "dashed")
    abline(v = v_h,lty = "dashed")
    points(x_val[ff],y_val[ff],cex = 0.5,pch = 19,col = cls_color[ff])
    dev.off()
  }

  
  png(sprintf("%s/%s_vs_%s_%s_%s_gated.png",fig_dir,xlab,ylab,clone,sort_date_nm))
  plot(x_val[cls_f],y_val[cls_f],pch = 19,cex = 0.5,ylab = ylab,xlab = xlab,col = gate_color,main = paste(clone, sort_date," chimera embryos",sep = " "),
       xlim =xlim,ylim = ylim)
  abline(v = v_l,lty = "dashed")
  abline(v = v_h,lty = "dashed")
  legend(x = "topright",pch = 19,col = c("blue","coral3","gray70"),legend = c(tag_l,tag_h,"unclear"))
  dev.off()
  
  return(cls_type)
}


apply_horizontal_gate = function(cls_type,h_l,h_h,tag_l,tag_h,x_val,y_val,cls_f,clone,sort_date,cls_color,xlab = "GFP",ylab = "tdTomato",xlim = c(-0.5,12),
                                 ylim = c(1,11)) {
  
  f = !(cls_color %in% c("#F6BFCB","#7F6874","#1A1A1A","#989898"))
  ff = !f & (cls %in% cls_f)
  
  f_low = y_val[cls_f] < h_l
  f_high = y_val[cls_f] > h_h
  
  cls_type[cls_f[f_low]] = tag_l
  cls_type[cls_f[f_high]] = tag_h
  
  gate_color = rep("gray70",length(cls_f))
  gate_color[f_low] = "blue"
  gate_color[f_high] = "coral3"
  
  sort_date_nm = gsub("/","",sort_date)
  
  if(0) {
    png(sprintf("%s/%s_vs_%s_%s_%s_gated.png",fig_dir,xlab,ylab,clone,sort_date_nm))
    plot(x_val[cls_f],y_val[cls_f],pch = 19,cex = 0.5,ylab = ylab,xlab = xlab,col = cls_color[cls_f],main = paste(clone, sort_date," chimera embryos",sep = " "),
         xlim =xlim,ylim = ylim)
    abline(h = h_l,lty = "dashed")
    abline(h = h_h,lty = "dashed")
    points(x_val[ff],y_val[ff],cex = 0.5,pch = 19,col = cls_color[ff])
    dev.off()
  }

  
  png(sprintf("%s/%s_vs_%s_%s_%s_gated.png",fig_dir,xlab,ylab,clone,sort_date_nm))
  plot(x_val[cls_f],y_val[cls_f],pch = 19,cex = 0.5,ylab = ylab,xlab = xlab,col = gate_color,main = paste(clone, sort_date," chimera embryos",sep = " "),
       xlim =xlim,ylim = ylim)
  abline(h = h_l,lty = "dashed")
  abline(h = h_h,lty = "dashed")
  legend(x = "topright",pch = 19,col = c("blue","coral3","gray70"),legend = c(tag_l,tag_h,"unclear"))
  dev.off()
  
  return(cls_type)
}

apply_vertical_and_horizontal_gate = function(cls_type,h_l,h_h,v_l,v_h,tag_ll,tag_lh,tag_hl,x_val,y_val,cls_f,clone,sort_date,cls_color,xlab = "GFP",ylab = "tdTomato",xlim = c(-0.5,12),
                                              ylim = c(1,11)) {
  
  f = !(cls_color %in% c("#F6BFCB","#7F6874","#1A1A1A","#989898"))
  ff = !f & (cls %in% cls_f)
  
  f_ll = ( x_val[cls_f] < v_l ) & ( y_val[cls_f] < h_l )
  f_hl = ( x_val[cls_f] > v_h ) & ( y_val[cls_f] < h_l )
  f_lh = ( x_val[cls_f] < v_l ) & ( y_val[cls_f] > h_h )
  
  cls_type[cls_f[f_ll]] = tag_ll
  cls_type[cls_f[f_lh]] = tag_lh
  cls_type[cls_f[f_hl]] = tag_hl
  
  gate_color = rep("gray70",length(cls_f))
  gate_color[f_ll] = "blue"
  gate_color[f_lh] = "coral3"
  gate_color[f_hl] = "chartreuse3"
  
  sort_date_nm = gsub("/","",sort_date)
  
  if(0) {
    png(sprintf("%s/%s_vs_%s_%s_%s_gated.png",fig_dir,xlab,ylab,clone,sort_date_nm))
    plot(x_val[cls_f],y_val[cls_f],pch = 19,cex = 0.5,ylab = ylab,xlab = xlab,col = cls_color[cls_f],main = paste(clone, sort_date," chimera embryos",sep = " "),
         xlim =xlim,ylim = ylim)
    abline(v = v_l,lty = "dashed")
    abline(v = v_h,lty = "dashed")
    abline(h = h_l,lty = "dashed")
    abline(h = h_h,lty = "dashed")
    points(x_val[ff],y_val[ff],cex = 0.5,pch = 19,col = cls_color[ff])
    dev.off()
  }

  
  png(sprintf("%s/%s_vs_%s_%s_%s_gated.png",fig_dir,xlab,ylab,clone,sort_date_nm))
  plot(x_val[cls_f],y_val[cls_f],pch = 19,cex = 0.5,ylab = ylab,xlab = xlab,col = gate_color,main = paste(clone, sort_date," chimera embryos",sep = " "),
       xlim =xlim,ylim = ylim)
  abline(v = v_l,lty = "dashed")
  abline(v = v_h,lty = "dashed")
  abline(h = h_l,lty = "dashed")
  abline(h = h_h,lty = "dashed")
  legend(x = "topright",pch = 19,col = c("blue","coral3","chartreuse3","gray70"),legend = c(tag_ll,tag_lh,tag_hl,"unclear"))
  dev.off()
  
  return(cls_type)

}

experiment = "Chimera Assay"

mat_nm = "tko_chim"
mat = scdb_mat(mat_nm)

if(!dir.exists(sprintf("%s/gating",.scfigs_base))) {
  dir.create(sprintf("%s/gating",.scfigs_base))
}

if(!dir.exists(sprintf("%s/gating/%s",.scfigs_base,mat_nm))) {
  dir.create(sprintf("%s/gating/%s",.scfigs_base,mat_nm))
}
fig_dir = sprintf("%s/gating/%s",.scfigs_base,mat_nm)

cls_type = rep("unclear",length(colnames(mat@mat)))
names(cls_type) = colnames(mat@mat)

cls_clone = rep("unclear",length(colnames(mat@mat)))
names(cls_clone) = colnames(mat@mat)

# Gating of KO cells

cls = colnames(mat@mat)

cls_mcherry = cls[is.na(mat@cell_metadata[cls,"tdtomato_a"]) & !is.na(mat@cell_metadata[cls,"mcherry_a"])]
cls_tomato = cls[!is.na(mat@cell_metadata[cls,"tdtomato_a"])]


# color cells by wt9 atlas color

cls_color = proj_on_wt9_atlas(mat_query = mat)

f = !(cls_color %in% c("#F6BFCB","#7F6874","#1A1A1A","#989898"))

cls_color[f] = "gray80"

cls_color[cls_color == "#989898"] = "indianred3"

# Gating of KO cells

a_x = 1
b_x = 1000
a_y = 1
b_y = 1000

xlim = c(-0.5,12)
ylim = c(1,11)


gfp_tr = arsinh(mat@cell_metadata[cls,"gfp_a"], a= a_x,b = b_x)
tdtomato_tr = arsinh(mat@cell_metadata[cls,"tdtomato_a"],a = a_y,b = b_y)
mcherry_tr = arsinh(mat@cell_metadata[cls,"mcherry_a"],a = a_y,b = b_y)
ssc_a = mat@cell_metadata[cls,"ssc_a"]
fsc_a = mat@cell_metadata[cls,"fsc_a"]

col_dens = densCols(ssc_a,fsc_a)
names(gfp_tr) = cls
names(tdtomato_tr) = cls
names(mcherry_tr) = cls
names(ssc_a) = cls
names(fsc_a) = cls


png(sprintf("%s/GFP_vs_tdTomato_all_cls.png",fig_dir))
plot(gfp_tr,tdtomato_tr,pch = 19,cex = 0.2,ylab = "tdTomato",col = cls_color,xlab = "gfp",xlim = xlim,ylim =ylim)
legend(x = "topleft",pch = 19,col = c("#F6BFCB","#7F6874","#1A1A1A","indianred3","gray80"),legend = c("Visceral Endoderm","ExE Endoderm","Parietal Endoderm","ExE Ectoderm","Embryonic cell types"))
dev.off()


for (clone in c("TKO23","TKO26","TKO29")) {
  
  cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
  
  ff = !f & (cls %in% cls_f)
  
  png(sprintf("%s/GFP_vs_tdTomato_%s.png",fig_dir,clone))
  plot(gfp_tr[cls_f],tdtomato_tr[cls_f],pch = 19,cex = 0.5,ylab = "tdTomato",xlab = "GFP",col = cls_color[cls_f],main = paste(clone," chimera embryos",sep = " "),
       xlim =xlim,ylim = ylim)
  points(gfp_tr[ff],tdtomato_tr[ff],cex = 0.5,pch = 19,col = cls_color[ff])
  dev.off()
  
  sort_dates = unique(as.character(mat@cell_metadata[cls_f,"Sort.Date"]))
  
  if(clone == "TKO29") {
    sort_dates = setdiff(sort_dates,"06/11/2020")
  }
  
  for (sort_date in sort_dates) {
    
    cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]
    
    sort_date_nm = gsub("/","",sort_date)
    
    ff = !f & (cls %in% cls_ff)
    
    png(sprintf("%s/GFP_vs_tdTomato_%s_%s.png",fig_dir,clone,sort_date_nm))
    plot(gfp_tr[cls_ff],tdtomato_tr[cls_ff],pch = 19,cex = 0.5,ylab = "tdTomato",xlab = "GFP",col = cls_color[cls_ff],main = paste(clone, sort_date," chimera embryos",sep = " "),
         xlim =xlim,ylim = ylim)
    points(gfp_tr[ff],tdtomato_tr[ff],cex = 0.5,pch = 19,col = cls_color[ff])
    dev.off()
    
    
  }
  
}

clone = "TKO29"
sort_date = "06/11/2020"

cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_f = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

sort_date_nm = gsub("/","",sort_date)

ff = !f & (cls %in% cls_ff)

png(sprintf("%s/fsc_vs_tdTomato_%s_%s.png",fig_dir,clone,sort_date_nm))
plot(fsc_a[cls_ff],tdtomato_tr[cls_ff],pch = 19,cex = 0.5,ylab = "tdTomtato",xlab = "FSC_A",col = cls_color[cls_ff],main = paste(clone, sort_date," chimera embryos",sep = " "),
     ylim = ylim)
points(fsc_a[ff],tdtomato_tr[ff],cex = 0.5,pch = 19,col = cls_color[ff])
dev.off()

png(sprintf("%s/ssc_vs_tdTomato_%s_%s.png",fig_dir,clone,sort_date_nm))
plot(ssc_a[cls_ff],tdtomato_tr[cls_ff],pch = 19,cex = 0.5,ylab = "tdTomtato",xlab = "SSC_A",col = cls_color[cls_ff],main = paste(clone, sort_date," chimera embryos",sep = " "),
     ylim = ylim)
points(ssc_a[ff],tdtomato_tr[ff],cex = 0.5,pch = 19,col = cls_color[ff])
dev.off()

png(sprintf("%s/ssc_vs_fsc_%s_%s.png",fig_dir,clone,sort_date_nm))
plot(ssc_a[cls_ff],fsc_a[cls_ff],pch = 19,cex = 0.5,ylab = "FSC-A",xlab = "SSC_A",col = cls_color[cls_ff],
     main = paste(clone, sort_date," chimera embryos",sep = " "))
points(ssc_a[ff],fsc_a[ff],cex = 0.5,pch = 19,col = cls_color[ff])
dev.off()





# next follows gating of the cells

clone = "TKO23"
sort_date = "06/09/2020"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

v_l = 5
v_h = 6

cls_type = apply_vertical_gate(cls_type = cls_type,
                               v_l = v_l,v_h = v_h,
                               tag_l = "host",
                               tag_h = "KO",
                               x_val = gfp_tr,
                               y_val = tdtomato_tr,cls_f = cls_ff,
                               clone = clone,
                               sort_date = sort_date,
                               cls_color = cls_color)


clone = "TKO23"
sort_date = "23/10/2020"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

v_l = 7
v_h = 7.5

cls_type = apply_vertical_gate(cls_type = cls_type,
                               v_l = v_l,v_h = v_h,
                               tag_l = "host",
                               tag_h = "KO",
                               x_val = gfp_tr,
                               y_val = tdtomato_tr,cls_f = cls_ff,
                               clone = clone,
                               sort_date = sort_date,
                               cls_color = cls_color)

clone = "TKO26"
sort_date = "06/01/2020"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

v_l = 6.5
v_h = 8.3
h_l = 6.5
h_h = 7.6


cls_type = apply_vertical_and_horizontal_gate(cls_type = cls_type,
                               v_l = v_l,v_h = v_h,h_l = h_l,h_h = h_h,tag_ll = "host",tag_lh = "KO",tag_hl = "control",
                               x_val = gfp_tr,
                               y_val = tdtomato_tr,cls_f = cls_ff,
                               clone = clone,
                               sort_date = sort_date,
                               cls_color = cls_color)

clone = "TKO26"
sort_date = "08/11/2019"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

v_l = 6.5
v_h = 8.3
h_l = 7.8
h_h = 8.0


cls_type = apply_vertical_and_horizontal_gate(cls_type = cls_type,
                                              v_l = v_l,v_h = v_h,h_l = h_l,h_h = h_h,tag_ll = "host",tag_lh = "KO",tag_hl = "control",
                                              x_val = gfp_tr,
                                              y_val = tdtomato_tr,cls_f = cls_ff,
                                              clone = clone,
                                              sort_date = sort_date,
                                              cls_color = cls_color)



clone = "TKO26"
sort_date = "23/10/2020"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

h_l = 6.5
h_h = 7.3

cls_type = apply_horizontal_gate(cls_type = cls_type,
                               h_l = h_l,h_h = h_h,
                               tag_l = "host",
                               tag_h = "KO",
                               x_val = gfp_tr,
                               y_val = tdtomato_tr,cls_f = cls_ff,
                               clone = clone,
                               sort_date = sort_date,
                               cls_color = cls_color)

# ----------------------------------------------------------------------------------------
clone = "TKO29"
sort_date = "06/08/2019"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

v_l = 6.4
v_h = 8.0
h_l = 6.8
h_h = 7.5


cls_type = apply_vertical_and_horizontal_gate(cls_type = cls_type,
                                              v_l = v_l,v_h = v_h,h_l = h_l,h_h = h_h,
                                              tag_ll = "host",tag_lh = "KO",tag_hl = "control",
                                              x_val = gfp_tr,
                                              y_val = tdtomato_tr,cls_f = cls_ff,
                                              clone = clone,
                                              sort_date = sort_date,
                                              cls_color = cls_color)


clone = "TKO29"
sort_date = "06/11/2020"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]


h_l = 7
h_h = 8


cls_type = apply_horizontal_gate(cls_type = cls_type,
                                 h_l = h_l,h_h = h_h,
                                 tag_l = "host",tag_h = "KO",
                                 x_val = fsc_a,
                                 y_val = tdtomato_tr,cls_f = cls_ff,
                                 clone = clone,
                                 sort_date = sort_date,
                                 cls_color = cls_color,
                                 xlim = c(35000,180000),
                                 xlab = "FSC_A")

clone = "TKO29"
sort_date = "18/03/2019"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

v_l = 7
v_h = 9
h_l = 7.8
h_h = 8.4


cls_type = apply_vertical_and_horizontal_gate(cls_type = cls_type,
                                              v_l = v_l,v_h = v_h,h_l = h_l,h_h = h_h,
                                              tag_ll = "host",tag_lh = "KO",tag_hl = "control",
                                              x_val = gfp_tr,
                                              y_val = tdtomato_tr,cls_f = cls_ff,
                                              clone = clone,
                                              sort_date = sort_date,
                                              cls_color = cls_color)


clone = "TKO29"
sort_date = "22/03/2019"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

h_l = 7.9
h_h = 8.0
v_l = 6.8
v_h = 9.5


cls_type = apply_vertical_and_horizontal_gate(cls_type = cls_type,
                                              v_l = v_l,v_h = v_h,h_l = h_l,h_h = h_h,
                                              tag_ll = "host",tag_lh = "KO",tag_hl = "control",
                                              x_val = gfp_tr,
                                              y_val = tdtomato_tr,cls_f = cls_ff,
                                              clone = clone,
                                              sort_date = sort_date,
                                              cls_color = cls_color)

clone = "TKO29"
sort_date = "30/07/2019"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

h_l = 7.1
h_h = 7.8
v_l = 6.5
v_h = 8.2


cls_type = apply_vertical_and_horizontal_gate(cls_type = cls_type,
                                              v_l = v_l,v_h = v_h,h_l = h_l,h_h = h_h,
                                              tag_ll = "host",tag_lh = "KO",tag_hl = "control",
                                              x_val = gfp_tr,
                                              y_val = tdtomato_tr,cls_f = cls_ff,
                                              clone = clone,
                                              sort_date = sort_date,
                                              cls_color = cls_color)


clone = "TKO29"
sort_date = "25/03/2019"
cls_f = cls[grep(clone,mat@cell_metadata[cls,"embryo"])]
cls_ff = cls_f[mat@cell_metadata[cls_f,"Sort.Date"] == sort_date]

h_l = 7.9
h_h = 8.7
v_l = 6.5
v_h = 9


cls_type = apply_vertical_and_horizontal_gate(cls_type = cls_type,
                                              v_l = v_l,v_h = v_h,h_l = h_l,h_h = h_h,
                                              tag_ll = "host",tag_lh = "KO",tag_hl = "control",
                                              x_val = gfp_tr,
                                              y_val = tdtomato_tr,cls_f = cls_ff,
                                              clone = clone,
                                              sort_date = sort_date,
                                              cls_color = cls_color)


# remove embryos TKO29_10 TKO29_11 TKO29_12 TKO29_20 TKO29_21
# We originally forgot to gate and include them.
# For consistency with downstream analysis we set they there gating to unclear

ig_embryos = c("TKO29_10","TKO29_11","TKO29_12","TKO29_20","TKO29_21")
f = mat@cell_metadata[colnames(mat@mat),"embryo"] %in% ig_embryos
cls_type[f] = "unclear"

if(1) {
  
  df_gating = data.frame(cell = colnames(mat@mat),cell_type = cls_type,stringsAsFactors = F)
  mat = scdb_mat(mat_nm)
  md = mat@cell_metadata
  
  f = colnames(mat@cell_metadata) %in% c("cell_type","clone_type","cell_genotype")
  md = md[,!f]
  
  
  md = left_join(md,df_gating,by = "cell")
  rownames(md) = md$cell
  mat@cell_metadata = md
  
  mat = add_clone_type_information_to_scmat(mat = mat)
  
  scdb_add_mat(id = mat_nm,mat = mat)
}


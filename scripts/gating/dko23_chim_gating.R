
arsinh = function(x,a = 1,b = 1) {
  
  return(log(x + sqrt(x^2/a^2 + b/a)))
}

calc_a_b_from_x_y = function(x1,y1,x2,y2) {
  
  b = (y2 - y1)/(x2 - x1)
  a = (y1*x2 - x1*y2)/(x2 - x1) 
  return(c(a,b))
}

dko23_gating = function(plot_pdf = F) {


mat_nm = "dko23_chim"
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
names(cls_type) = colnames(mat@mat)

cls = colnames(mat@mat)

# color cells by wt9 atlas color

cls_color = proj_on_wt9_atlas(mat_query = mat)

f = !(cls_color %in% c("#F6BFCB","#7F6874","#1A1A1A","#989898"))

cls_color[f] = "gray80"

cls_color[cls_color == "#989898"] = "indianred3"

# parameters needed for the transformation
a_x = 1
b_x = 1000
a_y = 1
b_y = 1000




gfp_tr = arsinh(mat@cell_metadata[cls,"gfp_a"], a= a_x,b = b_x)
ssc_a = mat@cell_metadata[cls,"ssc_a"]
fsc_a = mat@cell_metadata[cls,"fsc_a"]

names(gfp_tr) = cls
names(ssc_a) = cls
names(fsc_a) = cls

lim_gfp = c(min(gfp_tr),max(gfp_tr))
lim_fsc = c(min(fsc_a),max(fsc_a))
lim_ssc = c(min(ssc_a),max(ssc_a))

png(sprintf("%s/gfp_vs_fsc_all.png",fig_dir))
plot(gfp_tr,fsc_a,pch = 19,cex = 0.5,xlab = "GFP",ylab = "FSC-A",col = cls_color,main = paste("All DKO23 chimera embryos",sep = " "),xlim = lim_gfp,ylim = lim_fsc)
points(gfp_tr[!f],fsc_a[!f],cex = 0.5,pch = 19,col = cls_color[!f])
dev.off()

png(sprintf("%s/gfp_vs_ssc_all.png",fig_dir))
plot(gfp_tr,ssc_a,pch = 19,cex = 0.5,xlab = "GFP",ylab = "SSC-A",col = cls_color,main = paste("All DKO23 chimera embryos",sep = " "),xlim = lim_gfp,ylim = lim_ssc)
points(gfp_tr[!f],ssc_a[!f],cex = 0.5,pch = 19,col = cls_color[!f])
dev.off()

# next gating of the cells

a_l = -400000
b_l =   80000
a_h = -500000
b_h =   80000


f_host = a_l + b_l * gfp_tr < fsc_a
f_ko   = a_h + b_h * gfp_tr > fsc_a
cls_type[f_host] = "host"
cls_type[f_ko] = "DKO23"

cls_clone[f_host]= "host"
cls_ko = cls[f_ko]
cls_clone[f_ko] = "DKO23_5"

gate_color = rep("gray70",length(cls))
gate_color[f_host] = "blue"
gate_color[f_ko] = "coral3"

if(plot_pdf) {

  pdf(sprintf("%s/gfp_vs_fsc_all_gated.pdf",fig_dir),useDingbats = F)
plot(gfp_tr,fsc_a,pch = 19,cex = 0.5,xlab = "GFP",ylab = "FSC-A",col = gate_color,main = paste("All DKO23 chimera embryos",sep = " "),xlim = lim_gfp,ylim = lim_fsc)
abline(a = a_l,b = b_l,lty = "dashed")
abline(a = a_h,b = b_h,lty = "dashed")
legend(x = "topright",pch = 19,col = c("blue","coral3","gray70"),legend = c("Host","DKO23","unclear"))
dev.off()

pdf(sprintf("%s/gfp_vs_fsc_all.pdf",fig_dir),useDingbats = F)
plot(gfp_tr,fsc_a,pch = 19,cex = 0.5,xlab = "GFP",ylab = "FSC-A",col = cls_color,main = paste("All DKO23 chimera embryos",sep = " "),xlim = lim_gfp,ylim = lim_fsc)
points(gfp_tr[!f],fsc_a[!f],cex = 0.5,pch = 19,col = cls_color[!f])
abline(a = a_l,b = b_l,lty = "dashed")
abline(a = a_h,b = b_h,lty = "dashed")
dev.off()


} else {

  png(sprintf("%s/gfp_vs_fsc_all_gated.png",fig_dir))
plot(gfp_tr,fsc_a,pch = 19,cex = 0.5,xlab = "GFP",ylab = "FSC-A",col = gate_color,main = paste("All DKO23 chimera embryos",sep = " "),xlim = lim_gfp,ylim = lim_fsc)
abline(a = a_l,b = b_l,lty = "dashed")
abline(a = a_h,b = b_h,lty = "dashed")
legend(x = "topright",pch = 19,col = c("blue","coral3","gray70"),legend = c("Host","DKO23","unclear"))
dev.off()

png(sprintf("%s/gfp_vs_fsc_all.png",fig_dir))
plot(gfp_tr,fsc_a,pch = 19,cex = 0.5,xlab = "GFP",ylab = "FSC-A",col = cls_color,main = paste("All DKO23 chimera embryos",sep = " "),xlim = lim_gfp,ylim = lim_fsc)
points(gfp_tr[!f],fsc_a[!f],cex = 0.5,pch = 19,col = cls_color[!f])
abline(a = a_l,b = b_l,lty = "dashed")
abline(a = a_h,b = b_h,lty = "dashed")
dev.off()


}



if(1) {
  df_gating = data.frame(cell = colnames(mat@mat),
                         cell_type = cls_type,
                         clone_type = cls_clone,stringsAsFactors = F)
  mat = scdb_mat(mat_nm)
  md = mat@cell_metadata
  
  f = colnames(mat@cell_metadata) %in% c("cell_type","clone_type","cell_genotype")
  md = md[,!f]
  
  md = left_join(md,df_gating,by = "cell")
  rownames(md) = md$cell
  mat@cell_metadata = md
  scdb_add_mat(id = mat_nm,mat = mat)
}






}

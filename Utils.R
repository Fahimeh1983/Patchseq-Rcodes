######################################################################################################
### Author: Fahimeh Baftizadeh #######################################################################
### Date: 10/3/2020            #######################################################################
######################################################################################################

Get_size_of_core_I1_cells <- function(results) {
  core_I1_size <- as.data.frame(results %>% 
                                  filter(Tree_call=="Core" | Tree_call=="I1") %>% 
                                  group_by(Tree_first_cl, Tree_call) %>%
                                  summarise(size = n()))
  colnames(core_I1_size) <- c("Cluster_label", " Cell_identity", "Size")
  core_I1_size
}

Get_size_of_I2_confusions <- function(results) {
  I2_size <- as.data.frame(results %>% 
                             filter(Tree_call=="I2") %>%
                             group_by(Tree_first_cl, Tree_second_cl) %>%
                             summarise(size = n()))
  
  I2_size$Tree_first_cl <- factor(I2_size$Tree_first_cl, levels = select.cl)
  I2_size$Tree_second_cl <- factor(I2_size$Tree_second_cl, levels = select.cl)
  
  I2_size <- dcast(I2_size, formula = Tree_first_cl ~ Tree_second_cl, value.var = "size", drop = FALSE)
  I2_size[is.na(I2_size)] <- 0
  rownames(I2_size) <- I2_size$Tree_first_cl
  I2_size = I2_size[select.cl, select.cl]
  
  df <- data.frame(matrix(0, ncol= length(select.cl), nrow = length(select.cl)))
  
  for (i in 1:length(select.cl)) {
    for (j in 1:length(select.cl)){
      df[i,j] = ((I2_size[i,j] + I2_size[j,i]))
    }
  }
  df[lower.tri(df)] <- 0
  rownames(df) <- select.cl
  colnames(df) <- select.cl
  df[,"id"] <- rownames(df)
  I2_confusions <- melt(df) %>% filter(value != 0) 
  colnames(I2_confusions) <- c("First_cl", "second_cl", "I2_size")
  I2_confusions
}

Get_3_best_cor <- function(memb, ref.cl, cor){
  first = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 1)
  first = Get_matrix_member(first, cor)
  second = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 2)
  second = Get_matrix_member(second, cor)
  third = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 3)
  third = Get_matrix_member(third, cor)
  cbind(first, second, third)
}

Get_3_best_KL <- function(memb, ref.cl, KLdiv){
  first = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 1)
  first = Get_matrix_member(first, KLdiv)
  second = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 2)
  second = Get_matrix_member(second, KLdiv)
  third = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 3)
  third = Get_matrix_member(third, KLdiv)
  cbind(first, second, third)
}

Get_3_best_cl <- function(memb, ref.cl){
  first <- Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 1)
  second <- Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 2)
  third <- Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 3)
  cbind(first, second, third)
}

Get_3_best_bt <- function(memb, ref.cl){
  first <- Get_nth_predicted_bt(memb = memb, ref.cl = ref.cl, nth_bt = 1)
  second <- Get_nth_predicted_bt(memb = memb, ref.cl = ref.cl, nth_bt = 2)
  third <- Get_nth_predicted_bt(memb = memb, ref.cl = ref.cl, nth_bt = 3)
  cbind(first, second, third)
}

Compute_correlation_mat <- function(markers, cells, query.dat.norm, train.cl.med){
  cormat = cor(as.matrix(query.dat.norm[markers,cells]), train.cl.med[markers,])
  cormat
}

Compute_median_gene_expression <- function(anno_file, norm.dat, markers) {
  train.cl.med <- do.call("cbind", tapply(anno_file$sample_id, anno_file$cluster_label, function(x) rowMedians(norm.dat[markers, x, drop=F])))
  rownames(train.cl.med) <- markers
  train.cl.med
}

Get_nth_predicted_cl <- function(memb, ref.cl, nth_cl){
  sorted = t(apply(memb[ , ref.cl], 1, function(x) names(x)[order(x, decreasing = TRUE)]))
  sorted[, nth_cl]
}

Get_nth_predicted_bt <- function(memb, ref.cl, nth_bt){
  sorted = t(apply(memb[ ,ref.cl], 1, function(x) x[order(x, decreasing = TRUE)]))
  sorted[, nth_bt]
}

Read_anno_cellonly_region <- function(annopath, region){
  anno <- read_feather(annopath)
  region = region
  anno <- anno %>% filter(class_label %in% c("GABAergic", "Glutamatergic") & region_label ==region)
  anno <- as.data.frame(anno)
  rownames(anno) <- anno$sample_id
  anno
}

Read_patchseq_anno <- function(path) {
  patchseq_anno <- read_feather(path)
  patchseq_anno <- as.data.frame(patchseq_anno)
  rownames(patchseq_anno) <- patchseq_anno$sample_id
  patchseq_anno
}


Get_patchseq_norm <- function(cpm_rdata_path, samp_rdata_path){
  tmp <- load(cpm_rdata_path)
  query.dat <- cpmR
  tmp <- load(samp_rdata_path)
  keepcells = which(samp.dat$Region=="VISp" & samp.dat$Type=="patch_seq")
  samp.dat = samp.dat[c(keepcells, which(samp.dat$Region=="TCx"),which(samp.dat$Region=="FCx")),]   #FCx is for Brian.  Rat samples mapped in mouse
  query.dat = query.dat[,as.character(samp.dat$exp_component_name)]
  colnames(query.dat)=as.character(samp.dat$patched_cell_container)
  query.dat.norm = log2(as.matrix(query.dat+1))
  query.dat.norm
}

Get_sample_id <- function(anno, spec_id) {
  anno[anno$spec_id_label==spec_id,"sample_id"]
}
  
compute_mapping_probability <- function(memb, ref.cl, select.cells, select.cl) {
  library(reshape2)
  ref.cl <- as.data.frame(ref.cl[select.cells])
  colnames(ref.cl) <- c("cluster_label")
  tmp <- cbind.data.frame(memb[select.cells,], ref.cl)
  
  df <-do.call("rbind",tapply(row.names(tmp),tmp$cluster_label, function(x)colMeans(tmp[x,select.cl])))
  mapping_probability <- df[select.cl, select.cl] 
  mapping_probability
}

compute_KLdiv <- function(select.cl, select.cells, memb, mapping_probability){
  library("LaplacesDemon")
  memb <- memb[select.cells,select.cl]
  mapping_probability <- mapping_probability[select.cl, select.cl]
  KLdiv <- matrix(nrow = length(select.cells), ncol = length(select.cl))
  mapping_probability[mapping_probability==0] <- 0.0000001
  mapping_probability <- mapping_probability/rowSums(mapping_probability)
  memb[memb==0] <- 0.0000001
  memb <- memb/rowSums(memb)
  
  for (cell in 1:length(select.cells)) {
    for (cl in 1:length(select.cl)){
      P <- memb[cell,select.cl]
      Q <- mapping_probability[cl, select.cl]
      KLdiv[cell, cl] <- KLD(px= as.matrix(P),py =as.matrix(Q))$sum.KLD.px.py
    }
  }
  
  rownames(KLdiv) <- rownames(memb)
  colnames(KLdiv) <- as.character(select.cl)
  KLdiv
}

Renew_list <- function(ls, ref.df, label, new.label){
  rownames(ref.df) <- ref.df[, label]
  tmp <- ref.df[as.character(ls),new.label]
  names(tmp) <- names(ls)
  new.ls <- factor(tmp)
  new.ls
}

Get_matrix_member <- function(row_col_list, mat) {
  mat_members <- list()
  for (cell in names(row_col_list)) {
    row = cell
    col = row_col_list[[cell]]
    mat_members <- c(mat_members, mat[row, col])
  }
  names(mat_members) <- names(row_col_list)
  unlist(mat_members)
}

build_confusion_matrix <- function(mat1, mat2, colname1, colname2, cells){
  mat1 <- as.data.frame(mat1)
  mat2 <- as.data.frame(mat2)
  col1 <- mat1[cells, colname1]
  col2 <- mat2[cells, colname2]
  df <- cbind.data.frame(col1, col2, deparse.level = FALSE)
  colnames(df) <- c("cl_mat1", "cl_mat2")
  df <- df %>% group_by(cl_mat1, cl_mat2) %>%
    summarise(Value=n())
  colnames(df) <- c(as.character(colname1), as.character(colname2), "Value")
  as.data.frame(df)
  
}

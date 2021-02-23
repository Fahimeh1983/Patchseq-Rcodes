######################################################################################################
### Author: Fahimeh Baftizadeh #######################################################################
### Date: 10/23/2020            ######################################################################
######################################################################################################
# Helper functions for compute Kullback Leibler divergence and also to aggregate some of mapping results


#Computing Kullback Leibler divergence 
compute_KLdiv <- function(select.cl, select.cells, memb, mapping_probability){
  library("LaplacesDemon")
  memb <- memb[select.cells, select.cl]
  mapping_probability <- mapping_probability[select.cl, select.cl]
  KLdiv <- matrix(nrow = length(select.cells), ncol = length(select.cl))
  mapping_probability[mapping_probability==0] <- 0.0000001
  mapping_probability <- mapping_probability/rowSums(mapping_probability)
  memb[memb==0] <- 0.0000001
  memb <- memb/rowSums(memb)
  
  tt = 1
  for (cell in 1:length(select.cells)) {
    print(tt)
    for (cl in 1:length(select.cl)){
      P <- memb[cell,select.cl]
      Q <- mapping_probability[cl, select.cl]
      KLdiv[cell, cl] <- KLD(px= as.numeric(P),py =as.numeric(Q))$sum.KLD.px.py
    }
    tt= tt +1
  }
  
  rownames(KLdiv) <- rownames(memb)
  colnames(KLdiv) <- as.character(select.cl)
  KLdiv
}



Read_anno_cellonly_region <- function(annopath, region){
  anno <- read_feather(annopath)
  region = region
  anno <- anno %>% filter(class_label %in% c("GABAergic", "Glutamatergic") & region_label ==region)
  anno <- as.data.frame(anno)
  rownames(anno) <- anno$sample_id
  anno
}

Compute_median_gene_expression <- function(anno_file, norm.dat, markers) {
  train.cl.med <- do.call("cbind", tapply(anno_file$sample_id, anno_file$cluster_label, function(x) rowMedians(norm.dat[markers, x, drop=F])))
  rownames(train.cl.med) <- markers
  train.cl.med
}

Compute_correlation_mat <- function(markers, cells, query.dat.norm, train.cl.med){
  cormat = cor(as.matrix(query.dat.norm[markers,cells]), train.cl.med[markers,])
  cormat
}


Get_nth_predicted_cl <- function(memb, ref.cl, nth_cl){
  sorted = t(apply(memb[ , ref.cl], 1, function(x) names(x)[order(x, decreasing = TRUE)]))
  sorted[, nth_cl]
}

Get_nth_predicted_bt <- function(memb, ref.cl, nth_bt){
  sorted = t(apply(memb[ ,ref.cl], 1, function(x) x[order(x, decreasing = TRUE)]))
  sorted[, nth_bt]
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


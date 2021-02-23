######################################################################################################
### Author: Fahimeh Baftizadeh #######################################################################
### Date: 10/23/2020           #######################################################################
######################################################################################################
#This code is to get the correct information for patchseq data


Build_ref_dendclusterid_from_facs <- function(facs.anno.path, facs.dend.path, region) {
  facs.anno <- read_feather(facs.anno.path)
  facs.dend <- readRDS(facs.dend.path)
  facs.anno <- facs.anno[facs.anno$region_label == region,]
  facs.anno <- facs.anno[facs.anno$cluster_label %in% labels(facs.dend),]
  
  REF <-  as.data.frame(unique(facs.anno[,c("cluster_label", "cluster_id", "cluster_color", 
                                            "dendcluster_id", "dendcluster_label", "dendcluster_color")]))
  
  REF$dendcluster_label == REF$cluster_label
  REF$dendcluster_color == REF$cluster_color
  
  REF <- REF %>% arrange(dendcluster_id)
  sum(cbind.data.frame(REF[,c("dendcluster_id", "dendcluster_label")], labels(facs.dend))[,2] == 
        cbind.data.frame(REF[,c("dendcluster_id", "dendcluster_label")], labels(facs.dend))[,3])
  return(list("REF" = REF, "facs.anno" = facs.anno))
}



Assign_clusterlabel_clusterid_for_patchseq <- function(patchseq.path, facs.dend.path, ref_cluster_label_color){
  
  anno <- read_feather(paste0 (patchseq.path , "/anno.feather"))
  anno <- as.data.frame(anno)
  dend <- readRDS(facs.dend.path)
  ref_cluster_label_color <- ref_cluster_label_color
  
  #modifying cluster label and id and color
  if(!sum(anno$topLeaf_label %in% ref_cluster_label_color$cluster_label) == dim(anno)[1]){
    print("ERROR, cluster labels in FACS and Patch are not similar")} else {print("ALL GOOD")}
  anno$cluster_label <- anno$topLeaf_label
  anno$dendcluster_label <- anno$cluster_label
  
  anno <- anno[,setdiff(colnames(anno), "cluster_color")]
  anno  <- left_join(anno, ref_cluster_label_color[,c("cluster_color","cluster_label")], by = "cluster_label")
  anno$dendcluster_color <- anno$cluster_color
  
  anno <- anno[,setdiff(colnames(anno), "cluster_id")]
  anno  <- left_join(anno, ref_cluster_label_color[,c("cluster_label","cluster_id")], by = "cluster_label")
  anno  <- left_join(anno, ref_cluster_label_color[,c("dendcluster_id","dendcluster_label")], by = "dendcluster_label")
  rownames(anno) <- anno$sample_id
  return(anno)
}


Modify_layer_label_GABAcells <- function(GABAanno) {
  L1 <- c("VISp1")
  all_L1 <- unique(GABAanno$structure_label)[grep(c("VIS.*1"), unique(GABAanno$structure_label))]
  L1_ho <- setdiff(all_L1, L1)
  L23 <- c("VISp2/3")
  all_L23 <-  unique(GABAanno$structure_label)[grep(c("VIS.*2/3"), unique(GABAanno$structure_label))]
  L23_ho <- setdiff(all_L23, L23)
  L4 <- c("VISp4")
  all_L4 <- unique(GABAanno$structure_label)[grep(c("VIS.*4"), unique(GABAanno$structure_label))]
  L4_ho <- setdiff(all_L4, L4)
  L5 <- c("VISp5")
  all_L5 <- unique(GABAanno$structure_label)[grep(c("VIS.*5"), unique(GABAanno$structure_label))]
  L5_ho <- setdiff(all_L5, L5)
  L6a <- c("VISp6a")
  all_L6a <- unique(GABAanno$structure_label)[grep(c("VIS.*6a"), unique(GABAanno$structure_label))]
  L6a_ho <- setdiff(all_L6a, L6a) 
  L6b <- c("VISp6b")
  all_L6b <- c(unique(GABAanno$structure_label)[grep(c("VIS.*6b"), unique(GABAanno$structure_label))])
  L6b_ho <- setdiff(all_L6b, L6b)
  
  GABAanno <- GABAanno[GABAanno$structure_label %in% 
                         c(L1, L23, L4, L5, L6a, L6b, L1_ho, L23_ho, L4_ho, L5_ho, L6a_ho, L6b_ho),] %>%
    rownames_to_column("id") %>%
    mutate(Revisited_layer_label = case_when(structure_label %in% c(L1, L1_ho)  ~ "VISp1",
                                             structure_label %in% c(L23, L23_ho)  ~ "VISp2/3",
                                             structure_label %in% c(L4, L4_ho) ~ "VISp4", 
                                             structure_label %in% c(L5, L5_ho) ~ "VISp5",
                                             structure_label %in% c(L6a, L6a_ho, L6b, L6b_ho) ~ "VISp6")) %>% 
    mutate(Revisited_layer_color = case_when(structure_label %in% c(L1, L1_ho)  ~ "#3A1799",
                                             structure_label %in% c(L23, L23_ho) ~ "#5300FF", 
                                             structure_label %in% c(L4, L4_ho) ~ "#875CCC",
                                             structure_label %in% c(L5, L5_ho) ~ "#5D2E99",
                                             structure_label %in% c(L6a, L6a_ho, L6b, L6b_ho) ~ "#9326FF")) %>%
    mutate(Revisited_layer_id = case_when(structure_label %in% c(L1, L1_ho) ~ 75,
                                          structure_label %in% c(L23, L23_ho) ~ 76, 
                                          structure_label %in% c(L4, L4_ho) ~ 77,
                                          structure_label %in% c(L5, L5_ho) ~ 78,
                                          structure_label %in% c(L6a, L6a_ho, L6b, L6b_ho) ~ 79)) %>% 
    column_to_rownames("id")
  print("Done")
  GABAanno 
}





######################################################################################################
### Author: Fahimeh Baftizadeh #######################################################################
### Date: 10/3/2020            #######################################################################
######################################################################################################

# This code to map FACs data onto FACs, this is done to create the Ref mapping probability matrix

library(dendextend)
library(matrixStats)
library(feather)
library(dplyr)
library(scrattch.io)
library(scrattch.vis)
library(feather)
library(dendextend)
library(tibble)
library(Matrix)
library(scrattch.hicat)
source(paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/",
              "mouse_patchseq_script_repository/patchseq/patchseq.new.R"))
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Patchseq-Rcodes/Utils.R")

##########################################################################################
### Setting up the paths: ################################################################
##########################################################################################

work.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/EXC_patchseq_paper_2020/"

facs.anno.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/",
                        "mouse_V1_ALM_20180520/anno.feather")

facs.data.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/",
                        "mouse_V1_ALM_20180520/data.feather")


V1.ALM.dend.label.path = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/",
                                "patch_seq/star/mouse_patchseq_script_repository/V1_ALM/",
                                "dend.labeled.rda")

select.markers.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/",
                             "Taxonomies/AIT2.3.1/select.markers.rda")

cpm.dat.path1 = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/",
                       "Mouse/facs/R_Object/SparseMatrix/20181120_RSC-004-183_mouse_star2.0_cpm_sparse.Rdata")

cpm.dat.path2 = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/",
                       "Mouse/facs/R_Object/SparseMatrix/20191028_RSC-184-228_mouse_star2.0_cpm_sparse.Rdata")


##########################################################################################
### Reading inputs: ######################################################################
##########################################################################################

load(V1.ALM.dend.label.path)
plot(dend.labeled)

load(select.markers.path)
data <- read_feather(facs.data.path)
anno <- read_feather(facs.anno.path)
anno <- anno[anno$class_label %in% c("GABAergic", "Glutamatergic"),]
data <- data[data$sample_id %in% anno$sample_id,]
data <- data %>% column_to_rownames("sample_id")
norm.dat <- log2(data + 1)
norm.dat <- t(norm.dat)

dim(anno)
dim(norm.dat)

cl.df <- unique(as.data.frame(anno[,c("cluster_id", "cluster_label", 
                                     "cluster_color", "subclass_id", 
                                     "subclass_label", "class_id", 
                                     "class_label", "cl", "dendcluster_id")]))

dim(cl.df)
cl <- set_names(anno[, "cluster_label"]$cluster_label, anno[, "sample_id"]$sample_id)
cl <- factor(cl)


FACS.cells <- colnames(norm.dat)
FACS.memb <- as.data.frame(matrix(0,nrow = length(FACS.cells), ncol = length(labels(dend.labeled))))
rownames(FACS.memb) <- FACS.cells
colnames(FACS.memb) <- labels(dend.labeled)  
FACS.cells <- sample(FACS.cells) #randomize 
All.test.cells <- c()

for (i in 1:12){
  print(c("HERE WE ARE:",i))
  begin = (i-1) * 2000 + 1
  if ( i==12) {
    end = length(FACS.cells)
  } else {
    end = i *2000
  }
  print(c(begin, end))
  test.cells <- FACS.cells[begin : end]
  
  train.cells <- setdiff(FACS.cells, test.cells)
  cl.med = get_cl_means(norm.dat[, train.cells], cl[train.cells])
  rownames(cl.med)=rownames(norm.dat[, train.cells])
  
  memb <- map_dend_membership(dend.labeled, 
                              cl= cl[train.cells], 
                              cl.med = cl.med,
                              dat= norm.dat[, train.cells], 
                              map.dat = norm.dat[,test.cells], 
                              map.cells = test.cells, 
                              mc.cores=1, 
                              bs.num=100, 
                              p=0.7, 
                              low.th=0.15)
  
  mapped.cl <- colnames(memb)[colnames(memb) %in% labels(dend.labeled)]
  FACS.memb[test.cells, mapped.cl] <- FACS.memb[test.cells, mapped.cl] + memb[test.cells, mapped.cl]
  All.test.cells <- c(All.test.cells, test.cells)
}

save(FACS.memb, file= paste0(work.dir, "/V1_ALM_FACS_ON_FACS_REF_PROB.rda"))
select.cl = labels(dend.labeled)
Tree_mapping_probability = compute_mapping_probability(memb = FACS.memb, 
                                                       select.cells = FACS.cells,  
                                                       select.cl = select.cl, 
                                                       ref.cl= cl)



Tree_mapping_probability = Tree_mapping_probability[select.cl, select.cl]

save(Tree_mapping_probability, file = paste0(work.dir, "/V1_ALM_REF_PROB_MAT.rda"))


ggplot(data = melt(Tree_mapping_probability), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()+ theme(axis.text = element_text(size=7)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Clustering lable") + ylab("NN Mapping lables") +
  scale_fill_gradient(low = "white", high = "red")

